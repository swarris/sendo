'''
This script processes the interproscan XML output for GO term annotations and adds them to the neo4j database.
Please change hostname, username and password. 
Change levels for more or less details (line 89)
Input: interproscan.xml [proteinFiles.fasta]+ (python3 addAllCounts2DB interproscan.xml proteinFiles*.fasta)
Output: statistics on each protein file on level 1 and 2 GO terms  
'''

import sys
from collections import defaultdict
from neo4j.v1 import GraphDatabase, basic_auth

hostname = "wurnfs"
username = "neo4j"
password = "cytoscape"

driver = GraphDatabase.driver("bolt://{}".format(hostname), auth=basic_auth(username, password))
session = driver.session()


currentGene = ""
processed = set()
orgNames = set()
cazymes = defaultdict(list)

for f in sys.argv[2:]:
    proteinMapping = {}
    orgName = f.split("/")[-1].split(".")[0]
    orgNames.add(orgName)
    print(orgName)
    for l in open(f,"r"):
        if len(l) > 0 and l[0] == ">":
            proteinMapping[l[1:].split(" ")[0].strip()] = orgName
    #reset data
    
    session.run("create index on :CAzyme({})".format(orgName))
    for l in open(sys.argv[1], "r"):
        currentGene = l.split('\t')[0]
        l = l.split('\t')
        if currentGene in proteinMapping:
            # Query    Subject    E-value    Subject Start    Subject End    Query Start    Query End    Subject Covered fraction
            cazyme = l[1]
            session.run("merge (a:CAzyme {{name:'{cazyme}'}})".format(cazyme = cazyme))
            cazymes[cazyme].append(orgName)
            
for o in orgNames:
    session.run("match (a:CAzyme) set a.{orgName} = 0".format(orgName = o))

for c,org in cazymes.items():
    for o in org:
        session.run("match (a:CAzyme {{name:'{cazyme}'}}) SET a.{orgName} = a.{orgName} + 1".format(cazyme = c, orgName = o))
    
session.run("match (a:CAzyme) set a.ChytObl = a.MB42_proteins>0 or a.LEV6574_proteins > 0")
session.run("match (a:CAzyme) set a.ChytCult = a.Bd_JAM81_protein> 0 or a.Ppalu_CBS455_65> 0 or a.Bd_JEL423> 0 or a.Cconf_CBS675_73> 0 or a.Smicro_JEL517> 0 or a.Sp_BR117_protein > 0 or a.Gp_JEL478_protein > 0 or a.Hpoly_JEL142 > 0 or a.Phirt_CBS809_83 > 0")
session.run("match (a:CAzyme) set a.CtrlCult = a.SCE> 0 or a.NCR> 0 or a.CNE> 0")
session.run("match (a:CAzyme) set a.CtrlObl = a.UMA> 0 or a.PGR> 0 or a.MLR> 0")


session.run("match (a:CAzyme) set a.group ='{}'".format("unclassified"))
session.run("match (a:CAzyme) where a.ChytObl or a.ChytCult or a.CtrlCult or a.CtrlObl set a.group ='{}'".format("other"))
session.run("match (a:CAzyme) where a.ChytObl and a.ChytCult and a.CtrlCult and a.CtrlObl set a.group ='{}'".format("all"))
session.run("match (a:CAzyme) where not a.ChytObl and a.ChytCult and a.CtrlCult and not a.CtrlObl set a.group ='{}'".format("onlyCult"))
session.run("match (a:CAzyme) where a.ChytObl and not a.ChytCult and not a.CtrlCult and a.CtrlObl set a.group ='{}'".format("onlyObl"))
session.run("match (a:CAzyme) where not a.ChytObl and a.ChytCult and a.CtrlCult and a.CtrlObl set a.group ='{}'".format("noChytObl"))
session.run("match (a:CAzyme) where a.ChytObl and not a.ChytCult and not a.CtrlCult and not a.CtrlObl set a.group ='{}'".format("onlyChytObl"))
session.run("match (a:CAzyme) where a.ChytObl and a.ChytCult and not a.CtrlCult and not a.CtrlObl set a.group ='{}'".format("onlyChyt"))
session.run("match (a:CAzyme) where not a.ChytObl and not a.ChytCult and a.CtrlCult and a.CtrlObl set a.group ='{}'".format("higherFungi"))


session.close()
