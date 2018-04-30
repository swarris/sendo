"""
Sets per group of species statistics on the GO terms in a Neo4j database and reports these in a tsv file on level 1 and 2 of GO.
python3 goGrouping.py interproscan.xml *.fasta
"""
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

group = {'MB42_proteins':'ChytObl','LEV6574_proteins':'ChytObl',
         'Bd_JAM81_protein':'ChytCult', 'Ppalu_CBS455_65':'ChytCult', 'Bd_JEL423':'ChytCult', 'Cconf_CBS675_73':'ChytCult', 'Smicro_JEL517':'ChytCult', 'Sp_BR117_protein':'ChytCult','Gp_JEL478_protein':'ChytCult', 'Hpoly_JEL142':'ChytCult','Phirt_CBS809_83':'ChytCult',
         'SCE':'CtrlCult','NCR':'CtrlCult','CNE':'CtrlCult',
         'UMA':'CtrlObl','PGR':'CtrlObl','MLR':'CtrlObl'}

session.run("match (a:GOTerm) set a.ChytObl = 0, a.ChytCult = 0, a.CtrlCult = 0, a.CtrlObl = 0")
session.run("match (a:GOTerm) set a.allChytObl = 0, a.allChytCult = 0, a.allCtrlCult = 0, a.allCtrlObl = 0")


for f in sys.argv[2:]:
    proteinMapping = {}
    orgName = f.split("/")[-1].split(".")[0]
    print(orgName)
    for l in open(f,"r"):
        if len(l) > 0 and l[0] == ">":
            proteinMapping[l[1:].split(" ")[0].strip()] = orgName
    #reset data
    for l in open(sys.argv[1], "r"):
        if "<xref id=" in l:
            currentGene = l.split('"')[1]
        elif "<go-xref" in l and "<model" not in l:
            l = l.split('"')
            if len(l) >=8 and currentGene in proteinMapping:
                goID = l[5]
                goCat = l[1]
                goTerm = l[7]
                line = "\t".join([currentGene, goCat, goTerm, goID])
                print(goID)
                if line not in processed:
                    session.run("match (a:GOTerm {{id:'{id}'}}) SET a.{group} = a.{group} + 1".format(id=goID, group=group[orgName]))
                    processed.add(line)


for groupName in set(group.values()): 
    session.run("match(a:GOTerm) with collect(a) as allGo unwind allGo as goTerm match (goTerm)<-[r:ISA*]-(b:GOTerm) with collect(distinct b) as allB,goTerm set goTerm.all{group} = goTerm.{group} + reduce(all{group} = 0, n IN allB| all{group} + n.{group})".format(group=groupName))
    print("Toplevel nodes")
    #biological_process    GO:0008150
    #cellular_component    GO:0005575
    #molecular_function    GO:0003674
    goTerms = ['GO:0008150','GO:0005575','GO:0003674']
    for level in  ["1","2"]:
        topLine = True
        resultFile = open(level + "_" + groupName + ".tsv", "w")
        for go in goTerms:
            print("Level {} on term {}".format(level, go))
            result = session.run("Match (a:GOTerm {id:'"+go+"'})<-[:ISA*"+ level +"]-(c) where c.all{group} is not null and c.all{group} > 0 with distinct c as n return n.all{group}, n.name, n.id, n.namespace order by n.all{group} desc".format(group=groupName))
            #print(result)
            for r in result:
                if topLine:
                    resultFile.write("\t".join(r) + "\n")
                    topLine = False
                resultFile.write("\t".join([str(x) for x in r.values()]) + "\n")
            

        resultFile.close()

    
