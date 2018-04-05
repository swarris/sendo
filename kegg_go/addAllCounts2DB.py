import sys
from collections import defaultdict
from neo4j.v1 import GraphDatabase, basic_auth

hostname = "wurnfs"
username = "neo4j"
password = "cytoscape"

driver = GraphDatabase.driver("bolt://{}".format(hostname), auth=basic_auth(username, password))
session = driver.session()

colorCodes = {"green":"#25DE01",
              "orange":"#F8B10D",
              "blue":"#0D3CFD",
              "purple":"#8F00FF",
              "lightblue":"#3090C7",
              "red":"#FF0000",
              "silver": "#C0C0C0"
              }


currentGene = ""
processed = set()
orgNames = set()
for f in sys.argv[2:]:
    proteinMapping = {}
    orgName = f.split("/")[-1].split(".")[0]
    orgNames.add(orgName)
    #print(orgName)
    for l in open(f,"r"):
        if len(l) > 0 and l[0] == ">":
            proteinMapping[l[1:].split(" ")[0].strip()] = orgName
    #reset data
    session.run("create index on :GOTerm({})".format(orgName))
    session.run("match (a:GOTerm) set a.{orgName} = 0".format(orgName = orgName))
    session.run("match (a:GOTerm) set a.all{orgName} = 0".format(orgName = orgName))
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
                    session.run("match (a:GOTerm {{id:'{goID}'}}) SET a.{orgName} = a.{orgName} + 1".format(goID = goID, orgName = orgName))
    session.run("match(a:GOTerm) with collect(a) as allGo unwind allGo as goTerm match (goTerm)<-[r:ISA*]-(b:GOTerm) with collect(distinct b) as allB,goTerm set goTerm.all{orgName} = goTerm.{orgName} + reduce(allUp = 0, n IN allB| allUp + n.{orgName})".format(orgName = orgName))


session.run("match (a:GOTerm) set a.ChytObl = a.MB42_proteins>0 or a.LEV6574_proteins > 0")
session.run("match (a:GOTerm) set a.ChytCult = a.Bd_JAM81_protein> 0 or a.Ppalu_CBS455_65> 0 or a.Bd_JEL423> 0 or a.Cconf_CBS675_73> 0 or a.Smicro_JEL517> 0 or a.Sp_BR117_protein > 0 or a.Gp_JEL478_protein > 0 or a.Hpoly_JEL142 > 0 or a.Phirt_CBS809_83 > 0")
session.run("match (a:GOTerm) set a.CtrlCult = a.SCE> 0 or a.NCR> 0 or a.CNE> 0")
session.run("match (a:GOTerm) set a.CtrlObl = a.UMA> 0 or a.PGR> 0 or a.MLR> 0")


session.run("match (a:GOTerm) set a.group ='{}'".format("unclassified"))
session.run("match (a:GOTerm) where a.ChytObl and a.ChytCult and a.CtrlCult and a.CtrlObl set a.group ='{}'".format("all"))
session.run("match (a:GOTerm) where not a.ChytObl and a.ChytCult and a.CtrlCult and not a.CtrlObl set a.group ='{}'".format("onlyCult"))
session.run("match (a:GOTerm) where a.ChytObl and not a.ChytCult and not a.CtrlCult and a.CtrlObl set a.group ='{}'".format("onlyObl"))
session.run("match (a:GOTerm) where not a.ChytObl and a.ChytCult and a.CtrlCult and a.CtrlObl set a.group ='{}'".format("noChytObl"))
session.run("match (a:GOTerm) where a.ChytObl and not a.ChytCult and not a.CtrlCult and not a.CtrlObl set a.group ='{}'".format("onlyChytObl"))
session.run("match (a:GOTerm) where a.ChytObl and a.ChytCult and not a.CtrlCult and not a.CtrlObl set a.group ='{}'".format("onlyChyt"))
session.run("match (a:GOTerm) where not a.ChytObl and not a.ChytCult and a.CtrlCult and a.CtrlObl set a.group ='{}'".format("higherFungi"))


session.run("match (a:GOTerm) set a.allChytObl = a.allMB42_proteins>0 or a.allLEV6574_proteins > 0")
session.run("match (a:GOTerm) set a.allChytCult = a.allBd_JAM81_protein> 0 or a.allPpalu_CBS455_65> 0 or a.allBd_JEL423> 0 or a.allCconf_CBS675_73> 0 or a.allSmicro_JEL517> 0 or a.allSp_BR117_protein > 0 or a.allGp_JEL478_protein > 0 or a.allHpoly_JEL142 > 0 or a.allPhirt_CBS809_83 > 0")
session.run("match (a:GOTerm) set a.allCtrlCult = a.allSCE> 0 or a.allNCR> 0 or a.allCNE> 0")
session.run("match (a:GOTerm) set a.allCtrlObl = a.allUMA> 0 or a.allPGR> 0 or a.allMLR> 0")

session.run("match (a:GOTerm) set a.allgroup ='{}'".format("unclassified"))
session.run("match (a:GOTerm) where a.allChytObl and a.allChytCult and a.allCtrlCult and a.allCtrlObl set a.allgroup ='{}'".format("all"))
session.run("match (a:GOTerm) where not a.allChytObl and a.allChytCult and a.allCtrlCult and not a.allCtrlObl set a.allgroup ='{}'".format("onlyCult"))
session.run("match (a:GOTerm) where a.allChytObl and not a.allChytCult and not a.allCtrlCult and a.allCtrlObl set a.allgroup ='{}'".format("onlyObl"))
session.run("match (a:GOTerm) where not a.allChytObl and a.allChytCult and a.allCtrlCult and a.allCtrlObl set a.allgroup ='{}'".format("noChytObl"))
session.run("match (a:GOTerm) where a.allChytObl and not a.allChytCult and not a.allCtrlCult and not a.allCtrlObl set a.allgroup ='{}'".format("onlyChytObl"))
session.run("match (a:GOTerm) where a.allChytObl and a.allChytCult and not a.allCtrlCult and not a.allCtrlObl set a.allgroup ='{}'".format("onlyChyt"))
session.run("match (a:GOTerm) where not a.allChytObl and not a.allChytCult and a.allCtrlCult and a.allCtrlObl set a.allgroup ='{}'".format("higherFungi"))

print("Toplevel nodes")
#biological_process    GO:0008150
#cellular_component    GO:0005575
#molecular_function    GO:0003674
goTerms = ['GO:0008150','GO:0005575','GO:0003674']
for level in  ["1","2"]:
    resultFile = open(level + "_complete.tsv", "w")
    for go in goTerms:      
        result = session.run("Match (a:GOTerm {id:'"+go+"'})<-[:ISA*"+ level +"]-(c) where exists(c.allgroup) and (c.allChytObl or c.allChytCult or c.allCtrlCult or c.allCtrlObl) return distinct c")
                #print(result)                                                                                                                                                                                          
        allKeys = None
        for c in result:
            if allKeys == None:
                allKeys = c['c'].properties.keys()
                resultFile.write("\t".join(allKeys) + "\n")
            resultFile.write("\t".join([str(c['c'].properties[x]) for x in allKeys]) + "\n")


    resultFile.close()



session.close()
