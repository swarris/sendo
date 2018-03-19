import sys
from collections import defaultdict
from neo4j.v1 import GraphDatabase, basic_auth

hostname = "wurnfs"
username = "neo4j"
password = "cytoscape"

driver = GraphDatabase.driver("bolt://{}".format(hostname), auth=basic_auth(username, password))
session = driver.session()

allLabels = ["allHpoly_JEL142",
  "allUMA",
  "allBd_JAM81_protein",
  "CtrlObl",
  "allSmicro_JEL517",
  "allLEV6574_proteins",
  "ChytCult",
  "allPGR",
  "MB42_proteins",
  "MLR",
  "allCNE",
  "allNCR",
  "allCtrlCult",
  "id",
  "allSp_BR117_protein",
  "allSCE",
  "group",
  "LEV6574_proteins",
  "allgroup",
  "ChytObl",
  "Bd_JEL423",
  "Gp_JEL478_protein",
  "allPpalu_CBS455_65",
  "allGp_JEL478_protein",
  "allBd_JEL423",
  "Bd_JAM81_protein",
  "UMA",
  "allPhirt_CBS809_83",
  "Smicro_JEL517",
  "allChytCult",
  "allMB42_proteins",
  "allCconf_CBS675_73",
  "Phirt_CBS809_83",
  "allChytObl",
  "allCtrlObl",
  "SCE",
  "allMLR",
  "namespace",
  "name",
  "CtrlCult",
  "Hpoly_JEL142",
  "Sp_BR117_protein",
  "PGR",
  "Ppalu_CBS455_65",
  "Cconf_CBS675_73",
  "CNE",
  "NCR",
  "sendo_specific",
  "sendo_absent",
  "chytrid_core",
  "allsendo_specific",
  "allsendo_absent",
  "allchytrid_core"]

queryLabels = [  "id",
  "name",
  "namespace",
  "allHpoly_JEL142",
  "allUMA",
  "allBd_JAM81_protein",
  "allSmicro_JEL517",
  "allLEV6574_proteins",
  "allPGR",
  "allCNE",
  "allNCR",
  "allSp_BR117_protein",
  "allSCE",
  "allPpalu_CBS455_65",
  "allGp_JEL478_protein",
  "allBd_JEL423",
  "allPhirt_CBS809_83",
  "allMB42_proteins",
  "allCconf_CBS675_73",
  "allMLR",
  "allgroup",
  "allCtrlCult",
  "allChytCult",
  "allChytObl",
  "allCtrlObl",
  "allsendo_specific",
  "allsendo_absent",
  "allchytrid_core"]

print("Toplevel nodes")
#biological_process    GO:0008150
#cellular_component    GO:0005575
#molecular_function    GO:0003674
goTerms = ['GO:0008150','GO:0005575','GO:0003674']
for level in  ["1","2"]:
    resultFile = open(level + "_complete.tsv", "w")
    for go in goTerms:      
        result = session.run("Match (a:GOTerm {id:'"+go+"'})<-[:ISA*"+ level +"]-(c) where exists(c.allgroup) and (c.allChytObl or c.allChytCult or c.allCtrlCult or c.allCtrlObl) return distinct c")
        allKeys = None
        for c in result:
            if allKeys == None:
                allKeys = queryLabels
                resultFile.write("\t".join(allKeys) + "\n")
            resultFile.write("\t".join([str(c['c'].properties[x]) for x in allKeys]) + "\n")


    resultFile.close()



session.close()
