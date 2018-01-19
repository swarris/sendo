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
for f in sys.argv[2:]:
    proteinMapping = {}
    orgName = f.split("/")[-1].split(".")[0]
    print(orgName)
    for l in open(f,"r"):
        if len(l) > 0 and l[0] == ">":
            proteinMapping[l[1:].split(" ")[0].strip()] = orgName
    #reset data
    session.run("match (a:GOTerm) set a.fold = 0, a.up = 0, a.down = 0, a.neutral = 0, a.allUp = 0, a.allDown=0, a.allNeutral = 0, a.averageFold = 0, a.allFold = 0")
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
                    session.run("match (a:GOTerm {{id:'{}'}}) SET a.up = a.up + 1".format(goID))



    session.run("match(a:GOTerm) with collect(a) as allGo unwind allGo as goTerm match (goTerm)<-[r:ISA*]-(b:GOTerm) with collect(distinct b) as allB,goTerm set goTerm.allUp = reduce(allUp = 0, n IN allB| allUp + n.up)")

    
    print("Toplevel nodes")
    #biological_process    GO:0008150
    #cellular_component    GO:0005575
    #molecular_function    GO:0003674
    goTerms = ['GO:0008150','GO:0005575','GO:0003674']
    for level in  ["1","2"]:
        topLine = True
        resultFile = open(level + "_" + orgName + ".tsv", "w")
        for go in goTerms:
            print("Level {} on term {}".format(level, go))
            result = session.run("Match (a:GOTerm {id:'"+go+"'})<-[:ISA*"+ level +"]-(c) where c.up is not null and c.up > 0 with distinct c as n return n.up, n.name, n.id, n.namespace order by n.up desc")
            #print(result)
            for r in result:
                if topLine:
                    resultFile.write("\t".join(r) + "\n")
                    topLine = False
                resultFile.write("\t".join([str(x) for x in r.values()]) + "\n")
            

        resultFile.close()
