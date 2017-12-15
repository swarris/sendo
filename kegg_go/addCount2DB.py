import sys
from collections import defaultdict
from neo4j.v1 import GraphDatabase, basic_auth

hostname = "phi-host.ab.wurnet.nl"
username = "neo4j"
password = "Neo4J"

driver = GraphDatabase.driver("bolt://{}".format(hostname), auth=basic_auth(username, password))
session = driver.session()

currentGene = ""
processed = set()
for f in sys.argv[1:]:
    session.run("match (a:GOTerm) set a.fold = 0, a.up = 0, a.down = 0, a.neutral = 0, a.allUp = 0, a.allDown=0, a.allNeutral = 0, a.averageFold = 0, a.allFold = 0")
    for l in open(f, "r"):
        if "<xref id=" in l:
            currentGene = l.split('"')[1]
        elif "<go-xref" in l and "<model" not in l:
            l = l.split('"')
            if len(l) >=8:
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
        resultFile = open(level + "_" + f, "w")
        for go in goTerms:
            result = session.run("Match (a:GOTerm {id:'"+go+"'})<-[:ISA*"+ level +"]-(c) where c.up is not null and c.up > 0 with distinct c as n return n.up, n.name, n.id, n.namespace order by n.up desc")
            #print(result)
            for r in result:
                if topLine:
                    resultFile.write("\t".join(r) + "\n")
                    topLine = False
                resultFile.write("\t".join([str(x) for x in r.values()]) + "\n")
            

        resultFile.close()
