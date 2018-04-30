"""
Adds the GO basic file to the Neo4j database.
python3 goTerms.py go-basic.obo
"""
import sys
from neo4j.v1 import GraphDatabase, basic_auth
from Bio import SeqIO
from pprint import pprint
import re

hostname = "wurnfs"
username = "neo4j"
password = "cytoscape"

driver = GraphDatabase.driver("bolt://{}".format(hostname), auth=basic_auth(username, password))
session = driver.session()

class GO:
    def __init__(self,file):
        self.is_a = []
        try :
            l = file.next().strip()
            #print(l)
            while ( "[Term]" not in l):
                l = l.split(":")
                if len(l) > 2:
                    l[1] = ':'.join(l[1:])
                if len(l) > 1:
                    if l[0] in ["id", "name", "namespace"]:
                        setattr(self, l[0], l[1].strip())
                    elif l[0] == "is_a":
                        self.is_a.append(l[1].split("!")[0].strip())
                l = file.next().strip()
                #print(l)
                
        except Exception as e:
            raise e
        
    def store(self,session):
        if hasattr(self, 'id') and "GO" in self.id:
            newName= re.sub('[^0-9a-zA-Z]+', '_', self.name)
            cytoscapeLabel = re.sub('_', ' ', newName)
            session.run("create (a:GOTerm{{id:'{}', name:'{}', namespace:'{}', cytoscape:'{}'}})".format(self.id, newName, self.namespace, cytoscapeLabel))
            #print("create (a:GOTerm{{id:'{}', name:'{}', namespace:'{}'}})".format(self.id, self.name.replace("'","^"), self.namespace))

go = []
goFile = iter(open(sys.argv[1],"r").readlines())
try :
    
    l = goFile.next().strip()
    print(l)
    while (l != "[Term]"):
        l = goFile.next().strip()
    while True:
        goTerm = GO(goFile)
        go.append(goTerm)
        goTerm.store(session)
        
except Exception as e:
    print(e)
    print("Done")
    
for t in go:
    for link in t.is_a:
        #print("match (a:GOTerm), (b:GOTerm) where a.id = '{}' and b.id = '{}' create unique (a)-[r:ISA]->(b)".format(t.id, link))
        session.run("match (a:GOTerm), (b:GOTerm) where a.id = '{}' and b.id = '{}' create unique (a)-[r:ISA]->(b)".format(t.id, link))



# MATCH p=()-[r:Edge10x*]->() with p,length(p) as lP return p order by lP DESC limit 25
        
