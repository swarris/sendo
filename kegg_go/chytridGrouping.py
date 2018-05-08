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


session.run("match (a:GOTerm) set a.chytrids ='{}'".format("unclassified"))
session.run("match (a:GOTerm) where not a.ChytObl and a.ChytCult set a.chytrids ='only culturable'")
session.run("match (a:GOTerm) where a.ChytObl and not a.ChytCult set a.chytrids ='only obligate biotrophic'")
session.run("match (a:GOTerm) where a.ChytObl and a.ChytCult set a.chytrids ='both'")

session.run("match (a:enzyme) set a.chytrids ='{}'".format("unclassified"))
session.run("match (a:enzyme) where not a.ChytObl and a.ChytCult set a.chytrids ='only culturable'")
session.run("match (a:enzyme) where a.ChytObl and not a.ChytCult set a.chytrids ='only obligate biotrophic'")
session.run("match (a:enzyme) where a.ChytObl and a.ChytCult set a.chytrids ='both'")

session.close()
