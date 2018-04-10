from Bio import SeqIO
from Bio.KEGG.REST import *
from Bio.KEGG.KGML import KGML_parser
from Bio.Graphics.KGML_vis import KGMLCanvas
from Bio.Graphics.ColorSpiral import ColorSpiral
from collections import defaultdict
import sys
import os
import re

from neo4j.v1 import GraphDatabase, basic_auth

hostname = "wurnfs"
username = "neo4j"
password = "cytoscape"

driver = GraphDatabase.driver("bolt://{}".format(hostname), auth=basic_auth(username, password))
session = driver.session()

pattern = re.compile("[\A\s]GO: ([0-9]*)")

result = session.run("match(a:enzyme) return a.name as name")
for r in result:
    print(r["name"])
    for ecInfo in kegg_get(r["name"]):
        ecInfoLabel = ecInfo[:12]
        if "ORTHOLOGY" in ecInfoLabel:
            foundOrtho = True
            #KOToEC[ecInfo[12:18]].append(ec)
            for ko in kegg_get(ecInfo[12:18].strip()):
                match = pattern.findall(ko) 
                if len(match) > 0:
                    session.run("match (a:enzyme{{name: '{name}'}}), (g:GOTerm{{id: 'GO:{go}'}}) merge (a)-[r:xreference]->(g)".format(name=r["name"], go=match[0]))
session.close()


        
        
    
        
