"""
Links enzymes to compounds and GO terms.
python3 linkKEGG2GO.py 
"""

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

pattern = re.compile("[\A\s]GO: ([0-9\s]+)")
compound = re.compile("\[CPD:(C[0-9]*)\]")
result = session.run("match(a:enzyme) return a.name as name")
for r in result:
    print(r["name"])
    productFound = False
    substrateFound = False
    for ecInfo in kegg_get(r["name"]):
        ecInfoLabel = ecInfo[:12]
        if "ORTHOLOGY" in ecInfoLabel:
            foundOrtho = True
            #KOToEC[ecInfo[12:18]].append(ec)
            for ko in kegg_get(ecInfo[12:18].strip()):
                match = pattern.findall(ko) 
                if len(match) > 0:
                    match = ko.strip().split(" ")[1:]
                    for m in match:
                        print(m)
                        session.run("match (a:enzyme{{name: '{name}'}}), (g:GOTerm{{id: 'GO:{go}'}}) merge (a)-[r:crossConnect]->(g)".format(name=r["name"], go=m.strip()))
        if "PRODUCT" in ecInfoLabel or productFound:
            productFound = True
            match = compound.findall(ecInfo) 
            if len(match) > 0:
                for m in match:
                    session.run("match (a:enzyme{{name: '{name}'}}), (c:compound{{name: 'cpd:{compound}'}}) merge (a)-[r:produces]->(c)".format(name=r["name"], compound=m))
            else:
                productFound = False
            if ';' not in ecInfo:
                productFound = False
                    
        if "SUBSTRATE" in ecInfoLabel or substrateFound:
            substrateFound = True
            match = compound.findall(ecInfo) 
            if len(match) > 0:
                for m in match:
                    session.run("match (a:enzyme{{name: '{name}'}}), (c:compound{{name: 'cpd:{compound}'}}) merge (a)<-[r:usedby]-(c)".format(name=r["name"], compound=m))
            else:
                substrateFound = False
            if ';' not in ecInfo:
                substrateFound = False
session.close()


        
        
    
        
