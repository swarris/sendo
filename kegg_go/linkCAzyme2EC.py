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
import urllib.request


from neo4j.v1 import GraphDatabase, basic_auth

hostname = "wurnfs"
username = "neo4j"
password = "cytoscape"

driver = GraphDatabase.driver("bolt://{}".format(hostname), auth=basic_auth(username, password))
session = driver.session()


pattern = re.compile(".*http://www\.enzyme-database\.org/query\.php\?ec=(.*?)\"")
ecPattern = re.compile("\[(EC:.*?)\]")
result = session.run("match(a:CAzyme) return a.name as name")
for r in result:
    print(r["name"])
    """
    genes = kegg_find("genes", r["name"]).readlines()
    
    if len(genes) > 0 and len(genes[0].strip()) > 0:
        genes = genes[0].strip()
        print(genes)
        for ecInfo in kegg_get(genes):
            print(ecPattern.findall(ecInfo))
    else:
    """
    try:
        genesHTML = urllib.request.urlopen("http://www.cazy.org/{}_all.html".format(r["name"])).read().decode('utf-8')
        match = pattern.findall(genesHTML)
        if len(match) >0:
            m = match[0]
            print("match (a:CAzyme{{name:'{name}'}}), (e:enzyme) where e.name =~ 'ec:{ec}' merge (a)-[:crossConnect]->(e)".format(name=r["name"], ec=m.replace(".","\\\.").replace("*", ".*")))
            session.run("match (a:CAzyme{{name:'{name}'}}), (e:enzyme) where e.name =~ 'ec:{ec}' merge (a)-[:crossConnect]->(e)".format(name=r["name"], ec=m.replace(".","\\\.").replace("*", ".*")))
    except:
        print("Could not find name")
session.close()


        
        
    
        
