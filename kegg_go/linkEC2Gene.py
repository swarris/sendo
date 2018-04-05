from Bio import SeqIO
from Bio.KEGG.REST import *
from Bio.KEGG.KGML import KGML_parser
from Bio.Graphics.KGML_vis import KGMLCanvas
from Bio.Graphics.ColorSpiral import ColorSpiral
from collections import defaultdict
import sys
import os

from pprint import pprint

"""
python3 drawPathways.py kegg.txt ecFile *.faa'
kegg.txt: gene ids (<xref id=) + kegg ids
ecFile: sheet with EC numbers to be found
*.faa': protein files for IDs
"""

ecToGene = defaultdict(list)
ecToGeneOrg = {}
ecFilter = set()

stats = {}
orgList = []

currentGene = ""
# open kegg file other
# expand org list with files and extract IDs:

proteinMapping = {}
proteinSet = {}
proteinOut = open(sys.argv[2]+".faa", "w")

for f in sys.argv[3:]:
    orgName = f.split("/")[-1].split(".")[0]
    orgList.append(orgName)
    for l in SeqIO.parse(open(f,"r"), "fasta"):
        proteinMapping[l.id.split(" ")[0]] = orgName
        proteinSet[l.id.split(" ")[0]] = l

for l in open(sys.argv[2],"r"):
    l = l.strip().split()
    for e in l:
        if "EC:" in e:
            ecFilter.add(e)
            
        
for l in open(sys.argv[1],"r"):
    if "<xref id=" in l:
        currentGene = l.split('"')[1]
    else:
        keggID = l.split('"')[3].split("+")
        if len(keggID) > 1:
            if currentGene in proteinMapping: # skip those which are not found in the protein files
                org = proteinMapping[currentGene]
                for ec in keggID[1:]:
                    ec = "EC:" + ec
                    if ec in ecFilter:
                        if org not in ecToGeneOrg:
                            ecToGeneOrg[org] = defaultdict(list)
                        if currentGene not in ecToGeneOrg[org][ec]:
                            ecToGeneOrg[org][ec].append(currentGene)

for o in ecToGeneOrg:
    print(o)
    for e in ecToGeneOrg[o]:
        print("{}\t{}".format(e, "\t".join(ecToGeneOrg[o][e])))
        for g in ecToGeneOrg[o][e]:
            SeqIO.write(proteinSet[g], proteinOut, "fasta")
        


        
        
    
        
