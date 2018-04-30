"""
This script outputs a .csv per protein file with gene identifiers and GO terms found in Interproscan results.
python3 go2excel.py interproscan.xml proteins*.fasta
"""

from collections import defaultdict
import sys
import os
from _mysql import result

processed = set()



currentGene = ""


for f in sys.argv[2:]:
    proteinMapping = {}
    orgName = f.split("/")[-1].split(".")[0]
    print(orgName)
    resultFile = open("goterms_" + orgName + ".tsv", "w")
    resultFile.write("\t".join(["Gene", "Category", "Term", "ID"]) + "\n") 
    for l in open(f,"r"):
        if len(l) > 0 and l[0] == ">":
            proteinMapping[l[1:].split(" ")[0].strip()] = orgName
    for l in open(sys.argv[1], "r"):
        if "<xref id=" in l:
            currentGene = l.split('"')[1]
        elif "<go-xref" in l and "<model" not in l:
            l = l.split('"')
            if len(l) >=8 and currentGene in proteinMapping:
                goID = l[5].split(":")[1]
                goCat = l[1]
                goTerm = l[7]
                line = "\t".join([currentGene, goCat, goTerm, goID])
                if line not in processed:
                    resultFile.write(line+"\n")
                    processed.add(line)
