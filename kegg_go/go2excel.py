from collections import defaultdict
import sys
import os

"""
python3 go2excel.py go.csv
"""
processed = set()



currentGene = ""
print("\t".join(["Gene", "Category", "Term", "ID"]))


for l in open(sys.argv[2],"r"):
    if "<xref id=" in l:
        currentGene = l.split('"')[1]
    elif "<go-xref" in l and "<model" not in l:
        l = l.split('"')
        if len(l) >=8:
            goID = l[5].split(":")[1]
            goCat = l[1]
            goTerm = l[7]
            line = "\t".join([currentGene, goCat, goTerm, goID])
            if line not in processed:
                print(line)
                processed.add(line)
