from collections import defaultdict
import sys
import os

"""
python3 go2excel.py ../renaming/MB42_renamed_genes.csv ../renaming/LEV6574_renamed_genes.csv goMB42.csv goLEV.csv others.csv
"""
renamingMB42 = {}
renamingLEV = {}
processed = set()

# get ids MB42
for l in open(sys.argv[1],"r"):
    l = l.strip().split("\t")
    renamingMB42[l[0] + "MB42"] = l[1]

# get ids LEV
for l in open(sys.argv[2],"r"):
    l = l.strip().split("\t")
    renamingLEV[l[0] + "LEV"] = l[1]


currentGene = ""
print("\t".join(["Gene", "Category", "Term", "ID"]))

for l in open(sys.argv[3],"r"):
    if "<xref id=" in l:
        currentGene = l.split('"')[1] + "MB42"
    elif "<go-xref" in l:
        l = l.split('"')
        goID = l[5].split(":")[1]
        goCat = l[1]
        goTerm = l[7]
        if currentGene in renamingMB42:
            line = "\t".join([renamingMB42[currentGene], goCat, goTerm, goID])
            if line not in processed: # add only if valid gene after renaming
                print(line)
                processed.add(line)
            

for l in open(sys.argv[4],"r"):
    if "<xref id=" in l:
        currentGene = l.split('"')[1] + "LEV"
    elif "<go-xref" in l:
        l = l.split('"')
        if len(l) >=8:
            goID = l[5].split(":")[1]
            goCat = l[1]
            goTerm = l[7]

            if currentGene in renamingLEV:
                line = "\t".join([renamingLEV[currentGene], goCat, goTerm, goID])
                if line not in processed: # add only if valid gene after renaming
                    print(line)
                    processed.add(line)

for l in open(sys.argv[5],"r"):
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
