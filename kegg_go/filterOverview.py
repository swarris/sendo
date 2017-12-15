import sys

toKeep = set()
outFile = open(sys.argv[3],"w")

for i in open(sys.argv[1],"r"):
    for j in i.split(','):
        j = j.strip()
        if len(j) > 0:
            toKeep.add(j)


for i in open(sys.argv[2], "r"):
    j = i.strip().split("\t")
    if len(j) > 0 and j[0] in toKeep:
        outFile.write(i)
