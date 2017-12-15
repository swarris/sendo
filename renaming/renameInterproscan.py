import sys

renaming = {}

# get ids
for l in open(sys.argv[1],"r"):
    l = l.strip().split("\t")
    renaming[l[0]] = l[1]

# rename gff 
renamed = open(sys.argv[2] + ".renamed.gff3", "w")
for l in open(sys.argv[2], "r"):
    l = l.strip().split("\t")
    if len(l) > 0 and l[0] in renaming:
        currentTarget = "Target={}".format(l[0]) 
        renamedTarget = "Target={}".format(renaming[l[0]])
        l[-1] = l[-1].replace(currentTarget, renamedTarget)
        l[0] = renaming[l[0]]
    if l[0][0] != "#":
        renamed.write("\t".join(l) + "\n")
