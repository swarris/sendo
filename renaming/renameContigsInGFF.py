import sys

renaming = {}

# get ids
for l in open(sys.argv[1],"r"):
    l = l.strip().split("\t")
    renaming[l[0]] = l[1]

# rename first column:

for l in open(sys.argv[2], "r"):
    l = l.strip().split("\t")
    if len(l) > 0 and l[0] in renaming:
        l[0] = renaming[l[0]]
    print("\t".join(l))
