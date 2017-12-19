import sys

renaming = {}

# get ids
for l in open(sys.argv[1],"r"):
    l = l.strip().split("\t")
    renaming[l[0]] = l[1]

# rename first column:

renamed = open(sys.argv[2] + ".renamed.faa", "w")
for l in open(sys.argv[2], "r"):
    if len(l) > 0 and l[0] == ">":
        l = l.strip()[1:]
        if len(l) > 0 and l in renaming:
            l= renaming[l]
            renamed.write(">{}\n".format(l))
    else:
        renamed.write(l)