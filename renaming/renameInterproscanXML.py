import sys

renaming = {}

# get ids
for l in open(sys.argv[1],"r"):
    l = l.strip().split("\t")
    renaming[l[0]] = l[1]

# rename first column:

renamed = open(sys.argv[2] + ".renamed.xml", "w")
for l in open(sys.argv[2], "r"):
    if '<xref id=' in l:
        l = l.strip().split("\"")
        if len(l) > 0 and l[1] in renaming:
            l[1] = renaming[l[1]]
        renamed.write("\"".join(l) + "\n")
    else:
        renamed.write(l)
