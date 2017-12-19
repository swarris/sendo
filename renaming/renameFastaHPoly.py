import sys

result = open(sys.argv[2], "w")

for l in open(sys.argv[1], "r"):
    if l[0] == ">":
        label = l[1:].split("_")[2].strip()
        result.write(">{}\n".format(label))
        print("{}\t{}".format(l[1:].strip(), label))
    else:
        result.write(l)
