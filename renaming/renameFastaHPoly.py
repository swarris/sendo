import sys

result = open(sys.argv[2], "w")
index = 1
for l in open(sys.argv[1], "r"):
    if l[0] == ">":
        label = "__{:04}__".format(index) + l[1:].split("_")[2].strip()
        result.write(">{}\n".format(label))
        print("{}\t{}".format(l[1:].strip(), label))
        index += 1
    else:
        result.write(l)
