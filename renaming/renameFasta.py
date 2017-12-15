import sys

index = 1
label = sys.argv[2]
result = open(sys.argv[3], "w")

for l in open(sys.argv[1], "r"):
    if l[0] == ">":
        result.write(">{}{:04}\n".format(label, index))
        print("{}\t{}{:04}".format(l[1:].strip(), label, index))
        index += 1
    else:
        result.write(l)
