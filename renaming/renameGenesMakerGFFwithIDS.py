import sys
import re

ids = {}
currentMatch = ""
ids = sys.argv[1]
gffOut = open(sys.argv[3], "w")
pattern = re.compile("=.*?gene.*?\.[0-9]+")

renaming = {}

# get ids
for l in open(ids,"r"):
    l = l.strip().split("\t")
    renaming[l[0]] = l[1]


for l in open(sys.argv[2], "r"):
	match = pattern.findall(l)
	if not match == None and len(match) > 0:
		match = match[0][1:] # skip '='
		if match in renaming:
                        gffOut.write(pattern.sub("=" + renaming[match], l))			
                else:
                        gffOut.write(l)	
        else:
                gffOut.write(l)	

