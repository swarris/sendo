import sys
import re

ids = {}
currentMatch = ""
count = 1
prefix = sys.argv[2]
gffOut = open(sys.argv[3], "w")
pattern = re.compile("=.*?gene.*?\.[0-9]+")

for l in open(sys.argv[1], "r"):
	match = pattern.findall(l)
	if not match == None and len(match) > 0:
		match = match[0][1:] # skip '='
		if not match == currentMatch:
			currentMatch = match
			ids[match] = prefix + "{num:05d}".format(num=count)
			print("{}\t{}".format(match, ids[match]))
			count += 1
			
	
	if match == None or len(match) == 0:
		gffOut.write(l)	
	else:
		gffOut.write(pattern.sub("=" + ids[match], l))
