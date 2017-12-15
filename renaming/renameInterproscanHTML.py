import sys
import shutil

renaming = {}

# get ids
for l in open(sys.argv[1],"r"):
    l = l.strip().split("\t")
    renaming[l[0]] = l[1]


for i in renaming.keys():
    try:
        html = open(i+".html","r")
        renamed = open(renaming[i] + ".html", "w")

        for l in html:
            if 'h2 class="strapline">' or  '- InterPro' in l:
                l = l.replace(i, renaming[i])
            renamed.write(l)
        renamed.close()
        html.close()
        shutil.move(i+".html", "/tmp/"+i+".html")
    except Exception as e:
        print("Gene {} not annotated: {}".format(i,e))
