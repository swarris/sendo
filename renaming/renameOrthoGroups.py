import sys

renaming = {}

# get ids
for f in sys.argv[2:]:
    for l in open(f,"r"):
        l = l.strip().split("\t")
        renaming[l[0]] = l[1]

# rename first column:

renamed = open(sys.argv[1] + ".renamed.tsv", "w")
for l in open(sys.argv[1], "r"):
    if len(l) > 0:
        l = l.strip().split("\t")
        for geneList in xrange(0, len(l)):
            gene = [x.strip() for x in l[geneList].split(",")]
            for i in xrange(0,len(gene)):
                if '.1' == gene[i][-2:]:
                    geneTmp = gene[i].split('.')[0]
                    if geneTmp in renaming:
                        gene[i] = renaming[geneTmp] 
                elif gene[i] in renaming:
                    gene[i] = renaming[gene[i]]
            l[geneList] = ','.join(filter(lambda x: 'gene' not in x, gene))
        l = "\t".join(l)
        renamed.write("{}\n".format(l))
    else:
        renamed.write(l)