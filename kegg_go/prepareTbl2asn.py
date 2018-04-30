"""
Annotates a species for tbl2asn use prior to submission.
python3 prepareTbl2asn.py keggIDs.csv output/goIDs.csv interproscan.tsv species.faa
"""
from Bio import SeqIO
from Bio.KEGG.REST import *
from Bio.KEGG.KGML import KGML_parser
from Bio.Graphics.KGML_vis import KGMLCanvas
from Bio.Graphics.ColorSpiral import ColorSpiral
from collections import defaultdict
import sys
import os
import re

proteinMapping = {}

ecToGene = defaultdict(list)
ecToGeneOrg = {}
pfam = defaultdict(set)
goTerms = defaultdict(set)
tsvINfo = {}
shortName = sys.argv[5]
geneNumber = re.compile(".*_g([0-9]+)")

f = sys.argv[4]
orgName = f.split("/")[-1].split(".")[0]
print(orgName)
for l in open(f,"r"):
    if len(l) > 0 and l[0] == ">":
        proteinMapping[l[1:].split(" ")[0].strip()] = orgName


for l in open(sys.argv[1],"r"):
    if "<xref id=" in l:
        currentGene = l.split('"')[1]
    else:
        keggID = l.split('"')[3].split("+")
        if len(keggID) > 1:
            if currentGene in proteinMapping: # skip those which are not found in the protein files
                org = proteinMapping[currentGene]
                for ec in keggID[1:]:
                    if org not in ecToGeneOrg:
                        ecToGeneOrg[org] = defaultdict(list)
                    if ec not in ecToGeneOrg[org][currentGene]:
                        ecToGeneOrg[org][currentGene].append(ec)

for l in open(sys.argv[2], "r"):
    if "<xref id=" in l:
        currentGene = l.split('"')[1]
    elif "<go-xref" in l and "<model" not in l:
        l = l.split('"')
        if len(l) >=8 and currentGene in proteinMapping:
            goID = l[5]
            goCat = l[1]
            goTerm = l[7]
            goTerms[currentGene].add(goID)


for l in open(sys.argv[3], "r"):
    l = l.strip().split('\t')
    currentGene = l[0]
    if currentGene in proteinMapping and l[3] == "Pfam":
        pfam[currentGene].add(l[4])

                        
for o in ecToGeneOrg:
    outfile = open("annotation_{}.tsv".format(o), "w")
    print(o)
    for gene in ecToGeneOrg[o]:
        for pf in pfam[gene]:
            outfile.write("{}\tDbxref\tPFAM:{}\n".format(gene, pf))
        currentNames = set()
        if len(ecToGeneOrg[o][gene]) == 1:
            # only proteins with single EC are processed.
            ec = ecToGeneOrg[o][gene][0]
            for ecInfo in kegg_get("ec:{}".format(ec)):
                #print("{}".format(ecInfo.strip()))
                ecInfoLabel = ecInfo[:12]
                if "NAME" in ecInfoLabel[0:4]:
                    outfile.write("{}\t{}\t{}\n".format(gene, "product", ecInfo[12:].strip().strip(";")))
                elif "SCE: " in ecInfo:
                    # get name
                    for geneName in re.findall("\((.*?)\)", ecInfo):
                        currentNames.add(geneName)
            for go in goTerms[gene]:
                outfile.write("{}\t{}\t{}\n".format(gene, "Ontology_term", go))
        if len(currentNames) == 1:
            outfile.write("{}\t{}\t{}\n".format(gene, "name", currentNames.pop()))
        else:
            shortNameGene =geneNumber.findall(gene)
            if len(shortNameGene) > 0:
                outfile.write("{}\t{}\t{}\n".format(gene, "name", "{}{}".format(shortName, str(int(shortNameGene[0])))))
            
    outfile.close()
            
        
