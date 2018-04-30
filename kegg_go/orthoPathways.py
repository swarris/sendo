"""
Processes the Interproscan output for KEGG pathway information. Gets all relevant pathways from KEGG and creates for this orthogroup a TSV with statistics.
python3 orthoPathways.py orthoGroups.tsv Interproscan.xml *.fasta
"""

from Bio import SeqIO
from Bio.KEGG.REST import *
from Bio.KEGG.KGML import KGML_parser
from Bio.Graphics.KGML_vis import KGMLCanvas
from Bio.Graphics.ColorSpiral import ColorSpiral
from collections import defaultdict
import sys
import os

from pprint import pprint

orthoGroup = sys.argv[1]

kegg = defaultdict(list)

ecToGene = defaultdict(list)
ecToGeneOrg = {}
KOToEC = defaultdict(list)
KOToGene = defaultdict(list)

stats = {}
orgList = []

currentGene = ""
# open kegg file other
# expand org list with files and extract IDs:

orthoData = ''.join(open(orthoGroup, "r").readlines())
proteinMapping = {}
for f in sys.argv[3:]:
    orgName = f.split("/")[-1].split(".")[0]
    print(orgName)
    orgList.append(orgName)
    for l in open(f,"r"):
        if len(l) > 0 and l[0] == ">":
            proteinName = l[1:].split(" ")[0].strip()
            if proteinName in orthoData:
                proteinMapping[proteinName] = orgName
    
for l in open(sys.argv[2],"r"):
    if "<xref id=" in l:
        currentGene = l.split('"')[1]
    else:
        keggID = l.split('"')[3].split("+")
        if len(keggID) > 1:
            if currentGene in proteinMapping: # skip those which are not found in the protein files
                kegg[keggID[0]].extend([x for x in keggID[1:] if x not in kegg[keggID[0]]])
                org = proteinMapping[currentGene]
                for ec in keggID[1:]:
                    if org not in ecToGeneOrg:
                        ecToGeneOrg[org] = defaultdict(list)
                    if currentGene not in ecToGeneOrg[org][ec]:
                        ecToGeneOrg[org][ec].append(currentGene)

# process all found kegg pathways
for k in kegg:
    print("Processing: {}".format(k))
    stats[k] = defaultdict(int)
    processedIDs = set()
    # load current pathway
    pathway = KGML_parser.read(kegg_get("ko{}".format(k), "kgml"))
    
    # get information on EC numbers in kegg pathway
    for ec in kegg[k]:
        print(" EC: {}".format(ec))
        if True:
            foundOrtho = False
            # query KEGG
            for ecInfo in kegg_get("ec:{}".format(ec)):
                ecInfoLabel = ecInfo[:12]
                if "ORTHOLOGY" in ecInfoLabel:
                    foundOrtho = True
                    KOToEC[ecInfo[12:18]].append(ec)
#                    KOToGene[ecInfo[12:18]].extend(ecToGene[ec])
                else:
                    foundOrtho =  foundOrtho and len(ecInfoLabel.strip()) == 0
                    if foundOrtho:
#                        KOToGene[ecInfo[12:18]].extend(ecToGene[ec])
                        KOToEC[ecInfo[12:18]].append(ec)
    for element in pathway.orthologs:
        for graphic in element.graphics:
            #graphic.name = "myEC"
            
            if ":" not in graphic.name:
                graphicName =  graphic.name.split(" ")
            else:
                graphicName =  [x.split(":")[1] for x in graphic.name.split(" ")]
            graphicName = set([x.replace(".","") for x in graphicName])    
            
            if len(processedIDs & graphicName) == 0:
                stats[k]["all"] += 1 


            if len(list(filter(lambda x : x in KOToEC, graphicName))) > 0:
                orgOverview = defaultdict(list)
                orgCount = 0
                for o in ecToGeneOrg.keys():
                    orgOverview[o] = []
                    
                    for gN in graphicName:
                        orgOverview[o].extend([ec for ec in KOToEC[gN] if ec in ecToGeneOrg[o]])
                    if len(orgOverview[o]) > 0:
                        if len(processedIDs & graphicName) == 0:
                            stats[k][o] += 1

            processedIDs |= graphicName
    

statsF = open("stats_{}".format(orthoGroup), "w")
statsF.write("KEGG\tTotal\t" + "\t".join(ecToGeneOrg.keys()) + "\n")
for k in stats.keys():
    statsF.write("{}\t{}\t".format(k, stats[k]["all"]))
    statsF.write("\t".join([str(stats[k][o]) for o in ecToGeneOrg.keys()]))
    statsF.write("\n")
statsF.close()

