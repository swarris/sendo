from Bio import SeqIO
from Bio.KEGG.REST import *
from Bio.KEGG.KGML import KGML_parser
from Bio.Graphics.KGML_vis import KGMLCanvas
from Bio.Graphics.ColorSpiral import ColorSpiral
from collections import defaultdict
from enum import Enum
import sys
import os
import re

from pprint import pprint



"""
python3 drawPathways.py keggmap kegg.txt control.faa
keggmap: KEGG map identifier
kegg.txt: gene ids (<xref id=) + kegg ids
control.faa: protein file for IDs
"""

colorCodes = {"in":"#39818e",
              "out":"#c81837",
              "new": "#6ABEFF",
              "notRefIn":"#fcd549",
              "notRefOut": "#f77524"}

pattern = re.compile("\[EC:(.*?)\]")

keggmap = sys.argv[1]
kegg = defaultdict(list)

ecToGene = defaultdict(list)
ecToGeneOrg = {}
ecToRefGene = defaultdict(list)
KOToEC = defaultdict(list)
KOToGene = defaultdict(list)

stats = {}
orgList = []

class statKeys(Enum):
    complete="Total (pos + neg)"
    positives="Positives"
    negatives= "Negatives"
    TP="True Positives"
    TN="True negatives"
    FP="False positives"
    FN="False negatives"
    noEC="No EC found"



currentGene = ""

def ECclasses(ecNumber):
    dots= ecNumber.split(".")
    classes = []
    if len(dots) == 3:
        classes = [dots[0] + "._._._",
                   dots[0] + "." + dots[1] + "._._",
                   dots[0] + "." + dots[1] + "." + dots[2] + "._"]  
    
    return classes

# expand org list with files and extract IDs:

proteinMapping = {}
for f in sys.argv[3:]:
    orgName = f.split("/")[-1].split(".")[0]
    print(orgName)
    orgList.append(orgName)
    for l in open(f,"r"):
        if len(l) > 0 and l[0] == ">":
            proteinMapping[l[1:].split(" ")[0].strip()] = orgName
    
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
                        for c in ECclasses(ec):
                            ecToGeneOrg[org][c].append(currentGene)
# process all found kegg pathways
for k in kegg:
    print("Processing: {}".format(k))
    stats[k] = defaultdict(int)
    processedIDs = set()
    # load current pathway
    isRef = True
    try:
        pathway = KGML_parser.read(kegg_get("{}{}".format(keggmap,k), "kgml"))
    except:
        print("{} not a reference pathway.".format(k))
        pathway = KGML_parser.read(kegg_get("ko{}".format(k), "kgml"))
        isRef = False
        
        
    canvas = KGMLCanvas(pathway, import_imagemap=True)
    if isRef:
        for element in list(pathway.genes):
            f = kegg_find(keggmap, element.name)
            results = f.readline() 
            g = results.split("\t")[0]#[len(keggmap)+1:]
            
            for l in kegg_get(g).readlines():
                if "ORTHOLOGY" in l:
                    match = pattern.findall(l) # EC number
                    ortho= l[12:18]
                    #print(match) 
                    if len(match) == 0:
                        print("{} has no EC".format(ortho))
                        eName = set([x.split(":")[1] for x in element.name.split(" ")])
                        if ortho not in processedIDs:
                            stats[k][statKeys.noEC] += 1
                            stats[k][statKeys.complete] +=1 
                            processedIDs.add(ortho)
                            processedIDs |= eName 
                        elif len(eName & processedIDs) == 0:
                            stats[k][statKeys.noEC] += 1
                            stats[k][statKeys.complete] +=1 
                            processedIDs |= eName 
                            
                    if not match == None and len(match) > 0:
                        ecNumbers = match[0].split(" ")
                        for n in element.name.split(" "):
                            KOToEC[n.split(":")[1]].extend(ecNumbers)
                        KOToEC[ortho].extend(ecNumbers)
                        #for ecNumber in ecNumbers:
                        #    KOToEC[element.name].extend(ECclasses(ecNumber))
                        #    KOToEC[ortho].extend(ECclasses(ecNumber))
                        
                        #print("Added {} to {}".format(match[0].split(" "),element.name.split(":")[1]))
                        #print(results)
                        #print(l)
    
    # get information on EC numbers in kegg pathway
    for ec in kegg[k]:
        #print(" EC: {}".format(ec))
        if True:
            foundOrtho = False
            # query KEGG
            for ecInfo in kegg_get("ec:{}".format(ec)):
                ecInfoLabel = ecInfo[:12]
                #if len(ecInfoLabel.strip()) > 0: 
                #    print(ecInfo.strip())
                      
                if "ORTHOLOGY" in ecInfoLabel:
                    foundOrtho = True
                    KOToEC[ecInfo[12:18]].append(ec)
#                    KOToGene[ecInfo[12:18]].extend(ecToGene[ec])
                else:
                    foundOrtho =  foundOrtho and len(ecInfoLabel.strip()) == 0
                    if foundOrtho:
#                        KOToGene[ecInfo[12:18]].extend(ecToGene[ec])
                        KOToEC[ecInfo[12:18]].append(ec)
    #print(KOToEC)
    if not os.path.exists("paths/{}_{}_{}.pdf".format(k, pathway.title.replace("/","_"), keggmap)):
        if isRef:
            elements = list(pathway.genes) + list(pathway.orthologs)
        else:
            elements = list(pathway.orthologs)
        for element in elements:
            #print(element.name)
        #for element in pathway.genes:
            for graphic in element.graphics:
                #graphic.name = "myEC"
                #print("  " + graphic.name)

                if ":" not in graphic.name:
                    graphicName =  graphic.name.split(" ") + [x.split(":")[1] for x in element.name.split(" ")]
                else:
                    graphicName =  [x.split(":")[1] for x in graphic.name.split(" ")] + [x.split(":")[1] for x in element.name.split(" ")]
                
                graphicName = set([x.replace(".","") for x in graphicName])
                isInRef = False    
                if graphic.bgcolor == "#BFFFBF" or graphic.bgcolor == "#bfffbf":
                    if len(processedIDs & graphicName) == 0:
                        stats[k][statKeys.positives] += 1
                    isInRef = True
                elif isRef and len(processedIDs & graphicName) == 0:
                    stats[k][statKeys.negatives] += 1

                if len(processedIDs & graphicName) == 0:
                    stats[k][statKeys.complete] += 1
                orgCount = 0
                for o in ecToGeneOrg.keys():
                    for gN in graphicName:
                        orgCount += len([ec for ec in KOToEC[gN] if ec in ecToGeneOrg[o]])

                if len(list(filter(lambda x : x in KOToEC, graphicName))) > 0:
                        
                    if orgCount > 0 and isRef:
                        graphic.bgcolor = colorCodes["in"]
                        if len(processedIDs & graphicName) == 0:
                            if isInRef:
                                stats[k][statKeys.TP] += 1
                            else:
                                stats[k][statKeys.FP] += 1
                                graphic.bgcolor = colorCodes["new"]
                    elif isRef and orgCount == 0 and isInRef:
                        graphic.bgcolor = colorCodes["out"]
                        if len(processedIDs & graphicName) == 0 and isInRef:
                            stats[k][statKeys.FN] += 1
                    elif orgCount > 0 and not isRef:
                        graphic.bgcolor = colorCodes["notRefIn"]
                        if len(processedIDs & graphicName) == 0:
                            stats[k][statKeys.TP] +=1
                if orgCount == 0 and not isInRef and len(processedIDs & graphicName) == 0:
                    stats[k][statKeys.TN] += 1
                processedIDs |= graphicName
        canvas.draw("paths/{}_{}_{}.pdf".format(k, pathway.title.replace("/","_"), keggmap))
    else:
        print("paths/{}_{}_{}.pdf exists. Skipping.".format(k, pathway.title.replace("/","_"), keggmap))
    print(stats[k])

statsF = open("stats_{}.csv".format(keggmap), "w")
statsF.write("KEGG\t" + "\t".join(statKeys.__members__.keys()) + "\n")
for k in stats.keys():
    statsF.write("{}\t".format(k))
    statsF.write("\t".join([str(stats[k][s]) for s in statKeys]))
    statsF.write("\n")
statsF.close()

class Gene:
    def __init__(self):
        self.kegg = set()
        self.ec = set()
        self.ko = set()

        
        
    
        
