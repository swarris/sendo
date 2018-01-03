from Bio import SeqIO
from Bio.KEGG.REST import *
from Bio.KEGG.KGML import KGML_parser
from Bio.Graphics.KGML_vis import KGMLCanvas
from Bio.Graphics.ColorSpiral import ColorSpiral
from collections import defaultdict
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

currentGene = ""
# open kegg file other
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
# process all found kegg pathways
for k in kegg:
    print("Processing: {}".format(k))
    stats[k] = defaultdict(int)
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
                    if not match == None and len(match) > 0:
                        KOToEC[element.name].extend(match[0].split(" "))
                        KOToEC[ortho].extend(match[0].split(" "))
                        
                        #print("Added {} to {}".format(match[0].split(" "),element.name.split(":")[1]))
                        #print(results)
                        #print(l)
                    
                
    
    # get information on EC numbers in kegg pathway
    for ec in kegg[k]:
        print(" EC: {}".format(ec))
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
                stats[k]["all"] += 1 
                if graphic.name[:6] in KOToEC or element.name in KOToEC:
                    if graphic.name[:6] in KOToEC:
                        name = graphic.name[:6]
                    else:
                        name = element.name
                    #print(name)
                        
                    if graphic.bgcolor not in colorCodes.values():
                        defaultBG = graphic.bgcolor
                    orgOverview = defaultdict(list)
                    orgCount = 0
                    for o in ecToGeneOrg.keys():
                        orgOverview[o] = [ec for ec in KOToEC[name] if ec in ecToGeneOrg[o]]
                        if len(orgOverview[o]) > 0:
                            stats[k][o] += 1
                            orgCount += 1
                    if orgCount > 0 and isRef:
                        graphic.bgcolor = colorCodes["in"]
                    elif isRef and orgCount == 0:
                        graphic.bgcolor = colorCodes["out"]
                    elif orgCount > 0 and not isRef:
                        graphic.bgcolor = colorCodes["notRefIn"]
                    elif orgCount == 0 and not isRef:                        
                        graphic.bgcolor = colorCodes["notRefOut"]
                        
        canvas.draw("paths/{}_{}_{}.pdf".format(k, pathway.title.replace("/","_"), keggmap))
    else:
        print("paths/{}_{}_{}.pdf exists. Skipping.".format(k, pathway.title.replace("/","_"), keggmap))
    

statsF = open("stats_{}.csv".format(keggmap), "w")
statsF.write("KEGG\tTotal\t" + "\t".join(ecToGeneOrg.keys()) + "\n")
for k in stats.keys():
    statsF.write("{}\t{}\t".format(k, stats[k]["all"]))
    statsF.write("\t".join([str(stats[k][o]) for o in ecToGeneOrg.keys()]))
    statsF.write("\n")
statsF.close()

class Gene:
    def __init__(self):
        self.kegg = set()
        self.ec = set()
        self.ko = set()
        

        
        
    
        
