from Bio import SeqIO
from Bio.KEGG.REST import *
from Bio.KEGG.KGML import KGML_parser
from Bio.Graphics.KGML_vis import KGMLCanvas
from Bio.Graphics.ColorSpiral import ColorSpiral
from collections import defaultdict
import sys
import os

from pprint import pprint

"""
python3 drawPathways.py ../renaming/MB42_renamed_genes.csv ../renaming/LEV6574_renamed_genes.csv kegg.txt keggLEV.txt keggOthers.txt *.faa
renaming: only include genes which are valid (which means they are in the renaming file)
keggOthers: 9 other species
*.faa: protein files for IDs
"""

colorCodes = {"inAll":"#25DE01",
              "notInSendo":"#F8B10D",
              "onlyInSendo":"#0D3CFD"}

renamingMB42 = {}
renamingLEV = {}

# get ids MB42
for l in open(sys.argv[1],"r"):
    l = l.strip().split("\t")
    renamingMB42[l[0] + "MB42"] = l[1]

# get ids LEV
for l in open(sys.argv[2],"r"):
    l = l.strip().split("\t")
    renamingLEV[l[0] + "LEV"] = l[1]

    
kegg = defaultdict(list)

ecToGene = defaultdict(list)
ecToGeneMB42 = defaultdict(list) 
ecToGeneLEV = defaultdict(list) 
ecToGeneOrg = {}
KOToEC = defaultdict(list)
KOToGene = defaultdict(list)

stats = {}
orgList = []

currentGene = ""
defaultBG = "#FFFFFF"
# open kegg file MB42 
orgList.append("mb42")
for l in open(sys.argv[3],"r"):
    if "<xref id=" in l:
        currentGene = l.split('"')[1] + "MB42"
    else:
        keggID = l.split('"')[3].split("+")
        if len(keggID) > 1 and currentGene in renamingMB42: # add only if valid gene after renaming
            kegg[keggID[0]].extend([x for x in keggID[1:] if x not in kegg[keggID[0]]])
            for ec in keggID[1:]:
                #print(ec)
                if currentGene not in ecToGeneMB42[ec]:
                    ecToGeneMB42[ec].append(currentGene)
                    ecToGene[ec].append(currentGene)
currentGene = ""
# open kegg file LEV 
orgList.append("lev")
for l in open(sys.argv[4],"r"):
    if "<xref id=" in l:
        currentGene = l.split('"')[1] + "LEV"
    else:
        keggID = l.split('"')[3].split("+")
        if len(keggID) > 1 and currentGene in renamingLEV: # add only if valid gene after renaming
            kegg[keggID[0]].extend([x for x in keggID[1:] if x not in kegg[keggID[0]]])
            for ec in keggID[1:]:
                if currentGene not in ecToGeneLEV[ec]:
                    ecToGeneLEV[ec].append(currentGene)
                if currentGene not in ecToGene[ec]:
                    ecToGene[ec].append(currentGene)

currentGene = ""
# open kegg file other
# expand org list with files and extract IDs:
proteinMapping = {}
for f in sys.argv[6:]:
    orgName = f.split("/")[-1].split(".")[0]
    print(orgName)
    orgList.append(orgName)
    for l in open(f,"r"):
        if len(l) > 0 and l[0] == ">":
            proteinMapping[l[1:].split(" ")[0].strip()] = orgName
    
for l in open(sys.argv[5],"r"):
    if "<xref id=" in l:
        currentGene = l.split('"')[1]
    else:
        keggID = l.split('"')[3].split("+")
        if len(keggID) > 1:
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
    pathway = KGML_parser.read(kegg_get("ko{}".format(k), "kgml"))
    canvas = KGMLCanvas(pathway, import_imagemap=True)
    
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
    # check all species:
    if not os.path.exists("paths/{}_{}_ALL.pdf".format(k, pathway.title.replace("/","_"))):
        for element in pathway.orthologs:
            for graphic in element.graphics:
                #graphic.name = "myEC"
                stats[k]["all"] += 1 
                if graphic.name[:6] in KOToEC:
                    if graphic.bgcolor not in colorCodes.values():
                        defaultBG = graphic.bgcolor
                    mb42 = [ec for ec in KOToEC[graphic.name[:6]] if ec in ecToGeneMB42]
                    lev = [ec for ec in KOToEC[graphic.name[:6]] if ec in ecToGeneLEV]
                    orgOverview = defaultdict(list)
                    orgCount = 0
                    for o in ecToGeneOrg.keys():
                        orgOverview[o] = [ec for ec in KOToEC[graphic.name[:6]] if ec in ecToGeneOrg[o]]
                        if len(orgOverview[o]) > 0:
                            stats[k][o] += 1
                            orgCount += 1
                            
                    if len(mb42) > 0:
                        stats[k]["mb42"] += 1
                    if len(lev) > 0:
                        stats[k]["lev"] += 1

                    #colorCodes = {"inAll":"#25DE01", "notInSendo":"#F8B10D", "onlyInSendo":"#0D3CFD"}

                    if len(mb42) > 0 or len(lev) > 0: # sendo found
                        if len(mb42) > 0 and len(lev) > 0:
                            graphic.name = "S {}; C {}".format("2", orgCount)
                        elif len(mb42) > 0:
                            graphic.name = "S {}; C {}".format("1", orgCount)
                        else:
                            graphic.name = "S {}; C {}".format("1", orgCount)
                        if orgCount > 0:
                            graphic.bgcolor = colorCodes["inAll"]
                        else:
                            graphic.bgcolor = colorCodes["onlyInSendo"]
                            
                    else: # no sendo
                        graphic.name = "S {}; C {}".format(0, orgCount) 
                        graphic.bgcolor = colorCodes["notInSendo"]
                        
        canvas.draw("paths/{}_{}_ALL.pdf".format(k, pathway.title.replace("/","_")))
    else:
        print("paths/{}_{}_ALL.pdf exists. Skipping.".format(k, pathway.title.replace("/","_")))
"""
    # check both species:
    if not os.path.exists("paths/{}_{}_BOTH.pdf".format(k, pathway.title.replace("/","_"))):
        for element in pathway.orthologs:
            for graphic in element.graphics:
                #graphic.name = "myEC"
                stats[k][0] += 1 
                if graphic.name[:6] in KOToEC:
                    mb42 = [ec for ec in KOToEC[graphic.name[:6]] if ec in ecToGeneMB42]
                    lev = [ec for ec in KOToEC[graphic.name[:6]] if ec in ecToGeneLEV]
                    
                    if len(mb42) == 0:
                        graphic.bgcolor = "#FFFF00"
                    elif len(lev) == 0:
                        graphic.bgcolor = "#FF00FF"
                    else:
                        graphic.bgcolor = "#00FF00"
                else:
                    graphic.bgcolor = defaultBG
                    
        canvas.draw("paths/{}_{}_BOTH.pdf".format(k, pathway.title.replace("/","_")))
    else:
        print("paths/{}_{}_BOTH.pdf exists. Skipping.".format(k, pathway.title.replace("/","_")))

    # MB42
    if not os.path.exists("paths/{}_{}_MB42.pdf".format(k, pathway.title.replace("/","_"))):
        for element in pathway.orthologs:
            for graphic in element.graphics:
                #graphic.name = "myEC"
                if graphic.name[:6] in KOToEC:
                    mb42 = [ec for ec in KOToEC[graphic.name[:6]] if ec in ecToGeneMB42]
                    if len(mb42) > 0:
                        graphic.bgcolor = "#FF00FF"
                    else:
                        graphic.bgcolor = defaultBG
                else:
                    graphic.bgcolor = defaultBG
        canvas.draw("paths/{}_{}_MB42.pdf".format(k, pathway.title.replace("/","_")))
    else:
        print("paths/{}_{}_MB42.pdf exists. Skipping.".format(k, pathway.title.replace("/","_")))

    # LEV
    if not os.path.exists("paths/{}_{}_LEV6574.pdf".format(k, pathway.title.replace("/","_"))):
        for element in pathway.orthologs:
            for graphic in element.graphics:
                #graphic.name = "myEC"
                if graphic.name[:6] in KOToEC:
                    lev = [ec for ec in KOToEC[graphic.name[:6]] if ec in ecToGeneLEV]

                    if len(lev) > 0:
                        graphic.bgcolor = "#FFFF00"
                    else:
                        graphic.bgcolor = defaultBG
                else:
                    graphic.bgcolor = defaultBG
        canvas.draw("paths/{}_{}_LEV6574.pdf".format(k, pathway.title.replace("/","_")))
    else:
        print("paths/{}_{}_LEV6574.pdf exists. Skipping.".format(k, pathway.title.replace("/","_")))
"""

statsF = open("stats.csv", "w")
statsF.write("KEGG\tTotal\tMB52\tLEV\t" + "\t".join(ecToGeneOrg.keys()) + "\n")
for k in stats.keys():
    statsF.write("{}\t{}\t{}\t{}\t".format(k, stats[k]["all"], stats[k]["mb42"], stats[k]["lev"]))
    statsF.write("\t".join([str(stats[k][o]) for o in ecToGeneOrg.keys()]))
    statsF.write("\n")
statsF.close()

class Gene:
    def __init__(self):
        self.kegg = set()
        self.ec = set()
        self.ko = set()
        
for o in orgList:
    overviewF = open("overview_{}.csv".format(o), "w")
    selectedEcToGene = None
    geneData = defaultdict(Gene)
    if o == "lev":
        selectedEcToGene = ecToGeneLEV
    elif o == "mb42":
        selectedEcToGene = ecToGeneMB42
    else:
        selectedEcToGene = ecToGeneOrg[o]
    for k in kegg:
        for ec in kegg[k]:
            for gene in selectedEcToGene[ec]:
                if o == "lev" and gene in renamingLEV:
                    geneName = renamingLEV[gene]
                elif o == "mb42" and gene in renamingMB42:
                    geneName = renamingMB42[gene]
                else:
                    geneName = gene
                    
                geneData[geneName].kegg.add(k)
                geneData[geneName].ec.add(ec)
                for ko in KOToEC:
                    if ec in KOToEC[ko]:
                        geneData[geneName].ko.add(ko)
    overviewF.write("KEGG\n")
    
    for gene in geneData:
        for ko in geneData[gene].ko:
            for ec in geneData[gene].ec:
                for k in geneData[gene].kegg:
                    overviewF.write("{}\n".format('\t'.join([gene, ko, ec, k])))

        
        
    
        
