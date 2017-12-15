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
python3 drawPathways.py ../renaming/MB42_renamed_genes.csv ../renaming/LEV6574_renamed_genes.csv kegg.txt keggLEV.txt keggOthers.txt
renaming: only include genes which are valid (which means they are in the renaming file)
keggOthers: id will be split with '_' and then first 4 characters
"""

colorCodes = {"9-":"#FD0D27",
              "8-":"#FA0D41",
              "7-":"#F70D5B",
              "6-":"#F50D75",
              "5-":"#F20D8F",
              "4-":"#EF0DA9",
              "3-":"#ED0DC3",
              "2-":"#EA0DDD",
              "1-":"#E80DF8",
              "all":"#0D3CFD",
              "8+":"#25DE01",
              "7+":"#43D803",
              "6+":"#61D104",
              "5+":"#7FCB06",
              "4+":"#9DC408",
              "3+":"#BBBE09",
              "2+":"#D9B70B",
              "1+":"#F8B10D"}

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

currentGene = ""
defaultBG = "#FFFFFF"
# open kegg file MB42 
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
for l in open(sys.argv[5],"r"):
    if "<xref id=" in l:
        currentGene = l.split('"')[1]
    else:
        keggID = l.split('"')[3].split("+")
        if len(keggID) > 1:
            kegg[keggID[0]].extend([x for x in keggID[1:] if x not in kegg[keggID[0]]])
            org = currentGene.split("_")[0][:4]
            if "KXS" in org:
                org = "KXS"
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

                    #colorCodes = {"9-sendo": '#FD0D27',"#1-sendo": '#E80DF8', "all" : "#0D3CFD", "8+sendo": "#07E500", "1+sendo":"#F8B10D"}
                    if len(mb42) > 0 and len(lev) > 0: # sendo found
                        if orgCount == 9:
                            graphic.bgcolor = colorCodes["all"]
                            graphic.name = "all"
                        elif orgCount > 0:
                            graphic.bgcolor = colorCodes[str(orgCount) + "+"]
                            graphic.name = str(orgCount) + "+"
                        else:
                            graphic.bgcolor = "#000000"
                            graphic.fgcolor = "#FFFFFF"
                            graphic.name = "only"
                            
                    else: # no sendo
                        if orgCount > 0:
                            graphic.bgcolor = colorCodes[str(orgCount)+"-"]
                            graphic.name = str(orgCount) + "-"
                        
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

    
        
