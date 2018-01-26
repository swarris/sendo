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
python3 drawPathways.py kegg.txt *.faa'
kegg.txt: gene ids (<xref id=) + kegg ids
*.faa': protein files for IDs
"""

colorCodes = {"green":"#25DE01",
              "orange":"#F8B10D",
              "blue":"#0D3CFD",
              "purple":"#8F00FF",
              "lightblue":"#3090C7",
              "red":"#FF0000",
              "silver": "#C0C0C0"
              }



group = {'ChytObl':['MB42_proteins','LEV6574_proteins'],
         'ChytCult':['Bd_JAM81_protein', 'Ppalu_CBS455_65', 'Bd_JEL423', 'Cconf_CBS675_73', 'Smicro_JEL517', 'Sp_BR117_protein','Gp_JEL478_protein', 'Hpoly_JEL142','Phirt_CBS809_83'],
         'CtrlCult':['SCE','NCR','CNE'],
         'CtrlObl':['UMA','PGR','MLR']}



"""

+ means "1 or more", - is "0"


Pile    Color        ChytObl        ChytCult    CtrlCult    CtrlObl        

1    Green        +        +        
    Orange        -        +
    Blue        +        -

2    Green        +        +        +
    Orange        -        +        +
    Red        -        -        +

3    Green        +        +                +
    Orange        -        +                -
    Blue        +        -                +

4    Green        +        +        +        +
    Orange        -        +        +        -
    Blue        +        -        -        +
    Purple        -        +        +        +
    Light Blue    +        -        -        -


"""

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

proteinMapping = {}
for f in sys.argv[2:]:
    orgName = f.split("/")[-1].split(".")[0]
    print(orgName)
    orgList.append(orgName)
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
    for g in ["grouping1","grouping2","grouping3","grouping4"]:
        for element in pathway.orthologs:
            for graphic in element.graphics:
                #reset color:
                graphic.bgcolor = "#FFFFFF"
                if hasattr(graphic,"original"):
                    graphic.name = graphic.original
                else:
                    graphic.original = graphic.name
                    
                if ":" not in graphic.name:
                    graphicName =  graphic.name.split(" ")
                else:
                    graphicName =  [x.split(":")[1] for x in graphic.name.split(" ")]
                graphicName = set([x.replace(".","") for x in graphicName])    
                

                if len(list(filter(lambda x : x in KOToEC, graphicName))) > 0:
                    # check group
                    groupOverview = defaultdict(list)
                    for gN in graphicName:
                        for org in group['ChytObl']:
                            if len([ec for ec in KOToEC[gN] if ec in ecToGeneOrg[org]]) > 0:
                                groupOverview['ChytObl'].extend([org])
                        for org in group['ChytCult']:
                            if len([ec for ec in KOToEC[gN] if ec in ecToGeneOrg[org]]) > 0:
                                groupOverview['ChytCult'].extend([org])
                        for org in group['CtrlObl']:
                            if len([ec for ec in KOToEC[gN] if ec in ecToGeneOrg[org]]) > 0:
                                groupOverview['CtrlObl'].extend([org])
                        for org in group['CtrlCult']:
                            if len([ec for ec in KOToEC[gN] if ec in ecToGeneOrg[org]]) > 0:
                                groupOverview['CtrlCult'].extend([org])

                    if g == "grouping1":
                        count = "{}, {}".format(len(groupOverview['ChytObl']),len(groupOverview['ChytCult']))
                        if len(groupOverview['ChytObl']) > 0 and len(groupOverview['ChytCult']) > 0:  
                            graphic.bgcolor = colorCodes["green"]
                        elif len(groupOverview['ChytObl']) == 0 and len(groupOverview['ChytCult']) > 0:  
                            graphic.bgcolor = colorCodes["orange"]
                        elif len(groupOverview['ChytObl']) > 0 and len(groupOverview['ChytCult']) == 0:  
                            graphic.bgcolor = colorCodes["blue"]
                    if g == "grouping2":
                        count = "{}, {}, {}".format(len(groupOverview['ChytObl']), len(groupOverview['ChytCult']), len(groupOverview['CtrlCult']))
                        if len(groupOverview['ChytObl']) > 0 and len(groupOverview['ChytCult']) > 0 and len(groupOverview['CtrlCult']) > 0:  
                            graphic.bgcolor = colorCodes["green"]
                        elif len(groupOverview['ChytObl']) == 0 and len(groupOverview['ChytCult']) > 0 and len(groupOverview['CtrlCult']) > 0:  
                            graphic.bgcolor = colorCodes["orange"]
                        elif len(groupOverview['ChytObl']) == 0 and len(groupOverview['ChytCult']) == 0 and len(groupOverview['CtrlCult']) > 0:  
                            graphic.bgcolor = colorCodes["red"]
                    if g == "grouping3":
                        count = "{}, {}, {}".format(len(groupOverview['ChytObl']),len(groupOverview['ChytCult']), len(groupOverview['CtrlObl']))
                        if len(groupOverview['ChytObl']) > 0 and len(groupOverview['ChytCult']) > 0 and len(groupOverview['CtrlObl']) > 0:  
                            graphic.bgcolor = colorCodes["green"]
                        elif len(groupOverview['ChytObl']) == 0 and len(groupOverview['ChytCult']) > 0 and len(groupOverview['CtrlObl']) == 0:  
                            graphic.bgcolor = colorCodes["orange"]
                        elif len(groupOverview['ChytObl']) > 0 and len(groupOverview['ChytCult']) == 0 and len(groupOverview['CtrlCult']) > 0:  
                            graphic.bgcolor = colorCodes["blue"]
                    if g == "grouping4":
                        count = "{}, {}, {}, {}".format(len(groupOverview['ChytObl']), len(groupOverview['ChytCult']), len(groupOverview['CtrlCult']), len(groupOverview['CtrlObl']))
                        if len(groupOverview['ChytObl']) > 0 and len(groupOverview['ChytCult']) > 0 and len(groupOverview['CtrlCult']) > 0 and len(groupOverview['CtrlObl']) > 0:  
                            graphic.bgcolor = colorCodes["green"]
                        elif len(groupOverview['ChytObl']) == 0 and len(groupOverview['ChytCult']) > 0 and len(groupOverview['CtrlCult']) > 0 and len(groupOverview['CtrlObl']) == 0:  
                            graphic.bgcolor = colorCodes["orange"]
                        elif len(groupOverview['ChytObl']) > 0 and len(groupOverview['ChytCult']) == 0 and len(groupOverview['CtrlCult']) == 0 and len(groupOverview['CtrlObl']) > 0:  
                            graphic.bgcolor = colorCodes["blue"]
                        elif len(groupOverview['ChytObl']) == 0 and len(groupOverview['ChytCult']) > 0 and len(groupOverview['CtrlCult']) > 0 and len(groupOverview['CtrlObl']) > 0:  
                            graphic.bgcolor = colorCodes["purple"]
                        elif len(groupOverview['ChytObl']) > 0 and len(groupOverview['ChytCult']) == 0 and len(groupOverview['CtrlCult']) == 0 and len(groupOverview['CtrlObl']) == 0:  
                            graphic.bgcolor = colorCodes["lightblue"]
                        else:
                            graphic.bgcolor = colorCodes["silver"]
                    graphic.name = count
        if g == "grouping1":
            groupName =  "ChytObl_ChytCult"
        if g == "grouping2":
            groupName =  "ChytObl_ChytCult_CtrlCult"
        if g == "grouping3":
            groupName =  "ChytObl_ChytCult_CtrlObj"
        if g == "grouping4":
            groupName =  "ChytObl_ChytCult_CtrlCult_CtrlObj"
        
        canvas.draw("paths_compare/{}_{}_{}.pdf".format(k, pathway.title.replace("/","_"), groupName))


        
        
    
        
