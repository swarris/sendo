"""
Basic script to query KEGG for information on an KO number
python3 getKO.py KO01234
"""
from Bio import SeqIO
from Bio.KEGG.REST import *
from Bio.KEGG.KGML import KGML_parser
from Bio.Graphics.KGML_vis import KGMLCanvas
from Bio.Graphics.ColorSpiral import ColorSpiral
from collections import defaultdict
import sys
import os



for ecInfo in kegg_get("ko:{}".format(sys.argv[1])):
    print(ecInfo.strip())
