'''
Script for splitting a protein file in multiple files for batch processing on SLURM
Input: single protein fasta file (python3 splitFile.py protein.fasta)
'''
import sys
from Bio import SeqIO
import os

fileName = sys.argv[1]

splitAt = 100
index = 0
count = 0

if not os.path.exists("input/{}/{}".format(fileName,index)):
    os.makedirs("input/{}/{}".format(fileName,index))

if not os.path.exists("output/{}/{}".format(fileName,index)):
    os.makedirs("output/{}/{}".format(fileName,index))

outFile = open("input/{}/{}/{}.{}".format(fileName,index, fileName, index), "w")
for i in SeqIO.parse(open(fileName, "r"),"fasta"):
    if count == splitAt:
        count = 0
        index += 1
        outFile.close()
        if not os.path.exists("input/{}/{}".format(fileName,index)):
            os.makedirs("input/{}/{}".format(fileName,index))
        if not os.path.exists("output/{}/{}".format(fileName,index)):
            os.makedirs("output/{}/{}".format(fileName,index))
        outFile = open("input/{}/{}/{}.{}".format(fileName,index, fileName, index), "w")
    count += 1
    SeqIO.write(i,outFile,"fasta")

print("sbatch -a 0-{} run_interproscan.sh {}".format(index, fileName))
