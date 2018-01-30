import sys

fileName = sys.argv[1]

count = 0
outFilePE1 = open(fileName + ".PE_1.fastq", "w")
outFilePE2 = open(fileName + ".PE_2.fastq", "w")
outFile = outFilePE1
PE = "/1"
readCount = 0
for i in open(fileName, "r"):
    
    if count > 0 and count % 4 == 0:
        if PE == "/1":
            outFile = outFilePE2
            PE = "/2"
        else:
            outFile = outFilePE1
            PE = "/1"
    else:
        if PE == "/1":
            outFile = outFilePE1
        else:
            outFile = outFilePE2

    count += 1
    if "No name" in i:
        if PE == "/1":
                readCount += 1
        outFile.write("@Read_" + str(readCount) + PE + "\n")
    else:
        outFile.write(i)
