import sys
import os
import difflib
import pandas as pd
import gzip
import zipfile
import itertools
from Bio import SeqIO
import gc

path1 = (sys.argv[1])
path2 = (sys.argv[2])
bpath = (sys.argv[3])
outdir = (sys.argv[4])

def hamming(str1, str2):
    return sum(map(str.__ne__, str1, str2))

inf = pd.read_table(bpath, names = ["sample", "barcodes","adapter"])
names = inf['sample'].tolist()
barcodes = inf['barcodes'].tolist()
nadaptor = len(barcodes)

def fastq_spliter(path1, path2, nadaptor, names, outdir = outdir ):
    if path1.endswith(".gz"):
        file1 = gzip.open(path1, 'rt', encoding='utf-8')
        file2 = gzip.open(path2, 'rt', encoding='utf-8')
    else:
        file1 = open(path1, 'rt', encoding='utf-8')
        file2 = open(path2, 'rt', encoding='utf-8')

    print("Start reading the file...")
    a1 = file1.readline()
    # Create list outputs
    R1store = [[] for x in range(nadaptor)]
    R2store = [[] for x in range(nadaptor)]
    # Create a list with the adaptor sequences
    seq = barcodes
    count=0
    readN=0
    while a1 != '':
        # Read File 1
        b1 = file1.readline()
        c1 = file1.readline()
        c1 = "+"
        d1 = file1.readline()
        # Read File2
        a2 = file2.readline()
        b2 = file2.readline()
        c2 = file2.readline()
        c2 = "+"
        d2 = file2.readline()

        for i in range(len(R1store)):
            ###change by yuan_jianwen@gibh.ac.cn
            #ratio = difflib.SequenceMatcher(None,seq[i],b2[3:10]).ratio()
            mismatch = hamming(seq[i],b2[3:10])
            #Ratio >0.75 allows 1 mismatch, 0.99 equal sequencesq
            if mismatch < 1:
                barcode1 = b2[:10]
                R1store[i].append(a1.strip() + "+" + barcode1[:3] + "/" + barcode1[3:])
                R1store[i].append(b1.strip())
                R1store[i].append(c1)
                R1store[i].append(d1.strip())

            # And delete the barcode from the sequence
            ### change 20221108 yuan_jianwen@gibh.ac.cn
            #    R2store[i].append(b2[10:].strip())
            #    R2store[i].append(c2)
            #    R2store[i].append(d2[10:].strip())
            #    break
                R2store[i].append(a2.strip() + "+" + barcode1[:3] + "/" + barcode1[3:])
                R2store[i].append(b2.strip())
                R2store[i].append(c2)
                R2store[i].append(d2.strip())
                break
        # Continue reading
        a1 = file1.readline()
        count += 1
        if count % 5000000 == 0:
            print("treated reads ==" + str(count))

        readN += 1
        if readN >= 10**7:
            print("saving reads and cleaning space...")
            for i in range(len(R1store)):
                thefile = open(outdir + "/" + names[i] +"_R1"+ ".fastq", 'a')
                for item1 in R1store[i]:
                    thefile.writelines(item1 + '\n')
                    R1store[i] = []
                thefile.close()

            for i in range(len(R2store)):
                thefile = open(outdir + "/" + names[i] +"_R2"+".fastq", 'a')
                for item2 in R2store[i]:
                    thefile.writelines(item2 + '\n')
                    R2store[i] = []
                thefile.close()
            readN = 0
            gc.collect()


        # Save in files
    print("final saving reads...")
    for i in range(len(R1store)):
        thefile = open(outdir + "/" + names[i] +"_R1"+ ".fastq", 'a')
        for item1 in R1store[i]:
            thefile.writelines(item1 + '\n')
        thefile.close()

    for i in range(len(R2store)):
        thefile = open(outdir + "/" + names[i] +"_R2"+".fastq", 'a')
        for item2 in R2store[i]:
            thefile.writelines(item2 + '\n')
    thefile.close()

fastq_spliter(path1,path2,nadaptor,names, outdir)
