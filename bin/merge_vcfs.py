#! /usr/bin/env python

""" Creates an updated VCF file from the outputs from hapl;otype caller and Readbacked Phasing, to ensure genotyping fields are reported """
import sys
import re
from string import join

input1 = sys.argv[1]
input2 = sys.argv[2]

(geno,genvals,pos) = ([],[],[])

fh1 = open(input1,'r')
for lines in fh1:
    if lines.startswith("#"):
       continue
    fields = lines.rstrip("\r\n").split("\t")
    pos.append(fields[1])
    geno.append(fields[8])
    genvals.append(fields[9])
fh1.close()

fh2 = open(input2,'r')
for lines in fh2:
    lined = lines.rstrip("\r\n")
    fields = lines.rstrip("\r\n").split("\t")
    if lined.startswith("#"):
       print lined
    elif "AD" not in fields[8]:
       ind = pos.index(fields[1])
       fields[8] = geno[ind]
       fields[9] = genvals[ind]
       record = "\t".join(fields)
       print record
    else:
       print lined
fh2.close()
     
