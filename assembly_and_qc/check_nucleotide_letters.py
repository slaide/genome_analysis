#!/usr/bin/env python3

import sys
import os
import re
    
# check that my code does what it should
if False:
    thing="ATCGACTACCATACTACGACTACACCATACTATACACTCATNNXXNXN"
    print(re.sub('[ACGT\r\n]','',thing))

with open("dna_assembly_combined/out_as/combined.contigs.fa") as file:
    lines=file.readlines()
    # every second line, beginning at the second line
    lines=lines[1::2]

    # this should be a line containing a nucleotide sequence, which it is
    if False:
        print(lines[0])

    print("number of lines: {}",len(lines))

    weirdlines=0

    for i,line in enumerate(lines):
        line=re.sub('[ACGT\r\n]','',line)
        if(len(line)>0):
            print(line)
            weirdlines+=1

        if i%1000==0:
            print("done {}/{}".format(i,len(lines)))

    print("lines containing letters other than ACTG: {}".format(weirdlines))