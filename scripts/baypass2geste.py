#!/usr/bin/env python

# modified from: https://github.com/TheDBStern/poolseq_utils/blob/master/baypass2bayescan.py (26.07.2021)
# added poolsize as mandatory argument and use it to get effective size (see Feder et al., 2012)
# both genobaypass and poolsize files are produced by poolfstats's pooldata2genobaypass

import glob, os
import argparse
import numpy

parser = argparse.ArgumentParser(description='Script to convert a file from BayPass format to BayeScan format (created by poolfstat)')
parser.add_argument('-i', dest = 'input', type = str, required=True,  help = 'input baypass format file')
parser.add_argument('-o', dest = 'output', type = str, required=True,  help = 'name of output')
parser.add_argument('-p', dest = 'pops', type = int, required=True,  help = 'number of pools')
parser.add_argument('-s', dest = 'size', type = str, required=True,  help = 'haploid size of pools (poolsize file from poolfstat for BayPass)')

args = parser.parse_args()

## func to determine length of file
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1
flen = file_len(args.input)

#generate header lines
head = open('header.txt','w')
head.write('[loci]=%s\n\n[populations]=%s\n\n'%(str(flen),str(args.pops)))

#create tmp file for each pop
for i in range(1,args.pops+1):
    output = open('tmp'+str(i),'w')
    output.write('[pop]=%s\n'%(i))

with open(args.input,'rU') as f1, open(args.size,'rU') as f2:
    poolsizes = f2.readline()
    lcount = 1
    for line in f1:
        pcount = 1
        for pop in range(0,args.pops*2,2):
            c1 = int(line.split(' ')[pop]) # get first allele count for each pop
            c2 = int(line.split(' ')[pop+1]) # second allele count
            tot = c1 + c2
            poolsize = int(poolsizes.split(' ')[pop/2]) # read haploid poolsize
            if tot == 0 : # set this sample's site to 0 in case it was 0
                c1scaled = 0
                c2scaled = 0
            else : # scale read count to allele numbers using haploid poolsize
                #c1scaled = round((c1 / float(tot)) * poolsize,0)
                #c2scaled = round((c2 / float(tot)) * poolsize,0)
                if c1 == 0 :
                    c1scaled = 0
                else :
                    c1scaled = round((c1 * poolsize - 1) / float(c1 + poolsize)) # Feder et al., 2012
                if c2 == 0 :
                    c2scaled = 0
                else :
                    c2scaled = round((c2 * poolsize - 1) / float(c2 + poolsize)) # Feder et al., 2012
            totscaled = c1scaled + c2scaled # scaled total
            output = open('tmp'+str(pcount),'a')
            #output.write('%s %s 2 %s %s\n'%(lcount,tot,c1,c2))
            output.write('{:>7} {:>3.0f} 2 {:>3.0f} {:>3.0f}\n'.format(lcount,totscaled,c1scaled,c2scaled))
            pcount +=1
        lcount +=1 
