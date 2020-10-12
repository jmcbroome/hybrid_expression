#!/usr/bin/env python

#cousin to vcf_to_snpsplit, starting from an earlier point that hopefully bypasses any bullshit bcftools induces in its attempt to 'phase' obvious haplotypes

import sys

print('\t'.join(['ID', 'Chr', 'Position', 'SNP value', 'Ref/SNP']))
for entry in sys.stdin:
    if entry[0] == '#':
        continue
    spent = entry.strip().split()
    id = 'asm_aln'
    chr = spent[0]
    pos = spent[1]
    ref = spent[2]
    val = spent[5]
    if spent[4] in 'ACGT': #should be one of these.
        alt = spent[4]
        refsnp = ref + '/' + alt
        print("\t".join([id, chr, pos, val, refsnp]))