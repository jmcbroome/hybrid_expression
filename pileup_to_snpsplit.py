#!/usr/bin/env python

#cousin to vcf_to_snpsplit, starting from an earlier point that hopefully bypasses any bullshit bcftools induces in its attempt to 'phase' obvious haplotypes

import sys

print('\t'.join(['ID', 'Chr', 'Position', 'SNP value', 'Ref/SNP']))
for entry in sys.stdin:
    if entry[0] == '#':
        continue
    spent = entry.strip().split()
    if spent[3] != '1':
        continue #ignore regions with multiple representations in the alternative
    id = 'asm_aln'
    chr = spent[0]
    pos = spent[1]
    ref = spent[2]
    val = spent[5]
    if spent[4] in 'ACGT': #should be one of these. Only keep areas with depth 1 for simplicity (1-to-1 correspondence)
        alt = spent[4]
        refsnp = ref + '/' + alt
        print("\t".join([id, chr, pos, val, refsnp]))