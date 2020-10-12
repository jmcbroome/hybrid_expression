#!/usr/bin/env python3

#the purpose of this script is to take long read alignment data in bam format, parsed with pysam and using the fetch function
#then it will add a tag to each read, ps:Z:g1 or g2 or undecidable

#import
import argparse
import pysam

#define functions/classes

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbose', action = 'store_true', help = "Use to print status updates")
    parser.add_argument('-p', '--phaser', help = 'Path to a phase information file, formatted after what snpsplit uses')
    parser.add_argument('-g', '--genome_names', nargs = 2, help = 'Names to use for tagging in reference to genome one or two. Default g1, g2', default = ['g1', 'g2'])
    parser.add_argument('-b', '--bam', help = 'path to a BAM alignment file.')
    parser.add_argument('-c', '--conflict', type = float, help = 'Set to a float threshold to be the minimum that the majority genome needs to be for the decision. Default .9', default = .9)
    parser.add_argument('-u', '--update', choices = ['ref', 'alt', 'none'], help = 'Update mode- returns an updated phasing info file where either the reference or the alternative has been updated to fully match the input bam (removes conflict). DO NOT USE WITH HYBRID INPUT- SINGLE GENOME SAMPLES ONLY.', default = 'none')
    args = parser.parse_args()
    return args

def read_phaser(phase_path):
    phaser = {}
    with open(phase_path) as inf:
        for entry in inf:
            if entry[0] == 'I':
                continue
            _, chro, pos, _, refalt = entry.strip().split()
            phaser[(chro, int(pos))] = refalt.split('/')
    return phaser    

def update_phaser(phase_d, target, bampath, thresh = .9):
    print('\t'.join(['ID', 'Chr', 'Position', 'SNP value', 'Ref/SNP']))
    aln = pysam.AlignmentFile(bampath)
    for locd, ra in phase_d.items():
        chro, loc = locd
        ref, alt = ra
        try:
            for col in aln.pileup(chro, loc, loc+1, truncate = True): #iterator with only one entry- the current phaser location
                bc = {b:0 for b in 'ACGT'}
                for read in col.pileups:
                    if not read.is_refskip and not read.is_del:
                        nb = read.alignment.query_sequence[read.query_position]
                        bc[nb] += 1
                #check if the bc has high concurrence above the threshold
                ub, c = max(bc.items(), key = lambda x:x[1])
                if c > sum(bc.values()) * thresh:
                    if target == 'ref':
                        nk = ub + "/" + alt
                    elif target == 'alt':
                        nk = ref + "/" + ub
                    else:
                        print("Invalid genome targeting")
                        raise ValueError
                else:
                    #stick with the defaults if it doesn't pass the conflict threshold
                    nk = ref + "/" + alt
                print('\t'.join(['asm_aln', chro, str(loc), '~', nk]))
        except ValueError:
            #no coverage at this location, return the base
            nk = ref + "/" + alt
            print('\t'.join(['asm_aln', chro, str(loc), '~', nk]))

def tag_bam(bampath, phaser, named, thresh = .9):
    aln = pysam.AlignmentFile(bampath)
    print(aln.header)
    for r in aln:
        if r.is_unmapped or r.seq == None:
            continue
        g1 = 0
        g2 = 0
        locv = r.get_aligned_pairs()
        for rl, gl in locv:
            if rl != None and gl != None:
                pk = (r.reference_name, gl)
                if pk in phaser:
                    ref, alt = phaser[pk]
                    talt = r.seq[rl-1]
                    if talt == ref: #genome 1 is the first one to appear in the phaser file (ref/alt)
                        g1 += 1
                    elif talt == alt:
                        g2 += 1
        if g1 > (g1+g2) * thresh:
            r.set_tag(tag = 'ps', value_type = 'Z', value = named['g1'])
        elif g2 > (g1+g2) * thresh:
            r.set_tag(tag = 'ps', value_type = 'Z', value = named['g2'])
        elif g1 == 0 and g2 == 0:
            r.set_tag(tag = 'ps', value_type = 'Z', value = 'noinfo')
        else:
            r.set_tag(tag = 'ps', value_type = 'Z', value = 'conflict')
        print(r.to_string())

def main():
    args = argparser()
    ndict = {'g1':args.genome_names[0], 'g2':args.genome_names[1]}
    phaser_d = read_phaser(args.phaser)
    if args.update == 'none':
        tag_bam(args.bam, phaser_d, ndict, args.conflict)
    elif args.update == 'ref':
        update_phaser(phaser_d, 'ref', args.bam, args.conflict)
    elif args.update == 'alt':
        update_phaser(phaser_d, 'alt', args.bam, args.conflict)

if __name__ == "__main__":
    main()