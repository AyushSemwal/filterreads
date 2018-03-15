import os
import sys
import itertools
import argparse
from filterreads.fastq import Fastq

def kmerize(sequence, klen):
    '''K-mer creator'''
    mask = (1 << (klen*2)) -1
    kmer = 0
    bit_counter = 1
    seen = set()
    for nuc in sequence:
        kmer =  kmer << 2         #left shift k-kmer
        if nuc == 'A':            #add the value of the character using bitwise OR
            kmer = kmer | 0x00
        elif nuc == 'C':
            kmer = kmer | 0x01
        elif nuc == 'G':
            kmer = kmer | 0x02
        elif nuc == 'T':
            kmer = kmer | 0x03
        else:
            bit_counter = 0
        if bit_counter == klen:   #if length equals k-mer length, store k-mer
            frag = kmer & mask
            if frag in seen:
                continue
            else:
                seen.add(frag)
                yield frag
        else:
            bit_counter += 1


parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('-f', '--fastq_path', type=str, help='Path to Fastq file')
parser.add_argument('-o', '--output_path', type=str, help='Path to output file')
parser.add_argument('-k', '--klen', type=int, help='K-mer length')
args = parser.parse_args()

out_file = open(args.output_path,  'w', buffering=1048576000)
out_file.write('{0}\t{1}\n'.format("read name", '\t'.join([''.join(val) for val in itertools.product('ACGT', repeat=args.klen)])))
fastq_reader = Fastq(args.fastq_path, args.output_path, 'phred33')
for reads in fastq_reader.read():
    kmer_dict = {val:0 for val in range(0,4**args.klen)}
    kmer_reader = kmerize(reads.seq, args.klen)
    for kmers in kmer_reader:
        kmer_dict[kmers] = 1
    out_file.write('{0}\t{1}\n'.format(reads.header, '\t'.join([str(val) for val in kmer_dict.values()])))
out_file.close()
