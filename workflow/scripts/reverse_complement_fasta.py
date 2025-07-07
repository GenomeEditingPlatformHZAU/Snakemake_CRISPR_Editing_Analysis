#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse

# DNA碱基互补表
complement = str.maketrans("ACGTacgt", "TGCAtgca")

def reverse_complement_fasta(infile, outfile):
    with open(infile, 'r') as fin, open(outfile, 'w') as fout:
        header = None
        sequence = []

        for line in fin:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    # 输出上一条记录
                    rc_seq = "".join(sequence)[::-1].translate(complement)
                    fout.write(f"{header}\t{rc_seq}\n")
                header = line[1:]  # 去掉开头的 ">"
                sequence = []
            else:
                sequence.append(line)

        # 最后一条记录
        if header is not None:
            rc_seq = "".join(sequence)[::-1].translate(complement)
            fout.write(f"{header}\t{rc_seq}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Reverse complement FASTA and output tab-separated header and sequence.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Output file (TSV format)")

    args = parser.parse_args()
    reverse_complement_fasta(args.input, args.output)
