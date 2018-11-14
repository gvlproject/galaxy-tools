#!/usr/bin/env python

import sys
import argparse

def main():

    VERSION = 0.1

    parser = argparse.ArgumentParser(description="Calculates the means of rows of a numeric tab delimited table.")
    parser.add_argument("-i", "--infile", help="Input tab delimted numeric file")
    parser.add_argument("-o", "--outfile", help="Output list of means")
    parser.add_argument("-s", "--skip", default=0, help="Number of header rows to skip", type=int)
    parser.add_argument("--version", action='store_true')

    args = parser.parse_args()

    if args.version:
        print("means_by_row.py version: %.1f" % VERSION)
        return

    infile = str(args.infile)
    outfile = str(args.outfile)
    skip = int(args.skip)

    lines = open( infile, 'r').readlines()

    averages = []
    linenum = 0

    for line in lines:
        if int(skip) > 0:
            skip -= 1
            continue
        line.rstrip()
        sum = 0.0
        count = 0.0
        units = line.split('\t')
        for unit in units:
            try:
                unit = float(unit)
                sum += unit
                count += 1
            except:
                print("WARNING: Non numeric found at line %f.1" % linenum)
                continue
        avg = 0
        if count > 0:
            avg = sum/count
        averages.append(avg)
        linenum += 1

    fout = open(outfile, 'w')
    fout.write("Sample_Average\n")

    for avg in averages:
        fout.write(str('%.3f' % avg)+ '\n')

    return

if __name__ == "__main__": main()
