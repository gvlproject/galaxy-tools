#!/usr/bin/env python


import sys	 
import traceback 


def main():

    infile = sys.argv[1]
    outfile = sys.argv[2]

    lines = open( infile, 'r').readlines()

    averages = []
    linenum = 0

    for line in lines:
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
