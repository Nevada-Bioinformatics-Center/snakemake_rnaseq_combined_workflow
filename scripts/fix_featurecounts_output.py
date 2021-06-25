import pprint
import gzip
import argparse
#import pandas as pd
import numpy as np
"""
Quick Binning script for Nathan Ernster's project
By: Hans Vasquez-Gross
"""

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--featurecounts_file", help="featurecounts file")
parser.add_argument("-c", "--column", type=int, default=9, help="Column to fix")
args = parser.parse_args()


##Main Program
if __name__ == "__main__":
    if(args.featurecounts_file):
        featurecounts_fh = open(args.featurecounts_file, "r")
    else:
        raise Exception('Fastq file opened failed!')


    with open(args.featurecounts_file) as in_handle:
        for line in in_handle:
            line = line.rstrip()
            if line.startswith("#"):
                next
            elif line.startswith("Gene"):
                print(line+"\tuniq_column")
            else:
                col_dat = line.split("\t")[(args.column - 1)]
                if col_dat == "NA":
                    print(line+"\tNA")
                else:
                    col_list = list(np.unique(col_dat.split(";")))
                    if "NA" in col_list:
                        col_list.remove("NA")
                    print(line+"\t" + ";".join(col_list))

