#!/usr/bin/env python3

import sys
import os
import gzip
import getopt


def main():
    params = parseArgs()

    if params.popmap:
        popmap = parsePopmap(params.popmap)
        popID = numberPops(popmap)

    # open VCF file (gzip-aware)
    if params.vcf.endswith(".gz"):
        vcf_open = lambda f: gzip.open(f, "rt")
    else:
        vcf_open = lambda f: open(f, "r")

    with vcf_open(params.vcf) as vcf:
        last_header = None
        max_ploidy = 0
        sample_names = None
        data = []

        for line in vcf:
            line = line.strip()
            if not line:
                continue

            if line.startswith("#"):
                last_header = line.split("\t")
                continue

            if not sample_names:
                sample_names = last_header[9:]
            fields = line.split("\t")

            # check sample ploidy
            p = max(getPloidies(fields[9:]))
            if p > max_ploidy:
                max_ploidy = p

            data.append(fields[9:])

    # write STRUCTURE output
    with open(params.out, "w") as f:
        for i, sample in enumerate(sample_names):
            base = [sample]
            if params.popmap:
                if sample in popmap:
                    base.append(str(popID[popmap[sample]]))
                else:
                    print(f"Sample {sample} not in popmap. Skipping.")
                    continue

            if params.extracols:
                base.extend([""] * params.extracols)

            olines = [base.copy() for _ in range(max_ploidy)]

            for loc in data:
                gt_field = loc[i].split(":")[0]
                sample_loc = gt_field.replace("|", "/").split("/")

                # Replace all missing alleles "." with "-9"
                sample_loc = [a if a != "." else "-9" for a in sample_loc]

                # Pad or truncate to match max ploidy
                if len(sample_loc) < max_ploidy:
                    sample_loc += ["-9"] * (max_ploidy - len(sample_loc))
                elif len(sample_loc) > max_ploidy:
                    sample_loc = sample_loc[:max_ploidy]

                # handle depth filtering if required
                if params.minIndD is not None:
                    try:
                        dp = int(loc[i].split(":")[-1])
                        if dp < params.minIndD:
                            sample_loc = ["-9"] * len(sample_loc)
                    except (IndexError, ValueError):
                        sample_loc = ["-9"] * len(sample_loc)

                filled = 0
                for index, allele in enumerate(sample_loc):
                    olines[index].append(str(allele))
                    filled += 1

                # fill missing alleles if any
                for missing_index in range(filled, max_ploidy):
                    olines[missing_index].append("-9")

            for o in olines:
                f.write("\t".join(o) + "\n")


def numberPops(popmap):
    popID = dict()
    index = 1
    for pop in popmap.values():
        if pop not in popID:
            popID[pop] = index
            index += 1
    return popID


def parsePopmap(popmap):
    ret = dict()
    if os.path.exists(popmap):
        with open(popmap, "r") as fh:
            for line in fh:
                line = line.strip()
                if line:
                    stuff = line.split()
                    if len(stuff) < 2:
                        continue
                    ret[stuff[0]] = stuff[1]
        return ret
    else:
        raise FileNotFoundError(f"File {popmap} not found!")


def getPloidies(samples):
    return [len((s.split(":")[0]).split("/")) for s in samples]


class parseArgs:
    def __init__(self):
        try:
            options, remainder = getopt.getopt(sys.argv[1:], "hv:o:p:x:d:", ["help"])
        except getopt.GetoptError as err:
            print(err)
            self.display_help("\nExiting due to argument parsing error.")

        self.vcf = None
        self.out = "polyrad.str"
        self.popmap = None
        self.extracols = 0
        self.locnames = False
        self.minIndD = None

        for o, a in options:
            if o in ("-h", "--help"):
                self.display_help("Help menu requested.")

        for opt, arg_raw in options:
            arg = arg_raw.strip()
            opt = opt.replace("-", "")
            if opt == "v":
                self.vcf = arg
            elif opt == "o":
                self.out = arg
            elif opt == "p":
                self.popmap = arg
            elif opt == "x":
                self.extracols = int(arg)
            elif opt == "d":
                self.minIndD = int(arg)

        if not self.vcf:
            self.display_help("Need a VCF file.")

    def display_help(self, message=None):
        if message:
            print(f"\n{message}")
        print(
            """
polyVCFtoStructure.py

Author: Tyler K Chafin, University of Arkansas
Contact: tkchafin@uark.edu

Description: Converts polyRAD-formatted VCF to STRUCTURE format (.str)

Required:
  -v  : VCF input file formatted by polyRAD

Optional:
  -p  : popmap file with pop labels
  -x  : number of extra (blank) columns to add
  -d  : minimum individual DP to keep genotype [default=keep all]
  -o  : output file name [default=polyrad.str]
"""
        )
        sys.exit()


# Entry point
if __name__ == "__main__":
    main()
