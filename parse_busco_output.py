#!/usr/bin/env python3
"""Parse BUSCO summary data."""
import argparse
import os
import sys

__author__ = "Bogdan Kirilenko, 2021"

DESCRIPTION = "Parse BUSCO output and convert it to machine-readable tsv format"
INPUT_HELP = "BUSCO output directory OR summary file"
SHOW_HEADER_HELP = (
    "Output table header. May be helpful for batch processing, such as:\n"
    "Job 1: parse_busco_output.py [busco directory 1] --show_header\n"
    "Jobs 2-N: parse_busco_output.py [busco directory 2-N]"
)
SHORT_SUMM_SUFFIX = "short_summary."

COMPLETE_BUSCO_LINE = "Complete BUSCOs (C)"
COMPLETE_BUSCO_SC_LINE = "Complete and single-copy BUSCOs (S)"
COMPLETE_BUSCO_D_LINE = "Complete and duplicated BUSCOs (D)"
FRAGMENTED_BUSCO_LINE = "Fragmented BUSCOs (F)"
MISSING_BUSCO_LINE = "Missing BUSCOs (M)"

HEADER = "proj\tC\tS\tD\tF\tM\n"


def parse_args():
    """Argument parser."""
    app = argparse.ArgumentParser(description=DESCRIPTION)
    app.add_argument("input", help=INPUT_HELP)
    app.add_argument(
        "--show_header",
        "-s",
        action="store_true",
        dest="show_header",
        help=SHOW_HEADER_HELP,
    )
    app.add_argument(
        "--transpose",
        "-t",
        action="store_true",
        dest="transpose",
        help="Vertical results representation",
    )
    if len(sys.argv) < 2:
        app.print_help()
        sys.exit(0)
    args = app.parse_args()
    return args


def get_file_by_suffix(flist, suffix):
    files_start_with_suffix = [x for x in flist if x.startswith(suffix)]
    if len(files_start_with_suffix) == 0:
        return None
    else:
        return files_start_with_suffix[0]


def get_lines_that_contain(lines, to_find):
    lines_that_contain = [x for x in lines if to_find in x]
    if len(lines_that_contain) == 0:
        return None
    else:
        return lines_that_contain[0]


def get_busco_file(input_arg):
    if os.path.isdir(input_arg):
        files_in_dir = os.listdir(input_arg)
        short_summ_file = get_file_by_suffix(files_in_dir, SHORT_SUMM_SUFFIX)
        if short_summ_file is None:
            sys.exit(f"Error! Cannot find {SHORT_SUMM_SUFFIX}* file in {input_arg}!\n")
        return os.path.join(input_arg, short_summ_file)
    elif os.path.isfile(input_arg):
        return input_arg
    else:
        sys.exit(f"Error! {input_arg} is neither a BUSCO dir nor a file!\n")


def parse_busco_lines(busco_lines):
    """Extract busco numbers."""
    line_suffixes_to_find = [
        COMPLETE_BUSCO_LINE,
        COMPLETE_BUSCO_SC_LINE,
        COMPLETE_BUSCO_D_LINE,
        FRAGMENTED_BUSCO_LINE,
        MISSING_BUSCO_LINE,
    ]
    lines_detected = [
        get_lines_that_contain(busco_lines, s) for s in line_suffixes_to_find
    ]
    if any(x is None for x in lines_detected):
        sys.exit(f"Error! Corrupted busco file: {busco_lines}\n")
    numbers = tuple(int(x.split("\t")[1]) for x in lines_detected)
    return numbers


def get_proj_name(busco_res_filename):
    """Extract project name from summary filename."""
    basename_dot_split = os.path.basename(busco_res_filename).split(".")
    project_id = ".".join(basename_dot_split[3:-1])
    return project_id


def main():
    """Entry point."""
    args = parse_args()
    file_to_parse = get_busco_file(args.input)
    with open(file_to_parse, "r") as f:
        busco_lines = [x.rstrip() for x in f]
    busco_results = parse_busco_lines(busco_lines)
    proj_id = get_proj_name(file_to_parse)
    if args.transpose:  # show vertically and quit
        nums_ids = ["C", "S", "D", "F", "M"]
        print(f"{proj_id}\tclass")
        for i, n in zip(nums_ids, busco_results):
            print(f"{i}\t{n}")
        sys.exit(0)
    # write results row
    nums_line = "\t".join(map(str, busco_results))
    results_line = f"{proj_id}\t{nums_line}\n"
    if args.show_header:
        sys.stdout.write(HEADER)
    sys.stdout.write(results_line)


if __name__ == "__main__":
    main()
