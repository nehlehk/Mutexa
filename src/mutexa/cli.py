import argparse
from .extract_mutations import run_extract
from .merge_meta import run_merge
from .rollingAvePlots import run_plots
from .check_format import check_mutfile_format, check_metafile_format

def main():

    parser = argparse.ArgumentParser(
        prog="mutexa",
        description="MuTEXA - Mutation Explorer & Analysis"
    )

    parser.add_argument("-u", "--mut", required=True)
    parser.add_argument("-a", "--align", required=True)
    parser.add_argument("-m", "--meta", required=True)
    parser.add_argument("-c", "--cat", required=True)

    parser.add_argument("-s", "--start")
    parser.add_argument("-e", "--end")
    parser.add_argument("-p", "--prefix", default="output")
    parser.add_argument("-t", "--thresh", default=0)
    parser.add_argument("-d", "--days", default=14)
    parser.add_argument("-r", "--ref")

    args = parser.parse_args()

    print("MuTEXA v1.0")

    check_mutfile_format(args.mut)
    check_metafile_format(args.meta)

    run_extract(
        mut_file=args.mut,
        align_file=args.align,
        ref_genome=args.ref,
        prefix=args.prefix
    )

    run_merge(
        start=args.start,
        end=args.end,
        meta=args.meta,
        threshold=args.thresh,
        mut=args.mut,
        prefix=args.prefix,
        days=args.days,
        category=args.cat
    )

    run_plots(
        start=args.start,
        end=args.end,
        threshold=args.thresh,
        mut=args.mut,
        prefix=args.prefix,
        category=args.cat,
        days=args.days
    )