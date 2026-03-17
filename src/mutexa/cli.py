import argparse
from .extract_mutations import run_extract
from .merge_meta import run_merge
from .rollingAvePlots import run_plots
from .check_format import check_mutfile_format, check_metafile_format



GREEN = "\033[0;32m"
RED = "\033[0;31m"
NC = "\033[0m"


def print_usage():
    print()
    print(f"{GREEN}MuTEXA v1.0.0{NC}")
    print(f"{GREEN}Mutation EXplorer & Analyzer{NC}")
    print()
    print("Usage:")
    print("  mutexa -u <mut_file> -a <align_file> -m <meta_file> -c <category> -s <start_date> "
          "[-e <end_date>] [-p <prefix>] [-t <threshold>] [-d <days>] [-r <ref_genome>]")
    print()
    print("Required arguments:")
    print("  -u, --mut      Mutations file (CSV)")
    print("  -a, --align    Alignment file (FASTA/VCF)")
    print("  -m, --meta     Metadata file (CSV)")
    print("  -c, --cat      Category: country / continent / lineage")
    print("  -s, --start    Start date (YYYY-MM-DD)")
    print()
    print("Optional arguments:")
    print("  -e, --end      End date (default: today)")
    print("  -p, --prefix   Output prefix (default: output)")
    print("  -t, --thresh   Frequency threshold (default: 0)")
    print("  -d, --days     Sliding window days (default: 14)")
    print("  -r, --ref      Reference genome FASTA")
    print("  -h, --help     Show this help message")
    print()


def main():

    parser = argparse.ArgumentParser(
    prog="mutexa",
    description="MuTEXA - Mutation Explorer & Analysis",
    add_help=False
    )

    # parser = argparse.ArgumentParser(
    #     prog="mutexa",
    #     description="MuTEXA - Mutation Explorer & Analysis"
    # )

    parser.add_argument("-h", "--help", action="store_true")
    parser.add_argument("-u", "--mut")
    parser.add_argument("-a", "--align")
    parser.add_argument("-m", "--meta")
    parser.add_argument("-c", "--cat")

    parser.add_argument("-s", "--start")
    parser.add_argument("-e", "--end")
    parser.add_argument("-p", "--prefix", default="output")
    parser.add_argument("-t", "--thresh", default=0)
    parser.add_argument("-d", "--days", default=14)
    parser.add_argument("-r", "--ref")

    args = parser.parse_args()

    if args.help or len(vars(args)) == 0:
        print_usage()
        exit(0)

    
    # manual validation
    missing = []
    if not args.mut:
        missing.append("--mut")
    if not args.align:
        missing.append("--align")
    if not args.meta:
        missing.append("--meta")
    if not args.cat:
        missing.append("--cat")

    if missing:
        print(f"{RED}Error: Missing required arguments: {', '.join(missing)}{NC}")
        print_usage()
        exit(1)

    # print("MuTEXA v1.0")
    print(f"{GREEN}MuTEXA v1.0{NC}")

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