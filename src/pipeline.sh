#!/bin/bash


# set color variables
GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m'

usage() {
    echo ""
    echo ""
    echo ""
    echo -e "${GREEN} __  __ _    _ _______    _______ _____ ____  _   _   _______ _____            _____ _  _______ _   _  _____${NC}"
    echo -e "${GREEN}|  \/  | |  | |__   __|/\|__   __|_   _/ __ \| \ | | |__   __|  __ \     /\   / ____| |/ /_   _| \ | |/ ____|${NC}"
    echo -e "${GREEN}| \  / | |  | |  | |  /  \  | |    | || |  | |  \| |    | |  | |__) |   /  \ | |    | ' /  | | |  \| | |  __ ${NC}"
    echo -e "${GREEN}| |\/| | |  | |  | | / /\ \ | |    | || |  | | .   |    | |  |  _  /   / /\ \| |    |  <   | | | .   | | |_ |${NC}"
    echo -e "${GREEN}| |  | | |__| |  | |/ ____ \| |   _| || |__| | |\  |    | |  | | \ \  / ____ \ |____| . \ _| |_| |\  | |__| |${NC}"
    echo -e "${GREEN}|_|  |_|\____/   |_/_/    \_\_|  |_____\____/|_| \_|    |_|  |_|  \_\/_/    \_\_____|_|\_\_____|_| \_|\_____|${NC}"
    echo ""
    echo ""
    echo ""

    echo "Usage: $0 [-h] -u <mut_file> -a <align_file> -s <start_date> -e <end_date> -p <prefix>  -m <meta_file> -t <threshold> -d <days> -c <Categorizing> [-r <ref_genome>]"
    echo ""
    echo "Run the pipeline to extract mutations, merge metadata, and create plots."
    echo ""
    echo "Required arguments:"
    echo "  -u, --mut      Path to mutations file in CSV format (essential)."
    echo "  -a, --align    Path to alignment file in FASTA format (essential)."
    echo "  -m, --meta     Path to metadata file in CSV format (essential)."
    echo "  -c, --cat      Categorizing class. Based on the metadata file columns it can be: country/continent/lineage (essential)."
    echo "  -p, --prefix   Prefix to save outputs. Default value is 'output' "
    echo "  -s, --start    Start date for filtering metadata (YYYY-MM-DD)"
    echo "  -e, --end      End date for filtering metadata (YYYY-MM-DD). Default value is the date of today"
    echo "  -t, --thresh   Threshold for filtering haplotypes by frequency. Default value is 0 (Zero)"
    echo "  -d, --days     Sliding window days for ..... Default value is 14 "
    echo ""
    echo "Optional arguments:"
    echo "  -h, --help     Show this help message and exit."
    echo "  -r, --ref      Path to reference genome file in FASTA format (optional). If not provided, the first sequence in the alignment file will be used as the reference genome."
    exit 0
}



#start_date="2020-01-01"
end_date=$(date +"%Y-%m-%d")
threshold=0



if [ $# -eq 0 ]; then
    usage
fi

while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -h|--help)
    usage
    shift # past argument
    ;;
    -u|--mut)
    mut_file="$2"
    shift # past argument
    shift # past value
    ;;
    -a|--align)
    align_file="$2"
    shift # past argument
    shift # past value
    ;;
    -s|--start)
    start_date="$2"
    shift # past argument
    shift # past value
    ;;
    -e|--end)
    end_date="$2"
    shift # past argument
    shift # past value
    ;;
    -d|--days)
    days="$2"
    shift # past argument
    shift # past value
    ;;
    -c|--cat)
    cat="$2"
    shift # past argument
    shift # past value
    ;;
    -p|--prefix)
    prefix="$2"
    shift # past argument
    shift # past value
    ;;
    -m|--meta)
    meta_file="$2"
    shift # past argument
    shift # past value
    ;;
    -t|--thresh)
    threshold="$2"
    shift # past argument
    shift # past value
    ;;
    -r|--ref)
    ref_genome="$2"
    shift # past argument
    shift # past value
    ;;
    *)
    break
    ;;
esac
done

# Check if all required arguments are set
if [ -z "$mut_file" ] || [ -z "$align_file" ] || [ -z "$meta_file" ]; then
#    echo "Error: Missing required input file(s)."
    echo -e "${RED}Error: Missing required input file(s).${NC}"
    usage
fi

# Create outputs folder if it doesn't exist
if [ ! -d "outputs" ]; then
    mkdir outputs
fi

# set -e will cause the script to exit immediately if any command exits with a non-zero status
set -e

# Step 0: Create and activate virtual environment
python3 -m venv venv
source venv/bin/activate

# Step 1: Install dependencies
pip install --no-cache-dir -r requirements.txt

# Step 2: Run extract_mutations.py
python check_format.py -u "$mut_file"  -m "$meta_file"

# Step 3: Run extract_mutations.py
python extract_mutations.py -u "$mut_file" -a "$align_file" -r "$ref_genome" -p "$prefix"

# Step 4: Run merge_meta.py
python merge_meta.py -s "$start_date" -e "$end_date"  -m "$meta_file" -t "$threshold" -u "$mut_file" -p "$prefix" -d "$days" -c "$cat"

# Step 5: Run rollingAvePlots.py
python rollingAvePlots.py -s "$start_date" -e "$end_date"  -t "$threshold" -u "$mut_file" -p "$prefix" -c "$cat" -d "$days"

deactivate