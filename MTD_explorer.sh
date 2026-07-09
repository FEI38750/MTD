#!/bin/bash

# ------------------------------------------------------------
# Colors and text styles
# ------------------------------------------------------------

w=$(tput sgr0 2>/dev/null || true)
r=$(tput setaf 1 2>/dev/null || true)
g=$(tput setaf 2 2>/dev/null || true)
y=$(tput setaf 3 2>/dev/null || true)
p=$(tput setaf 5 2>/dev/null || true)

ital=$(tput sitm 2>/dev/null || printf '\033[3m')
noital=$(tput ritm 2>/dev/null || printf '\033[23m')

echo "${w}"

die() {
    echo "${r}[ERROR] $*${w}" >&2
    exit 1
}

require_file() {
    local f="$1"
    local label="${2:-file}"

    if [[ ! -s "$f" ]]; then
        echo "${r}[MISSING] $label: $f${w}" >&2
        exit 1
    fi
}

# ------------------------------------------------------------
# Paired FASTQ record-count validation
# ------------------------------------------------------------
# 1 = validate R1/R2 record counts
# 0 = skip validation
# ------------------------------------------------------------
VALIDATE_PAIRED_FASTQ="${VALIDATE_PAIRED_FASTQ:-1}"

FASTQ_RECORD_COUNT=0

count_fastq_records() {
    local fq="$1"
    local line_count=""

    if [[ ! -e "$fq" ]]; then
        echo "${r}[ERROR] FASTQ file not found: $fq${w}" >&2
        return 1
    fi

    # Empty Kraken2 outputs are valid when no reads belong
    # to a particular category.
    if [[ ! -s "$fq" ]]; then
        FASTQ_RECORD_COUNT=0
        return 0
    fi

    if [[ "$fq" == *.gz ]]; then
        if ! gzip -t -- "$fq"; then
            echo "${r}[ERROR] Invalid or corrupted gzip FASTQ: $fq${w}" >&2
            return 1
        fi

        line_count="$(
            gzip -cd -- "$fq" |
            awk 'END { print NR + 0 }'
        )"
    else
        line_count="$(
            awk 'END { print NR + 0 }' "$fq"
        )"
    fi

    if ! [[ "$line_count" =~ ^[0-9]+$ ]]; then
        echo "${r}[ERROR] Could not count FASTQ lines: $fq${w}" >&2
        return 1
    fi

    if (( line_count % 4 != 0 )); then
        echo "${r}[ERROR] FASTQ line count is not divisible by four.${w}" >&2
        echo "File: $fq" >&2
        echo "Lines detected: $line_count" >&2
        return 1
    fi

    FASTQ_RECORD_COUNT=$(( line_count / 4 ))
}


validate_fastq_pair() {
    local read1="$1"
    local read2="$2"
    local label="${3:-paired FASTQ}"

    local read1_count=0
    local read2_count=0

    if [[ "$VALIDATE_PAIRED_FASTQ" != "1" ]]; then
        return 0
    fi

    count_fastq_records "$read1" || \
        die "Could not validate R1 for $label: $read1"

    read1_count="$FASTQ_RECORD_COUNT"

    count_fastq_records "$read2" || \
        die "Could not validate R2 for $label: $read2"

    read2_count="$FASTQ_RECORD_COUNT"

    if [[ "$read1_count" -ne "$read2_count" ]]; then
        echo "${r}[ERROR] Paired FASTQ record counts do not match.${w}" >&2
        echo "Label: $label" >&2
        echo "R1: $read1" >&2
        echo "R1 records: $read1_count" >&2
        echo "R2: $read2" >&2
        echo "R2 records: $read2_count" >&2
        exit 1
    fi

    echo "${g}[OK] Paired FASTQ validation:${w} $label"
    echo "  Records per mate: $read1_count"
}

run_cmd() {
    echo "${g}[RUN]${w} $*"
    "$@" || die "Command failed: $*"
}

# ------------------------------------------------------------
# Default settings
# ------------------------------------------------------------

pdm="spearman"                         # HAllA method
length=35                              # fastp minimum read length
read_len=75                            # Bracken read length
threads="$(nproc)"                     # CPU threads
blast="hisat"                          # default host alignment method
no_trimm=0                             # default flag
metadata=""                            # optional metadata file
analysis_mode="auto"                   # auto, comparison, exploratory
NO_COMPARISON=0                        # set automatically later
SSGSEA_GMT="default"                  # default, auto, or path to custom GMT
SSGSEA_GMT_MODE="default_MSigDB_c2"
SSGSEA_PROTEIN_FASTA=""
SSGSEA_ANNOTATION_GFF=""
# Optional eggNOG-mapper database override.
# If empty, MTD resolves it from EGGNOG_DATA_DIR,
# offlineCachePath, or the legacy MTD/eggnog_db location.
EGGNOG_DB_DIR=""

# Exploratory-only taxonomic figures
# 1 = run taxonomic heatmap and stacked bar in exploratory mode
# 0 = skip
RUN_EXPLORATORY_TAXONOMIC_FIGURES="${RUN_EXPLORATORY_TAXONOMIC_FIGURES:-1}"

# Run exploratory/descriptive taxonomic figures even in comparison mode
# 1 = run also in comparison mode
# 0 = run only in exploratory/no-comparison mode
RUN_EXPLORATORY_FIGURES_IN_COMPARISON="${RUN_EXPLORATORY_FIGURES_IN_COMPARISON:-1}"

# Exploratory-only host expression QC
# 1 = run PCA, top variable genes heatmap, sample correlation and detected genes plot
# 0 = skip
RUN_EXPLORATORY_HOST_EXPRESSION_QC="${RUN_EXPLORATORY_HOST_EXPRESSION_QC:-1}"

# ------------------------------------------------------------
# Optional detected microbiome read extraction
# ------------------------------------------------------------
# Disabled by default because it can generate many files.
# When enabled, reads are extracted from the absolute detected microbiome ranking.
# These reads can be used later for sequence reconstruction, BLAST, or gene mapping.

EXTRACT_MICROBIOME_READS=0

# Number of top taxa from the absolute ranking to extract.
# 50 = top 50 taxa
# 0  = all detected taxa
EXTRACT_MICROBIOME_READS_TOP_N="${EXTRACT_MICROBIOME_READS_TOP_N:-50}"

# Minimum abundance required in a sample to extract reads for that taxon/sample.
# For absolute mode, this is the absolute abundance value from the Bracken-derived table.
EXTRACT_MICROBIOME_READS_MIN_ABUNDANCE="${EXTRACT_MICROBIOME_READS_MIN_ABUNDANCE:-0}"

# Include child taxa during KrakenTools extraction.
EXTRACT_MICROBIOME_READS_INCLUDE_CHILDREN="${EXTRACT_MICROBIOME_READS_INCLUDE_CHILDREN:-1}"

# Overwrite previous extracted files.
EXTRACT_MICROBIOME_READS_OVERWRITE="${EXTRACT_MICROBIOME_READS_OVERWRITE:-1}"

# Optional Kraken2 host-filtering DB override.
# If empty, DB_host is selected automatically from --hostid.
# If provided through --kraken-host-db, only Kraken2 host filtering uses this custom DB.
KRAKEN_HOST_DB=""

# Optional Kraken2 microbiome/target DB override.
# If empty, DB_micro uses the default MTD kraken2DB_micro.
# If provided through --kraken-micro-db, Kraken2/Bracken non-host classification uses this custom DB.
KRAKEN_MICRO_DB=""

# Kraken2 defaults for host filtering
# Rationale:
#   Host filtering should be sensitive enough to remove host reads,
#   but not completely permissive.
KRAKEN_HOST_CONF="0.05"
KRAKEN_HOST_MIN_HIT_GROUPS="3"

# Kraken2/Bracken defaults for microbiome classification
# Rationale:
#   Microbiome classification is slightly more conservative to reduce
#   false-positive microbial calls after host removal.
KRAKEN_MICRO_CONF="0.10"
KRAKEN_MICRO_MIN_HIT_GROUPS="3"
BRACKEN_THRESHOLD="10"

# Optional cache folder containing already compressed FASTQ files.
# Used only to speed up --no-trim mode.
# If empty, --no-trim uses FASTQ files from the samplesheet directory.
CUSTOM_PATH=""

# ------------------------------------------------------------
# Sequencing library layout
# ------------------------------------------------------------
# READ_LAYOUT_MODE is the mode requested by the user:
#   auto = detect from FASTQ filenames
#   se   = force single-end
#   pe   = force paired-end
#
# READ_LAYOUT will contain the effective layout used by the
# pipeline after input detection and validation.
READ_LAYOUT_MODE="auto"
READ_LAYOUT=""

show_help() {
cat << EOF
Usage:
  bash $(basename "$0") [options]

Required:
  -i, --input FILE                         Path to samplesheet.csv
  -o, --output DIR                         Output directory
  -h, --hostid TAXID                       Host species taxon ID used for annotation/downstream host analysis

Optional:
      --analysis-mode MODE                 Analysis mode: auto, comparison, exploratory
                                            auto: detects whether the samplesheet has >=2 groups
                                            comparison: requires experimental groups and runs DEG-dependent steps
                                            exploratory: skips DEG-dependent steps and generates non-comparison outputs
                                            Default: ${analysis_mode}

      --exploratory, --no-comparison        Alias for --analysis-mode exploratory outputs
                                            Default: ${analysis_mode}
  -m, --metadata FILE                      Metadata CSV file
  -p, --pdm METHOD                         HAllA metric: spearman, pearson, mi, nmi, xicor, dcor
                                           Default: ${pdm}
  -l, --trim-length INT                    Minimum read length required by fastp
                                           Default: ${length}
  -r, --bracken-read-len INT               Bracken read length
                                           Default: ${read_len}
      --threads INT                        Number of CPU threads
                                           Default: nproc = ${threads}
      --ssgsea-gmt default|custom_eggnog|FILE    ssGSEA GMT mode or custom GMT file.
                                          default: use default MSigDB/human GMT.
                                          custom_eggnog: build custom eggNOG/GO GMT from user-provided
                                          protein FASTA and GFF/GTF annotation.
                                          FILE: use an existing GMT file.

      --ssgsea-protein-fasta FILE          Required when --ssgsea-gmt custom_eggnog.
                                           Protein FASTA/proteome matching the host reference.

      --ssgsea-annotation-gff FILE         Required when --ssgsea-gmt custom_eggnog.
                                           GFF3/GTF annotation matching the host reference.
    --eggnog-db DIR
        Optional eggNOG-mapper database directory.
        Resolution order:
          1. --eggnog-db
          2. EGGNOG_DATA_DIR environment variable
          3. offlineCachePath/eggNOG/emapperdb-5.0.2
          4. MTD/eggnog_db
Host processing:
  --read-layout MODE
                                           Sequencing library layout: auto, se, or pe.
                                           auto: detect the layout from the FASTQ files.
                                           se: force single-end input.
                                           pe: force paired-end input.
                                           All samples in one run must use the same layout.
                                           Default: ${READ_LAYOUT_MODE}

  -b, --blast                              Use Magic-BLAST instead of HISAT2
  -t, --no-trim                            Skip fastp trimming
      --custom-raw-path DIR                Folder containing already compressed FASTQ files.
                                           Forces --no-trim mode automatically.
                                           Default: ${CUSTOM_PATH}

Kraken2 host filtering:
      --kraken-host-db DIR                 Optional custom Kraken2 host-filtering database.
                                           If not provided, DB_host is selected automatically from --hostid.
                                           This does NOT change GTF, BLAST/HISAT2, featureCounts, or host DEG resources.

      --kraken-host-confidence FLOAT       Kraken2 --confidence for host-filtering step
                                           Default: ${KRAKEN_HOST_CONF}

      --kraken-host-min-hit-groups INT     Kraken2 --minimum-hit-groups for host-filtering step
                                           Default: ${KRAKEN_HOST_MIN_HIT_GROUPS}

Kraken2 microbiome classification:
      --kraken-micro-db DIR                Optional custom Kraken2 microbiome/target database.
                                           If not provided, DB_micro uses:
                                           \$MTDIR/kraken2DB_micro
                                           Use this for targeted databases such as Trematoda.

      --kraken-micro-confidence FLOAT      Kraken2 --confidence for microbiome step
                                           Default: ${KRAKEN_MICRO_CONF}
      --kraken-micro-min-hit-groups INT    Kraken2 --minimum-hit-groups for microbiome step
                                           Default: ${KRAKEN_MICRO_MIN_HIT_GROUPS}

Bracken:
      --bracken-threshold INT              Bracken -t minimum read threshold
                                           Default: ${BRACKEN_THRESHOLD}
Detected microbiome read extraction:
      --extract-microbiome-reads           Extract Kraken-classified reads for detected microbiome taxa.
                                           Disabled by default.
                                           Uses the absolute detected microbiome ranking.

      --extract-microbiome-reads-top-n INT Number of top taxa to extract from the absolute ranking.
                                           Use 0 to extract all detected taxa.
                                           Default: ${EXTRACT_MICROBIOME_READS_TOP_N}

      --extract-microbiome-reads-min-abundance NUMERIC
                                           Minimum abundance required in a sample to extract reads.
                                           Default: ${EXTRACT_MICROBIOME_READS_MIN_ABUNDANCE}
Other:
      --help                               Show this help message

Examples:
  Default host DB from --hostid:
    bash $(basename "$0") \\
      --input samplesheet.csv \\
      --output MTD_results_Myotis_auto \\
      --hostid 59463 \\
      --blast \\
      --no-trim

  Custom Kraken2 host DB, while keeping --hostid for annotation:
    bash $(basename "$0") \\
      --input samplesheet.csv \\
      --output MTD_results_Carollia_Myotis \\
      --hostid 59463 \\
      --blast \\
      --no-trim \\
      --kraken-host-db /home/me/MTD/kraken2DB_Carollia_Myotis/ \\
      --kraken-micro-db /home/me/MTD/Kraken2DB_trematoda/ \\
      --kraken-host-confidence 0.05 \\
      --kraken-host-min-hit-groups 3 \\
      --kraken-micro-confidence 0.10 \\
      --kraken-micro-min-hit-groups 3 \\
      --bracken-threshold 10

  Custom eggNOG/GO ssGSEA GMT from protein FASTA and GFF/GTF:
    bash $(basename "$0") \\
      --input samplesheet.csv \\
      --output MTD_results_custom_ssGSEA \\
      --hostid 6526 \\
      --blast \\
      --no-trim \\
      --ssgsea-gmt custom_eggnog \\
      --ssgsea-protein-fasta /path/to/proteins.pep.all.fa.gz \\
      --ssgsea-annotation-gff /path/to/annotation.gff3.gz
EOF
}

# Save original command line before parsing arguments
ORIGINAL_COMMAND="$(printf '%q ' "$0" "$@")"
LAUNCH_DIR="$(pwd)"
SCRIPT_NAME="$(basename "$0")"
# ------------------------------------------------------------
# Parse command-line options
# ------------------------------------------------------------

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input)
            inputdr="$2"
            shift 2
            ;;
        --input=*)
            inputdr="${1#*=}"
            shift
            ;;

        -o|--output)
            outputdr="$2"
            shift 2
            ;;
        --output=*)
            outputdr="${1#*=}"
            shift
            ;;
        --analysis-mode)
            analysis_mode="$2"
            shift 2
            ;;
        --analysis-mode=*)
            analysis_mode="${1#*=}"
            shift
            ;;
        --no-comparison|--exploratory)
            analysis_mode="exploratory"
            shift
            ;;
        -h|--hostid)
            hostid="$2"
            shift 2
            ;;
        --hostid=*)
            hostid="${1#*=}"
            shift
            ;;
        -m|--metadata)
            metadata="$2"
            shift 2
            ;;
        --metadata=*)
            metadata="${1#*=}"
            shift
            ;;
        -p|--pdm)
            pdm="$2"
            shift 2
            ;;
        --pdm=*)
            pdm="${1#*=}"
            shift
            ;;
        -l|--trim-length)
            length="$2"
            shift 2
            ;;
        --trim-length=*)
            length="${1#*=}"
            shift
            ;;
        -r|--bracken-read-len)
            read_len="$2"
            shift 2
            ;;
        --bracken-read-len=*)
            read_len="${1#*=}"
            shift
            ;;
        --threads)
            threads="$2"
            shift 2
            ;;
        --threads=*)
            threads="${1#*=}"
            shift
            ;;
        --ssgsea-gmt)
            SSGSEA_GMT="$2"
            shift 2
            ;;
        --ssgsea-gmt=*)
            SSGSEA_GMT="${1#*=}"
            shift
            ;;
        --ssgsea-protein-fasta)
            SSGSEA_PROTEIN_FASTA="$2"
            shift 2
            ;;
        --ssgsea-protein-fasta=*)
            SSGSEA_PROTEIN_FASTA="${1#*=}"
            shift
            ;;
        --ssgsea-annotation-gff)
            SSGSEA_ANNOTATION_GFF="$2"
            shift 2
            ;;
        --ssgsea-annotation-gff=*)
            SSGSEA_ANNOTATION_GFF="${1#*=}"
            shift
            ;;
        --eggnog-db)
            if [[ $# -lt 2 || -z "${2:-}" ]]; then
                die "--eggnog-db requires a directory argument"
            fi
            EGGNOG_DB_DIR="$2"
            shift 2
            ;;
        --eggnog-db=*)
            EGGNOG_DB_DIR="${1#*=}"
            if [[ -z "$EGGNOG_DB_DIR" ]]; then
                die "--eggnog-db requires a non-empty directory"
            fi
            shift
            ;;
        --read-layout)
            if [[ $# -lt 2 || -z "${2:-}" ]]; then
                die "--read-layout requires one argument: auto, se, or pe"
            fi
            READ_LAYOUT_MODE="$2"
            shift 2
            ;;
        --read-layout=*)
            READ_LAYOUT_MODE="${1#*=}"
            if [[ -z "$READ_LAYOUT_MODE" ]]; then
                die "--read-layout requires one argument: auto, se, or pe"
            fi
            shift
            ;;
        -b|--blast)
            blast="blast"
            shift
            ;;
        -t|--no-trim)
            no_trimm=1
            shift
            ;;
        --custom-raw-path)
            CUSTOM_PATH="$2"
            shift 2
            ;;
        --custom-raw-path=*)
            CUSTOM_PATH="${1#*=}"
            shift
            ;;
        --kraken-host-db|--host-kraken-db)
            KRAKEN_HOST_DB="$2"
            shift 2
            ;;
        --kraken-host-db=*|--host-kraken-db=*)
            KRAKEN_HOST_DB="${1#*=}"
            shift
            ;;
        --kraken-micro-db|--micro-kraken-db)
            KRAKEN_MICRO_DB="$2"
            shift 2
            ;;
        --kraken-micro-db=*|--micro-kraken-db=*)
            KRAKEN_MICRO_DB="${1#*=}"
            shift
            ;;
        --kraken-host-confidence|--kraken-host-conf)
            KRAKEN_HOST_CONF="$2"
            shift 2
            ;;
        --kraken-host-confidence=*|--kraken-host-conf=*)
            KRAKEN_HOST_CONF="${1#*=}"
            shift
            ;;

        --kraken-host-min-hit-groups)
            KRAKEN_HOST_MIN_HIT_GROUPS="$2"
            shift 2
            ;;
        --kraken-host-min-hit-groups=*)
            KRAKEN_HOST_MIN_HIT_GROUPS="${1#*=}"
            shift
            ;;
        --kraken-micro-confidence|--kraken-micro-conf)
            KRAKEN_MICRO_CONF="$2"
            shift 2
            ;;
        --kraken-micro-confidence=*|--kraken-micro-conf=*)
            KRAKEN_MICRO_CONF="${1#*=}"
            shift
            ;;
        --kraken-micro-min-hit-groups)
            KRAKEN_MICRO_MIN_HIT_GROUPS="$2"
            shift 2
            ;;
        --kraken-micro-min-hit-groups=*)
            KRAKEN_MICRO_MIN_HIT_GROUPS="${1#*=}"
            shift
            ;;
        --bracken-threshold)
            BRACKEN_THRESHOLD="$2"
            shift 2
            ;;
        --bracken-threshold=*)
            BRACKEN_THRESHOLD="${1#*=}"
            shift
            ;;
        --extract-microbiome-reads)
            EXTRACT_MICROBIOME_READS=1
            shift
            ;;

        --extract-microbiome-reads-top-n)
            EXTRACT_MICROBIOME_READS_TOP_N="$2"
            shift 2
            ;;

        --extract-microbiome-reads-top-n=*)
            EXTRACT_MICROBIOME_READS_TOP_N="${1#*=}"
            shift
            ;;

        --extract-microbiome-reads-min-abundance)
            EXTRACT_MICROBIOME_READS_MIN_ABUNDANCE="$2"
            shift 2
            ;;

        --extract-microbiome-reads-min-abundance=*)
            EXTRACT_MICROBIOME_READS_MIN_ABUNDANCE="${1#*=}"
            shift
            ;;
        --help)
            show_help
            exit 0
            ;;
        *)
            echo "${r}[ERROR] Unknown option: $1${w}" >&2
            echo
            show_help
            exit 1
            ;;
    esac
done

# ------------------------------------------------------------
# Custom raw path behavior
# ------------------------------------------------------------
# --custom-raw-path is a cache of already compressed FASTQ files.
# It is intended only for --no-trim mode.
# If provided, force no_trimm=1 and copy .gz files directly.

if [[ -n "${CUSTOM_PATH:-}" ]]; then
    if [[ ! -d "$CUSTOM_PATH" ]]; then
        die "--custom-raw-path was provided but directory was not found: $CUSTOM_PATH"
    fi

    if [[ "$no_trimm" != "1" ]]; then
        echo "${y}[WARNING] --custom-raw-path was provided, so --no-trim mode will be enabled automatically.${w}"
        no_trimm=1
    fi
fi

# ------------------------------------------------------------
# Required argument checks
# ------------------------------------------------------------

if [[ -z "${inputdr:-}" ]]; then
    die "Missing required argument: -i or --input samplesheet.csv"
fi

if [[ -z "${outputdr:-}" ]]; then
    die "Missing required argument: -o or --output output_directory"
fi

if [[ -z "${hostid:-}" ]]; then
    die "Missing required argument: -h or --hostid TAXID"
fi
if [[ "$analysis_mode" != "auto" && "$analysis_mode" != "comparison" && "$analysis_mode" != "exploratory" ]]; then
    die "--analysis-mode must be one of: auto, comparison, exploratory. Got: $analysis_mode"
fi

# Normalize the requested library layout to lowercase.
READ_LAYOUT_MODE="${READ_LAYOUT_MODE,,}"

case "$READ_LAYOUT_MODE" in
    auto)
        # The effective layout will be determined after locating
        # the FASTQ files.
        READ_LAYOUT=""
        ;;

    se|pe)
        # The user explicitly selected the layout.
        READ_LAYOUT="$READ_LAYOUT_MODE"
        ;;

    *)
        die "--read-layout must be one of: auto, se, pe.
Got: $READ_LAYOUT_MODE"
        ;;
esac
# ------------------------------------------------------------
# Basic value validation
# ------------------------------------------------------------

if ! [[ "$threads" =~ ^[0-9]+$ ]] || [[ "$threads" -lt 1 ]]; then
    die "--threads must be a positive integer. Got: $threads"
fi

if ! [[ "$length" =~ ^[0-9]+$ ]] || [[ "$length" -lt 1 ]]; then
    die "--trim-length must be a positive integer. Got: $length"
fi

if ! [[ "$read_len" =~ ^[0-9]+$ ]] || [[ "$read_len" -lt 1 ]]; then
    die "--bracken-read-len must be a positive integer. Got: $read_len"
fi

# ------------------------------------------------------------
# Kraken2/Bracken parameter validation
# ------------------------------------------------------------

if ! [[ "$KRAKEN_HOST_MIN_HIT_GROUPS" =~ ^[0-9]+$ ]] || [[ "$KRAKEN_HOST_MIN_HIT_GROUPS" -lt 1 ]]; then
    die "--kraken-host-min-hit-groups must be a positive integer. Got: $KRAKEN_HOST_MIN_HIT_GROUPS"
fi

if ! [[ "$KRAKEN_MICRO_MIN_HIT_GROUPS" =~ ^[0-9]+$ ]] || [[ "$KRAKEN_MICRO_MIN_HIT_GROUPS" -lt 1 ]]; then
    die "--kraken-micro-min-hit-groups must be a positive integer. Got: $KRAKEN_MICRO_MIN_HIT_GROUPS"
fi

if ! [[ "$BRACKEN_THRESHOLD" =~ ^[0-9]+$ ]] || [[ "$BRACKEN_THRESHOLD" -lt 0 ]]; then
    die "--bracken-threshold must be an integer >= 0. Got: $BRACKEN_THRESHOLD"
fi

if ! awk -v x="$KRAKEN_HOST_CONF" 'BEGIN { exit !(x >= 0 && x <= 1) }'; then
    die "--kraken-host-confidence must be between 0 and 1. Got: $KRAKEN_HOST_CONF"
fi

if ! awk -v x="$KRAKEN_MICRO_CONF" 'BEGIN { exit !(x >= 0 && x <= 1) }'; then
    die "--kraken-micro-confidence must be between 0 and 1. Got: $KRAKEN_MICRO_CONF"
fi

if [[ "$pdm" != "spearman" && "$pdm" != "pearson" && "$pdm" != "mi" && "$pdm" != "nmi" && "$pdm" != "xicor" && "$pdm" != "dcor" ]]; then
    die "--pdm must be one of: spearman, pearson, mi, nmi, xicor, dcor. Got: $pdm"
fi

if ! [[ "$EXTRACT_MICROBIOME_READS_TOP_N" =~ ^[0-9]+$ ]]; then
    die "--extract-microbiome-reads-top-n must be an integer >= 0. Got: $EXTRACT_MICROBIOME_READS_TOP_N"
fi

if ! awk -v x="$EXTRACT_MICROBIOME_READS_MIN_ABUNDANCE" 'BEGIN { exit !(x >= 0) }'; then
    die "--extract-microbiome-reads-min-abundance must be numeric and >= 0. Got: $EXTRACT_MICROBIOME_READS_MIN_ABUNDANCE"
fi

# ------------------------------------------------------------
# ssGSEA custom GMT early validation
# ------------------------------------------------------------

if [[ "$SSGSEA_GMT" == "auto" || "$SSGSEA_GMT" == "custom_eggnog" ]]; then
    if [[ -z "${SSGSEA_PROTEIN_FASTA:-}" ]]; then
        die "--ssgsea-gmt $SSGSEA_GMT requires --ssgsea-protein-fasta FILE"
    fi

    if [[ -z "${SSGSEA_ANNOTATION_GFF:-}" ]]; then
        die "--ssgsea-gmt $SSGSEA_GMT requires --ssgsea-annotation-gff FILE"
    fi

    require_file "$SSGSEA_PROTEIN_FASTA" "ssGSEA protein FASTA"
    require_file "$SSGSEA_ANNOTATION_GFF" "ssGSEA GFF/GTF annotation"
fi

if [[ "$SSGSEA_GMT" != "default" &&
      "$SSGSEA_GMT" != "auto" &&
      "$SSGSEA_GMT" != "custom_eggnog" ]]; then
    require_file "$SSGSEA_GMT" "Custom ssGSEA GMT"
fi

# ------------------------------------------------------------
# Output directory behavior
# ------------------------------------------------------------

if [[ -d "$outputdr" ]]; then
    echo "The output directory '$outputdr' already exists."
    read -p "Do you want to delete it and overwrite the files? (y/n): " answer

    if [[ "$answer" =~ ^[Yy]$ ]]; then
        rm -rf "$outputdr"
        echo "Directory deleted."
    else
        echo "Operation cancelled by the user. Exiting."
        exit 1
    fi
fi

# ------------------------------------------------------------
# MTD path and Conda environment
# ------------------------------------------------------------

MTDIR=$(dirname "$(readlink -f "$0")")
echo "MTD directory is $MTDIR"

# ------------------------------------------------------------
# eggNOG database path resolution
# ------------------------------------------------------------
resolve_eggnog_db_dir() {
    local cache_root=""
    local candidate=""
    local required_file=""

    echo "${g}[INFO] Resolving eggNOG-mapper database directory...${w}"

    # 1. Explicit command-line option
    if [[ -n "${EGGNOG_DB_DIR:-}" ]]; then
        candidate="$EGGNOG_DB_DIR"
        echo "[INFO] Source: --eggnog-db"

    # 2. Standard eggNOG environment variable
    elif [[ -n "${EGGNOG_DATA_DIR:-}" ]]; then
        candidate="$EGGNOG_DATA_DIR"
        echo "[INFO] Source: EGGNOG_DATA_DIR"

    # 3. Persistent MTD installation cache
    elif [[ -s "$MTDIR/offlineCachePath" ]]; then
        cache_root="$(
            head -n 1 "$MTDIR/offlineCachePath" |
            tr -d '\r' |
            sed 's/[[:space:]]*$//'
        )"

        if [[ -z "$cache_root" ]]; then
            die "offlineCachePath exists but is empty: $MTDIR/offlineCachePath"
        fi

        candidate="$cache_root/eggNOG/emapperdb-5.0.2"
        echo "[INFO] Source: $MTDIR/offlineCachePath"

    # 4. Backward-compatible legacy directory
    elif [[ -d "$MTDIR/eggnog_db" ]]; then
        candidate="$MTDIR/eggnog_db"
        echo "[INFO] Source: legacy MTD/eggnog_db"
    fi

    if [[ -z "$candidate" ]]; then
        die "Could not resolve the eggNOG database directory.

Provide one of:
  --eggnog-db DIR
  EGGNOG_DATA_DIR=/path/to/database
  $MTDIR/offlineCachePath
  $MTDIR/eggnog_db"
    fi

    # Expand an initial ~/ when supplied through the command line.
    if [[ "$candidate" == "~/"* ]]; then
        candidate="$HOME/${candidate#~/}"
    fi

    if [[ ! -d "$candidate" ]]; then
        die "eggNOG database directory was not found:
  $candidate"
    fi

    candidate="$(readlink -f "$candidate")"

    for required_file in \
        eggnog.db \
        eggnog_proteins.dmnd \
        eggnog.taxa.db \
        eggnog.taxa.db.traverse.pkl
    do
        if [[ ! -s "$candidate/$required_file" ]]; then
            die "Required eggNOG database file is missing or empty:
  $candidate/$required_file"
        fi

        echo "${g}[OK]${w} $required_file"
    done

    EGGNOG_DB_DIR="$candidate"
    export EGGNOG_DATA_DIR="$EGGNOG_DB_DIR"

    echo "${g}[INFO] eggNOG database directory:${w}"
    echo "  $EGGNOG_DB_DIR"
}

# ------------------------------------------------------------
# Custom R library for MTD-generated OrgDb packages
# Example: custom org.Bglabrata.eg.db with BGLAX IDs
# ------------------------------------------------------------

mkdir -p "$MTDIR/custom_R_libs"
export R_LIBS_USER="$MTDIR/custom_R_libs:${R_LIBS_USER:-}"

echo "${g}[INFO] R custom library:${w}"
echo "  $R_LIBS_USER"

# ------------------------------------------------------------
# ssGSEA GMT default path
# Actual GMT selection happens later, after host.gct is generated.
# ------------------------------------------------------------

SSGSEA_DEFAULT_GMT="$MTDIR/Tools/ssGSEA2.0/db/msigdb/c2.all.v7.5.1.symbols.gmt"

if [[ "$SSGSEA_GMT" == "default" ]]; then
    SSGSEA_GMT_MODE="default_MSigDB_c2_symbols"
elif [[ "$SSGSEA_GMT" == "auto" || "$SSGSEA_GMT" == "custom_eggnog" ]]; then
    SSGSEA_GMT_MODE="custom_eggNOG_GO_from_user_files_pending"
else
    SSGSEA_GMT_MODE="custom_path_pending"
fi

condapath=$(head -n 1 "$MTDIR/condaPath")

source "$condapath/etc/profile.d/conda.sh"
conda deactivate
conda activate MTD

# ------------------------------------------------------------
# Input paths
# ------------------------------------------------------------

samplesheet_file="$inputdr"

echo "${g}[INFO] Requested read layout:${w} $READ_LAYOUT_MODE"

if [[ -n "$READ_LAYOUT" ]]; then
    echo "${g}[INFO] Effective read layout:${w} $READ_LAYOUT"
else
    echo "${g}[INFO] Effective read layout:${w} pending automatic detection"
fi

if [[ ! -s "$samplesheet_file" ]]; then
    die "Samplesheet not found or empty: $samplesheet_file"
fi

inputdr=$(dirname "$samplesheet_file")

mkdir -p "$outputdr"
mkdir -p "$outputdr/temp"

cd "$outputdr/temp" || die "Could not enter temp directory: $outputdr/temp"

# ------------------------------------------------------------
# Step 0: Host database auto selection from --hostid
# ------------------------------------------------------------

if [[ "$hostid" == 9606 ]]; then
    DB_host="$MTDIR/kraken2DB_human"
    DB_hisat2="$MTDIR/hisat2_index_human/genome_tran"
    DB_blast="$MTDIR/human_blastdb/human_blastdb"
    gtf="$MTDIR/ref_human/Homo_sapiens.GRCh38.104.gtf.gz"

elif [[ "$hostid" == 9544 ]]; then
    DB_host="$MTDIR/kraken2DB_rhesus"
    DB_hisat2="$MTDIR/hisat2_index_rhesus/genome_tran"
    DB_blast="$MTDIR/rhesus_blastdb/rhesus_blastdb"
    gtf="$MTDIR/ref_rhesus/Macaca_mulatta.Mmul_10.104.gtf.gz"

elif [[ "$hostid" == 10090 ]]; then
    DB_host="$MTDIR/kraken2DB_mice"
    DB_hisat2="$MTDIR/hisat2_index_mouse/genome_tran"
    DB_blast="$MTDIR/mouse_blastdb/mouse_blastdb"
    gtf="$MTDIR/ref_mouse/Mus_musculus.GRCm39.104.gtf.gz"

elif [[ -d "$MTDIR/kraken2DB_${hostid}" ]]; then
    DB_host="$MTDIR/kraken2DB_${hostid}"
    DB_hisat2="$MTDIR/hisat2_index_${hostid}/genome_tran"
    DB_blast="$MTDIR/blastdb_${hostid}/blastdb_${hostid}"
    gtf="$MTDIR/ref_${hostid}/ref_${hostid}.gtf.gz"

else
    echo "${r}[ERROR] Host species is not supported for --hostid $hostid.${w}"
    echo "You can use bash Customized_host.sh for building the required resources."
    exit 1
fi

# ------------------------------------------------------------
# Optional Kraken2 host-filtering DB override
# ------------------------------------------------------------
# --hostid still controls downstream host annotation resources:
#   GTF, HISAT2/Magic-BLAST DB, featureCounts and host DEG.
#
# --kraken-host-db, when provided, overrides only DB_host for
# Kraken2 host read filtering.

DB_host_from_hostid="$DB_host"

if [[ -n "${KRAKEN_HOST_DB:-}" ]]; then
    DB_host="$KRAKEN_HOST_DB"
    KRAKEN_HOST_DB_MODE="custom_path_from_--kraken-host-db"
else
    DB_host="$DB_host_from_hostid"
    KRAKEN_HOST_DB_MODE="auto_from_--hostid"
fi

DB_micro_default="$MTDIR/kraken2DB_micro"

if [[ -n "${KRAKEN_MICRO_DB:-}" ]]; then
    DB_micro="$KRAKEN_MICRO_DB"
    KRAKEN_MICRO_DB_MODE="custom_path_from_--kraken-micro-db"
else
    DB_micro="$DB_micro_default"
    KRAKEN_MICRO_DB_MODE="default_MTD_kraken2DB_micro"
fi

# ------------------------------------------------------------
# Scientific name from HostSpecies.csv
# ------------------------------------------------------------

species_name=$(awk -F, -v taxid="$hostid" '
    NR > 1 && $1 == taxid {
        print $3
        exit
    }
' "$MTDIR/HostSpecies.csv")

if [[ -z "$species_name" ]]; then
    species_name="scientific name not found"
    hostid_display="${hostid} (${y}${species_name}${w})"
else
    hostid_display="${hostid} (${ital}${species_name}${noital})"
fi

# ------------------------------------------------------------
# Opening summary
# ------------------------------------------------------------

echo "${g}Selected pipeline parameters:${w}"
echo "  input samplesheet:              $samplesheet_file"
echo "  input directory:                $inputdr"
echo "  output directory:               $outputdr"
echo "  host annotation taxid:          $hostid_display"
echo "  host Kraken2 DB:                $DB_host"
echo "  host Kraken2 DB mode:           $KRAKEN_HOST_DB_MODE"
echo "  microbiome Kraken2 DB:          $DB_micro"
echo "  microbiome Kraken2 DB mode:     $KRAKEN_MICRO_DB_MODE"
echo "  metadata:                       ${metadata:-none}"
echo "  HAllA metric:                   $pdm"
echo "  trim length:                    $length"
echo "  Bracken read length:            $read_len"
echo "  threads:                        $threads"
echo "  host alignment mode:            $blast"
echo "  no_trim:                        $no_trimm"
echo "  custom raw cache path:          ${CUSTOM_PATH:-not provided}"
echo "  sequencing read layout request: $READ_LAYOUT_MODE"
echo "  automatic SRA download:         enabled for SRR/ERR/DRR when no local FASTQ exists"
echo "  SRA FASTQ conversion mode:      fasterq-dump default split-3"
echo "  Kraken host confidence:         $KRAKEN_HOST_CONF"
echo "  Kraken host min hit groups:     $KRAKEN_HOST_MIN_HIT_GROUPS"
echo "  Kraken micro confidence:        $KRAKEN_MICRO_CONF"
echo "  Kraken micro min hit groups:    $KRAKEN_MICRO_MIN_HIT_GROUPS"
echo "  Bracken threshold:              $BRACKEN_THRESHOLD"
echo "  extract microbiome reads:       $EXTRACT_MICROBIOME_READS"
echo "  extract microbiome reads mode:  absolute"
echo "  extract microbiome top N:       $EXTRACT_MICROBIOME_READS_TOP_N"
echo "  extract microbiome min abundance: $EXTRACT_MICROBIOME_READS_MIN_ABUNDANCE"
echo "  ssGSEA GMT:                     $SSGSEA_GMT"
echo "  ssGSEA GMT mode:                $SSGSEA_GMT_MODE"

echo "${g}Kraken2 microbiome parameters:${w}"
if [[ "$KRAKEN_MICRO_DB_MODE" == "custom_path_from_--kraken-micro-db" ]]; then
    echo "  --kraken-micro-db $DB_micro"
else
    echo "  --kraken-micro-db not provided; using default MTD DB_micro"
fi
echo "  --kraken-micro-confidence $KRAKEN_MICRO_CONF"
echo "  --kraken-micro-min-hit-groups $KRAKEN_MICRO_MIN_HIT_GROUPS"

echo "${g}Kraken2 host filtering parameters:${w}"
if [[ "$KRAKEN_HOST_DB_MODE" == "custom_path_from_--kraken-host-db" ]]; then
    echo "  --kraken-host-db $DB_host"
else
    echo "  --kraken-host-db not provided; using DB from --hostid"
fi
echo "  --kraken-host-confidence $KRAKEN_HOST_CONF"
echo "  --kraken-host-min-hit-groups $KRAKEN_HOST_MIN_HIT_GROUPS"

echo "${g}============================================"
echo "Selected host annotation species:${w} ${ital}${species_name}${noital}${g}"
echo "Annotation Taxon ID:${w} $hostid ${g}"
echo "Host Kraken2 filtering DB:${w} $DB_host ${g}"
echo "Microbiome/target Kraken2 DB:${w} $DB_micro ${g}"
echo "${g}============================================${w}"

# ------------------------------------------------------------
# Export methods and run parameters
# ------------------------------------------------------------

write_methods_log() {
    local methods_dir="$outputdr/methods"
    local methods_csv="$methods_dir/mtd_methods_run_parameters.csv"
    local bracken_dist="$DB_micro/database${read_len}mers.kmer_distrib"

    mkdir -p "$methods_dir"

    csv_escape() {
        local s="${1:-}"
        s="${s//$'\r'/}"
        s="${s//$'\n'/ }"
        s="${s//\"/\"\"}"
        printf '"%s"' "$s"
    }

    csv_row() {
        csv_escape "$1"; printf ","
        csv_escape "$2"; printf ","
        csv_escape "$3"; printf ","
        csv_escape "$4"; printf ","
        csv_escape "$5"; printf "\n"
    }

    get_tool_path() {
        command -v "$1" 2>/dev/null || printf "not_found"
    }

    get_tool_version() {
        local cmd="$1"
        eval "$cmd" 2>&1 | head -n 1 | sed 's/\r//g'
    }

    {
        csv_row "category" "program" "parameter" "value" "description"

        # Run metadata
        csv_row "Run metadata" "$SCRIPT_NAME" "run_datetime" "$(date -Is)" "Date and time when the run was started"
        csv_row "Run metadata" "$SCRIPT_NAME" "launch_directory" "$LAUNCH_DIR" "Directory from which the command was launched"
        csv_row "Run metadata" "$SCRIPT_NAME" "original_command" "$ORIGINAL_COMMAND" "Original command line used to start the run"
        csv_row "Run metadata" "$SCRIPT_NAME" "MTD_directory" "$MTDIR" "Path to the MTD installation directory"
        csv_row "Run metadata" "$SCRIPT_NAME" "conda_path" "$condapath" "Path to Conda installation used by the pipeline"

        # Input/output
        csv_row "Input and output" "$SCRIPT_NAME" "samplesheet_file" "$samplesheet_file" "Input samplesheet CSV"
        csv_row "Input and output" "$SCRIPT_NAME" "input_directory" "$inputdr" "Directory containing the original FASTQ files or samplesheet"
        csv_row "Input and output" "$SCRIPT_NAME" "output_directory" "$outputdr" "Main output directory"
        csv_row "Input and output" "$SCRIPT_NAME" "metadata_file" "${metadata:-none}" "Optional metadata file"

        # Public sequencing archive input
        csv_row \
            "Input and output" \
            "SRA Toolkit" \
            "supported_run_accession_prefixes" \
            "SRR;ERR;DRR" \
            "Run accessions automatically converted when no local FASTQ input exists"

        csv_row \
            "Input and output" \
            "prefetch" \
            "maximum_accession_size" \
            "999G" \
            "Maximum accession size accepted by automatic prefetch"

        csv_row \
            "Input and output" \
            "fasterq-dump" \
            "output_mode" \
            "default_split_3" \
            "Produces ACCESSION.fastq for unpaired reads or ACCESSION_1.fastq and ACCESSION_2.fastq for paired reads"

        csv_row \
            "Input and output" \
            "fasterq-dump" \
            "threads" \
            "$threads" \
            "Threads used to convert SRA accessions to FASTQ"
        # Study design
        csv_row "Study design" "$SCRIPT_NAME" "host_annotation_taxid" "$hostid" "Taxon ID used for host annotation and downstream host analyses"
        csv_row "Study design" "HostSpecies.csv" "host_annotation_species" "$species_name" "Scientific name associated with host annotation taxid"
        csv_row "Study design" "$SCRIPT_NAME" "HAllA_metric" "$pdm" "Association metric used by HAllA"

# Read preparation
        csv_row \
            "Read preparation" \
            "fastp / pigz / cp" \
            "no_trim" \
            "$no_trimm" \
            "If 1, fastp trimming is skipped"

        csv_row \
            "Read preparation" \
            "fastp" \
            "trim_length" \
            "$length" \
            "Minimum read length required by fastp"

        csv_row \
            "Read preparation" \
            "$SCRIPT_NAME" \
            "read_layout" \
            "${READ_LAYOUT:-pending_auto_detection}" \
            "Effective sequencing library layout"

        csv_row \
            "Read preparation" \
            "$SCRIPT_NAME" \
            "custom_raw_cache_path" \
            "${CUSTOM_PATH:-not provided}" \
            "Optional cache directory containing already compressed FASTQ files"

        csv_row \
            "Read preparation" \
            "$SCRIPT_NAME" \
            "prepared_fastq_pattern_se" \
            "Trimmed_SAMPLE.fq.gz" \
            "Normalized prepared FASTQ pattern for single-end input"

        csv_row \
            "Read preparation" \
            "$SCRIPT_NAME" \
            "prepared_fastq_pattern_pe" \
            "Trimmed_SAMPLE_R1.fq.gz + Trimmed_SAMPLE_R2.fq.gz" \
            "Normalized prepared FASTQ patterns for paired-end input"
        csv_row \
            "Read preparation" \
            "$SCRIPT_NAME" \
            "paired_fastq_record_validation" \
            "$VALIDATE_PAIRED_FASTQ" \
            "Checks that paired-end R1 and R2 contain the same number of FASTQ records"
        # Host filtering
        csv_row "Host filtering" "Kraken2" "host_kraken_db" "$DB_host" "Kraken2 database used for host read filtering"
        csv_row "Host filtering" "Kraken2" "host_kraken_db_mode" "$KRAKEN_HOST_DB_MODE" "Whether DB_host came from --hostid or --kraken-host-db"
        csv_row "Host filtering" "Kraken2" "host_kraken_db_from_hostid" "$DB_host_from_hostid" "Default Kraken2 host DB selected from --hostid before optional override"
        csv_row "Host filtering" "Kraken2" "confidence" "$KRAKEN_HOST_CONF" "Kraken2 --confidence used for host filtering"
        csv_row "Host filtering" "Kraken2" "minimum_hit_groups" "$KRAKEN_HOST_MIN_HIT_GROUPS" "Kraken2 --minimum-hit-groups used for host filtering"
        csv_row \
            "Host filtering" \
            "Kraken2" \
            "classified_out_se" \
            "SAMPLE_host.fq" \
            "Single-end reads classified as host"

        csv_row \
            "Host filtering" \
            "Kraken2" \
            "classified_out_pe" \
            "SAMPLE_host_1.fq + SAMPLE_host_2.fq" \
            "Paired-end reads classified as host"

        csv_row \
            "Host filtering" \
            "Kraken2" \
            "unclassified_out_se" \
            "SAMPLE_non-host_raw.fq" \
            "Single-end reads not classified as host"

        csv_row \
            "Host filtering" \
            "Kraken2" \
            "unclassified_out_pe" \
            "SAMPLE_non-host_raw_1.fq + SAMPLE_non-host_raw_2.fq" \
            "Paired-end reads not classified as host"
        # Microbiome classification
        csv_row "Microbiome classification" "Kraken2" "microbiome_kraken_db" "$DB_micro" "Kraken2 database used for microbiome classification"
        csv_row "Microbiome classification" "Kraken2" "microbiome_kraken_db_mode" "$KRAKEN_MICRO_DB_MODE" "Whether DB_micro came from default MTD DB or --kraken-micro-db"
        csv_row "Microbiome classification" "Kraken2" "microbiome_kraken_db_default" "$DB_micro_default" "Default MTD microbiome Kraken2 DB before optional override"
        csv_row "Microbiome classification" "Kraken2" "confidence" "$KRAKEN_MICRO_CONF" "Kraken2 --confidence used for microbiome classification"
        csv_row "Microbiome classification" "Kraken2" "minimum_hit_groups" "$KRAKEN_MICRO_MIN_HIT_GROUPS" "Kraken2 --minimum-hit-groups used for microbiome classification"
        csv_row "Microbiome classification" "Kraken2" "raw_report_pattern" "Report_non-host.raw_SAMPLE.txt" "Raw microbiome Kraken2 report before optional contaminant removal"
        csv_row "Microbiome classification" "Kraken2" "final_report_pattern" "Report_non-host_SAMPLE.txt" "Final microbiome Kraken2 report used by Bracken"

        # Contaminant removal
        csv_row "Contaminant removal" "KrakenTools extract_kraken_reads.py" "contaminant_list" "$MTDIR/conta_ls.txt" "Optional list of contaminant taxids to exclude"
        csv_row \
            "Contaminant removal" \
            "KrakenTools extract_kraken_reads.py" \
            "output_format" \
            "FASTQ (--fastq-output)" \
            "Preserves nucleotide qualities in decontaminated SE or PE reads"

        csv_row \
            "Contaminant removal" \
            "KrakenTools extract_kraken_reads.py" \
            "single_end_output" \
            "SAMPLE_non-host.fq" \
            "Final decontaminated single-end reads"

        csv_row \
            "Contaminant removal" \
            "KrakenTools extract_kraken_reads.py" \
            "paired_end_output" \
            "SAMPLE_non-host_1.fq + SAMPLE_non-host_2.fq" \
            "Final decontaminated paired-end reads"
        if [[ -s "$MTDIR/conta_ls.txt" ]]; then
            csv_row "Contaminant removal" "KrakenTools extract_kraken_reads.py" "contaminant_list_status" "present" "conta_ls.txt was found and may be used"
        else
            csv_row "Contaminant removal" "KrakenTools extract_kraken_reads.py" "contaminant_list_status" "not_found_or_empty" "No contaminant list was found or file was empty"
        fi

        # Bracken
        csv_row "Abundance estimation" "Bracken" "bracken_read_length" "$read_len" "Read length used by Bracken -r"
        csv_row "Abundance estimation" "Bracken" "bracken_threshold" "$BRACKEN_THRESHOLD" "Minimum read threshold used by Bracken -t"
        csv_row "Abundance estimation" "Bracken" "bracken_distribution_file" "$bracken_dist" "Expected Bracken read-length-specific kmer distribution file"
        if [[ -s "$bracken_dist" ]]; then
            csv_row "Abundance estimation" "Bracken" "bracken_distribution_file_status" "present" "Bracken distribution file exists"
        else
            csv_row "Abundance estimation" "Bracken" "bracken_distribution_file_status" "missing" "Bracken distribution file not found at run-start logging time"
        fi
        csv_row "Abundance estimation" "Bracken" "taxonomic_levels" "P;G;S" "Bracken is run at phylum, genus and species levels"
# Detected microbiome read extraction
        csv_row \
            "Detected microbiome read extraction" \
            "$SCRIPT_NAME" \
            "extract_microbiome_reads" \
            "$EXTRACT_MICROBIOME_READS" \
            "Whether Kraken-classified reads are extracted for detected microbiome taxa"

        csv_row \
            "Detected microbiome read extraction" \
            "$SCRIPT_NAME" \
            "ranking_mode" \
            "absolute" \
            "Read extraction uses the absolute detected microbiome ranking"

        csv_row \
            "Detected microbiome read extraction" \
            "$SCRIPT_NAME" \
            "top_n" \
            "$EXTRACT_MICROBIOME_READS_TOP_N" \
            "Number of top taxa extracted; 0 means all detected taxa"

        csv_row \
            "Detected microbiome read extraction" \
            "$SCRIPT_NAME" \
            "min_abundance" \
            "$EXTRACT_MICROBIOME_READS_MIN_ABUNDANCE" \
            "Minimum sample abundance required for extraction"

        csv_row \
            "Detected microbiome read extraction" \
            "$SCRIPT_NAME" \
            "include_children" \
            "$EXTRACT_MICROBIOME_READS_INCLUDE_CHILDREN" \
            "Whether child taxa are included during KrakenTools read extraction"

        csv_row \
            "Detected microbiome read extraction" \
            "MTD_extract_reads_by_detected_microbiome.py" \
            "read_layout" \
            "${READ_LAYOUT:-pending_auto_detection}" \
            "Effective SE or PE layout used for taxon-specific read extraction"

        csv_row \
            "Detected microbiome read extraction" \
            "KrakenTools extract_kraken_reads.py" \
            "single_end_output_pattern" \
            "SAMPLE.taxidTAXID.TAXON.fastq" \
            "Taxon-specific output for single-end input"

        csv_row \
            "Detected microbiome read extraction" \
            "KrakenTools extract_kraken_reads.py" \
            "paired_end_output_pattern" \
            "SAMPLE.taxidTAXID.TAXON.R1.fastq + SAMPLE.taxidTAXID.TAXON.R2.fastq" \
            "Taxon-specific paired FASTQ outputs"

        csv_row \
            "Detected microbiome read extraction" \
            "MTD_extract_reads_by_detected_microbiome.py" \
            "pair_validation" \
            "R1_and_R2_sequence_counts_must_match" \
            "Mismatched paired outputs are reported and excluded from combined paired files"

        # Host downstream analysis
        csv_row "Host downstream analysis" "$blast" "host_alignment_mode" "$blast" "Host read alignment mode: magicblast if blast, otherwise HISAT2"
        csv_row \
            "Host downstream analysis" \
            "$SCRIPT_NAME" \
            "read_layout" \
            "${READ_LAYOUT:-pending_auto_detection}" \
            "Effective sequencing layout used for host alignment"

        csv_row \
            "Host downstream analysis" \
            "Magic-BLAST" \
            "single_end_input" \
            "SAMPLE_host.fq" \
            "Single-end reads classified as host by Kraken2"

        csv_row \
            "Host downstream analysis" \
            "Magic-BLAST" \
            "paired_end_input" \
            "SAMPLE_host_1.fq + SAMPLE_host_2.fq" \
            "Paired-end reads classified as host by Kraken2"

        csv_row \
            "Host downstream analysis" \
            "HISAT2" \
            "single_end_input" \
            "Trimmed_SAMPLE.fq.gz" \
            "Prepared single-end FASTQ"

        csv_row \
            "Host downstream analysis" \
            "HISAT2" \
            "paired_end_input" \
            "Trimmed_SAMPLE_R1.fq.gz + Trimmed_SAMPLE_R2.fq.gz" \
            "Prepared paired-end FASTQ mates"

        csv_row \
            "Host downstream analysis" \
            "featureCounts" \
            "paired_end_mode" \
            "-p plus --countReadPairs when supported" \
            "Paired-end alignments are counted as fragments/read pairs"
        csv_row "Host downstream analysis" "Magic-BLAST" "blast_database" "$DB_blast" "Magic-BLAST database selected from --hostid"
        csv_row "Host downstream analysis" "HISAT2" "hisat2_database" "$DB_hisat2" "HISAT2 database selected from --hostid"
        csv_row "Host downstream analysis" "featureCounts" "gtf_file" "$gtf" "GTF annotation selected from --hostid"
        csv_row "Host downstream analysis" "featureCounts" "output" "$outputdr/host_counts.txt" "Host count matrix output"
        csv_row "Host downstream analysis" "ssGSEA" "gmt_file" "$SSGSEA_GMT" "GMT gene set database used by ssGSEA"
        csv_row "Host downstream analysis" "ssGSEA" "gmt_mode" "$SSGSEA_GMT_MODE" "Whether ssGSEA used default MSigDB or a custom GMT"
        csv_row "Host downstream analysis" "ssGSEA" "protein_fasta" "${SSGSEA_PROTEIN_FASTA:-none}" "Protein FASTA used to build custom ssGSEA GMT"
        csv_row "Host downstream analysis" "ssGSEA" "annotation_gff" "${SSGSEA_ANNOTATION_GFF:-none}" "GFF/GTF annotation used to build custom ssGSEA GMT"
        # HUMAnN
        csv_row \
            "Functional profiling" \
            "HUMAnN" \
            "read_layout" \
            "${READ_LAYOUT:-pending_auto_detection}" \
            "Effective sequencing layout used to prepare HUMAnN input"

        csv_row \
            "Functional profiling" \
            "HUMAnN" \
            "single_end_input" \
            "SAMPLE_non-host.fq -> HUMAnN_output/input/SAMPLE.fq" \
            "Single-end final non-host FASTQ is copied to the normalized HUMAnN input"

        csv_row \
            "Functional profiling" \
            "HUMAnN" \
            "paired_end_input" \
            "SAMPLE_non-host_1.fq + SAMPLE_non-host_2.fq -> HUMAnN_output/input/SAMPLE.fq" \
            "Paired-end final non-host FASTQ mates are concatenated for HUMAnN"

        csv_row \
            "Functional profiling" \
            "HUMAnN" \
            "paired_end_handling" \
            "concatenate_R1_then_R2" \
            "HUMAnN analyzes reads independently rather than as synchronized pairs"

        csv_row \
            "Functional profiling" \
            "HUMAnN" \
            "threads" \
            "$threads" \
            "Threads used by HUMAnN"

        csv_row \
            "Functional profiling" \
            "HUMAnN" \
            "renormalization_units" \
            "relab" \
            "HUMAnN output is renormalized to relative abundance"

        # HAllA
        csv_row "Association analysis" "HAllA" "metric" "$pdm" "Association metric used for microbiome-host associations"
        csv_row "Association analysis" "HAllA" "RUN_EXTRA_PEARSON" "${RUN_EXTRA_PEARSON:-1}" "Whether extra Pearson HAllA analysis is enabled later in the script"
        csv_row "Association analysis" "HAllA" "RUN_FULL_HALLAGRAM" "${RUN_FULL_HALLAGRAM:-0}" "Whether full hallagram plots are enabled"
        csv_row "Association analysis" "HAllA" "HALLA_DIAGNOSTIC" "${HALLA_DIAGNOSTIC:-0}" "Whether HAllA diagnostic plot is enabled"

        # Outputs
        csv_row "Key outputs" "Kraken2" "raw_global_composition" "$outputdr/kraken/kraken_global_read_composition_raw.tsv" "Global host/microbiome/unclassified composition before optional contaminant removal"
        csv_row "Key outputs" "Kraken2" "final_global_composition" "$outputdr/kraken/kraken_global_read_composition_final.tsv" "Global host/microbiome/unclassified composition after optional contaminant removal"
        csv_row "Key outputs" "Bracken" "species_table" "$outputdr/bracken_species_all" "Combined Bracken species table"
        csv_row "Key outputs" "HAllA" "microbiome_input" "$outputdr/halla/Microbiomes.txt" "Microbiome matrix used by HAllA"

        # Software paths
        csv_row "Software path" "kraken2" "path" "$(get_tool_path kraken2)" "Executable path"
        csv_row "Software path" "bracken" "path" "$(get_tool_path bracken)" "Executable path"
        csv_row "Software path" "fastp" "path" "$(get_tool_path fastp)" "Executable path"
        csv_row "Software path" "pigz" "path" "$(get_tool_path pigz)" "Executable path"
        csv_row "Software path" "python" "path" "$(get_tool_path python)" "Executable path"
        csv_row "Software path" "Rscript" "path" "$(get_tool_path Rscript)" "Executable path"
        csv_row "Software path" "humann" "path" "$(get_tool_path humann)" "Executable path"
        csv_row "Software path" "magicblast" "path" "$(get_tool_path magicblast)" "Executable path"
        csv_row "Software path" "hisat2" "path" "$(get_tool_path hisat2)" "Executable path"
        csv_row "Software path" "featureCounts" "path" "$(get_tool_path featureCounts)" "Executable path"
        csv_row "Software path" "samtools" "path" "$(get_tool_path samtools)" "Executable path"

        # Software versions, best effort
        csv_row "Software version" "kraken2" "version" "$(get_tool_version 'kraken2 --version')" "Best-effort version capture"
        csv_row "Software version" "bracken" "version" "$(get_tool_version 'bracken -v')" "Best-effort version capture"
        csv_row "Software version" "fastp" "version" "$(get_tool_version 'fastp --version')" "Best-effort version capture"
        csv_row "Software version" "python" "version" "$(get_tool_version 'python --version')" "Best-effort version capture"
        csv_row "Software version" "Rscript" "version" "$(get_tool_version 'Rscript --version')" "Best-effort version capture"
        csv_row "Software version" "humann" "version" "$(get_tool_version 'humann --version')" "Best-effort version capture"
        csv_row "Software version" "magicblast" "version" "$(get_tool_version 'magicblast -version')" "Best-effort version capture"
        csv_row "Software version" "hisat2" "version" "$(get_tool_version 'hisat2 --version')" "Best-effort version capture"
        csv_row "Software version" "featureCounts" "version" "$(get_tool_version 'featureCounts -v')" "Best-effort version capture"
        csv_row "Software version" "samtools" "version" "$(get_tool_version 'samtools --version')" "Best-effort version capture"

    } > "$methods_csv"

    echo "${g}[OK] Methods/run parameters exported to:${w}"
    echo "  $methods_csv"
}

write_methods_log

# ------------------------------------------------------------
# Exploratory-only taxonomic figures
# ------------------------------------------------------------
# Runs only when NO_COMPARISON=1.
# Uses raw combined Bracken species table:
#   $outputdr/temp/bracken_raw_results/bracken_species_all
# ------------------------------------------------------------

run_exploratory_taxonomic_figures() {
    if [[ "${NO_COMPARISON:-0}" != "1" && "${RUN_EXPLORATORY_FIGURES_IN_COMPARISON:-1}" != "1" ]]; then
        echo "${y}[INFO] Comparison mode detected and RUN_EXPLORATORY_FIGURES_IN_COMPARISON=0. Skipping exploratory taxonomic figures.${w}"
        return 0
    fi

    if [[ "${RUN_EXPLORATORY_TAXONOMIC_FIGURES:-1}" != "1" ]]; then
        echo "${y}[INFO] RUN_EXPLORATORY_TAXONOMIC_FIGURES=0. Skipping exploratory taxonomic figures.${w}"
        return 0
    fi

    echo "============================================================"
    echo "${g}[EXPLORATORY FIGURES]${w}"
    echo "============================================================"

    local fig_script_dir="$MTDIR/aux_scripts/exploratory"

    local pheatmap_script="$fig_script_dir/MTD.taxonomic_pheatmap.R"
    local stacked_script="$fig_script_dir/MTD.taxonomic_stacked_bar.R"
    local prevalence_script="$fig_script_dir/MTD_exploratory_prevalence_abundance.R"
    local detected_microbiome_script="$fig_script_dir/MTD_exploratory_detected_microbiome_rank.R"
    local core_script="$fig_script_dir/MTD_exploratory_core_microbiome.R"
    local species_overlap_script="$fig_script_dir/MTD_species_overlap_venn_euler.R"
    local detected_species_pie_script="$fig_script_dir/MTD_detected_species_pie_by_phylum.R"
    local alpha_script="$fig_script_dir/MTD_exploratory_alpha_diversity.R"
    local beta_script="$fig_script_dir/MTD_exploratory_beta_diversity.R"
    local read_qc_script="$fig_script_dir/MTD_exploratory_read_composition_qc.R"
    local matrix_qc_script="$fig_script_dir/MTD_exploratory_matrix_qc.R"
    local detected_microbiome_read_extract_script="$fig_script_dir/MTD_extract_reads_by_detected_microbiome.py"

    # Default rank is genus, but you can override it:
    # EXPLORATORY_TAX_RANK=species bash $SCRIPT_NAME ...
    local tax_rank="${EXPLORATORY_TAX_RANK:-species}"

    local tax_input="$outputdr/temp/bracken_raw_results/bracken_${tax_rank}_all"

    # Fallback in case the raw Bracken table has not been moved yet
    if [[ ! -s "$tax_input" && -s "$outputdr/bracken_${tax_rank}_all" ]]; then
        tax_input="$outputdr/bracken_${tax_rank}_all"
    fi

    if [[ ! -s "$tax_input" ]]; then
        echo "${y}[WARNING] Taxonomic figure input not found. Skipping figures.${w}"
        echo "Expected:"
        echo "  $outputdr/temp/bracken_raw_results/bracken_${tax_rank}_all"
        echo "Fallback checked:"
        echo "  $outputdr/bracken_${tax_rank}_all"
        return 0
    fi

    if [[ ! -s "$samplesheet_file" ]]; then
        echo "${y}[WARNING] Samplesheet not found. Skipping taxonomic figures.${w}"
        echo "Expected:"
        echo "  $samplesheet_file"
        return 0
    fi

    # New exploratory output structure
    local exploratory_dir="$outputdr/exploratory"
    local taxonomy_dir="$exploratory_dir/taxonomy"

    local heatmap_base="$taxonomy_dir/heatmap"
    local stacked_base="$taxonomy_dir/stacked_bar"
    local prevalence_base="$taxonomy_dir/prevalence_abundance"
    local detected_microbiome_base="$taxonomy_dir/detected_microbiome"
    local detected_microbiome_reads_base="$taxonomy_dir/detected_microbiome_extracted_reads"
    local core_base="$taxonomy_dir/core_microbiome"
    local species_overlap_base="$taxonomy_dir/species_overlap"
    local detected_species_pie_base="$taxonomy_dir/detected_microbiome_pie"
    local alpha_base="$taxonomy_dir/alpha_diversity"
    local beta_base="$taxonomy_dir/beta_diversity"
    local microbiome_matrix_qc_base="$taxonomy_dir/microbiome_abundance_qc"

    local pipeline_qc_dir="$exploratory_dir/pipeline_qc"
    local read_qc_base="$pipeline_qc_dir/read_composition"

    local log_dir="$taxonomy_dir/logs"

mkdir -p "$heatmap_base" "$stacked_base" "$prevalence_base" "$detected_microbiome_base" "$detected_microbiome_reads_base" "$species_overlap_base" "$detected_species_pie_base" "$core_base" "$alpha_base" "$beta_base" "$microbiome_matrix_qc_base" "$pipeline_qc_dir" "$read_qc_base" "$log_dir"

    echo "[INFO] Taxonomic input:"
    echo "  $tax_input"
    echo "[INFO] Samplesheet:"
    echo "  $samplesheet_file"
    echo "[INFO] Taxonomic rank:"
    echo "  $tax_rank"
    echo "[INFO] Output directory:"
    echo "  $taxonomy_dir"
    echo "[INFO] Figure scripts:"
    echo "  Heatmap:     $pheatmap_script"
    echo "  Stacked bar: $stacked_script"
    echo "  Prevalence:  $prevalence_script"
    echo "  Detected microbiome ranking: $detected_microbiome_script"
    echo "  Core microbiome: $core_script"
    echo "  Alpha diversity: $alpha_script"
    echo "  Beta diversity:  $beta_script"
    echo "  Matrix QC:       $matrix_qc_script"
    echo "  Read composition QC: $read_qc_script"
    conda deactivate
    conda activate R412

    # ------------------------------------------------------------
    # 1) Taxonomic heatmap
    # ------------------------------------------------------------

    if [[ -s "$pheatmap_script" ]]; then
        for tax_mode in absolute relative; do

            local heatmap_out="$heatmap_base/${tax_rank}_${tax_mode}"
            local heatmap_log="$log_dir/MTD.taxonomic_pheatmap.${tax_rank}.${tax_mode}.log"
            local heatmap_title=""

            if [[ "$tax_mode" == "absolute" ]]; then
                heatmap_title="Absolute abundance heatmap with hierarchical clustering"
            else
                heatmap_title="Relative abundance heatmap with hierarchical clustering"
            fi

            mkdir -p "$heatmap_out"

            echo "------------------------------------------------------------"
            echo "[RUN] Taxonomic heatmap"
            echo "Rank: $tax_rank"
            echo "Mode: $tax_mode"
            echo "Output: $heatmap_out"
            echo "------------------------------------------------------------"

            if ! Rscript "$pheatmap_script" \
                --input "$tax_input" \
                --samplesheet "$samplesheet_file" \
                --rank "$tax_rank" \
                --top 100 \
                --mode "$tax_mode" \
                --transform log10 \
                --cluster_samples yes \
                --cluster_taxa yes \
                --cluster_distance correlation \
                --cluster_method complete \
                --title "$heatmap_title" \
                --width 8 \
                --height 9 \
                --fontsize_row 8 \
                --fontsize_col 10 \
                --output_dir "$heatmap_out" \
                > "$heatmap_log" 2>&1
            then
                echo "${y}[WARNING] Taxonomic heatmap failed for mode: $tax_mode. Continuing pipeline.${w}"
                echo "Log:"
                echo "  $heatmap_log"
                tail -n 40 "$heatmap_log" || true
            else
                echo "${g}[OK] Taxonomic heatmap generated:${w}"
                echo "  $heatmap_out"
            fi
        done
    else
        echo "${y}[WARNING] Heatmap script not found. Skipping heatmap.${w}"
        echo "Expected:"
        echo "  $pheatmap_script"
    fi

    # ------------------------------------------------------------
    # 2) Taxonomic stacked bar
    # ------------------------------------------------------------

    if [[ -s "$stacked_script" ]]; then
        for tax_mode in absolute relative; do

            local stacked_out="$stacked_base/${tax_rank}_${tax_mode}"
            local stacked_log="$log_dir/MTD.taxonomic_stacked_bar.${tax_rank}.${tax_mode}.log"
            local stacked_title=""

            if [[ "$tax_mode" == "absolute" ]]; then
                stacked_title="Genus-level taxonomic profile - absolute abundance"
            else
                stacked_title="Genus-level taxonomic profile - relative abundance"
            fi

            # If you switch rank to species/phylum/etc, make title follow the rank
            stacked_title="${tax_rank^}-level taxonomic profile - ${tax_mode} abundance"

            mkdir -p "$stacked_out"

            echo "------------------------------------------------------------"
            echo "[RUN] Taxonomic stacked bar"
            echo "Rank: $tax_rank"
            echo "Mode: $tax_mode"
            echo "Output: $stacked_out"
            echo "------------------------------------------------------------"

            if ! Rscript "$stacked_script" \
                --input "$tax_input" \
                --samplesheet "$samplesheet_file" \
                --rank "$tax_rank" \
                --top 50 \
                --mode "$tax_mode" \
                --taxon_format scientific \
                --title "$stacked_title" \
                --width auto \
                --height auto \
                --bar_width 0.40 \
                --legend_position bottom \
                --legend_ncol auto \
                --legend_max_rows 8 \
                --legend_text_size 10 \
                --legend_title_size 11 \
                --legend_key_size 0.36 \
                --x_text_size 11 \
                --y_text_size 11 \
                --axis_title_size 12 \
                --plot_title_size 14 \
                --output_dir "$stacked_out" \
                > "$stacked_log" 2>&1
            then
                echo "${y}[WARNING] Taxonomic stacked bar failed for mode: $tax_mode. Continuing pipeline.${w}"
                echo "Log:"
                echo "  $stacked_log"
                tail -n 40 "$stacked_log" || true
            else
                echo "${g}[OK] Taxonomic stacked bar generated:${w}"
                echo "  $stacked_out"
            fi
        done
    else
        echo "${y}[WARNING] Stacked bar script not found. Skipping stacked bar.${w}"
        echo "Expected:"
        echo "  $stacked_script"
    fi

    # ------------------------------------------------------------
    # 3) Prevalence vs abundance
    # ------------------------------------------------------------

    if [[ -s "$prevalence_script" ]]; then
        for tax_mode in absolute relative; do

            local prevalence_out="$prevalence_base/${tax_rank}_${tax_mode}"
            local prevalence_log="$log_dir/MTD.prevalence_abundance.${tax_rank}.${tax_mode}.log"

            mkdir -p "$prevalence_out"

            echo "------------------------------------------------------------"
            echo "[RUN] Prevalence vs abundance"
            echo "Rank: $tax_rank"
            echo "Mode: $tax_mode"
            echo "Output: $prevalence_out"
            echo "------------------------------------------------------------"

            if ! Rscript "$prevalence_script" \
                --input "$tax_input" \
                --output "$prevalence_out" \
                --rank "$tax_rank" \
                --mode "$tax_mode" \
                --top_labels 20 \
                > "$prevalence_log" 2>&1
            then
                echo "${y}[WARNING] Prevalence vs abundance failed for mode: $tax_mode. Continuing pipeline.${w}"
                echo "Log:"
                echo "  $prevalence_log"
                tail -n 40 "$prevalence_log" || true
            else
                echo "${g}[OK] Prevalence vs abundance generated:${w}"
                echo "  $prevalence_out"
            fi
        done
    else
        echo "${y}[WARNING] Prevalence script not found. Skipping prevalence vs abundance.${w}"
        echo "Expected:"
        echo "  $prevalence_script"
    fi
    # ------------------------------------------------------------
    # 3B) Detected microbiome ranked tables
    # ------------------------------------------------------------

    if [[ -s "$detected_microbiome_script" ]]; then

        for detected_mode in absolute relative; do

            local detected_in="$prevalence_base/${tax_rank}_${detected_mode}/prevalence_vs_abundance_${tax_rank}_${detected_mode}.tsv"
            local detected_out="$detected_microbiome_base/${tax_rank}_${detected_mode}_with_samples_distribution"
            local detected_log="$log_dir/MTD.detected_microbiome_ranked.${tax_rank}.${detected_mode}.log"

            mkdir -p "$detected_out"

            echo "------------------------------------------------------------"
            echo "[RUN] Detected microbiome ranked tables"
            echo "Rank: $tax_rank"
            echo "Mode: $detected_mode"
            echo "Input: $detected_in"
            echo "Output: $detected_out"
            echo "------------------------------------------------------------"

            if [[ ! -s "$detected_in" ]]; then
                echo "${y}[WARNING] Prevalence ${detected_mode} table not found. Skipping detected microbiome ranking for this mode.${w}"
                echo "Expected:"
                echo "  $detected_in"

            elif ! Rscript "$detected_microbiome_script" \
                --input "$detected_in" \
                --abundance_input "$tax_input" \
                --samplesheet "$samplesheet_file" \
                --output "$detected_out" \
                --rank "$tax_rank" \
                --mode "$detected_mode" \
                --presence_threshold 0 \
                --w_prevalence 0.45 \
                --w_mean 0.30 \
                --w_max 0.15 \
                --w_total 0.10 \
                --core_threshold 75 \
                --plot_top 20 \
                > "$detected_log" 2>&1
            then
                echo "${y}[WARNING] Detected microbiome ranking failed for mode: $detected_mode. Continuing pipeline.${w}"
                echo "Log:"
                echo "  $detected_log"
                tail -n 40 "$detected_log" || true

            else
                echo "${g}[OK] Detected microbiome ranked tables generated:${w}"
                echo "  $detected_out"
            fi

        done

    else
        echo "${y}[WARNING] Detected microbiome ranking script not found. Skipping ranking tables.${w}"
        echo "Expected:"
        echo "  $detected_microbiome_script"
    fi

    # ------------------------------------------------------------
    # 3B.1) Extract reads for detected microbiome taxa
    # ------------------------------------------------------------

    if [[ "${EXTRACT_MICROBIOME_READS:-0}" == "1" ]]; then

        local extract_mode="absolute"
        local extract_ranked_table="$detected_microbiome_base/${tax_rank}_${extract_mode}_with_samples_distribution/detected_microbiome_${tax_rank}_ranked_with_samples.tsv"
        local extract_out="$detected_microbiome_reads_base/${tax_rank}_${extract_mode}_extracted_reads"
        local extract_log="$log_dir/MTD.detected_microbiome_read_extraction.${tax_rank}.${extract_mode}.log"

        mkdir -p "$extract_out"

        echo "------------------------------------------------------------"
        echo "[RUN] Extract reads for detected microbiome taxa"
        echo "Rank: $tax_rank"
        echo "Mode: $extract_mode"
        echo "Ranked table: $extract_ranked_table"
        echo "Bracken table: $tax_input"
        echo "Output: $extract_out"
        echo "Top N: $EXTRACT_MICROBIOME_READS_TOP_N"
        echo "Minimum abundance: $EXTRACT_MICROBIOME_READS_MIN_ABUNDANCE"
        echo "------------------------------------------------------------"

        if [[ ! -s "$detected_microbiome_read_extract_script" ]]; then
            echo "${y}[WARNING] Detected microbiome read extraction script not found. Skipping read extraction.${w}"
            echo "Expected:"
            echo "  $detected_microbiome_read_extract_script"

        elif [[ ! -s "$MTDIR/Tools/KrakenTools/extract_kraken_reads.py" ]]; then
            echo "${y}[WARNING] KrakenTools extract_kraken_reads.py not found. Skipping read extraction.${w}"
            echo "Expected:"
            echo "  $MTDIR/Tools/KrakenTools/extract_kraken_reads.py"

        elif [[ ! -s "$extract_ranked_table" ]]; then
            echo "${y}[WARNING] Absolute ranked-with-samples table not found. Skipping read extraction.${w}"
            echo "Expected:"
            echo "  $extract_ranked_table"

        elif [[ ! -s "$tax_input" ]]; then
            echo "${y}[WARNING] Bracken combined table not found. Skipping read extraction.${w}"
            echo "Expected:"
            echo "  $tax_input"

        else
            conda deactivate
            conda activate MTD

            export PYTHONNOUSERSITE=1
            unset PYTHONPATH
            unset PYTHONHOME

            extract_args=(
                "$detected_microbiome_read_extract_script"
                --ranked_table "$extract_ranked_table"
                --bracken_table "$tax_input"
                --samplesheet "$samplesheet_file"
                --temp_dir "$outputdr/temp"
                --read-layout "$READ_LAYOUT"
                --krakentools_script "$MTDIR/Tools/KrakenTools/extract_kraken_reads.py"
                --output_dir "$extract_out"
                --top_n "$EXTRACT_MICROBIOME_READS_TOP_N"
                --min_abundance "$EXTRACT_MICROBIOME_READS_MIN_ABUNDANCE"
            )

            if [[ "${EXTRACT_MICROBIOME_READS_INCLUDE_CHILDREN:-1}" == "1" ]]; then
                extract_args+=( --include_children )
            fi

            if [[ "${EXTRACT_MICROBIOME_READS_OVERWRITE:-1}" == "1" ]]; then
                extract_args+=( --overwrite )
            fi

            if ! python "${extract_args[@]}" > "$extract_log" 2>&1
            then
                echo "${y}[WARNING] Detected microbiome read extraction failed. Continuing pipeline.${w}"
                echo "Log:"
                echo "  $extract_log"
                tail -n 60 "$extract_log" || true
            else
                echo "${g}[OK] Detected microbiome reads extracted:${w}"
                echo "  $extract_out"

                if [[ -s "$extract_out/extraction_summary.tsv" ]]; then
                    echo "[INFO] Extraction status summary:"
                    awk -F'\t' '
                        NR > 1 {
                            status[$9]++
                        }
                        END {
                            for (s in status) {
                                print "  " status[s], s
                            }
                        }
                    ' "$extract_out/extraction_summary.tsv" | sort -nr || true
                fi
            fi

            conda deactivate
            conda activate R412
        fi

    else
        echo "${y}[INFO] --extract-microbiome-reads not provided. Skipping detected microbiome read extraction.${w}"
    fi

    # ------------------------------------------------------------
    # 3C) Detected species pie chart by phylum
    # ------------------------------------------------------------

    if [[ "$tax_rank" == "species" && -s "$detected_species_pie_script" ]]; then

        local pie_ranked_file="$detected_microbiome_base/${tax_rank}_relative_with_samples_distribution/detected_microbiome_${tax_rank}_ranked_by_importance.tsv"
        local pie_taxonomy_file=""
        local pie_out="$detected_species_pie_base/${tax_rank}_by_phylum"
        local pie_log="$log_dir/MTD.detected_species_pie_by_phylum.${tax_rank}.log"
        local pie_category_level="${PIE_CATEGORY_LEVEL:-auto}"
        local pie_top_n="${PIE_TOP_N:-auto}"
        local pie_auto_min_categories="${PIE_AUTO_MIN_CATEGORIES:-4}"
        local pie_auto_max_categories="${PIE_AUTO_MAX_CATEGORIES:-18}"
        local pie_auto_target_categories="${PIE_AUTO_TARGET_CATEGORIES:-12}"

        mkdir -p "$pie_out"

        if [[ -s "$outputdr/Combined.mpa" ]]; then
            pie_taxonomy_file="$outputdr/Combined.mpa"
        elif [[ -s "$outputdr/graphlan/graphlan_input.clean.tsv" ]]; then
            pie_taxonomy_file="$outputdr/graphlan/graphlan_input.clean.tsv"
        fi

        echo "------------------------------------------------------------"
        echo "[RUN] Detected species pie chart by phylum"
        echo "Ranked file:   $pie_ranked_file"
        echo "Taxonomy file: $pie_taxonomy_file"
        echo "Output:        $pie_out"
        echo "------------------------------------------------------------"
        PIE_CATEGORY_LEVEL=auto
        if [[ ! -s "$pie_ranked_file" ]]; then
            echo "${y}[WARNING] Detected microbiome ranked table not found. Skipping species pie chart.${w}"
            echo "Expected:"
            echo "  $pie_ranked_file"
        elif [[ -z "$pie_taxonomy_file" || ! -s "$pie_taxonomy_file" ]]; then
            echo "${y}[WARNING] Taxonomy file for species-to-phylum mapping not found. Skipping species pie chart.${w}"
            echo "Checked:"
            echo "  $outputdr/Combined.mpa"
            echo "  $outputdr/graphlan/graphlan_input.clean.tsv"
        else
            if ! Rscript "$detected_species_pie_script" \
                  --ranked "$pie_ranked_file" \
                  --taxonomy "$pie_taxonomy_file" \
                  --category_level "$pie_category_level" \
                  --title auto \
                  --legend_title auto \
                  --panel_label "" \
                  --output_prefix "$pie_out/detected_species_by_${pie_category_level}_level" \
                  --top_n "$pie_top_n" \
                  --auto_min_categories "$pie_auto_min_categories" \
                  --auto_max_categories "$pie_auto_max_categories" \
                  --auto_target_categories "$pie_auto_target_categories"
                > "$pie_log" 2>&1
            then
                echo "${y}[WARNING] Detected species pie chart failed. Continuing pipeline.${w}"
                echo "Log:"
                echo "  $pie_log"
                tail -n 40 "$pie_log" || true
            else
                echo "${g}[OK] Detected species pie chart generated:${w}"
                echo "  $pie_out"
            fi
        fi

    else
        echo "${y}[INFO] Species pie chart by phylum skipped.${w}"
        echo "Reason: tax_rank=$tax_rank or script not found:"
        echo "  $detected_species_pie_script"
    fi


    # ------------------------------------------------------------
    # 3D) Species overlap Venn/Euler between two groups
    # ------------------------------------------------------------

    if [[ "$tax_rank" == "species" && -s "$species_overlap_script" ]]; then

        local species_overlap_out="$species_overlap_base"
        local species_overlap_log="$log_dir/MTD.species_overlap_venn_euler.${tax_rank}.log"

        mkdir -p "$species_overlap_out"

        # Detect group names from samplesheet column 2.
        local overlap_groups
        overlap_groups=$(awk -F',' '
            NR > 1 {
                g=$2
                gsub(/^[ \t\r\n]+|[ \t\r\n]+$/, "", g)
                if (g != "") groups[g]=1
            }
            END {
                for (g in groups) print g
            }
        ' "$samplesheet_file" | sort)

        local overlap_group_count
        overlap_group_count=$(printf "%s\n" "$overlap_groups" | awk 'NF > 0 {n++} END {print n+0}')

        if [[ "$overlap_group_count" -ne 2 ]]; then
            echo "${y}[INFO] Species Venn/Euler requires exactly 2 groups. Skipping overlap plot.${w}"
            echo "Detected groups:"
            printf "%s\n" "$overlap_groups"
        else
            local group1
            local group2

            group1=$(printf "%s\n" "$overlap_groups" | sed -n '1p')
            group2=$(printf "%s\n" "$overlap_groups" | sed -n '2p')

            echo "------------------------------------------------------------"
            echo "[RUN] Species overlap Venn/Euler"
            echo "Input:  $tax_input"
            echo "Group1: $group1"
            echo "Group2: $group2"
            echo "Output: $species_overlap_out"
            echo "------------------------------------------------------------"

            if ! conda run --no-capture-output -n MTD_fastp -- \
            Rscript "$species_overlap_script" \
                --input "$tax_input" \
                --samplesheet "$samplesheet_file" \
                --group1 "$group1" \
                --group2 "$group2" \
                --title "Total species detected" \
                --output_prefix "$species_overlap_out/${group1}_vs_${group2}_species_overlap" \
                --write_lists \
                > "$species_overlap_log" 2>&1
            then
                echo "${y}[WARNING] Species Venn/Euler failed. Continuing pipeline.${w}"
                echo "Log:"
                echo "  $species_overlap_log"
                tail -n 40 "$species_overlap_log" || true
            else
                echo "${g}[OK] Species Venn/Euler generated:${w}"
                echo "  $species_overlap_out"
            fi
        fi

    else
        echo "${y}[INFO] Species Venn/Euler skipped.${w}"
        echo "Reason: tax_rank=$tax_rank or script not found:"
        echo "  $species_overlap_script"
    fi
    # ------------------------------------------------------------
    # 4) Core microbiome
    # ------------------------------------------------------------

    if [[ -s "$core_script" ]]; then
        for tax_mode in absolute relative; do

            local core_out="$core_base/${tax_rank}_${tax_mode}"
            local core_log="$log_dir/MTD.core_microbiome.${tax_rank}.${tax_mode}.log"

            mkdir -p "$core_out"

            echo "------------------------------------------------------------"
            echo "[RUN] Core microbiome"
            echo "Rank: $tax_rank"
            echo "Mode: $tax_mode"
            echo "Output: $core_out"
            echo "------------------------------------------------------------"

            if ! Rscript "$core_script" \
                --input "$tax_input" \
                --output "$core_out" \
                --rank "$tax_rank" \
                --mode "$tax_mode" \
                --thresholds 25,50,75,90,100 \
                --core_threshold 75 \
                > "$core_log" 2>&1
            then
                echo "${y}[WARNING] Core microbiome failed for mode: $tax_mode. Continuing pipeline.${w}"
                echo "Log:"
                echo "  $core_log"
                tail -n 40 "$core_log" || true
            else
                echo "${g}[OK] Core microbiome generated:${w}"
                echo "  $core_out"
            fi
        done
    else
        echo "${y}[WARNING] Core microbiome script not found. Skipping core microbiome.${w}"
        echo "Expected:"
        echo "  $core_script"
    fi

    # ------------------------------------------------------------
    # 5) Alpha diversity
    # ------------------------------------------------------------

    if [[ -s "$alpha_script" ]]; then

        local alpha_mode="relative"
        local alpha_out="$alpha_base/${tax_rank}_${alpha_mode}"
        local alpha_log="$log_dir/MTD.alpha_diversity.${tax_rank}.${alpha_mode}.log"

        mkdir -p "$alpha_out"

        echo "------------------------------------------------------------"
        echo "[RUN] Alpha diversity"
        echo "Rank: $tax_rank"
        echo "Mode: $alpha_mode"
        echo "Output: $alpha_out"
        echo "------------------------------------------------------------"

        if ! Rscript "$alpha_script" \
            --input "$tax_input" \
            --samplesheet "$samplesheet_file" \
            --output "$alpha_out" \
            --rank "$tax_rank" \
            --mode "$alpha_mode" \
            > "$alpha_log" 2>&1
        then
            echo "${y}[WARNING] Alpha diversity failed. Continuing pipeline.${w}"
            echo "Log:"
            echo "  $alpha_log"
            tail -n 40 "$alpha_log" || true
        else
            echo "${g}[OK] Alpha diversity generated:${w}"
            echo "  $alpha_out"
        fi

    else
        echo "${y}[WARNING] Alpha diversity script not found. Skipping alpha diversity.${w}"
        echo "Expected:"
        echo "  $alpha_script"
    fi

    # ------------------------------------------------------------
    # 6) Beta diversity / ordination
    # ------------------------------------------------------------

    if [[ -s "$beta_script" ]]; then

        local beta_mode="relative"
        local beta_out="$beta_base/${tax_rank}_${beta_mode}"
        local beta_log="$log_dir/MTD.beta_diversity.${tax_rank}.${beta_mode}.log"

        mkdir -p "$beta_out"

        echo "------------------------------------------------------------"
        echo "[RUN] Beta diversity / ordination"
        echo "Rank: $tax_rank"
        echo "Mode: $beta_mode"
        echo "Output: $beta_out"
        echo "------------------------------------------------------------"

        if ! Rscript "$beta_script" \
            --input "$tax_input" \
            --samplesheet "$samplesheet_file" \
            --output "$beta_out" \
            --rank "$tax_rank" \
            --mode "$beta_mode" \
            --distance bray \
            --transform none \
            --label_samples yes \
            --ellipse no \
            > "$beta_log" 2>&1
        then
            echo "${y}[WARNING] Beta diversity failed. Continuing pipeline.${w}"
            echo "Log:"
            echo "  $beta_log"
            tail -n 40 "$beta_log" || true
        else
            echo "${g}[OK] Beta diversity generated:${w}"
            echo "  $beta_out"
        fi

    else
        echo "${y}[WARNING] Beta diversity script not found. Skipping beta diversity.${w}"
        echo "Expected:"
        echo "  $beta_script"
    fi

    # ------------------------------------------------------------
    # 7) Pipeline QC: read composition
    # ------------------------------------------------------------

    if [[ -s "$read_qc_script" ]]; then

        local read_qc_input="$outputdr/kraken/kraken_global_read_composition_final.tsv"
        local read_qc_out="$read_qc_base/final"
        local read_qc_log="$log_dir/MTD.read_composition_qc.final.log"

        mkdir -p "$read_qc_out"

        echo "------------------------------------------------------------"
        echo "[RUN] Pipeline QC: read composition"
        echo "Input:  $read_qc_input"
        echo "Output: $read_qc_out"
        echo "------------------------------------------------------------"

        if [[ ! -s "$read_qc_input" ]]; then
            echo "${y}[WARNING] Read composition table not found. Skipping read composition QC.${w}"
            echo "Expected:"
            echo "  $read_qc_input"
        elif ! Rscript "$read_qc_script" \
            --input "$read_qc_input" \
            --output "$read_qc_out" \
            --label final \
            --outlier_method mad \
            --outlier_cutoff 3.5 \
            > "$read_qc_log" 2>&1
        then
            echo "${y}[WARNING] Read composition QC failed. Continuing pipeline.${w}"
            echo "Log:"
            echo "  $read_qc_log"
            tail -n 40 "$read_qc_log" || true
        else
            echo "${g}[OK] Read composition QC generated:${w}"
            echo "  $read_qc_out"
        fi

    else
        echo "${y}[WARNING] Read composition QC script not found. Skipping pipeline QC.${w}"
        echo "Expected:"
        echo "  $read_qc_script"
    fi

    # ------------------------------------------------------------
    # 7) Microbiome abundance matrix QC
    # ------------------------------------------------------------

    if [[ -s "$matrix_qc_script" ]]; then

        local microbiome_qc_mode="relative"
        local microbiome_qc_out="$microbiome_matrix_qc_base/${tax_rank}_${microbiome_qc_mode}"
        local microbiome_qc_log="$log_dir/MTD.microbiome_matrix_qc.${tax_rank}.${microbiome_qc_mode}.log"

        mkdir -p "$microbiome_qc_out"

        echo "------------------------------------------------------------"
        echo "[RUN] Microbiome abundance matrix QC"
        echo "Rank: $tax_rank"
        echo "Mode: $microbiome_qc_mode"
        echo "Output: $microbiome_qc_out"
        echo "------------------------------------------------------------"

        if ! Rscript "$matrix_qc_script" \
            --input "$tax_input" \
            --samplesheet "$samplesheet_file" \
            --output "$microbiome_qc_out" \
            --label "microbiome_${tax_rank}" \
            --matrix_type bracken \
            --normalization relative \
            --top_variable 50 \
            --pca_top_features 5000 \
            > "$microbiome_qc_log" 2>&1
        then
            echo "${y}[WARNING] Microbiome abundance matrix QC failed. Continuing pipeline.${w}"
            echo "Log:"
            echo "  $microbiome_qc_log"
            tail -n 40 "$microbiome_qc_log" || true
        else
            echo "${g}[OK] Microbiome abundance matrix QC generated:${w}"
            echo "  $microbiome_qc_out"
        fi

    else
        echo "${y}[WARNING] Matrix QC script not found. Skipping microbiome abundance QC.${w}"
        echo "Expected:"
        echo "  $matrix_qc_script"
    fi

    conda deactivate
    conda activate MTD

    echo "============================================================"
    echo "${g}[OK] Exploratory taxonomic figure step finished${w}"
    echo "Main output:"
    echo "  $taxonomy_dir"
    echo "============================================================"
}

# ------------------------------------------------------------
# Exploratory-only host expression matrix QC
# ------------------------------------------------------------
# Runs only when NO_COMPARISON=1.
#
# This is NOT differential expression analysis.
# It generates non-comparison host expression QC:
#   - PCA of host expression
#   - top variable genes heatmap
#   - sample correlation heatmap
#   - detected genes per sample
#
# Input:
#   $outputdr/Host_DEG/host_counts_featureCounts_matrix.txt
#
# Script:
#   $MTDIR/aux_scripts/exploratory/MTD_exploratory_matrix_qc.R
# ------------------------------------------------------------

run_exploratory_host_expression_qc() {
    if [[ "${NO_COMPARISON:-0}" != "1" ]]; then
        echo "${y}[INFO] Comparison mode detected. Skipping exploratory host expression QC.${w}"
        return 0
    fi

    if [[ "${RUN_EXPLORATORY_HOST_EXPRESSION_QC:-1}" != "1" ]]; then
        echo "${y}[INFO] RUN_EXPLORATORY_HOST_EXPRESSION_QC=0. Skipping exploratory host expression QC.${w}"
        return 0
    fi

    local host_matrix_qc_script="$MTDIR/aux_scripts/exploratory/MTD_exploratory_matrix_qc.R"
    local host_matrix_qc_input="$outputdr/Host_DEG/host_counts_featureCounts_matrix.txt"
    local host_matrix_qc_out="$outputdr/exploratory/host_expression/matrix_qc"
    local host_matrix_qc_log_dir="$outputdr/exploratory/host_expression/logs"
    local host_matrix_qc_log="$host_matrix_qc_log_dir/MTD.host_expression_matrix_qc.log"
    local host_counts_summary="$outputdr/host_counts.txt.summary"

    mkdir -p "$host_matrix_qc_out" "$host_matrix_qc_log_dir"

    echo "============================================================"
    echo "${g}[EXPLORATORY HOST EXPRESSION QC]${w}"
    echo "============================================================"
    echo "Input matrix:"
    echo "  $host_matrix_qc_input"
    echo "featureCounts summary:"
    echo "  $host_counts_summary"
    echo "Output:"
    echo "  $host_matrix_qc_out"
    echo "Script:"
    echo "  $host_matrix_qc_script"
    echo "============================================================"

    if [[ ! -s "$host_matrix_qc_script" ]]; then
        echo "${y}[WARNING] Host matrix QC script not found. Skipping host expression QC.${w}"
        echo "Expected:"
        echo "  $host_matrix_qc_script"
        return 0
    fi

    if [[ ! -s "$host_matrix_qc_input" ]]; then
        echo "${y}[WARNING] Host count matrix not found. Skipping host expression QC.${w}"
        echo "Expected:"
        echo "  $host_matrix_qc_input"
        return 0
    fi

    # ------------------------------------------------------------
    # Check whether featureCounts assigned reads to host features
    # ------------------------------------------------------------

    if [[ -s "$host_counts_summary" ]]; then
        local assigned_total

        assigned_total=$(awk -F'\t' '
            $1 == "Assigned" {
                total = 0
                for (i = 2; i <= NF; i++) {
                    if ($i ~ /^[0-9]+$/) {
                        total += $i
                    }
                }
                print total
                found = 1
            }
            END {
                if (found != 1) print "NA"
            }
        ' "$host_counts_summary")

        echo "[INFO] Total Assigned reads by featureCounts: ${assigned_total}"

        if [[ "$assigned_total" == "0" ]]; then
            echo "${y}[WARNING] featureCounts assigned zero reads to host genes.${w}"
            echo "${y}[WARNING] Skipping host expression matrix QC because PCA/heatmap require non-zero gene counts.${w}"
            echo "${y}[WARNING] Likely causes:${w}"
            echo "  1. FASTA/index contig names do not match GTF contig names"
            echo "  2. Host alignment produced very few usable alignments"
            echo "  3. featureCounts parameters do not match the annotation"
            echo
            echo "[INFO] The pipeline will continue normally."
            return 0
        fi
    else
        echo "${y}[WARNING] featureCounts summary not found. Running host QC anyway.${w}"
        echo "Expected:"
        echo "  $host_counts_summary"
    fi

    # ------------------------------------------------------------
    # Run host expression QC
    # ------------------------------------------------------------

    conda deactivate
    conda activate R412

    if ! Rscript "$host_matrix_qc_script" \
        --input "$host_matrix_qc_input" \
        --samplesheet "$samplesheet_file" \
        --output "$host_matrix_qc_out" \
        --label host_expression \
        --matrix_type featurecounts \
        --normalization logcpm \
        --top_variable 50 \
        --pca_top_features 5000 \
        > "$host_matrix_qc_log" 2>&1
    then
        echo "${y}[WARNING] Host expression matrix QC failed. Continuing pipeline.${w}"
        echo "Log:"
        echo "  $host_matrix_qc_log"
        tail -n 40 "$host_matrix_qc_log" || true
    else
        echo "${g}[OK] Host expression matrix QC generated:${w}"
        echo "  $host_matrix_qc_out"
        echo
        echo "[INFO] Expected outputs include:"
        echo "  host_expression_pca.png"
        echo "  host_expression_top_variable_heatmap.png"
        echo "  host_expression_sample_correlation_heatmap.png"
        echo "  host_expression_detected_features_per_sample.png"
    fi

    conda deactivate
    conda activate MTD

    echo "============================================================"
    echo "${g}[OK] Exploratory host expression QC step finished${w}"
    echo "============================================================"
}

cd "$outputdr/temp" || die "Could not enter temp directory: $outputdr/temp"
PIPELINE_TEMP_DIR="$(pwd)"
# ------------------------------------------------------------
# Resolve FASTQ input for one sample
# ------------------------------------------------------------
# Output format:
#   se<TAB>READ1<TAB>-
#   pe<TAB>READ1<TAB>READ2
#
# Return codes:
#   0 = FASTQ successfully resolved
#   1 = no supported FASTQ found
#   2 = ambiguous, incomplete, or conflicting input
# ------------------------------------------------------------
resolve_fastq_for_sample() {
    local sample="$1"
    local search_dir="$2"
    local compressed_only="${3:-0}"

    local path
    local name

    local -a se_candidates=()
    local -a r1_candidates=()
    local -a r2_candidates=()

    while IFS= read -r -d '' path; do
        name="$(basename "$path")"

        # Accept FASTQ/FQ, optionally gzip-compressed.
        # In custom-cache mode, accept compressed files only.
        case "$name" in
            *.fastq.gz|*.fq.gz)
                ;;

            *.fastq|*.fq)
                if [[ "$compressed_only" == "1" ]]; then
                    continue
                fi
                ;;

            *)
                continue
                ;;
        esac

        case "$name" in

            # ------------------------------------------------
            # Single-end names
            # ------------------------------------------------
            "${sample}.fastq"|\
            "${sample}.fq"|\
            "${sample}.fastq.gz"|\
            "${sample}.fq.gz"|\
            "Trimmed_${sample}.fastq"|\
            "Trimmed_${sample}.fq"|\
            "Trimmed_${sample}.fastq.gz"|\
            "Trimmed_${sample}.fq.gz")
                se_candidates+=("$path")
                ;;

            # ------------------------------------------------
            # Paired-end R1 names
            #
            # Supported examples:
            # sample_R1.fastq.gz
            # sample_R1_001.fastq.gz
            # sample_1.fastq.gz
            # sample_1_001.fastq.gz
            # sample_S1_L001_R1_001.fastq.gz
            # ------------------------------------------------
            "${sample}_R1.fastq"|\
            "${sample}_R1.fq"|\
            "${sample}_R1.fastq.gz"|\
            "${sample}_R1.fq.gz"|\
            "${sample}_R1_"*.fastq|\
            "${sample}_R1_"*.fq|\
            "${sample}_R1_"*.fastq.gz|\
            "${sample}_R1_"*.fq.gz|\
            "${sample}_1.fastq"|\
            "${sample}_1.fq"|\
            "${sample}_1.fastq.gz"|\
            "${sample}_1.fq.gz"|\
            "${sample}_1_"*.fastq|\
            "${sample}_1_"*.fq|\
            "${sample}_1_"*.fastq.gz|\
            "${sample}_1_"*.fq.gz|\
            "${sample}_"*"_R1.fastq"|\
            "${sample}_"*"_R1.fq"|\
            "${sample}_"*"_R1.fastq.gz"|\
            "${sample}_"*"_R1.fq.gz"|\
            "${sample}_"*"_R1_"*.fastq|\
            "${sample}_"*"_R1_"*.fq|\
            "${sample}_"*"_R1_"*.fastq.gz|\
            "${sample}_"*"_R1_"*.fq.gz|\
            "Trimmed_${sample}_R1.fastq"|\
            "Trimmed_${sample}_R1.fq"|\
            "Trimmed_${sample}_R1.fastq.gz"|\
            "Trimmed_${sample}_R1.fq.gz"|\
            "Trimmed_${sample}_R1_"*.fastq|\
            "Trimmed_${sample}_R1_"*.fq|\
            "Trimmed_${sample}_R1_"*.fastq.gz|\
            "Trimmed_${sample}_R1_"*.fq.gz|\
            "Trimmed_${sample}_1.fastq"|\
            "Trimmed_${sample}_1.fq"|\
            "Trimmed_${sample}_1.fastq.gz"|\
            "Trimmed_${sample}_1.fq.gz")
                r1_candidates+=("$path")
                ;;

            # ------------------------------------------------
            # Paired-end R2 names
            # ------------------------------------------------
            "${sample}_R2.fastq"|\
            "${sample}_R2.fq"|\
            "${sample}_R2.fastq.gz"|\
            "${sample}_R2.fq.gz"|\
            "${sample}_R2_"*.fastq|\
            "${sample}_R2_"*.fq|\
            "${sample}_R2_"*.fastq.gz|\
            "${sample}_R2_"*.fq.gz|\
            "${sample}_2.fastq"|\
            "${sample}_2.fq"|\
            "${sample}_2.fastq.gz"|\
            "${sample}_2.fq.gz"|\
            "${sample}_2_"*.fastq|\
            "${sample}_2_"*.fq|\
            "${sample}_2_"*.fastq.gz|\
            "${sample}_2_"*.fq.gz|\
            "${sample}_"*"_R2.fastq"|\
            "${sample}_"*"_R2.fq"|\
            "${sample}_"*"_R2.fastq.gz"|\
            "${sample}_"*"_R2.fq.gz"|\
            "${sample}_"*"_R2_"*.fastq|\
            "${sample}_"*"_R2_"*.fq|\
            "${sample}_"*"_R2_"*.fastq.gz|\
            "${sample}_"*"_R2_"*.fq.gz|\
            "Trimmed_${sample}_R2.fastq"|\
            "Trimmed_${sample}_R2.fq"|\
            "Trimmed_${sample}_R2.fastq.gz"|\
            "Trimmed_${sample}_R2.fq.gz"|\
            "Trimmed_${sample}_R2_"*.fastq|\
            "Trimmed_${sample}_R2_"*.fq|\
            "Trimmed_${sample}_R2_"*.fastq.gz|\
            "Trimmed_${sample}_R2_"*.fq.gz|\
            "Trimmed_${sample}_2.fastq"|\
            "Trimmed_${sample}_2.fq"|\
            "Trimmed_${sample}_2.fastq.gz"|\
            "Trimmed_${sample}_2.fq.gz")
                r2_candidates+=("$path")
                ;;
        esac

    done < <(
        find "$search_dir" \
            -maxdepth 1 \
            -type f \
            -print0 |
        sort -z
    )

    # Any R1 or R2 candidate means that paired-end input is expected.
    if (( ${#r1_candidates[@]} > 0 || ${#r2_candidates[@]} > 0 )); then

        if (( ${#se_candidates[@]} > 0 )); then
            echo "${r}[ERROR] Both single-end and paired-end files were found for sample:${w} $sample" >&2
            echo "Single-end candidates:" >&2
            printf '  %s\n' "${se_candidates[@]}" >&2
            echo "R1 candidates:" >&2
            printf '  %s\n' "${r1_candidates[@]}" >&2
            echo "R2 candidates:" >&2
            printf '  %s\n' "${r2_candidates[@]}" >&2
            return 2
        fi

        if (( ${#r1_candidates[@]} != 1 || ${#r2_candidates[@]} != 1 )); then
            echo "${r}[ERROR] Could not resolve exactly one R1/R2 pair for sample:${w} $sample" >&2

            echo "R1 candidates found: ${#r1_candidates[@]}" >&2
            if (( ${#r1_candidates[@]} > 0 )); then
                printf '  %s\n' "${r1_candidates[@]}" >&2
            fi

            echo "R2 candidates found: ${#r2_candidates[@]}" >&2
            if (( ${#r2_candidates[@]} > 0 )); then
                printf '  %s\n' "${r2_candidates[@]}" >&2
            fi

            return 2
        fi

        printf 'pe\t%s\t%s\n' \
            "${r1_candidates[0]}" \
            "${r2_candidates[0]}"

        return 0
    fi

    # No paired-end candidates: check for a unique SE file.
    if (( ${#se_candidates[@]} == 1 )); then
        printf 'se\t%s\t-\n' "${se_candidates[0]}"
        return 0
    fi

    if (( ${#se_candidates[@]} > 1 )); then
        echo "${r}[ERROR] More than one single-end FASTQ was found for sample:${w} $sample" >&2
        printf '  %s\n' "${se_candidates[@]}" >&2
        return 2
    fi

    return 1
}


# ------------------------------------------------------------
# Resolve FASTQ from the appropriate input source
# ------------------------------------------------------------
# Search order:
#   1. --custom-raw-path, when provided
#   2. samplesheet directory
#   3. pipeline temp directory, for files downloaded from SRA
# ------------------------------------------------------------
resolve_fastq_input() {
    local sample="$1"
    local resolved=""
    local status=0

    if [[ -n "${CUSTOM_PATH:-}" ]]; then
        resolve_fastq_for_sample "$sample" "$CUSTOM_PATH" 1
        return $?
    fi

    resolved="$(resolve_fastq_for_sample "$sample" "$inputdr" 0)"
    status=$?

    if [[ "$status" -eq 0 ]]; then
        printf '%s\n' "$resolved"
        return 0
    fi

    # An ambiguous result must not be hidden by searching elsewhere.
    if [[ "$status" -eq 2 ]]; then
        return 2
    fi

    # SRA Toolkit files are generated inside the temporary directory.
    resolve_fastq_for_sample "$sample" "$PIPELINE_TEMP_DIR" 0
}

# ------------------------------------------------------------
# Download and convert one SRA run when no local FASTQ exists
# ------------------------------------------------------------
# Supported run accession prefixes:
#   SRR = NCBI SRA
#   ERR = European Nucleotide Archive
#   DRR = DDBJ Sequence Read Archive
#
# fasterq-dump is intentionally used in its default split-3
# mode:
#
#   SE:
#     ACCESSION.fastq
#
#   PE:
#     ACCESSION_1.fastq
#     ACCESSION_2.fastq
#
# Mixed paired and unpaired output from one accession is not
# accepted automatically because the MTD run requires one
# consistent library layout per sample.
# ------------------------------------------------------------
download_sra_fastq_if_needed() {
    local accession="$1"

    local resolved=""
    local status=0
    local detected_layout=""
    local detected_read1=""
    local detected_read2=""

    local sra_download_dir="$inputdr/$accession"
    local sra_scratch_dir="$PIPELINE_TEMP_DIR/sra_fasterq_tmp/$accession"

    echo "============================================================"
    echo "[SRA] Checking accession: $accession"
    echo "FASTQ output directory:"
    echo "  $inputdr"
    echo "============================================================"

    # First check whether a complete SE file or PE pair already
    # exists locally.
    resolved="$(resolve_fastq_for_sample "$accession" "$inputdr" 0)"
    status=$?

    if [[ "$status" -eq 0 ]]; then
        IFS=$'\t' read -r \
            detected_layout \
            detected_read1 \
            detected_read2 \
            <<< "$resolved"

        echo "${g}[OK] Local FASTQ already exists for SRA accession:${w} $accession"
        echo "Layout: $detected_layout"
        echo "Read 1: $detected_read1"

        if [[ "$detected_layout" == "pe" ]]; then
            echo "Read 2: $detected_read2"
        fi

        return 0
    fi

    # Status 2 means that local files exist but are ambiguous,
    # incomplete, or contain conflicting SE/PE inputs.
    if [[ "$status" -eq 2 ]]; then
        die "Ambiguous or incomplete local FASTQ input was found for SRA accession: $accession"
    fi

    echo "${y}[INFO] No complete local FASTQ found for accession:${w} $accession"
    echo "[INFO] Starting SRA download and FASTQ conversion."

    mkdir -p "$sra_scratch_dir"

    echo
    echo "============================================================"
    echo "[SRA PREFETCH] Accession: $accession"
    echo "Working directory: $inputdr"
    echo "Maximum accession size: 999G"
    echo "============================================================"

    if ! (
        cd "$inputdr" || exit 1

        prefetch \
            -X 999G \
            "$accession"
    ); then
        die "SRA prefetch failed for accession: $accession"
    fi

    echo
    echo "============================================================"
    echo "[FASTERQ-DUMP] Accession: $accession"
    echo "Output directory: $inputdr"
    echo "Scratch directory: $sra_scratch_dir"
    echo "Threads: $threads"
    echo "Output mode: default split-3"
    echo "============================================================"

    if ! fasterq-dump \
        "$sra_download_dir" \
        --outdir "$inputdr" \
        --temp "$sra_scratch_dir" \
        --threads "$threads" \
        --progress
    then
        die "fasterq-dump failed for accession: $accession"
    fi

    # Validate the files using the same resolver used by the
    # rest of the pipeline.
    resolved="$(resolve_fastq_for_sample "$accession" "$inputdr" 0)"
    status=$?

    if [[ "$status" -ne 0 ]]; then
        echo "${r}[ERROR] SRA conversion did not produce one valid SE input or one valid PE pair.${w}"
        echo "Accession:"
        echo "  $accession"
        echo "FASTQ directory:"
        echo "  $inputdr"
        echo
        echo "Files generated for this accession:"

        find "$inputdr" \
            -maxdepth 1 \
            -type f \
            \( \
                -name "${accession}.fastq" -o \
                -name "${accession}.fq" -o \
                -name "${accession}_1.fastq" -o \
                -name "${accession}_2.fastq" -o \
                -name "${accession}_1.fq" -o \
                -name "${accession}_2.fq" \
            \) \
            -printf '  %f\n' \
            2>/dev/null || true

        die "SRA FASTQ validation failed for accession: $accession"
    fi

    IFS=$'\t' read -r \
        detected_layout \
        detected_read1 \
        detected_read2 \
        <<< "$resolved"

    echo
    echo "${g}============================================================"
    echo "[OK] SRA FASTQ conversion completed"
    echo "Accession:${w} $accession"
    echo "${g}Layout:${w} $detected_layout"
    echo "${g}Read 1:${w} $detected_read1"

    if [[ "$detected_layout" == "pe" ]]; then
        echo "${g}Read 2:${w} $detected_read2"
    fi

    echo "${g}============================================================${w}"

    # The prefetched accession directory and fasterq scratch
    # files are no longer required after successful conversion.
    rm -rf -- \
        "$sra_download_dir" \
        "$sra_scratch_dir"
}

# ------------------------------------------------------------
# Load one sample from the validated FASTQ manifest
# ------------------------------------------------------------
# Global variables populated:
#   INPUT_LAYOUT
#   INPUT_READ1
#   INPUT_READ2
# ------------------------------------------------------------
load_fastq_manifest_row() {
    local sample="$1"
    local row=""
    local match_count=0

    match_count="$(
        awk -F $'\t' -v sample="$sample" '
            NR > 1 && $1 == sample {
                count++
            }
            END {
                print count + 0
            }
        ' "$FASTQ_INPUT_MANIFEST"
    )"

    if [[ "$match_count" -ne 1 ]]; then
        die "Expected exactly one FASTQ manifest entry for sample '$sample', but found $match_count."
    fi

    row="$(
        awk -F $'\t' -v sample="$sample" '
            NR > 1 && $1 == sample {
                print
                exit
            }
        ' "$FASTQ_INPUT_MANIFEST"
    )"

    IFS=$'\t' read -r \
        MANIFEST_SAMPLE \
        INPUT_LAYOUT \
        INPUT_READ1 \
        INPUT_READ2 \
        <<< "$row"

    if [[ "$MANIFEST_SAMPLE" != "$sample" ]]; then
        die "FASTQ manifest sample mismatch. Expected '$sample', got '$MANIFEST_SAMPLE'."
    fi

    if [[ "$INPUT_LAYOUT" != "se" && "$INPUT_LAYOUT" != "pe" ]]; then
        die "Invalid FASTQ layout in manifest for sample '$sample': $INPUT_LAYOUT"
    fi

    require_file "$INPUT_READ1" "FASTQ read 1 for sample $sample"

    if [[ "$INPUT_LAYOUT" == "pe" ]]; then
        if [[ -z "$INPUT_READ2" || "$INPUT_READ2" == "-" ]]; then
            die "Paired-end sample '$sample' has no read 2 in the FASTQ manifest."
        fi

        require_file "$INPUT_READ2" "FASTQ read 2 for sample $sample"
    else
        INPUT_READ2="-"
    fi

    if [[ "$INPUT_LAYOUT" != "$READ_LAYOUT" ]]; then
        die "FASTQ manifest layout mismatch for sample '$sample'. Run layout: $READ_LAYOUT; sample layout: $INPUT_LAYOUT"
    fi
}
# ------------------------------------------------------------
# Extract sample names from samplesheet.csv
# ------------------------------------------------------------

lsn=$(awk -F',' '
NR > 1 {
    s=$1
    gsub(/^[ \t\r\n]+|[ \t\r\n]+$/, "", s)
    if (s != "") print s
}
' "$samplesheet_file" | sort -u | tr "\n" " " | sed 's/[[:space:]]*$//')

if [[ -z "$lsn" ]]; then
    die "No sample names were found in samplesheet: $samplesheet_file"
fi

echo "[INFO] Samples detected from samplesheet:"
echo "  $lsn"

# ------------------------------------------------------------
# Detect and validate FASTQ layout for all samples
# ------------------------------------------------------------
# ------------------------------------------------------------
# Download SRA run accessions that do not already have FASTQ
# files in the input directory
# ------------------------------------------------------------

if [[ -z "${CUSTOM_PATH:-}" ]]; then
    for i in $lsn; do
        case "$i" in
            SRR[0-9]*|ERR[0-9]*|DRR[0-9]*)
                download_sra_fastq_if_needed "$i"
                ;;

            *)
                # Ordinary local sample name; no SRA download.
                ;;
        esac
    done
else
    echo "${y}[INFO] --custom-raw-path is active. Automatic SRA download is disabled.${w}"
fi

FASTQ_INPUT_MANIFEST="$PIPELINE_TEMP_DIR/fastq_input_manifest.tsv"

printf 'sample\tlayout\tread1\tread2\n' > "$FASTQ_INPUT_MANIFEST"

detected_run_layout=""
fastq_errors=0

echo "${g}============================================"
echo "FASTQ input detection"
echo "Requested read layout:${w} $READ_LAYOUT_MODE"
echo "${g}============================================${w}"

for i in $lsn; do
    resolved="$(resolve_fastq_input "$i")"
    resolve_status=$?

    if [[ "$resolve_status" -ne 0 ]]; then
        if [[ "$resolve_status" -eq 1 ]]; then
            echo "${r}[ERROR] No supported FASTQ input was found for sample:${w} $i" >&2

            if [[ -n "${CUSTOM_PATH:-}" ]]; then
                echo "Search directory:" >&2
                echo "  $CUSTOM_PATH" >&2
            else
                echo "Search directories:" >&2
                echo "  $inputdr" >&2
                echo "  $PIPELINE_TEMP_DIR" >&2
            fi
        fi

        fastq_errors=1
        continue
    fi

    IFS=$'\t' read -r sample_layout read1 read2 <<< "$resolved"

    if [[ -z "$sample_layout" || -z "$read1" ]]; then
        echo "${r}[ERROR] Invalid FASTQ resolver output for sample:${w} $i" >&2
        fastq_errors=1
        continue
    fi

    # All samples in one run must use the same layout.
    if [[ -z "$detected_run_layout" ]]; then
        detected_run_layout="$sample_layout"
    elif [[ "$sample_layout" != "$detected_run_layout" ]]; then
        echo "${r}[ERROR] Mixed SE and PE samples were detected in the same run.${w}" >&2
        echo "Sample: $i" >&2
        echo "Sample layout: $sample_layout" >&2
        echo "Previously detected run layout: $detected_run_layout" >&2
        fastq_errors=1
    fi

    # Validate an explicitly requested layout.
    if [[ "$READ_LAYOUT_MODE" != "auto" && "$sample_layout" != "$READ_LAYOUT_MODE" ]]; then
        echo "${r}[ERROR] FASTQ layout does not match --read-layout for sample:${w} $i" >&2
        echo "Requested: $READ_LAYOUT_MODE" >&2
        echo "Detected:  $sample_layout" >&2
        fastq_errors=1
    fi

    printf '%s\t%s\t%s\t%s\n' \
        "$i" \
        "$sample_layout" \
        "$read1" \
        "$read2" \
        >> "$FASTQ_INPUT_MANIFEST"

    echo
    echo "[FASTQ] Sample: $i"
    echo "Layout: $sample_layout"
    echo "Read 1: $read1"

    if [[ "$sample_layout" == "pe" ]]; then
        echo "Read 2: $read2"
    fi
done

if [[ "$fastq_errors" == "1" ]]; then
    echo
    die "FASTQ input validation failed."
fi

if [[ -z "$detected_run_layout" ]]; then
    die "Could not determine the sequencing read layout."
fi

if [[ "$READ_LAYOUT_MODE" == "auto" ]]; then
    READ_LAYOUT="$detected_run_layout"
else
    READ_LAYOUT="$READ_LAYOUT_MODE"
fi

echo
echo "${g}============================================"
echo "[OK] FASTQ input validation completed"
echo "Effective read layout:${w} $READ_LAYOUT"
echo "FASTQ input manifest:"
echo "  $FASTQ_INPUT_MANIFEST"
echo "${g}============================================${w}"

column -s $'\t' -t "$FASTQ_INPUT_MANIFEST" 2>/dev/null || \
    cat "$FASTQ_INPUT_MANIFEST"

# Rewrite the methods log now that the effective read layout
# has been detected and validated.
write_methods_log

# ------------------------------------------------------------
# Study design summary
# ------------------------------------------------------------

echo "${g}============================================"
echo "Main study design:${w}"

awk -F',' '
NR > 1 {
    groups[$2]++
}
END {
    for (g in groups) {
        printf "Group: %s - Number of samples: %d\n", g, groups[g]
    }
}' "$samplesheet_file"

echo "${g}============================================${w}"

# ------------------------------------------------------------
# Detect whether experimental comparison is possible
# ------------------------------------------------------------

GROUP_COUNT=$(awk -F',' '
NR > 1 {
    g=$2
    gsub(/^[ \t\r\n]+|[ \t\r\n]+$/, "", g)
    if (g != "") groups[g]=1
}
END {
    print length(groups)
}' "$samplesheet_file")

if [[ -z "$GROUP_COUNT" ]]; then
    GROUP_COUNT=0
fi

if [[ "$analysis_mode" == "auto" ]]; then
    if [[ "$GROUP_COUNT" -lt 2 ]]; then
        NO_COMPARISON=1
        analysis_mode_resolved="exploratory"
    else
        NO_COMPARISON=0
        analysis_mode_resolved="comparison"
    fi
elif [[ "$analysis_mode" == "exploratory" ]]; then
    NO_COMPARISON=1
    analysis_mode_resolved="exploratory"
else
    NO_COMPARISON=0
    analysis_mode_resolved="comparison"

    if [[ "$GROUP_COUNT" -lt 2 ]]; then
        echo "${r}[ERROR] --analysis-mode comparison was requested, but samplesheet has fewer than 2 groups.${w}"
        echo "Detected groups: $GROUP_COUNT"
        echo "Use --analysis-mode exploratory for single-group projects."
        exit 1
    fi
fi

echo "${g}============================================"
echo "Analysis mode:${w} $analysis_mode_resolved"
echo "Groups detected in samplesheet: $GROUP_COUNT"

if [[ "$NO_COMPARISON" == "1" ]]; then
    echo "${y}[INFO] No experimental comparison will be performed.${w}"
    echo "${y}[INFO] DEG/DESeq2-dependent steps will be skipped or treated as optional.${w}"
fi

echo "${g}============================================${w}"
# ------------------------------------------------------------
# Optional metadata summary
# ------------------------------------------------------------

if [[ -n "$metadata" ]]; then
    if [[ ! -s "$metadata" ]]; then
        die "Metadata file not found or empty: $metadata"
    fi

    echo "============================================"

    header=$(head -n 1 "$metadata")
    IFS=',' read -ra columns <<< "$header"

    for ((i=3; i<=${#columns[@]}; i++)); do
        col="${columns[$i-1]}"

        echo "Metadata column: $col,"
        echo "Meta-groups:"

        awk -v col_index="$i" -F',' '
        NR > 1 {
            values[$col_index]++
        }
        END {
            for (value in values) {
                printf "  %s: %d\n", value, values[value]
            }
        }' "$metadata"
    done

    echo "============================================"
    echo ""
fi

echo "${g}MTD running  progress:"
echo ">>                  [10%]"

echo "Raw reads preparation${w}"

# ------------------------------------------------------------
# Prepare FASTQ files for the pipeline
# ------------------------------------------------------------
# The host Kraken2 step consumes normalized FASTQ files recorded
# in PREPARED_FASTQ_MANIFEST:
#
#   SE:
#     Trimmed_${sample}.fq.gz
#
#   PE:
#     Trimmed_${sample}_R1.fq.gz
#     Trimmed_${sample}_R2.fq.gz
#
# Modes:
#   1) Normal mode:
#        no --no-trim and no --custom-raw-path
#        -> run fastp using FASTQ files from samplesheet directory
#
#   2) No-trim mode without custom cache:
#        --no-trim
#        -> copy .gz files or compress uncompressed FASTQ files from samplesheet directory
#
#   3) No-trim mode with custom cache:
#        --custom-raw-path DIR
#        -> copy already compressed FASTQ files from DIR
#        -> no fastp, no pigz compression
# ------------------------------------------------------------
total_cores=$(nproc)

if [[ "$total_cores" -le 4 ]]; then
    threads_per_job=1
elif [[ "$total_cores" -le 8 ]]; then
    threads_per_job=2
elif [[ "$total_cores" -le 16 ]]; then
    threads_per_job=4
else
    threads_per_job=10
fi

max_jobs=$(( total_cores / threads_per_job ))

if [[ "$max_jobs" -lt 1 ]]; then
    max_jobs=1
fi

# ------------------------------------------------------------
# Prepared FASTQ manifest
# ------------------------------------------------------------
# This manifest records the normalized FASTQ filenames created
# by the read-preparation step. Downstream tools must use these
# paths instead of searching for FASTQ files again.
# ------------------------------------------------------------
PREPARED_FASTQ_MANIFEST="$PIPELINE_TEMP_DIR/prepared_fastq_manifest.tsv"

printf 'sample\tlayout\tread1\tread2\n' > "$PREPARED_FASTQ_MANIFEST"


# ------------------------------------------------------------
# Copy or compress a FASTQ without trimming
# ------------------------------------------------------------
prepare_fastq_without_trimming() {
    local input_fq="$1"
    local output_fq="$2"
    local label="$3"

    require_file "$input_fq" "$label input FASTQ"

    if [[ "$input_fq" == *.gz ]]; then
        cp -- "$input_fq" "$output_fq" || \
            die "Failed to copy $label FASTQ: $input_fq"
    else
        pigz \
            -p "$threads_per_job" \
            -c "$input_fq" \
            > "$output_fq" || \
            die "Failed to compress $label FASTQ: $input_fq"
    fi

    require_file "$output_fq" "$label prepared FASTQ"
}


# ------------------------------------------------------------
# Register prepared FASTQ files
# ------------------------------------------------------------
register_prepared_fastq() {
    local sample="$1"
    local layout="$2"
    local read1="$3"
    local read2="${4:--}"

    require_file "$read1" "Prepared read 1 for sample $sample"

    if [[ "$layout" == "pe" ]]; then
        if [[ -z "$read2" || "$read2" == "-" ]]; then
            die "Prepared paired-end sample '$sample' has no read 2."
        fi

        require_file "$read2" "Prepared read 2 for sample $sample"

        validate_fastq_pair \
            "$read1" \
            "$read2" \
            "prepared FASTQ pair for sample $sample"
    else
        read2="-"
    fi

    printf '%s\t%s\t%s\t%s\n' \
        "$sample" \
        "$layout" \
        "$read1" \
        "$read2" \
        >> "$PREPARED_FASTQ_MANIFEST"
}

# ------------------------------------------------------------
# Load one sample from the prepared FASTQ manifest
# ------------------------------------------------------------
# Global variables populated:
#   PREPARED_LAYOUT
#   PREPARED_READ1
#   PREPARED_READ2
# ------------------------------------------------------------
load_prepared_fastq_manifest_row() {
    local sample="$1"
    local row=""
    local match_count=0

    match_count="$(
        awk -F $'\t' -v sample="$sample" '
            NR > 1 && $1 == sample {
                count++
            }
            END {
                print count + 0
            }
        ' "$PREPARED_FASTQ_MANIFEST"
    )"

    if [[ "$match_count" -ne 1 ]]; then
        die "Expected exactly one prepared FASTQ manifest entry for sample '$sample', but found $match_count."
    fi

    row="$(
        awk -F $'\t' -v sample="$sample" '
            NR > 1 && $1 == sample {
                print
                exit
            }
        ' "$PREPARED_FASTQ_MANIFEST"
    )"

    IFS=$'\t' read -r \
        PREPARED_SAMPLE \
        PREPARED_LAYOUT \
        PREPARED_READ1 \
        PREPARED_READ2 \
        <<< "$row"

    if [[ "$PREPARED_SAMPLE" != "$sample" ]]; then
        die "Prepared FASTQ manifest sample mismatch. Expected '$sample', got '$PREPARED_SAMPLE'."
    fi

    if [[ "$PREPARED_LAYOUT" != "se" && "$PREPARED_LAYOUT" != "pe" ]]; then
        die "Invalid prepared FASTQ layout for sample '$sample': $PREPARED_LAYOUT"
    fi

    require_file \
        "$PREPARED_READ1" \
        "Prepared read 1 for sample $sample"

    if [[ "$PREPARED_LAYOUT" == "pe" ]]; then
        if [[ -z "$PREPARED_READ2" || "$PREPARED_READ2" == "-" ]]; then
            die "Prepared paired-end sample '$sample' has no read 2."
        fi

        require_file \
            "$PREPARED_READ2" \
            "Prepared read 2 for sample $sample"
    else
        PREPARED_READ2="-"
    fi

    if [[ "$PREPARED_LAYOUT" != "$READ_LAYOUT" ]]; then
        die "Prepared FASTQ layout mismatch for sample '$sample'. Run layout: $READ_LAYOUT; sample layout: $PREPARED_LAYOUT"
    fi
}

# ------------------------------------------------------------
# Mode 1: custom cache path provided
# ------------------------------------------------------------

if [[ -n "${CUSTOM_PATH:-}" ]]; then
    echo "${y}[INFO] Using --custom-raw-path cache mode.${w}"
    echo "This mode expects already compressed FASTQ files."
    echo "No fastp and no pigz compression will be performed."
    echo "Custom FASTQ cache:"
    echo "  $CUSTOM_PATH"
    echo

        for i in $lsn; do
        load_fastq_manifest_row "$i"

        echo "============================================================"
        echo "[COPY CACHED FASTQ] Sample: $i"
        echo "Layout: $INPUT_LAYOUT"

        if [[ "$INPUT_LAYOUT" == "pe" ]]; then
            out_fq1="$PIPELINE_TEMP_DIR/Trimmed_${i}_R1.fq.gz"
            out_fq2="$PIPELINE_TEMP_DIR/Trimmed_${i}_R2.fq.gz"

            echo "Input R1:  $INPUT_READ1"
            echo "Input R2:  $INPUT_READ2"
            echo "Output R1: $out_fq1"
            echo "Output R2: $out_fq2"
            echo "============================================================"

            prepare_fastq_without_trimming \
                "$INPUT_READ1" \
                "$out_fq1" \
                "$i R1"

            prepare_fastq_without_trimming \
                "$INPUT_READ2" \
                "$out_fq2" \
                "$i R2"

            register_prepared_fastq \
                "$i" \
                "pe" \
                "$out_fq1" \
                "$out_fq2"
        else
            out_fq="$PIPELINE_TEMP_DIR/Trimmed_${i}.fq.gz"

            echo "Input:  $INPUT_READ1"
            echo "Output: $out_fq"
            echo "============================================================"

            prepare_fastq_without_trimming \
                "$INPUT_READ1" \
                "$out_fq" \
                "$i SE"

            register_prepared_fastq \
                "$i" \
                "se" \
                "$out_fq" \
                "-"
        fi
    done

# ------------------------------------------------------------
# Mode 2: no-trim without custom cache
# ------------------------------------------------------------

elif [[ "$no_trimm" == "1" ]]; then
    echo "${y}[INFO] --no-trim was declared.${w}"
    echo "Using FASTQ files from samplesheet directory:"
    echo "  $inputdr"
    echo "Compressed files will be copied; uncompressed files will be compressed with pigz."
    echo

        for i in $lsn; do
        load_fastq_manifest_row "$i"

        echo "============================================================"
        echo "[NO TRIM] Sample: $i"
        echo "Layout: $INPUT_LAYOUT"

        if [[ "$INPUT_LAYOUT" == "pe" ]]; then
            out_fq1="$PIPELINE_TEMP_DIR/Trimmed_${i}_R1.fq.gz"
            out_fq2="$PIPELINE_TEMP_DIR/Trimmed_${i}_R2.fq.gz"

            echo "Input R1:  $INPUT_READ1"
            echo "Input R2:  $INPUT_READ2"
            echo "Output R1: $out_fq1"
            echo "Output R2: $out_fq2"
            echo "============================================================"

            prepare_fastq_without_trimming \
                "$INPUT_READ1" \
                "$out_fq1" \
                "$i R1"

            prepare_fastq_without_trimming \
                "$INPUT_READ2" \
                "$out_fq2" \
                "$i R2"

            register_prepared_fastq \
                "$i" \
                "pe" \
                "$out_fq1" \
                "$out_fq2"
        else
            out_fq="$PIPELINE_TEMP_DIR/Trimmed_${i}.fq.gz"

            echo "Input:  $INPUT_READ1"
            echo "Output: $out_fq"
            echo "============================================================"

            prepare_fastq_without_trimming \
                "$INPUT_READ1" \
                "$out_fq" \
                "$i SE"

            register_prepared_fastq \
                "$i" \
                "se" \
                "$out_fq" \
                "-"
        fi
    done

# ------------------------------------------------------------
# Mode 3: normal fastp trimming
# ------------------------------------------------------------

else
    echo "[INFO] Running fastp trimming using FASTQ files from samplesheet directory:"
    echo "  $inputdr"
    echo

    conda deactivate
    conda activate MTD_fastp

    mkdir -p "$outputdr/fastp"

    for i in $lsn; do
        load_fastq_manifest_row "$i"

        fastp_report_base="Trimmed_${i}"
        fastp_html="$outputdr/fastp/${fastp_report_base}.fastp.html"
        fastp_json="$outputdr/fastp/${fastp_report_base}.fastp.json"

        echo "============================================================"
        echo "[FASTP] Sample: $i"
        echo "Layout: $INPUT_LAYOUT"
        echo "Threads: $threads"
        echo "Minimum length: $length"

        # Common fastp arguments for both SE and PE.
        fastp_common_args=(
            --trim_poly_x
            --qualified_quality_phred 15
            --unqualified_percent_limit 40
            --n_base_limit 5
            --cut_front
            --cut_front_window_size 1
            --cut_front_mean_quality 5
            --cut_tail
            --cut_tail_window_size 1
            --cut_tail_mean_quality 5
            --length_required "$length"
            --thread "$threads"
            --html "$fastp_html"
            --json "$fastp_json"
        )

        if [[ "$INPUT_LAYOUT" == "pe" ]]; then
            out_fq1="$PIPELINE_TEMP_DIR/Trimmed_${i}_R1.fq.gz"
            out_fq2="$PIPELINE_TEMP_DIR/Trimmed_${i}_R2.fq.gz"

            echo "Input R1:  $INPUT_READ1"
            echo "Input R2:  $INPUT_READ2"
            echo "Output R1: $out_fq1"
            echo "Output R2: $out_fq2"
            echo "HTML:      $fastp_html"
            echo "JSON:      $fastp_json"
            echo "============================================================"

            if ! fastp \
                "${fastp_common_args[@]}" \
                -i "$INPUT_READ1" \
                -I "$INPUT_READ2" \
                -o "$out_fq1" \
                -O "$out_fq2"
            then
                die "fastp failed for paired-end sample: $i"
            fi

            require_file "$out_fq1" "fastp R1 output for sample $i"
            require_file "$out_fq2" "fastp R2 output for sample $i"

            register_prepared_fastq \
                "$i" \
                "pe" \
                "$out_fq1" \
                "$out_fq2"
        else
            out_fq="$PIPELINE_TEMP_DIR/Trimmed_${i}.fq.gz"

            echo "Input:  $INPUT_READ1"
            echo "Output: $out_fq"
            echo "HTML:   $fastp_html"
            echo "JSON:   $fastp_json"
            echo "============================================================"

            if ! fastp \
                "${fastp_common_args[@]}" \
                -i "$INPUT_READ1" \
                -o "$out_fq"
            then
                die "fastp failed for single-end sample: $i"
            fi

            require_file "$out_fq" "fastp SE output for sample $i"

            register_prepared_fastq \
                "$i" \
                "se" \
                "$out_fq" \
                "-"
        fi

        require_file "$fastp_html" "fastp HTML report for sample $i"
        require_file "$fastp_json" "fastp JSON report for sample $i"
    done
fi
#$MTDIR/MTD_scripts/data_trimming.sh 
conda deactivate
conda activate MTD

echo
echo "${g}============================================"
echo "[OK] FASTQ preparation completed"
echo "Read layout:${w} $READ_LAYOUT"
echo "Prepared FASTQ manifest:"
echo "  $PREPARED_FASTQ_MANIFEST"
echo "${g}============================================${w}"

column -s $'\t' -t "$PREPARED_FASTQ_MANIFEST" 2>/dev/null || \
    cat "$PREPARED_FASTQ_MANIFEST"

echo "${g}MTD running  progress:"
echo '>>>>                [20%]'
echo "Reads classification by kraken2; 1st step for host ${w}"
echo "Host DB: $DB_host"

if [[ ! -d "$DB_host" ]]; then
    echo "[ERROR] Host Kraken2 DB folder not found:"
    echo "$DB_host"
    exit 1
fi

if [[ ! -s "$DB_host/hash.k2d" || ! -s "$DB_host/opts.k2d" || ! -s "$DB_host/taxo.k2d" ]]; then
    echo "[ERROR] Host Kraken2 DB appears incomplete."
    echo "Expected files:"
    echo "  $DB_host/hash.k2d"
    echo "  $DB_host/opts.k2d"
    echo "  $DB_host/taxo.k2d"
    exit 1
fi

summary_file="kraken_host_summary.tsv"

echo -e "sample\thost_classified_reads\thost_classified_pct\thost_unclassified_reads\thost_unclassified_pct\tcount_unit" \
    > "$summary_file"

# Threshold for warning about low host classification.
HOST_LOW_WARN=50

for i in $lsn; do
    load_prepared_fastq_manifest_row "$i"

    report="Report_host_${i}.txt"
    kraken_output="Report_host_${i}.kraken"

    if [[ "$PREPARED_LAYOUT" == "pe" ]]; then
        count_unit="read_pairs"

        host_classified_pattern="${i}_host#.fq"
        host_unclassified_pattern="${i}_non-host_raw#.fq"

        host_read1="${i}_host_1.fq"
        host_read2="${i}_host_2.fq"

        nonhost_read1="${i}_non-host_raw_1.fq"
        nonhost_read2="${i}_non-host_raw_2.fq"
    else
        count_unit="reads"

        host_classified_pattern="${i}_host.fq"
        host_unclassified_pattern="${i}_non-host_raw.fq"

        host_read1="${i}_host.fq"
        host_read2="-"

        nonhost_read1="${i}_non-host_raw.fq"
        nonhost_read2="-"
    fi

    echo "============================================================"
    echo "[HOST] Sample: $i"
    echo "Layout: $PREPARED_LAYOUT"
    echo "Count unit: $count_unit"
    echo "Host Kraken2 DB: $DB_host"
    echo "Kraken2 host confidence: $KRAKEN_HOST_CONF"
    echo "Kraken2 host minimum hit groups: $KRAKEN_HOST_MIN_HIT_GROUPS"

    if [[ "$PREPARED_LAYOUT" == "pe" ]]; then
        echo "Input R1: $PREPARED_READ1"
        echo "Input R2: $PREPARED_READ2"
        echo "Classified output R1: $host_read1"
        echo "Classified output R2: $host_read2"
        echo "Unclassified output R1: $nonhost_read1"
        echo "Unclassified output R2: $nonhost_read2"
    else
        echo "Input: $PREPARED_READ1"
        echo "Classified output: $host_read1"
        echo "Unclassified output: $nonhost_read1"
    fi

    echo "============================================================"

    kraken_host_args=(
        --db "$DB_host"
        --use-names
        --confidence "$KRAKEN_HOST_CONF"
        --minimum-hit-groups "$KRAKEN_HOST_MIN_HIT_GROUPS"
        --report "$report"
        --output "$kraken_output"
        --threads "$threads"
        --gzip-compressed
    )

    if [[ "$PREPARED_LAYOUT" == "pe" ]]; then
        kraken_host_args+=(
            --paired
            --classified-out "$host_classified_pattern"
            --unclassified-out "$host_unclassified_pattern"
            "$PREPARED_READ1"
            "$PREPARED_READ2"
        )
    else
        kraken_host_args+=(
            --classified-out "$host_classified_pattern"
            --unclassified-out "$host_unclassified_pattern"
            "$PREPARED_READ1"
        )
    fi

    if ! kraken2 "${kraken_host_args[@]}"; then
        die "Kraken2 host filtering failed for sample: $i"
    fi

    require_file \
        "$report" \
        "Kraken2 host report for sample $i"

    require_file \
        "$kraken_output" \
        "Kraken2 host classification output for sample $i"

    # Kraken2 may legitimately create an empty classified or
    # unclassified FASTQ when no sequences belong to that category.
    # Therefore, check existence here rather than non-zero size.
    if [[ "$PREPARED_LAYOUT" == "pe" ]]; then
        for output_file in \
            "$host_read1" \
            "$host_read2" \
            "$nonhost_read1" \
            "$nonhost_read2"
        do
            if [[ ! -e "$output_file" ]]; then
                die "Expected paired-end Kraken2 output was not created for sample '$i': $output_file"
            fi
        done
    else
        for output_file in \
            "$host_read1" \
            "$nonhost_read1"
        do
            if [[ ! -e "$output_file" ]]; then
                die "Expected single-end Kraken2 output was not created for sample '$i': $output_file"
            fi
        done
    fi
    if [[ "$PREPARED_LAYOUT" == "pe" ]]; then
        validate_fastq_pair \
            "$nonhost_read1" \
            "$nonhost_read2" \
            "Kraken2 non-host pair for sample $i"

        # Host-classified reads are used by Magic-BLAST.
        if [[ "$blast" == "blast" ]]; then
            validate_fastq_pair \
                "$host_read1" \
                "$host_read2" \
                "Kraken2 host-classified pair for sample $i"
        fi
    fi
    host_unclassified_pct="$(
        awk '$4 == "U" {
            print $1
            exit
        }' "$report"
    )"

    host_unclassified_reads="$(
        awk '$4 == "U" {
            print $2
            exit
        }' "$report"
    )"

    host_classified_pct="$(
        awk '$4 == "R" && $5 == 1 {
            print $1
            exit
        }' "$report"
    )"

    host_classified_reads="$(
        awk '$4 == "R" && $5 == 1 {
            print $2
            exit
        }' "$report"
    )"

    # A Kraken2 report may omit a category when its count is zero.
if [[ -z "$host_unclassified_pct" ]]; then
    host_unclassified_pct="0.00"
fi

if [[ -z "$host_unclassified_reads" ]]; then
    host_unclassified_reads="0"
fi

if [[ -z "$host_classified_reads" ]]; then
    host_classified_reads="0"
fi

    # Percentage fallback in case the root percentage is absent.
    if [[ -z "$host_classified_pct" ]]; then
        host_classified_pct="$(
            awk -v u="$host_unclassified_pct" '
                BEGIN {
                    printf "%.2f", 100 - u
                }
            '
        )"
    fi

    echo
    echo "[RESULT] Sample: $i"
    echo "  Classified as host: ${host_classified_pct}% (${host_classified_reads} ${count_unit})"
    echo "  Unclassified:       ${host_unclassified_pct}% (${host_unclassified_reads} ${count_unit})"

    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$i" \
        "$host_classified_reads" \
        "$host_classified_pct" \
        "$host_unclassified_reads" \
        "$host_unclassified_pct" \
        "$count_unit" \
        >> "$summary_file"

    echo
    echo "[HOST] Main taxa with >=1%:"

    awk '
        $4 != "U" &&
        !($4 == "R" && $5 == 1) &&
        $1 >= 1 {
            name = $6

            for (j = 7; j <= NF; j++) {
                name = name " " $j
            }

            printf "    %7s%%  %12s units  rank=%-4s taxid=%-10s %s\n",
                   $1, $2, $4, $5, name
        }
    ' "$report" | head -n 20

    if awk \
        -v p="$host_classified_pct" \
        -v t="$HOST_LOW_WARN" \
        'BEGIN {
            exit !(p < t)
        }'
    then
        echo
        echo "  [WARNING] Low host classification for sample $i"
        echo "  Host classified: ${host_classified_pct}%"
        echo "  This sample may contain more microbial reads, contamination,"
        echo "  or lower host RNA content."
    fi

    echo
done

echo "============================================================"
echo "[OK] Host Kraken2 summary saved to:"
echo "$summary_file"
echo "============================================================"
column -t "$summary_file"

echo "${g}MTD running  progress:"
echo '>>>>>               [25%]'

# ------------------------------------------------------------
# Reads classification by Kraken2; 2nd step for non-host reads
# Microbiome classification using reads not classified as host
# ------------------------------------------------------------

echo "Reads classification by kraken2; 2nd step for non-host reads ${w}"
echo "Microbiome DB: $DB_micro"

if [[ ! -d "$DB_micro" ]]; then
    echo "[ERROR] Microbiome Kraken2 DB folder not found:"
    echo "$DB_micro"
    exit 1
fi

if [[ ! -s "$DB_micro/hash.k2d" || ! -s "$DB_micro/opts.k2d" || ! -s "$DB_micro/taxo.k2d" ]]; then
    echo "[ERROR] Microbiome Kraken2 DB appears incomplete."
    echo "Expected files:"
    echo "  $DB_micro/hash.k2d"
    echo "  $DB_micro/opts.k2d"
    echo "  $DB_micro/taxo.k2d"
    exit 1
fi

micro_summary="kraken_nonhost_raw_summary.tsv"

printf 'sample\tmicro_classified_reads\tmicro_classified_pct\tmicro_unclassified_reads\tmicro_unclassified_pct\tcount_unit\n' \
    > "$micro_summary"

# Threshold for warning about high microbial classification
# among the reads that were not classified as host.
MICRO_HIGH_WARN=20

for i in $lsn; do

    report="Report_non-host.raw_${i}.txt"
    kraken_output="Report_non-host_raw_${i}.kraken"

    if [[ "$READ_LAYOUT" == "pe" ]]; then
        count_unit="read_pairs"

        raw_nonhost_read1="${i}_non-host_raw_1.fq"
        raw_nonhost_read2="${i}_non-host_raw_2.fq"

        raw_classified_pattern="${i}_raw_cseqs#.fq"
        raw_unclassified_pattern="${i}_raw_ucseqs#.fq"

        raw_classified_read1="${i}_raw_cseqs_1.fq"
        raw_classified_read2="${i}_raw_cseqs_2.fq"

        raw_unclassified_read1="${i}_raw_ucseqs_1.fq"
        raw_unclassified_read2="${i}_raw_ucseqs_2.fq"
    else
        count_unit="reads"

        raw_nonhost_read1="${i}_non-host_raw.fq"
        raw_nonhost_read2="-"

        raw_classified_pattern="${i}_raw_cseqs.fq"
        raw_unclassified_pattern="${i}_raw_ucseqs.fq"

        raw_classified_read1="${i}_raw_cseqs.fq"
        raw_classified_read2="-"

        raw_unclassified_read1="${i}_raw_ucseqs.fq"
        raw_unclassified_read2="-"
    fi

    echo "============================================================"
    echo "[MICRO RAW] Sample: $i"
    echo "Layout: $READ_LAYOUT"
    echo "Count unit: $count_unit"
    echo "Microbiome Kraken2 DB: $DB_micro"
    echo "Kraken2 confidence: $KRAKEN_MICRO_CONF"
    echo "Kraken2 minimum hit groups: $KRAKEN_MICRO_MIN_HIT_GROUPS"

    if [[ "$READ_LAYOUT" == "pe" ]]; then
        echo "Input R1: $raw_nonhost_read1"
        echo "Input R2: $raw_nonhost_read2"
        echo "Classified output R1: $raw_classified_read1"
        echo "Classified output R2: $raw_classified_read2"
        echo "Unclassified output R1: $raw_unclassified_read1"
        echo "Unclassified output R2: $raw_unclassified_read2"
    else
        echo "Input: $raw_nonhost_read1"
        echo "Classified output: $raw_classified_read1"
        echo "Unclassified output: $raw_unclassified_read1"
    fi

    echo "============================================================"

    # The host-filtering outputs are uncompressed FASTQ files.
    if [[ "$READ_LAYOUT" == "pe" ]]; then
        if [[ ! -s "$raw_nonhost_read1" ||
              ! -s "$raw_nonhost_read2" ]]; then
            echo "${r}[ERROR] Missing or empty paired non-host input from host-filtering step:${w}"
            echo "  R1: $raw_nonhost_read1"
            echo "  R2: $raw_nonhost_read2"
            exit 1
        fi
    else
        if [[ ! -s "$raw_nonhost_read1" ]]; then
            echo "${r}[ERROR] Missing or empty non-host input from host-filtering step:${w}"
            echo "  $raw_nonhost_read1"
            exit 1
        fi
    fi

    kraken_micro_args=(
        --db "$DB_micro"
        --use-names
        --confidence "$KRAKEN_MICRO_CONF"
        --minimum-hit-groups "$KRAKEN_MICRO_MIN_HIT_GROUPS"
        --report "$report"
        --output "$kraken_output"
        --threads "$threads"
    )

    if [[ "$READ_LAYOUT" == "pe" ]]; then
        kraken_micro_args+=(
            --paired
            --classified-out "$raw_classified_pattern"
            --unclassified-out "$raw_unclassified_pattern"
            "$raw_nonhost_read1"
            "$raw_nonhost_read2"
        )
    else
        kraken_micro_args+=(
            --classified-out "$raw_classified_pattern"
            --unclassified-out "$raw_unclassified_pattern"
            "$raw_nonhost_read1"
        )
    fi

    if ! kraken2 "${kraken_micro_args[@]}"; then
        die "Kraken2 raw microbiome classification failed for sample: $i"
    fi

    require_file \
        "$report" \
        "Raw microbiome Kraken2 report for sample $i"

    require_file \
        "$kraken_output" \
        "Raw microbiome Kraken2 classification output for sample $i"

    # Classified or unclassified FASTQ files may legitimately be empty.
    # Verify their existence rather than requiring a non-zero size.
    if [[ "$READ_LAYOUT" == "pe" ]]; then
        for output_file in \
            "$raw_classified_read1" \
            "$raw_classified_read2" \
            "$raw_unclassified_read1" \
            "$raw_unclassified_read2"
        do
            if [[ ! -e "$output_file" ]]; then
                die "Expected paired-end raw microbiome output was not created for sample '$i': $output_file"
            fi
        done
    else
        for output_file in \
            "$raw_classified_read1" \
            "$raw_unclassified_read1"
        do
            if [[ ! -e "$output_file" ]]; then
                die "Expected single-end raw microbiome output was not created for sample '$i': $output_file"
            fi
        done
    fi

    micro_unclassified_pct="$(
        awk '$4 == "U" {
            print $1
            exit
        }' "$report"
    )"

    micro_unclassified_reads="$(
        awk '$4 == "U" {
            print $2
            exit
        }' "$report"
    )"

    micro_classified_pct="$(
        awk '$4 == "R" && $5 == 1 {
            print $1
            exit
        }' "$report"
    )"

    micro_classified_reads="$(
        awk '$4 == "R" && $5 == 1 {
            print $2
            exit
        }' "$report"
    )"

    # Kraken2 may omit a zero-count category from the report.
    if [[ -z "$micro_unclassified_pct" ]]; then
        micro_unclassified_pct="0.00"
    fi

    if [[ -z "$micro_unclassified_reads" ]]; then
        micro_unclassified_reads="0"
    fi

    if [[ -z "$micro_classified_reads" ]]; then
        micro_classified_reads="0"
    fi

    if [[ -z "$micro_classified_pct" ]]; then
        micro_classified_pct="$(
            awk -v u="$micro_unclassified_pct" '
                BEGIN {
                    printf "%.2f", 100 - u
                }
            '
        )"
    fi

    echo
    echo "[RESULT] Sample: $i"
    echo "  Classified in DB_micro:   ${micro_classified_pct}% (${micro_classified_reads} ${count_unit})"
    echo "  Unclassified in DB_micro: ${micro_unclassified_pct}% (${micro_unclassified_reads} ${count_unit})"

    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$i" \
        "$micro_classified_reads" \
        "$micro_classified_pct" \
        "$micro_unclassified_reads" \
        "$micro_unclassified_pct" \
        "$count_unit" \
        >> "$micro_summary"

    if awk \
        -v p="$micro_classified_pct" \
        -v t="$MICRO_HIGH_WARN" \
        'BEGIN {
            exit !(p >= t)
        }'
    then
        echo
        echo "  [WARNING] High DB_micro classification for sample $i"
        echo "  Microbial classification here is ${micro_classified_pct}%"
        echo "  of the NON-HOST sequence units, not of the original total."

        echo
        echo "  Top taxa with >=1% in report:"

        awk '
            $4 != "U" &&
            !($4 == "R" && $5 == 1) &&
            $1 >= 1 {
                name = $6

                for (j = 7; j <= NF; j++) {
                    name = name " " $j
                }

                printf "    %7s%%  %12s units  rank=%-4s taxid=%-10s %s\n",
                       $1, $2, $4, $5, name
            }
        ' "$report" | head -n 20
    fi

    echo
done

echo "============================================================"
echo "[OK] Non-host raw Kraken2 summary saved to:"
echo "$micro_summary"
echo "============================================================"
column -s $'\t' -t "$micro_summary"

echo "${g}MTD running  progress:"
echo '>>>>>>              [30%]'

# ------------------------------------------------------------
# Global read composition summary
# Percentages are calculated relative to the original total reads
# ------------------------------------------------------------

echo "============================================================"
echo "[SUMMARY] Creating global Kraken read composition table"
echo "============================================================"

host_summary="kraken_host_summary.tsv"
micro_summary="kraken_nonhost_raw_summary.tsv"
out_summary="kraken_global_read_composition_raw.tsv"

if [[ ! -s "$host_summary" ]]; then
    echo "[ERROR] Missing file: $host_summary"
    exit 1
fi

if [[ ! -s "$micro_summary" ]]; then
    echo "[ERROR] Missing file: $micro_summary"
    exit 1
fi

echo -e "sample\ttotal_reads\thost\tmicrobiome\tunclassified\tcheck_pct_sum" > "$out_summary"

awk '
BEGIN {
    FS=OFS="\t"
}

NR==FNR {
    if (FNR == 1) next

    sample=$1

    host_reads[sample]=$2
    host_unclassified_reads[sample]=$4

    total_reads[sample]=$2 + $4

    next
}

FNR > 1 {
    sample=$1

    micro_reads=$2
    micro_unclassified_reads=$4

    if (!(sample in total_reads)) {
        print "[WARNING] Sample found in micro summary but not in host summary: " sample > "/dev/stderr"
        next
    }

    total=total_reads[sample]
    host=host_reads[sample]
    micro=micro_reads
    unclassified=micro_unclassified_reads

    host_pct=(host/total)*100
    micro_pct=(micro/total)*100
    unclassified_pct=(unclassified/total)*100

    check_sum=host_pct + micro_pct + unclassified_pct

    host_label=sprintf("%d (%.2f%%)", host, host_pct)
    micro_label=sprintf("%d (%.2f%%)", micro, micro_pct)
    unclassified_label=sprintf("%d (%.2f%%)", unclassified, unclassified_pct)

    printf "%s\t%d\t%s\t%s\t%s\t%.2f\n", \
        sample, total, host_label, micro_label, unclassified_label, check_sum
}
' "$host_summary" "$micro_summary" >> "$out_summary"

echo
echo "[OK] Global read composition saved to:"
echo "$out_summary"
echo

column -s $'\t' -t "$out_summary"

echo
if [[ "$READ_LAYOUT" == "pe" ]]; then
    composition_count_unit="read pairs"
else
    composition_count_unit="reads"
fi

echo "[INFO] RAW interpretation:"
echo "  count unit     = $composition_count_unit"
echo "  total_reads    = original sequence units after trimming/compression step"
echo "  host           = reads classified as host in Kraken2 host step"
echo "  microbiome     = reads classified by DB_micro after host removal"
echo "  unclassified   = reads not classified as host and not classified by DB_micro"
echo "  check_pct_sum  = should be close to 100.00"
echo "============================================================"

echo "${g}MTD running  progress:"
echo '>>>>>>              [30%]'
mkdir -p "$outputdr/kraken"
mv kraken_global_read_composition_raw.tsv kraken_host_summary.tsv kraken_nonhost_raw_summary.tsv "$outputdr/kraken/"

echo "Decontamination step${w}"

source "$condapath/etc/profile.d/conda.sh"
conda activate MTD

export PYTHONNOUSERSITE=1
unset PYTHONPATH
unset PYTHONHOME

# ------------------------------------------------------------
# Optional contaminant removal before final microbiome reports
# ------------------------------------------------------------
# If conta_ls.txt exists and contains taxids, remove those taxa from
# the raw non-host reads and re-run Kraken2 to create final reports.
#
# If conta_ls.txt does not exist, or is empty, use the raw non-host
# Kraken2 outputs as the final non-host outputs.
#
# This guarantees that Bracken always receives:
#   Report_non-host_${sample}.txt
# ------------------------------------------------------------

# ------------------------------------------------------------
# Resolve host-classified FASTQ filenames
# ------------------------------------------------------------
# These files were generated by the first Kraken2 host-
# filtering step.
#
# Global variables populated:
#   HOST_KRAKEN_R1
#   HOST_KRAKEN_R2
# ------------------------------------------------------------
set_host_kraken_fastq_paths() {
    local sample="$1"

    if [[ "$READ_LAYOUT" == "pe" ]]; then
        HOST_KRAKEN_R1="$PIPELINE_TEMP_DIR/${sample}_host_1.fq"
        HOST_KRAKEN_R2="$PIPELINE_TEMP_DIR/${sample}_host_2.fq"
    else
        HOST_KRAKEN_R1="$PIPELINE_TEMP_DIR/${sample}_host.fq"
        HOST_KRAKEN_R2="-"
    fi
}

# ------------------------------------------------------------
# Resolve raw and final microbiome FASTQ filenames
# ------------------------------------------------------------
# Global variables populated for the requested sample.
# ------------------------------------------------------------
set_microbiome_fastq_paths() {
    local sample="$1"

    RAW_KRAKEN_REPORT="Report_non-host.raw_${sample}.txt"
    RAW_KRAKEN_OUTPUT="Report_non-host_raw_${sample}.kraken"

    FINAL_KRAKEN_REPORT="Report_non-host_${sample}.txt"
    FINAL_KRAKEN_OUTPUT="Report_non-host_${sample}.kraken"

    if [[ "$READ_LAYOUT" == "pe" ]]; then
        RAW_NONHOST_R1="${sample}_non-host_raw_1.fq"
        RAW_NONHOST_R2="${sample}_non-host_raw_2.fq"

        FINAL_NONHOST_R1="${sample}_non-host_1.fq"
        FINAL_NONHOST_R2="${sample}_non-host_2.fq"

        RAW_CLASSIFIED_R1="${sample}_raw_cseqs_1.fq"
        RAW_CLASSIFIED_R2="${sample}_raw_cseqs_2.fq"

        RAW_UNCLASSIFIED_R1="${sample}_raw_ucseqs_1.fq"
        RAW_UNCLASSIFIED_R2="${sample}_raw_ucseqs_2.fq"

        FINAL_CLASSIFIED_R1="${sample}_cseqs_1.fq"
        FINAL_CLASSIFIED_R2="${sample}_cseqs_2.fq"

        FINAL_UNCLASSIFIED_R1="${sample}_ucseqs_1.fq"
        FINAL_UNCLASSIFIED_R2="${sample}_ucseqs_2.fq"
    else
        RAW_NONHOST_R1="${sample}_non-host_raw.fq"
        RAW_NONHOST_R2="-"

        FINAL_NONHOST_R1="${sample}_non-host.fq"
        FINAL_NONHOST_R2="-"

        RAW_CLASSIFIED_R1="${sample}_raw_cseqs.fq"
        RAW_CLASSIFIED_R2="-"

        RAW_UNCLASSIFIED_R1="${sample}_raw_ucseqs.fq"
        RAW_UNCLASSIFIED_R2="-"

        FINAL_CLASSIFIED_R1="${sample}_cseqs.fq"
        FINAL_CLASSIFIED_R2="-"

        FINAL_UNCLASSIFIED_R1="${sample}_ucseqs.fq"
        FINAL_UNCLASSIFIED_R2="-"
    fi
}

# ------------------------------------------------------------
# Use raw microbiome outputs as final outputs
# ------------------------------------------------------------
# Used when no valid contaminant taxids are available.
# ------------------------------------------------------------
copy_raw_microbiome_outputs_to_final() {
    local sample="$1"

    set_microbiome_fastq_paths "$sample"

    require_file \
        "$RAW_KRAKEN_REPORT" \
        "Raw microbiome Kraken2 report for sample $sample"

    require_file \
        "$RAW_KRAKEN_OUTPUT" \
        "Raw microbiome Kraken2 output for sample $sample"

    require_file \
        "$RAW_NONHOST_R1" \
        "Raw non-host read 1 for sample $sample"

    if [[ "$READ_LAYOUT" == "pe" ]]; then
        require_file \
            "$RAW_NONHOST_R2" \
            "Raw non-host read 2 for sample $sample"
    fi

    cp -- \
        "$RAW_KRAKEN_REPORT" \
        "$FINAL_KRAKEN_REPORT" || \
        die "Failed to copy final Kraken2 report for sample: $sample"

    cp -- \
        "$RAW_KRAKEN_OUTPUT" \
        "$FINAL_KRAKEN_OUTPUT" || \
        die "Failed to copy final Kraken2 output for sample: $sample"

    cp -- \
        "$RAW_NONHOST_R1" \
        "$FINAL_NONHOST_R1" || \
        die "Failed to copy final non-host read 1 for sample: $sample"

    # Classified and unclassified FASTQ outputs may legitimately
    # exist with zero bytes, so test existence instead of size.
    if [[ ! -e "$RAW_CLASSIFIED_R1" ]]; then
        die "Raw classified FASTQ was not created for sample '$sample': $RAW_CLASSIFIED_R1"
    fi

    if [[ ! -e "$RAW_UNCLASSIFIED_R1" ]]; then
        die "Raw unclassified FASTQ was not created for sample '$sample': $RAW_UNCLASSIFIED_R1"
    fi

    cp -- \
        "$RAW_CLASSIFIED_R1" \
        "$FINAL_CLASSIFIED_R1" || \
        die "Failed to copy classified FASTQ for sample: $sample"

    cp -- \
        "$RAW_UNCLASSIFIED_R1" \
        "$FINAL_UNCLASSIFIED_R1" || \
        die "Failed to copy unclassified FASTQ for sample: $sample"

    if [[ "$READ_LAYOUT" == "pe" ]]; then
        cp -- \
            "$RAW_NONHOST_R2" \
            "$FINAL_NONHOST_R2" || \
            die "Failed to copy final non-host read 2 for sample: $sample"

        if [[ ! -e "$RAW_CLASSIFIED_R2" ]]; then
            die "Raw classified R2 FASTQ was not created for sample '$sample': $RAW_CLASSIFIED_R2"
        fi

        if [[ ! -e "$RAW_UNCLASSIFIED_R2" ]]; then
            die "Raw unclassified R2 FASTQ was not created for sample '$sample': $RAW_UNCLASSIFIED_R2"
        fi

        cp -- \
            "$RAW_CLASSIFIED_R2" \
            "$FINAL_CLASSIFIED_R2" || \
            die "Failed to copy classified R2 FASTQ for sample: $sample"

        cp -- \
            "$RAW_UNCLASSIFIED_R2" \
            "$FINAL_UNCLASSIFIED_R2" || \
            die "Failed to copy unclassified R2 FASTQ for sample: $sample"
    fi
}

conta_file="$MTDIR/conta_ls.txt"

# Safety defaults, caso você ainda não tenha definido esses parâmetros antes
KRAKEN_MICRO_CONF="${KRAKEN_MICRO_CONF:-0.10}"
KRAKEN_MICRO_MIN_HIT_GROUPS="${KRAKEN_MICRO_MIN_HIT_GROUPS:-3}"

# ------------------------------------------------------------
# Detect valid contaminant taxids
# ------------------------------------------------------------
conta_ls=""
valid_contaminant_list=0

if [[ -s "$conta_file" ]]; then
    echo "[INFO] Contaminant list found:"
    echo "  $conta_file"

    # Taxids are expected in column 2, tab-separated.
    conta_ls="$(
        awk -F $'\t' '
            {
                gsub(/\r/, "", $2)

                if ($2 ~ /^[0-9]+$/) {
                    print $2
                }
            }
        ' "$conta_file" |
        sort -u |
        paste -sd ' ' -
    )"

    if [[ -n "$conta_ls" ]]; then
        valid_contaminant_list=1

        echo "[INFO] TaxIDs to exclude:"
        echo "  $conta_ls"
    else
        echo "${y}[WARNING] conta_ls.txt exists, but no valid taxids were found in column 2.${w}"
    fi
else
    echo "[INFO] No contaminant list found or file is empty:"
    echo "  $conta_file"
fi


# ------------------------------------------------------------
# Path A: remove contaminants and reclassify final reads
# ------------------------------------------------------------
if [[ "$valid_contaminant_list" == "1" ]]; then

    read -r -a contaminant_taxids <<< "$conta_ls"

    for i in $lsn; do
        set_microbiome_fastq_paths "$i"

        echo "============================================================"
        echo "[DECONTAMINATION] Sample: $i"
        echo "Layout: $READ_LAYOUT"
        echo "Input Kraken output: $RAW_KRAKEN_OUTPUT"
        echo "Input Kraken report: $RAW_KRAKEN_REPORT"

        if [[ "$READ_LAYOUT" == "pe" ]]; then
            echo "Input R1:  $RAW_NONHOST_R1"
            echo "Input R2:  $RAW_NONHOST_R2"
            echo "Output R1: $FINAL_NONHOST_R1"
            echo "Output R2: $FINAL_NONHOST_R2"
        else
            echo "Input:  $RAW_NONHOST_R1"
            echo "Output: $FINAL_NONHOST_R1"
        fi

        echo "============================================================"

        require_file \
            "$RAW_KRAKEN_OUTPUT" \
            "Raw microbiome Kraken2 output for sample $i"

        require_file \
            "$RAW_KRAKEN_REPORT" \
            "Raw microbiome Kraken2 report for sample $i"

        require_file \
            "$RAW_NONHOST_R1" \
            "Raw non-host read 1 for sample $i"

        extract_args=(
            "$MTDIR/Tools/KrakenTools/extract_kraken_reads.py"
            -k "$RAW_KRAKEN_OUTPUT"
            -s1 "$RAW_NONHOST_R1"
            -o "$FINAL_NONHOST_R1"
            -r "$RAW_KRAKEN_REPORT"
            --fastq-output
        )

        if [[ "$READ_LAYOUT" == "pe" ]]; then
            require_file \
                "$RAW_NONHOST_R2" \
                "Raw non-host read 2 for sample $i"

            extract_args+=(
                -s2 "$RAW_NONHOST_R2"
                -o2 "$FINAL_NONHOST_R2"
            )
        fi

        extract_args+=(
            --taxid "${contaminant_taxids[@]}"
            --exclude
            --include-children
        )

        if ! python "${extract_args[@]}"; then
            die "KrakenTools contaminant removal failed for sample: $i"
        fi

        require_file \
            "$FINAL_NONHOST_R1" \
            "Decontaminated non-host read 1 for sample $i"

                if [[ "$READ_LAYOUT" == "pe" ]]; then
            require_file \
                "$FINAL_NONHOST_R2" \
                "Decontaminated non-host read 2 for sample $i"

            validate_fastq_pair \
                "$FINAL_NONHOST_R1" \
                "$FINAL_NONHOST_R2" \
                "decontaminated non-host pair for sample $i"
        fi
    done

    echo "${g}MTD running  progress:"
    echo '>>>>>>>             [35%]'

    echo "Reads classification by Kraken2; 3rd step for decontaminated non-host reads to get final reports ${w}"

    for i in $lsn; do
        set_microbiome_fastq_paths "$i"

        if [[ "$READ_LAYOUT" == "pe" ]]; then
            final_classified_pattern="${i}_cseqs#.fq"
            final_unclassified_pattern="${i}_ucseqs#.fq"
        else
            final_classified_pattern="$FINAL_CLASSIFIED_R1"
            final_unclassified_pattern="$FINAL_UNCLASSIFIED_R1"
        fi

        echo "============================================================"
        echo "[MICRO FINAL] Sample: $i"
        echo "Layout: $READ_LAYOUT"
        echo "Kraken2 confidence: $KRAKEN_MICRO_CONF"
        echo "Kraken2 minimum hit groups: $KRAKEN_MICRO_MIN_HIT_GROUPS"

        if [[ "$READ_LAYOUT" == "pe" ]]; then
            echo "Input R1: $FINAL_NONHOST_R1"
            echo "Input R2: $FINAL_NONHOST_R2"
        else
            echo "Input: $FINAL_NONHOST_R1"
        fi

        echo "Final report: $FINAL_KRAKEN_REPORT"
        echo "============================================================"

        kraken_final_args=(
            --db "$DB_micro"
            --use-names
            --confidence "$KRAKEN_MICRO_CONF"
            --minimum-hit-groups "$KRAKEN_MICRO_MIN_HIT_GROUPS"
            --report "$FINAL_KRAKEN_REPORT"
            --output "$FINAL_KRAKEN_OUTPUT"
            --threads "$threads"
        )

        if [[ "$READ_LAYOUT" == "pe" ]]; then
            kraken_final_args+=(
                --paired
                --classified-out "$final_classified_pattern"
                --unclassified-out "$final_unclassified_pattern"
                "$FINAL_NONHOST_R1"
                "$FINAL_NONHOST_R2"
            )
        else
            kraken_final_args+=(
                --classified-out "$final_classified_pattern"
                --unclassified-out "$final_unclassified_pattern"
                "$FINAL_NONHOST_R1"
            )
        fi

        if ! kraken2 "${kraken_final_args[@]}"; then
            die "Final Kraken2 microbiome classification failed for sample: $i"
        fi

        require_file \
            "$FINAL_KRAKEN_REPORT" \
            "Final Kraken2 microbiome report for sample $i"

        require_file \
            "$FINAL_KRAKEN_OUTPUT" \
            "Final Kraken2 microbiome classification output for sample $i"

        if [[ "$READ_LAYOUT" == "pe" ]]; then
            for output_file in \
                "$FINAL_CLASSIFIED_R1" \
                "$FINAL_CLASSIFIED_R2" \
                "$FINAL_UNCLASSIFIED_R1" \
                "$FINAL_UNCLASSIFIED_R2"
            do
                if [[ ! -e "$output_file" ]]; then
                    die "Expected final paired-end Kraken2 output was not created for sample '$i': $output_file"
                fi
            done
        else
            for output_file in \
                "$FINAL_CLASSIFIED_R1" \
                "$FINAL_UNCLASSIFIED_R1"
            do
                if [[ ! -e "$output_file" ]]; then
                    die "Expected final single-end Kraken2 output was not created for sample '$i': $output_file"
                fi
            done
        fi
    done


# ------------------------------------------------------------
# Path B: no valid contaminant taxids
# ------------------------------------------------------------
else
    echo "[INFO] Skipping contaminant removal."
    echo "[INFO] Using raw non-host Kraken2 outputs as final outputs."

    for i in $lsn; do
        echo "============================================================"
        echo "[NO DECONTAMINATION] Sample: $i"
        echo "Layout: $READ_LAYOUT"
        echo "Using raw microbiome outputs as final outputs"
        echo "============================================================"

        copy_raw_microbiome_outputs_to_final "$i"
    done
fi

# ------------------------------------------------------------
# Final Kraken read composition summary
# Uses final non-host reports after optional contaminant removal
# ------------------------------------------------------------

echo "============================================================"
echo "[SUMMARY] Creating FINAL Kraken read composition table"
echo "============================================================"

final_micro_summary="$outputdr/kraken/kraken_nonhost_final_summary.tsv"
final_out_summary="$outputdr/kraken/kraken_global_read_composition_final.tsv"
host_summary="$outputdr/kraken/kraken_host_summary.tsv"

if [[ ! -s "$host_summary" ]]; then
    echo "${r}[ERROR] Missing host summary:${w}"
    echo "  $host_summary"
    exit 1
fi

printf 'sample\tmicro_classified_reads\tmicro_classified_pct\tmicro_unclassified_reads\tmicro_unclassified_pct\tcount_unit\n' \
    > "$final_micro_summary"

for i in $lsn; do
    report="Report_non-host_${i}.txt"
    if [[ "$READ_LAYOUT" == "pe" ]]; then
        count_unit="read_pairs"
    else
        count_unit="reads"
    fi

    echo "============================================================"
    echo "[FINAL MICRO SUMMARY] Sample: $i"
    echo "Final Kraken2 report: $report"
    echo "============================================================"

    if [[ ! -s "$report" ]]; then
        echo "${r}[ERROR] Missing final Kraken2 report:${w}"
        echo "  $report"
        exit 1
    fi

    micro_unclassified_pct=$(awk '$4=="U"{print $1; exit}' "$report")
    micro_unclassified_reads=$(awk '$4=="U"{print $2; exit}' "$report")

    micro_classified_pct=$(awk '$4=="R" && $5==1{print $1; exit}' "$report")
    micro_classified_reads=$(awk '$4=="R" && $5==1{print $2; exit}' "$report")

# Kraken2 may omit zero-count categories from the report.
    if [[ -z "$micro_unclassified_pct" ]]; then
        micro_unclassified_pct="0.00"
    fi

    if [[ -z "$micro_unclassified_reads" ]]; then
        micro_unclassified_reads="0"
    fi

    if [[ -z "$micro_classified_reads" ]]; then
        micro_classified_reads="0"
    fi

    if [[ -z "$micro_classified_pct" ]]; then
        micro_classified_pct="$(
            awk -v u="$micro_unclassified_pct" '
                BEGIN {
                    printf "%.2f", 100 - u
                }
            '
        )"
    fi

    echo "[RESULT] Final DB_micro classified:   ${micro_classified_pct}%  (${micro_classified_reads} ${count_unit})"
    echo "[RESULT] Final DB_micro unclassified: ${micro_unclassified_pct}%  (${micro_unclassified_reads} ${count_unit})"

        printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$i" \
        "$micro_classified_reads" \
        "$micro_classified_pct" \
        "$micro_unclassified_reads" \
        "$micro_unclassified_pct" \
        "$count_unit" \
        >> "$final_micro_summary"
done

echo
echo "[OK] Final non-host Kraken2 summary saved to:"
echo "$final_micro_summary"
echo
column -s $'\t' -t "$final_micro_summary"

echo -e "sample\ttotal_reads\thost\tmicrobiome\tunclassified\tcheck_pct_sum" > "$final_out_summary"

awk '
BEGIN {
    FS=OFS="\t"
}

NR==FNR {
    if (FNR == 1) next

    sample=$1

    host_reads[sample]=$2
    host_unclassified_reads[sample]=$4
    total_reads[sample]=$2 + $4

    next
}

FNR > 1 {
    sample=$1

    micro_reads=$2
    micro_unclassified_reads=$4

    if (!(sample in total_reads)) {
        print "[WARNING] Sample found in final micro summary but not in host summary: " sample > "/dev/stderr"
        next
    }

    total=total_reads[sample]
    host=host_reads[sample]
    micro=micro_reads
    unclassified=micro_unclassified_reads

    host_pct=(host/total)*100
    micro_pct=(micro/total)*100
    unclassified_pct=(unclassified/total)*100

    check_sum=host_pct + micro_pct + unclassified_pct

    host_label=sprintf("%d (%.2f%%)", host, host_pct)
    micro_label=sprintf("%d (%.2f%%)", micro, micro_pct)
    unclassified_label=sprintf("%d (%.2f%%)", unclassified, unclassified_pct)

    printf "%s\t%d\t%s\t%s\t%s\t%.2f\n", \
        sample, total, host_label, micro_label, unclassified_label, check_sum
}
' "$host_summary" "$final_micro_summary" >> "$final_out_summary"

echo
echo "[OK] FINAL global read composition saved to:"
echo "$final_out_summary"
echo
column -s $'\t' -t "$final_out_summary"

echo
if [[ "$READ_LAYOUT" == "pe" ]]; then
    composition_count_unit="read pairs"
else
    composition_count_unit="reads"
fi

echo "[INFO] FINAL interpretation:"
echo "  count unit     = $composition_count_unit"
echo "  total_reads    = original sequence units after trimming/compression step"
echo "  host           = reads classified as host in Kraken2 host step"
echo "  microbiome     = reads classified by DB_micro after optional contaminant removal"
echo "  unclassified   = reads not classified as host and not classified by final DB_micro step"
echo "  check_pct_sum  = should be close to 100.00"
echo "============================================================"

# ------------------------------------------------------------
# Save individual Kraken2 reports and outputs
# ------------------------------------------------------------

echo "============================================================"
echo "[SUMMARY] Saving individual Kraken2 reports to kraken folder"
echo "============================================================"

mkdir -p "$outputdr/kraken/reports_host"
mkdir -p "$outputdr/kraken/reports_micro_raw"
mkdir -p "$outputdr/kraken/reports_micro_final"

# Host reports
cp Report_host_*.txt Report_host_*.kraken "$outputdr/kraken/reports_host/" 2>/dev/null || true

# Raw microbiome reports
cp Report_non-host.raw_*.txt Report_non-host_raw_*.kraken "$outputdr/kraken/reports_micro_raw/" 2>/dev/null || true

# Final microbiome reports used by Bracken
cp Report_non-host_*.txt Report_non-host_*.kraken "$outputdr/kraken/reports_micro_final/" 2>/dev/null || true

echo "[OK] Kraken2 individual reports copied to:"
echo "  $outputdr/kraken/reports_host"
echo "  $outputdr/kraken/reports_micro_raw"
echo "  $outputdr/kraken/reports_micro_final"
echo "============================================================"

echo "${g}MTD running  progress:"
echo '>>>>>>>>            [40%]'

echo "Bracken analysis ${w}"

# ------------------------------------------------------------
# Check Bracken read-length distribution file
# ------------------------------------------------------------
# Bracken needs a read-length-specific distribution file.
# Example: if read_len=75, the database must contain:
#   database75mers.kmer_distrib
#
# This prevents running Bracken with a read length that was not built
# for the current Kraken2 microbiome database.

BRACKEN_DIST="$DB_micro/database${read_len}mers.kmer_distrib"

echo "Bracken DB: $DB_micro"
echo "Bracken read length: $read_len"
echo "Expected Bracken distribution file: $BRACKEN_DIST"

if [[ ! -s "$BRACKEN_DIST" ]]; then
    echo "${r}[ERROR] Bracken distribution file not found for read length ${read_len}.${w}"
    echo
    echo "Expected file:"
    echo "  $BRACKEN_DIST"
    echo
    echo "Available Bracken distribution files in DB_micro:"
    ls -lh "$DB_micro"/database*mers.kmer_distrib 2>/dev/null || echo "  None found."
    echo
    echo "You probably need to build the Bracken distribution file for read length ${read_len}."
    echo "Example command:"
    echo "  bracken-build -d \"$DB_micro\" -t \"$threads\" -k 35 -l \"$read_len\""
    echo
    exit 1
fi

echo "${g}[OK] Bracken distribution file found:${w} $BRACKEN_DIST"

# Bracken -t is the minimum read threshold, not threads.
# Keep it fixed so results do not change when CPU thread number changes.

echo "Bracken threshold: $BRACKEN_THRESHOLD"

for i in $lsn; do
    echo "============================================================"
    echo "[BRACKEN] Sample: $i"
    echo "Input report: Report_non-host_${i}.txt"
    echo "Read length: $read_len"
    echo "Threshold: $BRACKEN_THRESHOLD"
    echo "============================================================"

    if [[ ! -s "Report_non-host_${i}.txt" ]]; then
        echo "${r}[ERROR] Missing Kraken2 report for Bracken:${w}"
        echo "  Report_non-host_${i}.txt"
        exit 1
    fi

    bracken -d "$DB_micro" \
        -i "Report_non-host_${i}.txt" \
        -o "Report_$i.phylum.bracken" \
        -r "$read_len" \
        -l P \
        -t "$BRACKEN_THRESHOLD"

    bracken -d "$DB_micro" \
        -i "Report_non-host_${i}.txt" \
        -o "Report_$i.genus.bracken" \
        -r "$read_len" \
        -l G \
        -t "$BRACKEN_THRESHOLD"

    bracken -d "$DB_micro" \
        -i "Report_non-host_${i}.txt" \
        -o "Report_$i.species.bracken" \
        -r "$read_len" \
        -l S \
        -t "$BRACKEN_THRESHOLD"
done

echo "${g}MTD running  progress:"
echo '>>>>>>>>>           [45%]'
echo "combined .bracken files (table like) into a single outputdr for Deseq2 ${w}"
python $MTDIR/Tools/combine_bracken_outputs.py --files *.phylum.bracken -o $outputdr/bracken_phylum_all
python $MTDIR/Tools/combine_bracken_outputs.py --files *.genus.bracken -o $outputdr/bracken_genus_all
python $MTDIR/Tools/combine_bracken_outputs.py --files *.species.bracken -o $outputdr/bracken_species_all

echo "${g}Move _bracken report files (tree like) to a separate folder${w}"
mkdir -p Report_non-host_bracken_species_normalized
mv *_bracken_species.txt Report_non-host_bracken_species_normalized
cd Report_non-host_bracken_species_normalized

echo "${g}Trim the name of _bracken report files (tree like) to the sample name (eg. DJ01) ${w}"
for i in $lsn; do
    mv *${i}_* $i
done

echo "${g}Converted original _bracken report files (tree like) into .biom file for ANCOMBC and diversity analysis in phyloseq (R) etc. in DEG_Anno_Plot.R ${w}"
kraken-biom * -o $outputdr/temp/bracken_species_all0.biom --fmt json

echo "${g}Adjust bracken file (tree like) for additional visualization (.biom, .mpa, .krona) ${w}"

if [[ "$NO_COMPARISON" == "1" ]]; then
    echo "${y}[INFO] Exploratory mode: skipping DESeq2-based Bracken normalization.${w}"
    echo "${y}[INFO] Using original Bracken tree-like reports for BIOM, Krona, MPA and GraPhlAn.${w}"
else
    conda deactivate
    conda activate R412

    norm_args=(
        "$outputdr/bracken_species_all"
        "$samplesheet_file"
        "$outputdr/temp/Report_non-host_bracken_species_normalized"
    )

    if [[ -n "${metadata:-}" ]]; then
        norm_args+=( "$metadata" )
    fi

    if ! Rscript "$MTDIR/Normalization_afbr.R" "${norm_args[@]}"; then
        echo "${y}[WARNING] Normalization_afbr.R failed.${w}"
        echo "${y}[WARNING] Continuing with original Bracken tree-like reports for visualization.${w}"
    fi

    conda deactivate
    conda activate MTD
fi

echo "${g}MTD running  progress:"
echo '>>>>>>>>>>          [50%]'

echo "Converted adjusted _bracken report files for visualization: GraPhlAn, MPA, Krona ${w}"

# ------------------------------------------------------------
# Important:
# At this point the script should still be inside:
#   $outputdr/temp/Report_non-host_bracken_species_normalized
#
# These files are normalized Kraken/Bracken-style reports.
# We will:
#   1) keep BIOM generation for compatibility with the old pipeline
#   2) generate Krona files
#   3) generate MPA files
#   4) combine MPA files
#   5) use the combined MPA table as input for export2graphlan
#
# This avoids the fragile BIOM parser inside export2graphlan.py.
# ------------------------------------------------------------

cd "$outputdr/temp/Report_non-host_bracken_species_normalized" || exit 1

conda deactivate
conda activate MTD

# ------------------------------------------------------------
# Keep BIOM generation for compatibility with downstream outputs
# ------------------------------------------------------------

echo "${g}Creating BIOM file from normalized species reports${w}"

kraken-biom * -o "$outputdr/bracken_species_all.biom" --fmt json

if [[ ! -s "$outputdr/bracken_species_all.biom" ]]; then
    echo "${r}[ERROR] Failed to create BIOM:${w}"
    echo "  $outputdr/bracken_species_all.biom"
    exit 1
fi

echo "${g}[OK] BIOM saved:${w}"
echo "  $outputdr/bracken_species_all.biom"

# ------------------------------------------------------------
# Prepare output folders
# ------------------------------------------------------------

mkdir -p "$outputdr/graphlan"
mkdir -p "$outputdr/krona"
mkdir -p "$outputdr/temp/graphlan_mpa_input"

GRAPHLAN_DIR="$outputdr/graphlan"
GRAPHLAN_MPA_DIR="$outputdr/temp/graphlan_mpa_input"
GRAPHLAN_COMBINED_MPA="$GRAPHLAN_DIR/Combined_for_graphlan.mpa"
GRAPHLAN_INPUT="$GRAPHLAN_DIR/graphlan_input.clean.tsv"

rm -f "$GRAPHLAN_MPA_DIR"/*.mpa.txt
rm -f "$GRAPHLAN_COMBINED_MPA" "$GRAPHLAN_INPUT"

# ------------------------------------------------------------
# Generate Krona, Krona HTML, and MPA files from normalized reports
# ------------------------------------------------------------

echo "${g}Generating Krona, Krona HTML, and MPA files from normalized reports${w}"

mkdir -p "$outputdr/krona"
mkdir -p "$GRAPHLAN_MPA_DIR"

if ! command -v ktImportText >/dev/null 2>&1; then
    echo "${r}[ERROR] ktImportText was not found in PATH.${w}"
    echo "Please install/activate KronaTools first, for example:"
    echo "  conda install -c bioconda krona"
    exit 1
fi

krona_all_inputs=()

for i in $lsn; do
    if [[ ! -s "$i" ]]; then
        echo "${r}[ERROR] Missing normalized Bracken/Kraken report for sample:${w} $i"
        echo "Expected file:"
        echo "  $outputdr/temp/Report_non-host_bracken_species_normalized/$i"
        exit 1
    fi

    sample=$(basename "$i")

    krona_file="$outputdr/krona/${sample}-bracken.krona"
    krona_html="$outputdr/krona/${sample}-bracken.html"
    mpa_file="$GRAPHLAN_MPA_DIR/${sample}-bracken.mpa.txt"

    echo "  [Krona] $sample"
    python "$MTDIR/Tools/KrakenTools/kreport2krona.py" \
        -r "$i" \
        -o "$krona_file"

    if [[ ! -s "$krona_file" ]]; then
        echo "${r}[ERROR] Krona file was not generated or is empty:${w} $krona_file"
        exit 1
    fi

    echo "  [Krona HTML] $sample"
    ktImportText "$krona_file" \
        -o "$krona_html"

    if [[ ! -s "$krona_html" ]]; then
        echo "${r}[ERROR] Krona HTML was not generated or is empty:${w} $krona_html"
        exit 1
    fi

    krona_all_inputs+=("$krona_file")

    echo "  [MPA] $sample"
    python "$MTDIR/Tools/KrakenTools/kreport2mpa.py" \
        --display-header \
        -r "$i" \
        -o "$mpa_file"
done

# ------------------------------------------------------------
# Generate combined Krona HTML for all detected taxa across animals
# ------------------------------------------------------------

if [[ ${#krona_all_inputs[@]} -gt 0 ]]; then
    all_krona="$outputdr/krona/all_animals-bracken.krona"
    all_html="$outputdr/krona/all_animals-bracken.html"

    echo "${g}Generating combined Krona file for all animals${w}"
    cat "${krona_all_inputs[@]}" > "$all_krona"

    echo "${g}Generating combined Krona HTML for all animals${w}"
    ktImportText "$all_krona" \
        -o "$all_html"

    if [[ ! -s "$all_html" ]]; then
        echo "${r}[ERROR] Combined Krona HTML was not generated or is empty:${w} $all_html"
        exit 1
    fi

    echo "${g}Combined Krona HTML generated:${w}"
    echo "  $all_html"
fi
# ------------------------------------------------------------
# Combine MPA files
# ------------------------------------------------------------

echo "${g}Combining MPA files${w}"

python "$MTDIR/Tools/KrakenTools/combine_mpa.py" \
    -i "$GRAPHLAN_MPA_DIR"/*.mpa.txt \
    -o "$GRAPHLAN_COMBINED_MPA"

if [[ ! -s "$GRAPHLAN_COMBINED_MPA" ]]; then
    echo "${r}[ERROR] Combined MPA file was not created:${w}"
    echo "  $GRAPHLAN_COMBINED_MPA"
    exit 1
fi

# Preserve old expected output location/name
cp "$GRAPHLAN_COMBINED_MPA" "$outputdr/Combined.mpa"

echo "${g}[OK] Combined MPA saved:${w}"
echo "  $GRAPHLAN_COMBINED_MPA"
echo "  $outputdr/Combined.mpa"

# ------------------------------------------------------------
# Clean MPA table for export2graphlan
# ------------------------------------------------------------
# This avoids the BIOM parser entirely.
# It ensures:
#   - first column is always a taxonomy string
#   - abundance columns are numeric
#   - existing k__/p__/c__/... prefixes are preserved
#   - empty / NA / unknown rows are skipped
#   - spaces and problematic characters are cleaned
# ------------------------------------------------------------

echo "${g}Preparing clean MPA/LEfSe-like table for export2graphlan${w}"

python - "$GRAPHLAN_COMBINED_MPA" "$GRAPHLAN_INPUT" << 'PY'
import sys
import csv
import math
import re

inp = sys.argv[1]
out = sys.argv[2]

rank_prefixes = ["k__", "p__", "c__", "o__", "f__", "g__", "s__", "t__"]

def bad(x):
    if x is None:
        return True
    s = str(x).strip()
    return s == "" or s.lower() in ("nan", "na", "none", "null")

def clean_name(x):
    x = str(x).strip()
    x = x.replace("sp. ", "sp_")
    x = x.replace(" ", "_")
    x = x.replace("'", "")
    x = x.replace('"', "")
    x = x.replace("[", "")
    x = x.replace("]", "")
    x = x.replace("{", "")
    x = x.replace("}", "")
    x = x.replace("(", "")
    x = x.replace(")", "")
    x = x.replace("=", "_")
    x = x.replace(",", "_")
    x = re.sub(r"_+", "_", x)
    x = x.strip("|._")
    return x

def clean_taxonomy(x):
    if bad(x):
        return None

    x = str(x).strip()

    if x.lower() in ("root", "unclassified", "unknown"):
        return None

    if "|" in x:
        parts = x.split("|")
    elif ";" in x:
        parts = x.split(";")
    else:
        parts = [x]

    clean_parts = []

    for idx, part in enumerate(parts):
        part = str(part).strip()

        if bad(part):
            continue

        # Preserve existing prefixes such as k__, p__, c__, o__, f__, g__, s__
        m = re.match(r"^([a-zA-Z])__(.+)$", part)

        if m:
            rank = m.group(1).lower()
            name = clean_name(m.group(2))
            if name:
                clean_parts.append(f"{rank}__{name}")
        else:
            name = clean_name(part)
            if name:
                prefix = rank_prefixes[idx] if idx < len(rank_prefixes) else "x__"
                clean_parts.append(prefix + name)

    if not clean_parts:
        return None

    return "|".join(clean_parts)

def clean_number(x):
    if bad(x):
        return 0.0

    s = str(x).strip()
    s = s.replace("%", "")
    s = s.replace(",", "")
    s = s.replace('"', "")
    s = s.replace("'", "")

    try:
        v = float(s)
        if math.isnan(v) or math.isinf(v):
            return 0.0
        return v
    except Exception:
        return 0.0

with open(inp, "r", newline="") as f:
    reader = csv.reader(f, delimiter="\t")
    rows = [r for r in reader if r and any(str(c).strip() for c in r)]

if not rows:
    raise SystemExit("[ERROR] Combined MPA file is empty.")

header = rows[0]

if len(header) < 2:
    raise SystemExit("[ERROR] Combined MPA header has fewer than 2 columns.")

first_cell = header[0].strip().lower()

if first_cell.startswith("#") or "classification" in first_cell or "clade" in first_cell:
    data_rows = rows[1:]
    sample_names = [
        str(x).strip() if str(x).strip() else "sample_%d" % idx
        for idx, x in enumerate(header[1:], start=1)
    ]
else:
    data_rows = rows
    n_samples = len(rows[0]) - 1
    sample_names = ["sample_%d" % i for i in range(1, n_samples + 1)]

clean_rows = []
fixed_short_rows = 0
debug_seen = []
debug_zero_sum = 0
debug_no_tax = 0

for r in data_rows:
    if len(r) < 2:
        continue

    raw_tax = r[0]

    if str(raw_tax).startswith("#"):
        continue

    tax = clean_taxonomy(raw_tax)

    if tax is None:
        debug_no_tax += 1
        if len(debug_seen) < 10:
            debug_seen.append(("NO_TAX", raw_tax))
        continue

    raw_vals = r[1:]

    # Fix rare combine_mpa formatting issue like:
    # 0 0 0 020
    # where "020" should probably be "0" and "20".
    if len(raw_vals) == len(sample_names) - 1:
        last = str(raw_vals[-1]).strip()
        if re.match(r"^0[0-9]+$", last):
            raw_vals = raw_vals[:-1] + ["0", last[1:]]
            fixed_short_rows += 1

    vals = [clean_number(x) for x in raw_vals]

    if len(vals) < len(sample_names):
        vals = vals + [0.0] * (len(sample_names) - len(vals))
    elif len(vals) > len(sample_names):
        vals = vals[:len(sample_names)]

    if sum(vals) <= 0:
        debug_zero_sum += 1
        if len(debug_seen) < 10:
            debug_seen.append(("ZERO_SUM", raw_tax))
        continue

    clean_rows.append([tax] + vals)

if not clean_rows:
    print("[DEBUG] Total input rows:", len(rows))
    print("[DEBUG] Data rows checked:", len(data_rows))
    print("[DEBUG] Rows skipped because taxonomy could not be parsed:", debug_no_tax)
    print("[DEBUG] Rows skipped because abundance sum was zero:", debug_zero_sum)
    print("[DEBUG] First skipped examples:")
    for reason, example in debug_seen:
        print("  ", reason, repr(example))
    raise SystemExit("[ERROR] No valid taxonomic rows remained for GraPhlAn.")

with open(out, "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t", lineterminator="\n")
    writer.writerow(["ID"] + sample_names)
    for r in clean_rows:
        writer.writerow([r[0]] + ["%.10g" % x for x in r[1:]])

print("[OK] Clean GraPhlAn input created")
print("Input:", inp)
print("Output:", out)
print("Taxa rows:", len(clean_rows))
print("Samples:", len(sample_names))
print("Short rows fixed:", fixed_short_rows)
print("[INFO] First cleaned rows:")
for r in clean_rows[:5]:
    print("  ", r[0])
PY

if [[ ! -s "$GRAPHLAN_INPUT" ]]; then
    echo "${r}[ERROR] Clean GraPhlAn input was not created:${w}"
    echo "  $GRAPHLAN_INPUT"
    exit 1
fi

echo "${g}[OK] Clean GraPhlAn input:${w}"
echo "  $GRAPHLAN_INPUT"

# Remove possible Windows-style carriage returns from the clean GraPhlAn input
sed -i 's/\r$//' "$GRAPHLAN_INPUT"
echo
echo "${g}Preview of GraPhlAn input:${w}"
head -n 5 "$GRAPHLAN_INPUT"

# ------------------------------------------------------------
# GraPhlAn cladogram-like settings
# ------------------------------------------------------------
# You can override these only if GRAPHLAN_KEEP_ENV=1 is set.
# Example:
#   GRAPHLAN_KEEP_ENV=1 GRAPHLAN_TOP=20 bash $SCRIPT_NAME ...
# ------------------------------------------------------------

if [[ "${GRAPHLAN_KEEP_ENV:-0}" != "1" ]]; then
    unset GRAPHLAN_TOP
    unset GRAPHLAN_SIZE
    unset GRAPHLAN_DPI
    unset GRAPHLAN_MAX_CLADE_SIZE
    unset GRAPHLAN_LEVELS
    unset GRAPHLAN_EXTERNAL_LEVELS
    unset GRAPHLAN_BACKGROUND_LEVELS
    unset GRAPHLAN_LEAST_BIOMARKERS
    unset GRAPHLAN_ABUNDANCE_THRESHOLD
fi

GRAPHLAN_TOP="${GRAPHLAN_TOP:-80}"
GRAPHLAN_SIZE="${GRAPHLAN_SIZE:-13.0}"
GRAPHLAN_DPI="${GRAPHLAN_DPI:-600}"
GRAPHLAN_MAX_CLADE_SIZE="${GRAPHLAN_MAX_CLADE_SIZE:-300}"

# Taxonomic levels:
# 1 = kingdom
# 2 = phylum
# 3 = class
# 4 = order
# 5 = family
# 6 = genus
# 7 = species
GRAPHLAN_LEVELS="${GRAPHLAN_LEVELS:-2,3,4,5,6,7}"

# For species-level exploration, put the external labels mainly at species level.
GRAPHLAN_EXTERNAL_LEVELS="${GRAPHLAN_EXTERNAL_LEVELS:-7}"

# Keep broad background clades at phylum/class level.
GRAPHLAN_BACKGROUND_LEVELS="${GRAPHLAN_BACKGROUND_LEVELS:-2,3}"

GRAPHLAN_LEAST_BIOMARKERS="${GRAPHLAN_LEAST_BIOMARKERS:-5}"

# Lower threshold helps species-level labels appear.
GRAPHLAN_ABUNDANCE_THRESHOLD="${GRAPHLAN_ABUNDANCE_THRESHOLD:-0.001}"

echo "${g}GraPhlAn cladogram-like settings:${w}"
echo "  Input table:                $GRAPHLAN_INPUT"
echo "  Top taxa:                   $GRAPHLAN_TOP"
echo "  Figure size:                $GRAPHLAN_SIZE"
echo "  DPI:                        $GRAPHLAN_DPI"
echo "  Max clade size:             $GRAPHLAN_MAX_CLADE_SIZE"
echo "  Internal annotation levels: $GRAPHLAN_LEVELS"
echo "  External annotation levels: $GRAPHLAN_EXTERNAL_LEVELS"
echo "  Background levels:          $GRAPHLAN_BACKGROUND_LEVELS"
echo "  Least biomarkers:           $GRAPHLAN_LEAST_BIOMARKERS"
echo "  Abundance threshold:        $GRAPHLAN_ABUNDANCE_THRESHOLD"

# ------------------------------------------------------------
# Run export2graphlan using MPA/LEfSe-like TSV, not BIOM
# ------------------------------------------------------------

cd "$GRAPHLAN_DIR" || exit 1

rm -f annot.txt annot_original.txt corrected_annot.txt tree.txt outtree.txt
rm -f outimg.png outimg.pdf
rm -f outimg.cladogram_top*.png outimg.cladogram_top*.pdf

conda deactivate
conda activate py2

python "$MTDIR/Tools/export2graphlan/export2graphlan.py" \
    -i "$GRAPHLAN_INPUT" \
    -a annot.txt \
    -t tree.txt \
    --most_abundant "$GRAPHLAN_TOP" \
    --least_biomarkers "$GRAPHLAN_LEAST_BIOMARKERS" \
    --annotations "$GRAPHLAN_LEVELS" \
    --external_annotations "$GRAPHLAN_EXTERNAL_LEVELS" \
    --background_levels "$GRAPHLAN_BACKGROUND_LEVELS" \
    --internal_levels \
    --max_clade_size "$GRAPHLAN_MAX_CLADE_SIZE" \
    --abundance_threshold "$GRAPHLAN_ABUNDANCE_THRESHOLD"

if [[ ! -s tree.txt || ! -s annot.txt ]]; then
    echo "${r}[ERROR] export2graphlan did not create tree.txt/annot.txt.${w}"
    echo "Check input:"
    echo "  $GRAPHLAN_INPUT"
    exit 1
fi

conda deactivate
conda activate MTD

# ------------------------------------------------------------
# Continue microbiome DEG / annotation / diversity preprocessing
# ------------------------------------------------------------

cd "$outputdr/temp" || exit 1

echo "${g}DEG & Annotation & Plots & Diversity & Preprocess for Microbiome ${w}"

if [[ "$NO_COMPARISON" == "1" ]]; then
    echo "${y}[INFO] Exploratory mode: skipping microbiome DEG/DESeq2 analysis.${w}"
    echo "${y}[INFO] Microbiome profiling outputs were already generated: Bracken, BIOM, MPA, Krona and GraPhlAn.${w}"
else
    conda deactivate
    conda activate R412

    deg_args=(
        "$outputdr/bracken_species_all"
        "$samplesheet_file"
        "$hostid"
        "$MTDIR/HostSpecies.csv"
    )

    if [[ -n "${metadata:-}" ]]; then
        deg_args+=( "$metadata" )
    fi

    Rscript "$MTDIR/DEG_Anno_Plot.R" "${deg_args[@]}"

    conda deactivate
    conda activate MTD
fi

cd "$outputdr/temp" || exit 1

mkdir -p bracken_raw_results

# Save raw combined Bracken tables if they are still present
mv ../bracken_*_all bracken_raw_results/ 2>/dev/null || true

# ------------------------------------------------------------
# Exploratory-only taxonomic figures
# ------------------------------------------------------------
# This uses:
#   $outputdr/temp/bracken_raw_results/bracken_species_all
#
# It runs only when NO_COMPARISON=1.
# ------------------------------------------------------------

run_exploratory_taxonomic_figures

# ------------------------------------------------------------
# Annotate and render GraPhlAn
# ------------------------------------------------------------

cd "$GRAPHLAN_DIR" || exit 1

echo "${g}Applying a fix for both tree.txt and annot.txt${w}"

python "$MTDIR/Tools/graphlan/verify_and_correct_annotations.py" \
    tree.txt \
    annot.txt \
    corrected_annot.txt

if [[ ! -s corrected_annot.txt ]]; then
    echo "${r}[ERROR] corrected_annot.txt was not created.${w}"
    exit 1
fi

mv annot.txt annot_original.txt
mv corrected_annot.txt annot.txt

echo "${g}Adding cladogram-like GraPhlAn style settings${w}"

{
    printf "title\tMicrobiome cladogram - top %s taxa\n" "$GRAPHLAN_TOP"
    printf "title_font_size\t13\n"

    # Figure geometry
    printf "total_plotted_degrees\t320\n"
    printf "start_rotation\t90\n"

    # Legends
    printf "annotation_legend_font_size\t7\n"
    printf "class_legend_font_size\t7\n"

    # Tree style
    printf "branch_thickness\t0.75\n"
    printf "clade_marker_size\t30\n"
    printf "clade_marker_edge_width\t0.45\n"
} >> annot.txt

python "$MTDIR/Tools/graphlan/graphlan_annotate.py" \
    --annot annot.txt \
    tree.txt \
    outtree.txt

if [[ ! -s outtree.txt ]]; then
    echo "${r}[ERROR] outtree.txt was not created.${w}"
    exit 1
fi

python "$MTDIR/Tools/graphlan/graphlan.py" \
    --dpi "$GRAPHLAN_DPI" \
    --size "$GRAPHLAN_SIZE" \
    outtree.txt \
    "outimg.cladogram_species_top${GRAPHLAN_TOP}.png"

python "$MTDIR/Tools/graphlan/graphlan.py" \
    --size "$GRAPHLAN_SIZE" \
    outtree.txt \
    "outimg.cladogram_species_top${GRAPHLAN_TOP}.pdf"

if [[ ! -s "outimg.cladogram_species_top${GRAPHLAN_TOP}.png" ]]; then
    echo "${r}[ERROR] GraPhlAn PNG was not created.${w}"
    exit 1
fi

if [[ ! -s "outimg.cladogram_species_top${GRAPHLAN_TOP}.pdf" ]]; then
    echo "${r}[ERROR] GraPhlAn PDF was not created.${w}"
    exit 1
fi

# Keep old output names too, so the rest of the pipeline/user habits do not break
cp "outimg.cladogram_species_top${GRAPHLAN_TOP}.png" outimg.png
cp "outimg.cladogram_species_top${GRAPHLAN_TOP}.pdf" outimg.pdf

echo "${g}[OK] Cladogram-like GraPhlAn outputs:${w}"
echo "  $outputdr/graphlan/outimg.cladogram_top${GRAPHLAN_TOP}.png"
echo "  $outputdr/graphlan/outimg.cladogram_top${GRAPHLAN_TOP}.pdf"
echo "  $outputdr/graphlan/outimg.png"
echo "  $outputdr/graphlan/outimg.pdf"

cd "$outputdr/temp" || exit 1

echo "${g}MTD running  progress:"
echo '>>>>>>>>>>>         [55%]'

echo "HUMAnN3${w}"

# ------------------------------------------------------------
# Prepare universal HUMAnN inputs
# ------------------------------------------------------------
# HUMAnN consumes one FASTQ per sample.
#
# SE:
#   copy SAMPLE_non-host.fq to SAMPLE.fq
#
# PE:
#   concatenate SAMPLE_non-host_1.fq and
#   SAMPLE_non-host_2.fq into SAMPLE.fq
#
# HUMAnN treats the reads independently during functional
# profiling, so paired-end mates are concatenated rather than
# passed as a synchronized pair.
# ------------------------------------------------------------

HUMANN_WORK_DIR="$PIPELINE_TEMP_DIR/HUMAnN_output"
HUMANN_INPUT_DIR="$HUMANN_WORK_DIR/input"
HUMANN_RESULTS_DIR="$HUMANN_WORK_DIR/hmn_output"
HUMANN_INPUT_MANIFEST="$HUMANN_WORK_DIR/humann_input_manifest.tsv"

mkdir -p \
    "$HUMANN_WORK_DIR" \
    "$HUMANN_INPUT_DIR" \
    "$HUMANN_RESULTS_DIR"

printf 'sample\tlayout\tsource_read1\tsource_read2\thumann_input\tpreparation\n' \
    > "$HUMANN_INPUT_MANIFEST"

echo "============================================================"
echo "[HUMAnN] Preparing functional profiling inputs"
echo "Read layout: $READ_LAYOUT"
echo "Input directory:"
echo "  $HUMANN_INPUT_DIR"
echo "Results directory:"
echo "  $HUMANN_RESULTS_DIR"
echo "============================================================"

for i in $lsn; do
    set_microbiome_fastq_paths "$i"

    humann_input="$HUMANN_INPUT_DIR/${i}.fq"
    source_read1="$PIPELINE_TEMP_DIR/$FINAL_NONHOST_R1"

    rm -f -- "$humann_input"

    echo
    echo "------------------------------------------------------------"
    echo "[HUMAnN INPUT] Sample: $i"
    echo "Layout: $READ_LAYOUT"

    if [[ "$READ_LAYOUT" == "pe" ]]; then
        source_read2="$PIPELINE_TEMP_DIR/$FINAL_NONHOST_R2"

        require_file \
            "$source_read1" \
            "Final non-host R1 for HUMAnN sample $i"

        require_file \
            "$source_read2" \
            "Final non-host R2 for HUMAnN sample $i"

        echo "Source R1: $source_read1"
        echo "Source R2: $source_read2"
        echo "Output:    $humann_input"
        echo "Method:    concatenate R1 followed by R2"

        if ! cat -- \
            "$source_read1" \
            "$source_read2" \
            > "$humann_input"
        then
            die "Failed to concatenate paired-end HUMAnN input for sample: $i"
        fi

        humann_preparation="concatenated_R1_then_R2"
    else
        source_read2="-"

        require_file \
            "$source_read1" \
            "Final non-host SE input for HUMAnN sample $i"

        echo "Source: $source_read1"
        echo "Output: $humann_input"
        echo "Method: copy single-end FASTQ"

        if ! cp -- \
            "$source_read1" \
            "$humann_input"
        then
            die "Failed to copy single-end HUMAnN input for sample: $i"
        fi

        humann_preparation="copied_single_end"
    fi

    require_file \
        "$humann_input" \
        "Prepared HUMAnN FASTQ for sample $i"

    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$i" \
        "$READ_LAYOUT" \
        "$source_read1" \
        "$source_read2" \
        "$humann_input" \
        "$humann_preparation" \
        >> "$HUMANN_INPUT_MANIFEST"

    echo "${g}[OK] HUMAnN input prepared:${w}"
    echo "  $humann_input"
done

echo
echo "${g}============================================"
echo "[OK] HUMAnN input preparation completed"
echo "Manifest:"
echo "  $HUMANN_INPUT_MANIFEST"
echo "${g}============================================${w}"

column -s $'\t' -t "$HUMANN_INPUT_MANIFEST" 2>/dev/null || \
    cat "$HUMANN_INPUT_MANIFEST"


# ------------------------------------------------------------
# Run HUMAnN once per sample
# ------------------------------------------------------------

echo "${g}Run HUMAnN3${w}"

for i in $lsn; do
    humann_input="$HUMANN_INPUT_DIR/${i}.fq"

    require_file \
        "$humann_input" \
        "HUMAnN input for sample $i"

    echo "============================================================"
    echo "[HUMAnN] Sample: $i"
    echo "Layout source: $READ_LAYOUT"
    echo "Input: $humann_input"
    echo "Output: $HUMANN_RESULTS_DIR"
    echo "Threads: $threads"
    echo "============================================================"

    if ! humann \
        --input "$humann_input" \
        --output "$HUMANN_RESULTS_DIR" \
        --threads "$threads" \
        --verbose
    then
        die "HUMAnN failed for sample: $i"
    fi

    humann_genefamilies="$HUMANN_RESULTS_DIR/${i}_genefamilies.tsv"
    humann_pathabundance="$HUMANN_RESULTS_DIR/${i}_pathabundance.tsv"

    require_file \
        "$humann_genefamilies" \
        "HUMAnN gene families output for sample $i"

    require_file \
        "$humann_pathabundance" \
        "HUMAnN pathway abundance output for sample $i"

    echo "${g}[OK] HUMAnN completed for sample:${w} $i"
done

cd "$HUMANN_WORK_DIR" || \
    die "Could not enter HUMAnN working directory: $HUMANN_WORK_DIR"

echo "${g}>>>>>>>>>>>>        [60%]"

echo "Join all gene family and pathway abudance files${w}"
humann_join_tables -i hmn_output/ -o humann_pathabundance.tsv --file_name pathabundance
humann_join_tables -i hmn_output/ -o humann_genefamilies.tsv --file_name genefamilies

# #Normalizing RPKs to CPM
# humann_renorm_table --input humann_pathabundance.tsv --output humann_pathabundance_cpm.tsv --units cpm --update-snames
# humann_renorm_table --input humann_genefamilies.tsv --output humann_genefamilies_cpm.tsv --units cpm --update-snames

echo "${g}Normalizing RPKs to "relab" (relative abundance)${w}"
humann_renorm_table --input humann_pathabundance.tsv --output humann_pathabundance_relab.tsv --units relab --update-snames
humann_renorm_table --input humann_genefamilies.tsv --output humann_genefamilies_relab.tsv --units relab --update-snames

echo "${g}Generate stratified tables; This utility will split a table into two files (one stratified and one unstratified). ${w}"
humann_split_stratified_table --input humann_pathabundance_relab.tsv --output ./
humann_split_stratified_table --input humann_genefamilies_relab.tsv --output ./
    echo "${g}Stratify unnormalized table (for Deseq2)${w}"
    humann_split_stratified_table --input humann_pathabundance.tsv --output ./
    humann_split_stratified_table --input humann_genefamilies.tsv --output ./

echo "${g}Regroup gene familites table into KEGG orthologs and GO terms${w}"
humann_regroup_table --input humann_genefamilies_relab_stratified.tsv --groups uniref90_ko --output humann_genefamilies_relAbundance_kegg.tsv
humann_regroup_table --input humann_genefamilies_relab_stratified.tsv --groups uniref90_go --output humann_genefamilies_relAbundance_go.tsv
    echo "${g}Regroup unnormalized table (for Deseq2${w}"
    humann_regroup_table --input humann_genefamilies_stratified.tsv --groups uniref90_ko --output humann_genefamilies_Abundance_kegg.tsv
    humann_regroup_table --input humann_genefamilies_stratified.tsv --groups uniref90_go --output humann_genefamilies_Abundance_go.tsv

echo "${g}Translate KEGG and GO ID to human readable terms${w}"
conda deactivate
conda activate R412
#Rscript $MTDIR/humann_ID_translation.R \
Rscript $MTDIR/humann_ID_translation_adjusted.R $outputdr/temp/HUMAnN_output/humann_genefamilies_relAbundance_kegg.tsv $outputdr/temp/HUMAnN_output/humann_genefamilies_relAbundance_go.tsv $MTDIR
    # Tranlate unnormalized table (for Deseq2)
#    Rscript $MTDIR/humann_ID_translation.R \
Rscript $MTDIR/humann_ID_translation_adjusted.R $outputdr/temp/HUMAnN_output/humann_genefamilies_Abundance_kegg.tsv $outputdr/temp/HUMAnN_output/humann_genefamilies_Abundance_go.tsv $MTDIR
conda deactivate
conda activate MTD

#Cleaning up file structure
mkdir -p $outputdr/hmn_pathway_abundance_files
mkdir -p $outputdr/hmn_genefamily_abundance_files
mv *pathabundance* $outputdr/hmn_pathway_abundance_files/
mv *genefamilies* $outputdr/hmn_genefamily_abundance_files/

# #Translate KEGG and GO ID to human readable terms
# Rscript $MTDIR/humann_ID_translation.R $outputdr/hmn_genefamily_abundance_files/humann_genefamilies_Abundance_kegg.tsv \
#     $outputdr/hmn_genefamily_abundance_files/humann_genefamilies_Abundance_go.tsv

echo "${g}DEG & Annotation & Plots & Diversity & Preprocess${w}"
cd $outputdr/hmn_genefamily_abundance_files
conda deactivate
conda activate R412
if [[ "$NO_COMPARISON" == "1" ]]; then
    echo "${y}[INFO] Exploratory mode: skipping HUMAnN DEG analysis.${w}"
else
    Rscript "$MTDIR/DEG_Anno_Plot.R" \
        "$outputdr/hmn_genefamily_abundance_files/humann_genefamilies_Abundance_kegg_translated.tsv" \
        "$samplesheet_file"

    Rscript "$MTDIR/DEG_Anno_Plot.R" \
        "$outputdr/hmn_genefamily_abundance_files/humann_genefamilies_Abundance_go_translated.tsv" \
        "$samplesheet_file"
fi

conda deactivate
conda activate MTD

#humann_barplot
# humann_barplot --input $outputdr/hmn_pathway_abundance_files/humann_pathabundance_cpm_stratified.tsv \
#     --focal-metadatum Group --last-metadatum Group \
#     --focal-feature PWY-3781 \
#     --output $outputdr/hmn_pathway_abundance_files/humann_pathabundance_barplot.png
# humann_barplot --input $outputdr/hmn_genefamily_abundance_files/humann_genefamilies_cpm_stratified.tsv \
#     --output $outputdr/hmn_genefamily_abundance_files/humann_genefamilies_barplot.png

echo "${g}MTD running  progress:"
echo '>>>>>>>>>>>>>       [65%]'

echo "Starting to process the host reads...${w}"

# ------------------------------------------------------------
# Host alignment
# ------------------------------------------------------------
# Current MTD behavior is preserved:
#
#   Magic-BLAST:
#     aligns reads classified as host by Kraken2
#
#   HISAT2:
#     aligns the complete prepared FASTQ input
#
# PE handling:
#   Magic-BLAST -> -query R1 -query_mate R2
#   HISAT2      -> -1 R1 -2 R2
# ------------------------------------------------------------

cd "$PIPELINE_TEMP_DIR" || \
    die "Could not enter pipeline temporary directory: $PIPELINE_TEMP_DIR"

host_sam_files=()

if [[ "$blast" == "blast" ]]; then
    echo "${g}Magic-BLAST${w}"

    for i in $lsn; do
        set_host_kraken_fastq_paths "$i"

        sam_file="${i}.sam"
        host_sam_files+=("$sam_file")

        rm -f -- "$sam_file"

        echo "============================================================"
        echo "[MAGIC-BLAST] Sample: $i"
        echo "Layout: $READ_LAYOUT"
        echo "Database: $DB_blast"
        echo "Output SAM: $sam_file"
        echo "Threads: $threads"

        magicblast_args=(
            -query "$HOST_KRAKEN_R1"
            -db "$DB_blast"
            -infmt fastq
            -out "$sam_file"
            -num_threads "$threads"
        )

        if [[ "$READ_LAYOUT" == "pe" ]]; then
            require_file \
                "$HOST_KRAKEN_R1" \
                "Kraken2 host-classified R1 for Magic-BLAST sample $i"

            require_file \
                "$HOST_KRAKEN_R2" \
                "Kraken2 host-classified R2 for Magic-BLAST sample $i"

            echo "Input R1: $HOST_KRAKEN_R1"
            echo "Input R2: $HOST_KRAKEN_R2"

            magicblast_args+=(
                -query_mate "$HOST_KRAKEN_R2"
            )
        else
            require_file \
                "$HOST_KRAKEN_R1" \
                "Kraken2 host-classified SE FASTQ for Magic-BLAST sample $i"

            echo "Input: $HOST_KRAKEN_R1"
        fi

        echo "============================================================"

        if ! magicblast "${magicblast_args[@]}"; then
            die "Magic-BLAST failed for sample: $i"
        fi

        require_file \
            "$sam_file" \
            "Magic-BLAST SAM output for sample $i"

        echo "${g}[OK] Magic-BLAST completed for sample:${w} $i"
    done

else
    echo "${g}HISAT2 alignment${w}"

    for i in $lsn; do
        load_prepared_fastq_manifest_row "$i"

        sam_file="${i}.sam"
        hisat2_summary="${i}_hisat2_summary.txt"

        host_sam_files+=("$sam_file")

        rm -f -- "$sam_file" "$hisat2_summary"

        echo "============================================================"
        echo "[HISAT2] Sample: $i"
        echo "Layout: $PREPARED_LAYOUT"
        echo "Database: $DB_hisat2"
        echo "Output SAM: $sam_file"
        echo "Summary: $hisat2_summary"
        echo "Threads: $threads"

        hisat2_args=(
            -p "$threads"
            -q
            -x "$DB_hisat2"
            --summary-file "$hisat2_summary"
            -S "$sam_file"
        )

        if [[ "$PREPARED_LAYOUT" == "pe" ]]; then
            echo "Input R1: $PREPARED_READ1"
            echo "Input R2: $PREPARED_READ2"

            hisat2_args+=(
                -1 "$PREPARED_READ1"
                -2 "$PREPARED_READ2"
            )
        else
            echo "Input: $PREPARED_READ1"

            hisat2_args+=(
                -U "$PREPARED_READ1"
            )
        fi

        echo "============================================================"

        if ! hisat2 "${hisat2_args[@]}"; then
            die "HISAT2 alignment failed for sample: $i"
        fi

        require_file \
            "$sam_file" \
            "HISAT2 SAM output for sample $i"

        require_file \
            "$hisat2_summary" \
            "HISAT2 summary for sample $i"

        echo "${g}[OK] HISAT2 completed for sample:${w} $i"
    done
fi


# ------------------------------------------------------------
# Validate host SAM collection
# ------------------------------------------------------------

if [[ "${#host_sam_files[@]}" -eq 0 ]]; then
    die "No host SAM files were registered for featureCounts."
fi

for sam_file in "${host_sam_files[@]}"; do
    require_file \
        "$sam_file" \
        "Host SAM input for featureCounts"
done


# ------------------------------------------------------------
# featureCounts
# ------------------------------------------------------------
# SE:
#   count reads
#
# PE:
#   declare paired-end input with -p
#   explicitly count fragments/read pairs with --countReadPairs
#   when supported by the installed featureCounts version
# ------------------------------------------------------------

echo "${g}featureCounts${w}"

featurecounts_args=(
    -T "$threads"
    -a "$gtf"
    -o "$outputdr/host_counts.txt"
)

if [[ "$READ_LAYOUT" == "pe" ]]; then
    featurecounts_args+=(
        -p
    )

    # Newer featureCounts versions require --countReadPairs
    # explicitly. In older versions, -p already implied fragment
    # counting, and --countReadPairs may not exist.
    if featureCounts -h 2>&1 |
        grep -q -- '--countReadPairs'
    then
        featurecounts_args+=(
            --countReadPairs
        )

        echo "[INFO] featureCounts paired-end mode:"
        echo "  -p --countReadPairs"
    else
        echo "${y}[WARNING] Installed featureCounts does not advertise --countReadPairs.${w}"
        echo "[WARNING] Using legacy paired-end behavior with -p only."
    fi
else
    echo "[INFO] featureCounts single-end mode."
fi

featurecounts_args+=(
    "${host_sam_files[@]}"
)

echo "============================================================"
echo "[FEATURECOUNTS]"
echo "Layout: $READ_LAYOUT"
echo "Annotation: $gtf"
echo "Output: $outputdr/host_counts.txt"
echo "SAM inputs: ${#host_sam_files[@]}"
printf '  %s\n' "${host_sam_files[@]}"
echo "============================================================"

if ! featureCounts "${featurecounts_args[@]}"; then
    die "featureCounts failed for host alignments."
fi

require_file \
    "$outputdr/host_counts.txt" \
    "Host featureCounts matrix"

require_file \
    "$outputdr/host_counts.txt.summary" \
    "Host featureCounts assignment summary"

echo "${g}[OK] featureCounts completed.${w}"
echo "Count matrix:"
echo "  $outputdr/host_counts.txt"
echo "Assignment summary:"
echo "  $outputdr/host_counts.txt.summary"


for i in $lsn; do
    sam_file="${i}.sam"
    bam_file="${i}.bam"
    sorted_bam="${i}.sorted.bam"

    require_file \
        "$sam_file" \
        "SAM input for samtools sample $i"

    echo "============================================================"
    echo "[SAMTOOLS] Sample: $i"
    echo "Input SAM: $sam_file"
    echo "Output BAM: $bam_file"
    echo "Sorted BAM: $sorted_bam"
    echo "============================================================"

    if ! samtools view \
        -@ "$threads" \
        -bS "$sam_file" \
        -o "$bam_file"
    then
        die "samtools view failed for sample: $i"
    fi

    if ! samtools sort \
        -@ "$threads" \
        -o "$sorted_bam" \
        "$bam_file"
    then
        die "samtools sort failed for sample: $i"
    fi

    if ! samtools index \
        -@ "$threads" \
        "$sorted_bam"
    then
        die "samtools index failed for sample: $i"
    fi

    require_file \
        "$sorted_bam" \
        "Sorted BAM for sample $i"

    require_file \
        "${sorted_bam}.bai" \
        "BAM index for sample $i"

    rm -f -- "$bam_file"
done
#Comando abaixo [e o mesmo acima, mas em uma unica linha
#for i in $lsn; do samtools view -bS $i.sam > $i.bam -@ $threads && samtools sort $i.bam -o $i.sorted.bam -@ $threads && samtools index $i.sorted.bam -@ $threads; done

mkdir -p BAM
mv *.sorted.bam *.sorted.bam.bai BAM/

cd $outputdr
# trim the featureCounts output(host_counts.txt) for downstream analysis
echo "${g}Delete the first line/row of a file then trim the sample name${w}"
sed '1d; 2 s/\.sam//g' host_counts.txt > tmpfile; mv tmpfile host_counts.txt

echo "${g}DEG & Annotation & Plots & preprocess for host${w}"
conda deactivate
conda activate R412
cd $outputdr
echo "${w}"
echo $MTDIR
echo $outputdr
echo $inputdr
echo $hostid 
echo $metadata
#### BEGIN FUNCTION: prepare_gene_id_cache_from_gtf ####
prepare_gene_id_cache_from_gtf() {
    local outputdr="$1"
    local gtf_file="$2"

    local cache_dir="${outputdr}/Host_DEG"
    local cache_file="${cache_dir}/gene_ID_cache.csv"

    local rscript_bin
    rscript_bin="$(command -v Rscript)"

    echo "------------------------------------------------------------"
    echo "Preparing offline host gene annotation cache from GTF"
    echo "Rscript: $rscript_bin"
    echo "GTF: $gtf_file"
    echo "Cache: $cache_file"
    echo "------------------------------------------------------------"

    mkdir -p "$cache_dir"

    if [[ -s "$cache_file" ]]; then
        echo "[INFO] Existing gene_ID_cache.csv found. Validating..."

        "$rscript_bin" - <<RS
cache_file <- "$cache_file"

gene_ID <- read.csv(cache_file, header = TRUE, check.names = FALSE)

required <- c(
  "gene_name",
  "ensembl_gene_id",
  "chromosome_name",
  "start_position",
  "end_position",
  "strand",
  "gene_biotype",
  "description",
  "gene_length"
)

missing <- setdiff(required, names(gene_ID))

if (length(missing) > 0) {
  stop("Cache is missing columns: ", paste(missing, collapse = ", "))
}

if (nrow(gene_ID) == 0) {
  stop("Cache has zero rows.")
}

cat("[OK] Existing cache is valid. Genes:", nrow(gene_ID), "\\n")
RS

        if [[ "$?" -eq 0 ]]; then
            return 0
        else
            echo "[WARN] Existing cache is invalid. Rebuilding from GTF..."
            rm -f "$cache_file"
        fi
    fi

    if [[ ! -s "$gtf_file" ]]; then
        echo "[ERROR] GTF file not found or empty: $gtf_file"
        return 1
    fi

    if [[ -z "$rscript_bin" || ! -x "$rscript_bin" ]]; then
        echo "[ERROR] Could not find Rscript in PATH."
        return 1
    fi

    "$rscript_bin" - <<RS
gtf_file <- "$gtf_file"
cache_file <- "$cache_file"

cat("[INFO] Reading GTF:", gtf_file, "\\n")

open_gtf <- function(path) {
  if (grepl("\\\\.gz$", path)) {
    gzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }
}

extract_attr <- function(x, key) {
  pattern <- paste0(key, ' "([^"]+)"')
  hit <- regexpr(pattern, x, perl = TRUE)
  out <- rep(NA_character_, length(x))
  ok <- hit > 0
  out[ok] <- sub(pattern, "\\\\1", regmatches(x, hit)[ok], perl = TRUE)
  out
}

con <- open_gtf(gtf_file)
gtf <- read.delim(
  con,
  header = FALSE,
  sep = "\\t",
  comment.char = "#",
  quote = "",
  stringsAsFactors = FALSE
)
close(con)

if (ncol(gtf) < 9) {
  stop("GTF has fewer than 9 columns. Is this a valid GTF?")
}

names(gtf)[1:9] <- c(
  "seqname", "source", "feature", "start", "end",
  "score", "strand", "frame", "attribute"
)

genes <- gtf[gtf\$feature == "gene", , drop = FALSE]

if (nrow(genes) == 0) {
  cat("[WARN] No 'gene' features found. Falling back to exon-derived gene ranges.\\n")

  exons <- gtf[gtf\$feature == "exon", , drop = FALSE]

  if (nrow(exons) == 0) {
    stop("No gene or exon features found in GTF.")
  }

  exons\$gene_id <- extract_attr(exons\$attribute, "gene_id")
  exons\$gene_name <- extract_attr(exons\$attribute, "gene_name")
  exons\$gene_biotype <- extract_attr(exons\$attribute, "gene_biotype")

  if (all(is.na(exons\$gene_biotype))) {
    exons\$gene_biotype <- extract_attr(exons\$attribute, "gene_type")
  }

  exons <- exons[!is.na(exons\$gene_id), , drop = FALSE]

  genes <- aggregate(
    cbind(start, end) ~ gene_id + seqname + strand,
    data = exons,
    FUN = function(z) c(min = min(z), max = max(z))
  )

  genes\$start_position <- genes\$start[, "min"]
  genes\$end_position <- genes\$end[, "max"]
  genes\$start <- NULL
  genes\$end <- NULL

  meta <- exons[!duplicated(exons\$gene_id), c("gene_id", "gene_name", "gene_biotype"), drop = FALSE]
  genes <- merge(genes, meta, by = "gene_id", all.x = TRUE)

} else {
  genes\$gene_id <- extract_attr(genes\$attribute, "gene_id")
  genes\$gene_name <- extract_attr(genes\$attribute, "gene_name")
  genes\$gene_biotype <- extract_attr(genes\$attribute, "gene_biotype")

  if (all(is.na(genes\$gene_biotype))) {
    genes\$gene_biotype <- extract_attr(genes\$attribute, "gene_type")
  }

  genes\$start_position <- genes\$start
  genes\$end_position <- genes\$end
}

genes <- genes[!is.na(genes\$gene_id), , drop = FALSE]
genes <- genes[!duplicated(genes\$gene_id), , drop = FALSE]

genes\$gene_name[is.na(genes\$gene_name) | genes\$gene_name == ""] <- genes\$gene_id[is.na(genes\$gene_name) | genes\$gene_name == ""]
genes\$gene_biotype[is.na(genes\$gene_biotype) | genes\$gene_biotype == ""] <- "unknown"

gene_ID <- data.frame(
  gene_name = genes\$gene_name,
  ensembl_gene_id = genes\$gene_id,
  chromosome_name = genes\$seqname,
  start_position = genes\$start_position,
  end_position = genes\$end_position,
  strand = genes\$strand,
  gene_biotype = genes\$gene_biotype,
  description = genes\$gene_name,
  gene_length = abs(as.numeric(genes\$end_position) - as.numeric(genes\$start_position)) + 1,
  stringsAsFactors = FALSE
)

gene_ID <- gene_ID[!is.na(gene_ID\$ensembl_gene_id), , drop = FALSE]
gene_ID <- gene_ID[!duplicated(gene_ID\$ensembl_gene_id), , drop = FALSE]

write.csv(gene_ID, cache_file, row.names = FALSE, quote = TRUE)

cat("[OK] Offline cache created:", cache_file, "\\n")
cat("[OK] Genes:", nrow(gene_ID), "\\n")
cat("[INFO] First rows:\\n")
print(utils::head(gene_ID, 3))
RS

    local status=$?

    if [[ "$status" -ne 0 || ! -s "$cache_file" ]]; then
        echo "[ERROR] Failed to create gene_ID_cache.csv from GTF."
        return 1
    fi

    echo "[OK] gene_ID_cache.csv is ready."
    return 0
}
#### END FUNCTION: prepare_gene_id_cache_from_gtf ####

#### BEGIN CALL: host gene annotation cache update ####

prepare_gene_id_cache_from_gtf "$outputdr" "$gtf" || {
    echo "${r}[ERROR] Could not prepare gene_ID_cache.csv from local GTF.${w}"
    exit 1
}

# Opcional: só tenta BioMart depois que já existe um cache local seguro.
# Se a internet/DNS falhar, o pipeline continua usando o cache do GTF.
if declare -F update_host_gene_cache_online >/dev/null; then
    update_host_gene_cache_online "$outputdr" "$hostid" "$MTDIR" || true
else
    echo "[INFO] update_host_gene_cache_online function not defined; skipping online cache update."
fi

#### END CALL: host gene annotation cache update ####

if [[ "$NO_COMPARISON" == "1" ]]; then
    echo "${y}[INFO] Exploratory mode: skipping host DEG analysis.${w}"
    echo "${y}[INFO] Host count matrix was still generated by featureCounts:${w}"
    echo "  $outputdr/host_counts.txt"

    mkdir -p "$outputdr/Host_DEG"

    cp "$outputdr/host_counts.txt" "$outputdr/Host_DEG/host_counts_featureCounts_matrix.txt"

    run_exploratory_host_expression_qc

else
    host_deg_args=(
        "$outputdr/host_counts.txt"
        "$samplesheet_file"
        "$hostid"
        "$MTDIR/HostSpecies.csv"
    )

    if [[ -n "${metadata:-}" ]]; then
        host_deg_args+=( "$metadata" )
    fi


    Rscript "$MTDIR/DEG_Anno_Plot.R" "${host_deg_args[@]}"

    # ------------------------------------------------------------
    # Extra DEG volcano plots with EnhancedVolcano
    # Runs EV.volcano.R in a selected conda environment.
    # Default: base, because EnhancedVolcano is already installed there.
    # Outputs are written inside each comparison folder.
    # ------------------------------------------------------------

    EV_VOLCANO_SCRIPT="$MTDIR/aux_scripts/EV/EV.volcano.R"
    EV_VOLCANO_LABEL_TOP="${EV_VOLCANO_LABEL_TOP:-25}"
    EV_VOLCANO_ENV="${EV_VOLCANO_ENV:-base}"

    run_ev_volcano_for_deg_folder() {
        local deg_label="$1"
        local deg_dir_main="$2"
        local csv_pattern="$3"

        if [[ ! -s "$EV_VOLCANO_SCRIPT" ]]; then
            echo "${y}[WARNING] EV volcano script not found. Skipping ${deg_label} EV volcano plots.${w}"
            echo "Expected:"
            echo "  $EV_VOLCANO_SCRIPT"
            return 0
        fi

        echo "${g}Generating extra EV volcano plots for ${deg_label}${w}"
        echo "[EV VOLCANO] Conda environment:"
        echo "  $EV_VOLCANO_ENV"

        shopt -s nullglob
        local de_csv_files=( "$deg_dir_main"/*/$csv_pattern )
        shopt -u nullglob

        if [[ "${#de_csv_files[@]}" -eq 0 ]]; then
            echo "${y}[WARNING] No ${deg_label} DEG comparison CSV files found for EV volcano plots.${w}"
            echo "Expected pattern:"
            echo "  $deg_dir_main/*/$csv_pattern"
            return 0
        fi

        for de_csv in "${de_csv_files[@]}"; do
            local de_dir
            local de_file
            local volcano_log

            de_dir="$(dirname "$de_csv")"
            de_file="$(basename "$de_csv")"
            volcano_log="$de_dir/${de_file%.csv}.EV.volcano.log"

            echo "------------------------------------------------------------"
            echo "[EV VOLCANO] Dataset:"
            echo "  $deg_label"
            echo "[EV VOLCANO] Input:"
            echo "  $de_csv"
            echo "[EV VOLCANO] Output directory:"
            echo "  $de_dir"
            echo "[EV VOLCANO] label_top:"
            echo "  $EV_VOLCANO_LABEL_TOP"
            echo "[EV VOLCANO] conda env:"
            echo "  $EV_VOLCANO_ENV"
            echo "------------------------------------------------------------"

            if ! (
                conda activate "$EV_VOLCANO_ENV" || exit 1

                echo "[INFO] Rscript used for EV volcano:"
                which Rscript
                Rscript -e 'cat("[INFO] EnhancedVolcano available: ", requireNamespace("EnhancedVolcano", quietly=TRUE), "\n")'

                cd "$de_dir" || exit 1

                Rscript "$EV_VOLCANO_SCRIPT" \
                    --de_results "$de_file" \
                    --label_top "$EV_VOLCANO_LABEL_TOP"
            ) > "$volcano_log" 2>&1
            then
                echo "${y}[WARNING] EV volcano failed for:${w}"
                echo "  $de_csv"
                echo "Log:"
                echo "  $volcano_log"
                tail -n 40 "$volcano_log" || true
            else
                echo "${g}[OK] EV volcano generated in:${w}"
                echo "  $de_dir"
            fi
        done
    }

    run_ev_volcano_for_deg_folder \
        "host DEG" \
        "$outputdr/Host_DEG" \
        "host_counts_*.csv"

    run_ev_volcano_for_deg_folder \
        "non-host DEG" \
        "$outputdr/Nonhost_DEG" \
        "bracken_species_all_*.csv"

    require_file "$outputdr/Host_DEG/host_counts_TPM.csv" "Host TPM matrix generated by DEG_Anno_Plot.R"
fi

if [[ "$NO_COMPARISON" == "1" ]]; then
    echo "${y}[INFO] Exploratory mode: skipping ssGSEA and HAllA association stages.${w}"
    echo "${y}[INFO] Reason: no experimental comparison groups were detected.${w}"
    echo "${y}[INFO] The pipeline will finish after generating profiling outputs.${w}"

    conda deactivate
    conda activate MTD

    echo "${g}"
    echo 'MTD running  progress:'
    echo '>>>>>>>>>>>>>>>>>>>>[100%]'
    echo "MTD exploratory run is finished"
    echo -e "${w}"
    exit 0
fi

echo "${g}MTD running  progress:"
echo '>>>>>>>>>>>>>>>     [75%]'

echo "ssGSEA${w}"

require_file "$outputdr/Host_DEG/host_counts_TPM.csv" "Host TPM matrix"

run_cmd Rscript "$MTDIR/gct_making.R" \
    "$outputdr/Host_DEG/host_counts_TPM.csv" \
    "$samplesheet_file"

require_file "$outputdr/ssGSEA/host.gct" "ssGSEA input GCT"

echo "${g}[INFO] ssGSEA GMT file:${w}"
echo "  $SSGSEA_GMT"
echo "${g}[INFO] ssGSEA GMT mode:${w}"
echo "  $SSGSEA_GMT_MODE"
write_methods_log

# ------------------------------------------------------------
# Select or build ssGSEA GMT
# ------------------------------------------------------------

SSGSEA_DEFAULT_GMT="$MTDIR/Tools/ssGSEA2.0/db/msigdb/c2.all.v7.5.1.symbols.gmt"

if [[ "$SSGSEA_GMT" == "auto" || "$SSGSEA_GMT" == "custom_eggnog" ]]; then
    echo "${g}[INFO] Building custom ssGSEA GMT from user-provided protein FASTA and annotation.${w}"
    echo "[INFO] This mode does not download NCBI datasets by taxid."
    echo "[INFO] TaxID: $hostid"
    echo "[INFO] host.gct: $outputdr/ssGSEA/host.gct"

    if [[ -z "${SSGSEA_PROTEIN_FASTA:-}" ]]; then
        die "--ssgsea-gmt auto requires --ssgsea-protein-fasta FILE"
    fi

    if [[ -z "${SSGSEA_ANNOTATION_GFF:-}" ]]; then
        die "--ssgsea-gmt auto requires --ssgsea-annotation-gff FILE"
    fi

    require_file "$SSGSEA_PROTEIN_FASTA" "ssGSEA protein FASTA"
    require_file "$SSGSEA_ANNOTATION_GFF" "ssGSEA GFF/GTF annotation"

    resolve_eggnog_db_dir

    CUSTOM_GMT_DIR="$outputdr/ssGSEA/custom_gmt_user_files"
    mkdir -p "$CUSTOM_GMT_DIR"

    MTD_BIN="$condapath/envs/MTD/bin"

    if [[ ! -d "$MTD_BIN" ]]; then
        die "MTD conda env bin folder not found: $MTD_BIN"
    fi

    AUTO_GMT_ORIGINAL_PATH="$PATH"
    AUTO_GMT_PATH="$PATH:$MTD_BIN"

    echo "${g}[INFO] Custom GMT tool resolution:${w}"

    missing_auto_gmt_tools=0

    for tool in Rscript python3 emapper.py diamond; do
        tool_path=$(PATH="$AUTO_GMT_PATH" command -v "$tool" 2>/dev/null || true)

        if [[ -z "$tool_path" ]]; then
            echo "${r}[ERROR] Required tool not found for custom ssGSEA GMT:${w} $tool"
            missing_auto_gmt_tools=1
        else
            echo "${g}[OK]${w} $tool -> $tool_path"
        fi
    done

    if [[ "$missing_auto_gmt_tools" == "1" ]]; then
        echo
        echo "${y}[INFO] Install missing tools inside the MTD environment:${w}"
        echo "  conda activate MTD"
        echo "  conda install -c conda-forge -c bioconda eggnog-mapper diamond -y"
        exit 1
    fi

    export PATH="$AUTO_GMT_PATH"

    conda deactivate
    conda activate MTD

    run_cmd bash "$MTDIR/aux_scripts/ssGSEA/make_custom_ssgsea_gmt_from_taxid_auto.sh" \
        --taxid "$hostid" \
        --host-gct "$outputdr/ssGSEA/host.gct" \
        --outdir "$CUSTOM_GMT_DIR" \
        --threads "$threads" \
        --eggnog-db "$EGGNOG_DB_DIR" \
        --protein-fasta "$SSGSEA_PROTEIN_FASTA" \
        --annotation-gff "$SSGSEA_ANNOTATION_GFF" \
        --force

    SSGSEA_GMT="$CUSTOM_GMT_DIR/custom_taxid_${hostid}_eggNOG_GO.gmt"
    SSGSEA_GMT_MODE="custom_eggNOG_GO_from_user_files"

    require_file "$SSGSEA_GMT" "Custom ssGSEA GMT"

    export PATH="$AUTO_GMT_ORIGINAL_PATH"

    conda deactivate
    conda activate R412

elif [[ "$SSGSEA_GMT" == "default" ]]; then
    SSGSEA_GMT="$SSGSEA_DEFAULT_GMT"
    SSGSEA_GMT_MODE="default_MSigDB_c2_symbols"
    require_file "$SSGSEA_GMT" "Default ssGSEA GMT"

else
    SSGSEA_GMT_MODE="custom_path"
    require_file "$SSGSEA_GMT" "Custom ssGSEA GMT"
fi

echo "${g}[INFO] ssGSEA GMT file:${w}"
echo "  $SSGSEA_GMT"
echo "${g}[INFO] ssGSEA GMT mode:${w}"
echo "  $SSGSEA_GMT_MODE"

write_methods_log

run_cmd Rscript "$MTDIR/Tools/ssGSEA2.0/ssgsea-cli.R" \
    -i "$outputdr/ssGSEA/host.gct" \
    -o "$outputdr/ssGSEA/ssgsea-results" \
    -d "$SSGSEA_GMT" \
    -y "$MTDIR/Tools/ssGSEA2.0/config.yaml" \
    -u "$threads"

require_file "$outputdr/ssGSEA/ssgsea-results-scores.gct" "ssGSEA scores"

run_cmd Rscript "$MTDIR/for_halla.R" \
    "$outputdr/ssGSEA/ssgsea-results-scores.gct" \
    "$samplesheet_file" \
    $metadata

require_file "$outputdr/halla/Microbiomes.txt" "HAllA microbiome input"
require_file "$outputdr/halla/Host_gene.txt" "HAllA host gene input"
require_file "$outputdr/halla/Host_score.txt" "HAllA host pathway input"

# ------------------------------------------------------------
# ssGSEA visualization
# ------------------------------------------------------------

RUN_SSGSEA_PLOTS="${RUN_SSGSEA_PLOTS:-1}"

# Optional GO-name corrected ssGSEA plots
# 1 = generate an additional set of plots using readable GO names
# 0 = keep only the original GO-ID plots
RUN_SSGSEA_GO_NAME_PLOTS="${RUN_SSGSEA_GO_NAME_PLOTS:-1}"

# Optional online QuickGO fallback for GO IDs not found in local GO.db
# 1 = use QuickGO when internet is available
# 0 = offline mode only; GO.db + built-in/manual corrections
RUN_SSGSEA_QUICKGO="${RUN_SSGSEA_QUICKGO:-1}"

# Optional QuickGO secondary-ID/replacement search
# 1 = try to resolve old/secondary GO IDs after direct lookup fails
# 0 = skip this slower fallback
RUN_SSGSEA_QUICKGO_SEARCH="${RUN_SSGSEA_QUICKGO_SEARCH:-1}"

if [[ "$RUN_SSGSEA_PLOTS" == "1" ]]; then
    echo "${g}Generating ssGSEA plots${w}"

    SSGSEA_PLOT_SCRIPT="$MTDIR/aux_scripts/ssGSEA/plot_ssgsea_results.R"
    SSGSEA_GO_RESOLVER_SCRIPT="$MTDIR/aux_scripts/ssGSEA/resolve_ssgsea_go_terms.py"
    SSGSEA_MANUAL_GO_MAP="$MTDIR/aux_scripts/ssGSEA/go_replacement_manual_map.tsv"

    require_file "$SSGSEA_PLOT_SCRIPT" "ssGSEA plotting script"

    # ------------------------------------------------------------
    # 1) Original ssGSEA plots using original feature names / GO IDs
    # ------------------------------------------------------------

    run_cmd Rscript "$SSGSEA_PLOT_SCRIPT" \
        --scores "$outputdr/ssGSEA/ssgsea-results-scores.gct" \
        --samplesheet "$samplesheet_file" \
        --outdir "$outputdr/ssGSEA/plots" \
        --top-var 50 \
        --top-diff 20 \
        --pca-top-var 500

    require_file "$outputdr/ssGSEA/plots/ssGSEA_top_variable_heatmap.png" "ssGSEA top variable heatmap"
    require_file "$outputdr/ssGSEA/plots/ssGSEA_PCA_samples.png" "ssGSEA PCA plot"
    require_file "$outputdr/ssGSEA/plots/ssGSEA_differential_scores.tsv" "ssGSEA differential score table"

    echo "${g}[OK] ssGSEA plots saved to:${w}"
    echo "  $outputdr/ssGSEA/plots"

    # ------------------------------------------------------------
    # 2) Optional corrected GO-name ssGSEA plots
    # ------------------------------------------------------------

    if [[ "$RUN_SSGSEA_GO_NAME_PLOTS" == "1" ]]; then
        echo "${g}Generating corrected GO-name ssGSEA plots${w}"

        require_file "$SSGSEA_GO_RESOLVER_SCRIPT" "ssGSEA GO term resolver script"

        SSGSEA_GO_RESOLUTION_DIR="$outputdr/ssGSEA/go_term_resolution"
        SSGSEA_GO_NAME_DIR="$outputdr/ssGSEA/plots_GO_names"
        SSGSEA_GO_NAME_GCT="$outputdr/ssGSEA/ssgsea-results-scores_GO_names_corrected.gct"
        SSGSEA_GO_NAME_MAP="$SSGSEA_GO_RESOLUTION_DIR/ssGSEA_GO_ID_to_name_map.tsv"
        SSGSEA_QUICKGO_CACHE="$SSGSEA_GO_RESOLUTION_DIR/quickgo_cache.tsv"

        mkdir -p "$SSGSEA_GO_RESOLUTION_DIR" "$SSGSEA_GO_NAME_DIR"

        if [[ "$RUN_SSGSEA_QUICKGO" == "1" ]]; then
            SSGSEA_QUICKGO_ARG="yes"
        else
            SSGSEA_QUICKGO_ARG="no"
        fi

        if [[ "$RUN_SSGSEA_QUICKGO_SEARCH" == "1" ]]; then
            SSGSEA_QUICKGO_SEARCH_ARG="yes"
        else
            SSGSEA_QUICKGO_SEARCH_ARG="no"
        fi

        echo "[INFO] ssGSEA GO-name resolver settings:"
        echo "  Resolver script:        $SSGSEA_GO_RESOLVER_SCRIPT"
        echo "  Input GCT:              $outputdr/ssGSEA/ssgsea-results-scores.gct"
        echo "  Corrected GCT:          $SSGSEA_GO_NAME_GCT"
        echo "  GO map:                 $SSGSEA_GO_NAME_MAP"
        echo "  QuickGO enabled:        $SSGSEA_QUICKGO_ARG"
        echo "  QuickGO search enabled: $SSGSEA_QUICKGO_SEARCH_ARG"
        echo "  QuickGO cache:          $SSGSEA_QUICKGO_CACHE"

        if [[ -s "$SSGSEA_MANUAL_GO_MAP" ]]; then
            echo "  Manual GO map:          $SSGSEA_MANUAL_GO_MAP"

            run_cmd python3 "$SSGSEA_GO_RESOLVER_SCRIPT" \
                --scores "$outputdr/ssGSEA/ssgsea-results-scores.gct" \
                --out-gct "$SSGSEA_GO_NAME_GCT" \
                --out-map "$SSGSEA_GO_NAME_MAP" \
                --manual-map "$SSGSEA_MANUAL_GO_MAP" \
                --quickgo "$SSGSEA_QUICKGO_ARG" \
                --quickgo-search "$SSGSEA_QUICKGO_SEARCH_ARG" \
                --quickgo-cache "$SSGSEA_QUICKGO_CACHE"
        else
            echo "  Manual GO map:          not found; using built-in resolver fallbacks only"

            run_cmd python3 "$SSGSEA_GO_RESOLVER_SCRIPT" \
                --scores "$outputdr/ssGSEA/ssgsea-results-scores.gct" \
                --out-gct "$SSGSEA_GO_NAME_GCT" \
                --out-map "$SSGSEA_GO_NAME_MAP" \
                --quickgo "$SSGSEA_QUICKGO_ARG" \
                --quickgo-search "$SSGSEA_QUICKGO_SEARCH_ARG" \
                --quickgo-cache "$SSGSEA_QUICKGO_CACHE"
        fi

        require_file "$SSGSEA_GO_NAME_GCT" "Corrected GO-name ssGSEA GCT"
        require_file "$SSGSEA_GO_NAME_MAP" "ssGSEA GO ID to name resolution map"

        run_cmd Rscript "$SSGSEA_PLOT_SCRIPT" \
            --scores "$SSGSEA_GO_NAME_GCT" \
            --samplesheet "$samplesheet_file" \
            --outdir "$SSGSEA_GO_NAME_DIR" \
            --top-var 50 \
            --top-diff 20 \
            --pca-top-var 500

        require_file "$SSGSEA_GO_NAME_DIR/ssGSEA_top_variable_heatmap.png" "Corrected GO-name ssGSEA top variable heatmap"
        require_file "$SSGSEA_GO_NAME_DIR/ssGSEA_PCA_samples.png" "Corrected GO-name ssGSEA PCA plot"
        require_file "$SSGSEA_GO_NAME_DIR/ssGSEA_differential_scores.tsv" "Corrected GO-name ssGSEA differential score table"

        echo "${g}[OK] Corrected GO-name ssGSEA plots saved to:${w}"
        echo "  $SSGSEA_GO_NAME_DIR"
        echo "${g}[OK] GO ID resolution map saved to:${w}"
        echo "  $SSGSEA_GO_NAME_MAP"

        echo "[INFO] GO resolution summary:"
        awk -F'\t' '
        NR > 1 {
            total++
            source[$5]++
        }
        END {
            print "  Total rows: " total
            for (s in source) {
                print "  " s ": " source[s]
            }
            if (total > 0) {
                unresolved = source["unresolved"] + 0
                printf "  Resolved percent: %.2f%%\n", ((total - unresolved) / total) * 100
            }
        }' "$SSGSEA_GO_NAME_MAP" || true
    else
        echo "${y}[INFO] RUN_SSGSEA_GO_NAME_PLOTS=0; skipping corrected GO-name ssGSEA plots.${w}"
    fi

else
    echo "${y}[INFO] RUN_SSGSEA_PLOTS=0; skipping ssGSEA plots.${w}"
fi
echo "${g}MTD running  progress:"
echo '>>>>>>>>>>>>>>>>    [80%]'
echo "MTD DEG analyses are done. Starting microbiome x host association analyses..."

echo "halla: association analysis${w}"

conda deactivate
conda activate halla0820

export PYTHONNOUSERSITE=1
unset PYTHONPATH
unset PYTHONHOME
export MPLBACKEND=Agg
export PYTHONWARNINGS="ignore"

HALLA_THREADS="${threads:-$(nproc)}"
export OMP_NUM_THREADS="$HALLA_THREADS"
export OPENBLAS_NUM_THREADS="$HALLA_THREADS"
export MKL_NUM_THREADS="$HALLA_THREADS"
export NUMEXPR_NUM_THREADS="$HALLA_THREADS"

RUN_EXTRA_PEARSON="${RUN_EXTRA_PEARSON:-1}"
RUN_FULL_HALLAGRAM="${RUN_FULL_HALLAGRAM:-0}"
HALLA_DIAGNOSTIC="${HALLA_DIAGNOSTIC:-0}"

echo "[INFO] HAllA threads: $HALLA_THREADS"
echo "[INFO] RUN_EXTRA_PEARSON: $RUN_EXTRA_PEARSON"
echo "[INFO] RUN_FULL_HALLAGRAM: $RUN_FULL_HALLAGRAM"
echo "[INFO] HALLA_DIAGNOSTIC: $HALLA_DIAGNOSTIC"

run_halla_safe() {
    local xfile="$1"
    local yfile="$2"
    local outdir="$3"
    local metric="$4"
    local xlabel="$5"
    local ylabel="$6"
    local logfile="${outdir}.halla.log"

    echo "============================================================"
    echo "[HALLA] X: $xfile"
    echo "[HALLA] Y: $yfile"
    echo "[HALLA] Output: $outdir"
    echo "[HALLA] Metric: $metric"
    echo "[HALLA] Threads: $HALLA_THREADS"
    echo "============================================================"

    if [[ ! -s "$xfile" ]]; then
        echo "[WARNING] Missing HAllA X input: $xfile"
        return 0
    fi

    if [[ ! -s "$yfile" ]]; then
        echo "[WARNING] Missing HAllA Y input: $yfile"
        return 0
    fi

    mkdir -p "$(dirname "$outdir")"

    if [[ -d "$outdir" ]]; then
        local backup="${outdir}.previous_$(date +%Y%m%d_%H%M%S)"
        echo "[INFO] Existing HAllA output folder found. Moving to: $backup"
        mv "$outdir" "$backup"
    fi

    local halla_diag_args=""
if [[ "${HALLA_DIAGNOSTIC:-0}" == "1" ]]; then
    halla_diag_args="--diagnostic_plot"
fi

halla -x "$xfile" -y "$yfile" -o "$outdir" --x_dataset_label "$xlabel" --y_dataset_label "$ylabel" $halla_diag_args -m "$metric" --num_threads "$HALLA_THREADS" > "$logfile" 2>&1
    local status=$?

    if [[ "$status" -ne 0 ]]; then
        echo "[WARNING] HAllA exited with status $status for: $outdir"
        echo "[WARNING] This often happens during report/hallagram plotting after the statistics were computed."
        echo "[WARNING] Log saved at: $logfile"

        mkdir -p "$outdir"

        {
            echo "HAllA finished with warning/error status: $status"
            echo "This may be caused by report/hallagram plotting, especially MatplotlibDeprecationWarning."
            echo "The MTD pipeline continued instead of stopping."
            echo "Log file:"
            echo "$logfile"
            echo
            echo "Relevant log lines:"
            grep -E "Number of significant|significant clusters|Traceback|MatplotlibDeprecationWarning|ERROR|WARNING" "$logfile" | tail -n 100
        } > "$outdir/HAllA_finished_with_warning.txt"

        echo "[INFO] Relevant HAllA log lines:"
        grep -E "Number of significant|significant clusters|Traceback|MatplotlibDeprecationWarning|ERROR|WARNING" "$logfile" | tail -n 40

        return 0
    fi

    echo "[OK] HAllA completed successfully: $outdir"
    echo "[OK] Log saved at: $logfile"
    return 0
}

run_hallagram_safe() {
    local indir="$1"
    local outfile="$2"
    local block_num="$3"
    local xlabel="$4"
    local ylabel="$5"
    local cbar="$6"
    local logfile="${outfile}.log"

    echo "============================================================"
    echo "[HALLAGRAM] Input: $indir"
    echo "[HALLAGRAM] Output: $outfile"
    echo "[HALLAGRAM] block_num: $block_num"
    echo "============================================================"

    if [[ ! -d "$indir" ]]; then
        echo "[WARNING] HAllA folder does not exist; skipping hallagram: $indir"
        return 0
    fi

    mkdir -p "$(dirname "$outfile")"

    hallagram -i "$indir" --cbar_label "$cbar" --x_dataset_label "$xlabel" --y_dataset_label "$ylabel" --output "$outfile" --block_num "$block_num" > "$logfile" 2>&1
    local status=$?

    if [[ "$status" -ne 0 ]]; then
        echo "[WARNING] hallagram failed for: $outfile"
        echo "[WARNING] Log saved at: $logfile"
        grep -E "Traceback|Error|WARNING|MatplotlibDeprecationWarning" "$logfile" | tail -n 40
        return 0
    fi

    echo "[OK] hallagram saved: $outfile"
    return 0
}

run_python_plot_safe() {
    local label="$1"
    local cmd="$2"
    local logfile="$3"

    echo "============================================================"
    echo "[PYTHON] $label"
    echo "============================================================"

    eval "$cmd" > "$logfile" 2>&1
    local status=$?

    if [[ "$status" -ne 0 ]]; then
        echo "[WARNING] $label failed; continuing."
        echo "[WARNING] Log saved at: $logfile"
        tail -n 40 "$logfile"
        return 0
    fi

    echo "[OK] $label completed."
    return 0
}

if [[ "$pdm" == "spearman" ]]; then
    pdm_name='Pairwise Spearman'
elif [[ "$pdm" == "pearson" ]]; then
    pdm_name='Pairwise Pearson'
elif [[ "$pdm" == "mi" ]]; then
    pdm_name='mi'
elif [[ "$pdm" == "nmi" ]]; then
    pdm_name='nmi'
elif [[ "$pdm" == "xicor" ]]; then
    pdm_name='xicor'
elif [[ "$pdm" == "dcor" ]]; then
    pdm_name='dcor'
else
    pdm_name="$pdm"
fi

echo "${g}Analyzing microbiome x host_genes associations...${w}"

run_halla_safe "$outputdr/halla/Microbiomes.txt" "$outputdr/halla/Host_gene.txt" "$outputdr/halla/host_gene" "$pdm" "Microbiomes" "Host_gene"

run_python_plot_safe "PLS-DA microbiome x host_gene" "python $MTDIR/pls_da_analysis.py -x $outputdr/halla/Microbiomes.txt -y $outputdr/halla/Host_gene.txt -o $outputdr/halla/pls_da_results.pdf" "$outputdr/halla/pls_da_analysis.log"

run_python_plot_safe "k-means microbiome x host_gene" "python $MTDIR/kmeans_clustering.py -x $outputdr/halla/Microbiomes.txt -y $outputdr/halla/Host_gene.txt -o $outputdr/halla/kmeans_results.pdf -k 3" "$outputdr/halla/kmeans_clustering.log"

if [[ "$RUN_EXTRA_PEARSON" == "1" ]]; then
    echo "[INFO] RUN_EXTRA_PEARSON=1, running extra Pearson HAllA analysis."
    run_halla_safe "$outputdr/halla/Microbiomes.txt" "$outputdr/halla/Host_gene.txt" "$outputdr/halla/pearson" "pearson" "Microbiomes" "Host_gene"
else
    echo "[INFO] Skipping extra Pearson HAllA run."
    echo "[INFO] To enable it: RUN_EXTRA_PEARSON=1 bash $SCRIPT_NAME ..."
fi

run_hallagram_safe "$outputdr/halla/host_gene" "$outputdr/halla/host_gene/hallagram_Top5.pdf" 5 "Microbiomes" "Host_gene" "$pdm_name"
run_hallagram_safe "$outputdr/halla/host_gene" "$outputdr/halla/host_gene/hallagram_Top10.pdf" 10 "Microbiomes" "Host_gene" "$pdm_name"
run_hallagram_safe "$outputdr/halla/host_gene" "$outputdr/halla/host_gene/hallagram_Top25.pdf" 25 "Microbiomes" "Host_gene" "$pdm_name"
run_hallagram_safe "$outputdr/halla/host_gene" "$outputdr/halla/host_gene/hallagram_Top50.pdf" 50 "Microbiomes" "Host_gene" "$pdm_name"

if [[ "$RUN_FULL_HALLAGRAM" == "1" ]]; then
    run_hallagram_safe "$outputdr/halla/host_gene" "$outputdr/halla/host_gene/hallagram_all.pdf" -1 "Microbiomes" "Host_gene" "$pdm_name"
else
    echo "[INFO] Skipping full host_gene hallagram by default."
    echo "[INFO] To enable it: RUN_FULL_HALLAGRAM=1 bash $SCRIPT_NAME ..."
fi

echo "${g}"
echo 'MTD running  progress:'
echo '>>>>>>>>>>>>>>>>>>  [90%]'

echo 'Analyzing microbiome x host_pathways associations...'
echo "${w}"

run_halla_safe "$outputdr/halla/Microbiomes.txt" "$outputdr/halla/Host_score.txt" "$outputdr/halla/pathway" "$pdm" "Microbiomes" "Host_pathway"

run_hallagram_safe "$outputdr/halla/pathway" "$outputdr/halla/pathway_hallagram_Top5.pdf" 5 "Microbiomes" "Host_pathway" "$pdm_name"
run_hallagram_safe "$outputdr/halla/pathway" "$outputdr/halla/pathway_hallagram_Top10.pdf" 10 "Microbiomes" "Host_pathway" "$pdm_name"
run_hallagram_safe "$outputdr/halla/pathway" "$outputdr/halla/pathway_hallagram_Top25.pdf" 25 "Microbiomes" "Host_pathway" "$pdm_name"
run_hallagram_safe "$outputdr/halla/pathway" "$outputdr/halla/pathway_hallagram_Top50.pdf" 50 "Microbiomes" "Host_pathway" "$pdm_name"

if [[ "$RUN_FULL_HALLAGRAM" == "1" ]]; then
    run_hallagram_safe "$outputdr/halla/pathway" "$outputdr/halla/pathway_hallagram_all.pdf" -1 "Microbiomes" "Host_pathway" "$pdm_name"
else
    echo "[INFO] Skipping full pathway hallagram by default."
fi

conda deactivate
conda activate MTD

echo "${g}"
echo 'MTD running  progress:'
echo '>>>>>>>>>>>>>>>>>>>>[100%]'
echo "MTD running is finished"
echo -e "${w}"
