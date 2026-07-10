#!/usr/bin/env bash

# Automatically accept the Anaconda Terms of Service required by
# the default Conda channels during unattended MTD installation.
export CONDA_PLUGINS_AUTO_ACCEPT_TOS=yes

# ==============================================================================
# MTD installer
# ==============================================================================
# This version reorganizes the original installer into functions while preserving
# the execution order and the local Kraken2 helper-script workflow.
#
# Deliberately not using `set -e` here. The original installer continues after
# some non-critical commands fail, and enabling it would change existing behavior.
# ==============================================================================
#
# ------------------------------------------------------------------------------
# Defaults and global paths
# ------------------------------------------------------------------------------

kmer=""                     # Kraken2 --kmer-len, used only with --build
min_l=""                    # Kraken2 --minimizer-len, used only with --build
min_s=""                    # Kraken2 --minimizer-spaces, used only with --build
read_len=75                 # Bracken read length
threads="$(nproc)"
condapath="${HOME}/miniconda3"
offline_files_folder=""

dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"

KRAKEN_ENV_LIBEXEC=""
KRAKEN_PKG_LIBEXEC=""
kraken_build_opts=()
ORIGINAL_CHANNEL_PRIORITY=""
KRAKEN_TAXONOMY_CACHE=""

# ------------------------------------------------------------------------------
# Terminal formatting and messages
# ------------------------------------------------------------------------------

init_colors() {
    if command -v tput >/dev/null 2>&1 && [[ -t 1 ]]; then
        w="$(tput sgr0)"
        r="$(tput setaf 1)"
        g="$(tput setaf 2)"
        y="$(tput setaf 3)"
        p="$(tput setaf 5)"
    else
        w=""
        r=""
        g=""
        y=""
        p=""
    fi
}

print_rule() {
    printf '%s\n' "============================================================"
}

log_info() {
    printf '%s[INFO]%s %s\n' "$g" "$w" "$*"
}

log_ok() {
    printf '%s[OK]%s %s\n' "$g" "$w" "$*"
}

log_warning() {
    printf '%s[WARNING]%s %s\n' "$y" "$w" "$*" >&2
}

log_error() {
    printf '%s[ERROR]%s %s\n' "$r" "$w" "$*" >&2
}

show_progress() {
    local bar="$1"
    local percent="$2"
    local message="$3"

    echo "${g}"
    echo "MTD installation progress:"
    printf '%-20s[%s]\n' "$bar" "$percent"
    echo "$message"
    echo "${w}"
}

restore_conda_channel_priority() {
    if [[ -n "$ORIGINAL_CHANNEL_PRIORITY" ]] && command -v conda >/dev/null 2>&1; then
        conda config --set channel_priority "$ORIGINAL_CHANNEL_PRIORITY" >/dev/null 2>&1 || true
    fi
}

on_exit() {
    local exit_status=$?
    restore_conda_channel_priority

    if (( exit_status != 0 && exit_status != 130 )); then
        echo
        log_error "MTD installation failed with exit status $exit_status."
    fi
}

on_interrupt() {
    echo
    log_error "Installation stopped by user."
    exit 130
}

trap on_exit EXIT
trap on_interrupt INT TERM

# ------------------------------------------------------------------------------
# Command-line arguments and validation
# ------------------------------------------------------------------------------

usage() {
    cat <<USAGE

Usage:

  $0 -o <offline_files_folder> [options]

Required:

  -o PATH  Persistent installation cache; created and populated if absent

Optional:

  -p PATH  Miniconda installation path
           Default: \$HOME/miniconda3

  -k INT   Kraken2 k-mer length used during --build
  -m INT   Kraken2 minimizer length used during --build
  -s INT   Kraken2 minimizer spaces used during --build
  -r INT   Bracken read length (default: 75)
  -h       Show this help message

The installer automatically downloads and installs Miniconda.

WARNING:
If the selected Miniconda directory already exists, the installer will ask
for explicit confirmation before permanently deleting it.

USAGE
}

parse_arguments() {
    if [[ $# -eq 1 && "$1" == "-h" ]]; then
        usage
        exit 0
    fi

   if [[ $# -lt 2 ]]; then
        usage
        exit 1
    fi

    while getopts ":p:o:k:m:s:r:h" option; do
        case "$option" in
            p) condapath="$OPTARG" ;;
            o) offline_files_folder="$OPTARG" ;;
            k) kmer="$OPTARG" ;;
            m) min_l="$OPTARG" ;;
            s) min_s="$OPTARG" ;;
            r) read_len="$OPTARG" ;;
            h)
                usage
                exit 0
                ;;
            :)
                log_error "Option -$OPTARG requires a value."
                usage
                exit 1
                ;;
            \?)
                log_error "Unknown option: -$OPTARG"
                usage
                exit 1
                ;;
        esac
    done
}

validate_arguments() {
    if [[ -z "$offline_files_folder" ]]; then
        log_error "The persistent installation cache path is required."
        usage
        exit 1
    fi

    # Expand a literal leading ~, including when the user passed it quoted.
    if [[ "$condapath" == "~" ]]; then
        condapath="$HOME"
    elif [[ "$condapath" == "~/"* ]]; then
        condapath="$HOME/${condapath#~/}"
    fi

    if [[ "$offline_files_folder" == "~" ]]; then
        offline_files_folder="$HOME"
    elif [[ "$offline_files_folder" == "~/"* ]]; then
        offline_files_folder="$HOME/${offline_files_folder#~/}"
    fi

    # -m allows the final path to be resolved even before Miniconda exists.
    if ! condapath="$(readlink -m -- "$condapath")"; then
        log_error "Could not resolve the Miniconda installation path."
        exit 1
    fi

    if [[ -z "$condapath" || "$condapath" == "/" || "$condapath" == "$HOME" ]]; then
        log_error "Unsafe Miniconda installation path:"
        log_error "  $condapath"
        log_error "Refusing to continue because this path must never be removed."
        exit 1
    fi

    if [[ -e "$offline_files_folder" && ! -d "$offline_files_folder" ]]; then
        log_error "The cache path exists but is not a directory:"
        log_error "  $offline_files_folder"
        exit 1
    fi

    if ! mkdir -p "$offline_files_folder"; then
        log_error "Could not create the installation cache directory:"
        log_error "  $offline_files_folder"
        exit 1
    fi

    if ! offline_files_folder="$(readlink -f -- "$offline_files_folder")"; then
        log_error "Could not resolve the installation cache path."
        exit 1
    fi

    if [[ ! -w "$offline_files_folder" ]]; then
        log_error "Installation cache directory is not writable:"
        log_error "  $offline_files_folder"
        exit 1
    fi

    if ! [[ "$threads" =~ ^[1-9][0-9]*$ ]]; then
        log_error "Invalid CPU thread count: $threads"
        exit 1
    fi

    if ! [[ "$read_len" =~ ^[1-9][0-9]*$ ]]; then
        log_error "Invalid Bracken read length: $read_len"
        exit 1
    fi

    log_info "Miniconda installation path:"
    log_info "  $condapath"

    log_info "Persistent installation cache:"
    log_info "  $offline_files_folder"
}

configure_paths_and_options() {
    KRAKEN_ENV_LIBEXEC="$condapath/envs/MTD/libexec"
    KRAKEN_PKG_LIBEXEC="$condapath/pkgs/kraken2-2.1.2-pl5262h7d875b9_0/libexec"
    KRAKEN_TAXONOMY_CACHE="$offline_files_folder/Kraken2_taxonomy_cache"

    kraken_build_opts=()

    if [[ -n "$kmer" ]]; then
        kraken_build_opts+=(--kmer-len "$kmer")
    fi
    if [[ -n "$min_l" ]]; then
        kraken_build_opts+=(--minimizer-len "$min_l")
    fi
    if [[ -n "$min_s" ]]; then
        kraken_build_opts+=(--minimizer-spaces "$min_s")
    fi
}

# ------------------------------------------------------------------------------
# General helpers
# ------------------------------------------------------------------------------

ensure_sudo_credentials() {
    if (( EUID == 0 )); then
        return 0
    fi

    if ! command -v sudo >/dev/null 2>&1; then
        log_error "sudo is required to install system dependencies."
        return 127
    fi

    # Reuse an existing valid sudo authentication, when available.
    if sudo -n -v >/dev/null 2>&1; then
        return 0
    fi

    if [[ ! -t 0 ]]; then
        log_error "Interactive sudo authentication is required."
        log_error "Run the installer directly from an interactive terminal."
        return 1
    fi

    log_info "Administrator privileges are required to install system dependencies."

    if ! sudo -v; then
        log_error "Could not obtain administrator privileges."
        return 1
    fi
}

run_as_root() {
    if (( EUID == 0 )); then
        "$@"
        return $?
    fi

    ensure_sudo_credentials || return $?

    # Authentication was validated above, so this command must not prompt.
    sudo -n "$@"
}
confirm_miniconda_removal() {
    local confirmation=""

    print_rule
    log_warning "An existing Miniconda installation was found at:"
    log_warning "  $condapath"
    echo
    log_warning "Continuing will permanently delete this directory."
    log_warning "All Conda environments and packages stored inside it will be lost."
    log_warning "This operation cannot be undone."
    print_rule

    if [[ ! -t 0 ]]; then
        log_error "Interactive confirmation is required before deleting Miniconda."
        log_error "Run the installer directly from an interactive terminal."
        exit 1
    fi

    printf '%s' "${y}Remove the existing Miniconda installation? [y/N]: ${w}"
    read -r confirmation

    case "$confirmation" in
        y | Y)
            echo
            log_warning "Removal confirmed by the user."
            ;;
        *)
            echo
            log_warning "Miniconda removal was not confirmed."
            log_warning "Installation cancelled without deleting anything."
            exit 0
            ;;
    esac
}

detect_miniconda_architecture() {
    local machine_arch

    machine_arch="$(uname -m)"

    case "$machine_arch" in
        x86_64 | amd64)
            printf '%s\n' "x86_64"
            ;;

        aarch64 | arm64)
            printf '%s\n' "aarch64"
            ;;

        *)
            log_error "Unsupported CPU architecture for Miniconda:"
            log_error "  $machine_arch"
            return 1
            ;;
    esac
}

remove_existing_miniconda() {
    if [[ ! -e "$condapath" ]]; then
        return 0
    fi

    case "$condapath" in
        "" | "/" | "$HOME")
            log_error "Refusing to remove unsafe path:"
            log_error "  $condapath"
            exit 1
            ;;
    esac

    confirm_miniconda_removal

    log_info "Removing the previous Miniconda installation:"
    log_info "  $condapath"

    if ! rm -rf -- "$condapath"; then
        log_error "Could not remove the previous Miniconda installation:"
        log_error "  $condapath"
        exit 1
    fi

    if [[ -e "$condapath" ]]; then
        log_error "The previous Miniconda directory still exists:"
        log_error "  $condapath"
        exit 1
    fi

    log_ok "Previous Miniconda installation removed."
}

install_miniconda() {
    local miniconda_arch
    local miniconda_url
    local miniconda_installer

    miniconda_arch="$(detect_miniconda_architecture)" || exit 1

    miniconda_url="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-${miniconda_arch}.sh"
    miniconda_installer="${TMPDIR:-/tmp}/mtd-miniconda-installer-${USER:-user}-$$.sh"

    remove_existing_miniconda

    log_info "Downloading the latest Miniconda installer..."
    log_info "Architecture: $miniconda_arch"
    log_info "URL: $miniconda_url"

    rm -f -- "$miniconda_installer"

    if ! curl \
        --fail \
        --location \
        --show-error \
        --retry 5 \
        --retry-delay 10 \
        --connect-timeout 30 \
        --output "$miniconda_installer" \
        "$miniconda_url"; then

        rm -f -- "$miniconda_installer"

        log_error "Could not download the Miniconda installer."
        exit 1
    fi

    if [[ ! -s "$miniconda_installer" ]]; then
        rm -f -- "$miniconda_installer"

        log_error "The downloaded Miniconda installer is empty."
        exit 1
    fi

    log_ok "Miniconda installer downloaded successfully."

    log_info "Installing Miniconda in batch mode:"
    log_info "  $condapath"

    if ! bash "$miniconda_installer" -b -p "$condapath"; then
        rm -f -- "$miniconda_installer"

        log_error "Miniconda installation failed."
        exit 1
    fi

    rm -f -- "$miniconda_installer"

    if [[ ! -x "$condapath/bin/conda" ]]; then
        log_error "The Conda executable was not created:"
        log_error "  $condapath/bin/conda"
        exit 1
    fi

    if [[ ! -f "$condapath/etc/profile.d/conda.sh" ]]; then
        log_error "The Conda initialization script was not created:"
        log_error "  $condapath/etc/profile.d/conda.sh"
        exit 1
    fi

    # Remove cached shell command paths inherited from an older installation.
    hash -r

    log_ok "Miniconda installed and validated successfully."
    log_info "Installed Conda version:"
    "$condapath/bin/conda" --version
}


safe_conda_deactivate() {
    conda deactivate >/dev/null 2>&1 || true
}

run_required_command() {
    local description="$1"
    shift

    log_info "$description"
    "$@"
    local exit_status=$?

    if (( exit_status != 0 )); then
        log_error "$description failed with exit status $exit_status."
        printf '[ERROR] Command:' >&2
        printf ' %q' "$@" >&2
        printf '\n' >&2
        exit "$exit_status"
    fi
}

activate_required_env() {
    local env_name="$1"

    if ! conda activate "$env_name"; then
        log_error "Could not activate required Conda environment: $env_name"
        exit 1
    fi

    if [[ "${CONDA_DEFAULT_ENV:-}" != "$env_name" ]]; then
        log_error "Conda reported an unexpected active environment."
        log_error "Expected: $env_name"
        log_error "Active:   ${CONDA_DEFAULT_ENV:-none}"
        exit 1
    fi

    log_ok "Activated Conda environment: $env_name"
}

require_env_command() {
    local env_name="$1"
    local command_name="$2"

    if ! conda run -n "$env_name" bash -c "command -v '$command_name' >/dev/null 2>&1"; then
        log_error "Required command '$command_name' is missing from Conda environment '$env_name'."
        exit 1
    fi

    log_ok "Found '$command_name' in Conda environment '$env_name'."
}

copy_required_file() {
    local source_file="$1"
    local destination="$2"

    if [[ ! -f "$source_file" ]]; then
        log_error "Required file not found: $source_file"
        exit 1
    fi

    if ! cp -f "$source_file" "$destination"; then
        log_error "Could not copy required file:"
        log_error "  Source:      $source_file"
        log_error "  Destination: $destination"
        exit 1
    fi
}

run_required_script() {
    local script="$1"
    shift

    if [[ ! -f "$script" ]]; then
        log_error "Required script not found: $script"
        exit 1
    fi

    chmod +x "$script"
    if ! "$script" "$@"; then
        log_error "Required script failed: $script"
        exit 1
    fi
}

retry_until_success() {
    local description="$1"
    shift

    local attempt=1
    local retry_delay="${RETRY_INITIAL_DELAY:-20}"
    local max_retry_delay="${RETRY_MAX_DELAY:-300}"
    local exit_status=0

    while true; do
        print_rule
        echo "[RETRY] $description"
        echo "[RETRY] Attempt: $attempt"
        printf '[RETRY] Command:'
        printf ' %q' "$@"
        printf '\n'
        print_rule

        if "$@"; then
            print_rule
            log_ok "$description"
            print_rule
            return 0
        else
            exit_status=$?
        fi

        if (( exit_status == 126 || exit_status == 127 )); then
            print_rule
            log_error "Permanent command failure while running: $description"
            log_error "Exit status: $exit_status"
            log_error "The command is unavailable or cannot be executed; retrying will not help."
            print_rule
            return "$exit_status"
        fi

        print_rule
        log_warning "Command failed: $description"
        log_warning "Exit status: $exit_status"
        log_warning "Retrying in $retry_delay seconds. Press Ctrl+C to stop."
        print_rule

        sleep "$retry_delay"
        attempt=$((attempt + 1))

        if (( retry_delay < max_retry_delay )); then
            retry_delay=$((retry_delay * 2))
            if (( retry_delay > max_retry_delay )); then
                retry_delay=$max_retry_delay
            fi
        fi
    done
}

# ------------------------------------------------------------------------------
# Persistent installation cache
# ------------------------------------------------------------------------------

validate_downloaded_cache_file() {
    local downloaded_file="$1"
    local expected_name="$2"

    if [[ ! -s "$downloaded_file" ]]; then
        log_warning "Downloaded cache file is empty: $downloaded_file"
        return 1
    fi

    case "$expected_name" in
        *.tar.gz)
            tar -tzf "$downloaded_file" >/dev/null 2>&1
            ;;
        *.gz)
            gzip -t "$downloaded_file" >/dev/null 2>&1
            ;;
        *)
            return 0
            ;;
    esac
}

download_cache_file_once() {
    local description="$1"
    local url="$2"
    local destination="$3"
    local partial_file="${destination}.part"
    local exit_status=0

    if ! mkdir -p "$(dirname "$destination")"; then
        log_error "Could not create cache destination directory:"
        log_error "  $(dirname "$destination")"
        return 1
    fi

    log_info "$description"
    log_info "URL:         $url"
    log_info "Destination: $destination"

    if command -v curl >/dev/null 2>&1; then
        curl \
            --fail \
            --location \
            --connect-timeout 30 \
            --retry 5 \
            --retry-delay 10 \
            --continue-at - \
            --output "$partial_file" \
            "$url"
        exit_status=$?
    elif command -v wget >/dev/null 2>&1; then
        wget \
            --continue \
            --tries=5 \
            --timeout=60 \
            --read-timeout=60 \
            --output-document="$partial_file" \
            "$url"
        exit_status=$?
    else
        log_error "Neither curl nor wget is available for downloading cache files."
        return 127
    fi

    if (( exit_status != 0 )); then
        return "$exit_status"
    fi

    if ! validate_downloaded_cache_file "$partial_file" "$destination"; then
        log_warning "Downloaded file failed archive validation:"
        log_warning "  $destination"
        rm -f "$partial_file"
        return 1
    fi

    if ! mv -f "$partial_file" "$destination"; then
        log_error "Could not finalize cached file:"
        log_error "  $destination"
        return 1
    fi

    log_ok "Cached successfully: $destination"
}

ensure_cached_file() {
    local description="$1"
    local url="$2"
    local destination="$3"

    if [[ -s "$destination" ]]; then
        log_ok "Using cached file: $destination"
        return 0
    fi

    if [[ -e "$destination" ]]; then
        log_warning "Removing empty cache file: $destination"
        rm -f "$destination"
    fi

    if ! retry_until_success \
        "$description" \
        download_cache_file_once \
        "$description" \
        "$url" \
        "$destination"; then

        log_error "Could not prepare required cache file:"
        log_error "  $destination"
        exit 1
    fi
}

prepare_installation_cache() {
    log_info "Preparing persistent MTD installation cache:"
    log_info "  $offline_files_folder"

    if ! mkdir -p \
        "$offline_files_folder/Ref_genomes/MTD_virus" \
        "$offline_files_folder/HUMAnN" \
        "$offline_files_folder/Kraken2_taxonomy_cache" \
        "$offline_files_folder/Kraken2DB_micro/library/viral/all" \
        "$offline_files_folder/Kraken2DB_micro/library/bacteria/all" \
        "$offline_files_folder/Kraken2DB_micro/library/archaea/all" \
        "$offline_files_folder/Kraken2DB_micro/library/plasmid" \
        "$offline_files_folder/eggNOG/emapperdb-5.0.2" \
        "$offline_files_folder/Customized_hosts"
    then
        log_error "Could not create the persistent cache structure:"
        log_error "  $offline_files_folder"
        exit 1
    fi

    ensure_cached_file \
        "Virus-Host DB reference" \
        "https://github.com/patrick-douglas/MTD/releases/download/virushostdb-cache-2026.07.07/virushostdb.genomic.fna.gz" \
        "$offline_files_folder/Ref_genomes/MTD_virus/virushostdb.genomic.fna.gz"
    
    ensure_cached_file \
    "HUMAnN utility mapping database" \
    "https://huttenhower.sph.harvard.edu/humann_data/full_mapping_v201901b.tar.gz" \
    "$offline_files_folder/HUMAnN/full_mapping_v201901b.tar.gz"

    ensure_cached_file \
        "HUMAnN UniRef90 annotated database" \
        "https://huttenhower.sph.harvard.edu/humann_data/uniprot/uniref_annotated/uniref90_annotated_v201901b_full.tar.gz" \
        "$offline_files_folder/HUMAnN/uniref90_annotated_v201901b_full.tar.gz"

    ensure_cached_file \
        "HUMAnN ChocoPhlAn database" \
        "https://huttenhower.sph.harvard.edu/humann_data/chocophlan/full_chocophlan.v201901_v31.tar.gz" \
        "$offline_files_folder/HUMAnN/full_chocophlan.v201901_v31.tar.gz"

    log_ok "Persistent installation cache is ready."
}

# ------------------------------------------------------------------------------
# Kraken2 helpers
# ------------------------------------------------------------------------------

download_kraken2_library_until_success() {
    local database="$1"
    local library="$2"
    shift 2

    if ! retry_until_success \
        "Kraken2 library '$library' for database '$database'" \
        kraken2-build \
        "$@" \
        --download-library "$library" \
        --threads "$threads" \
        --db "$database"; then
        log_error "Kraken2 library download cannot continue: $library"
        exit 1
    fi
}

validate_kraken2_taxonomy_dir() {
    local taxonomy_dir="$1"
    local required_file

    for required_file in \
        names.dmp \
        nodes.dmp \
        nucl_gb.accession2taxid \
        nucl_wgs.accession2taxid
    do
        if [[ ! -s "$taxonomy_dir/$required_file" ]]; then
            log_warning "Missing or empty taxonomy file:"
            log_warning "  $taxonomy_dir/$required_file"
            return 1
        fi
    done

    if ! head -n 1 "$taxonomy_dir/nucl_gb.accession2taxid" |
        grep -q $'^accession\taccession.version\ttaxid'; then
        log_warning "Invalid header in nucl_gb.accession2taxid."
        return 1
    fi

    if ! head -n 1 "$taxonomy_dir/nucl_wgs.accession2taxid" |
        grep -q $'^accession\taccession.version\ttaxid'; then
        log_warning "Invalid header in nucl_wgs.accession2taxid."
        return 1
    fi

    if ! head -n 10 "$taxonomy_dir/names.dmp" | grep -q $'\t|\t'; then
        log_warning "names.dmp does not appear to have the expected format."
        return 1
    fi

    if ! head -n 10 "$taxonomy_dir/nodes.dmp" | grep -q $'\t|\t'; then
        log_warning "nodes.dmp does not appear to have the expected format."
        return 1
    fi

    return 0
}

download_shared_kraken2_taxonomy_once() {
    local downloader="$dir/kraken2-build-download-taxonomy"
    local taxonomy_dir="$KRAKEN_TAXONOMY_CACHE/taxonomy"
    local complete_marker="$KRAKEN_TAXONOMY_CACHE/.mtd_taxonomy_complete"

    if [[ ! -f "$downloader" ]]; then
        log_error "Kraken2 taxonomy downloader not found:"
        log_error "  $downloader"
        return 127
    fi

    if ! chmod +x "$downloader"; then
        log_error "Could not make taxonomy downloader executable:"
        log_error "  $downloader"
        return 126
    fi

    log_info "Removing any incomplete taxonomy cache before download..."

    rm -rf "$taxonomy_dir"
    rm -f "$complete_marker"

    if ! mkdir -p "$KRAKEN_TAXONOMY_CACHE"; then
        log_error "Could not create taxonomy cache directory:"
        log_error "  $KRAKEN_TAXONOMY_CACHE"
        return 1
    fi

    log_info "Downloading the shared NCBI taxonomy cache..."

    if ! "$downloader" \
        --download-taxonomy \
        --threads "$threads" \
        --db "$KRAKEN_TAXONOMY_CACHE"; then

        log_warning "Shared taxonomy download failed."
        log_warning "Removing incomplete files before the next attempt."

        rm -rf "$taxonomy_dir"
        rm -f "$complete_marker"
        return 1
    fi

    if ! validate_kraken2_taxonomy_dir "$taxonomy_dir"; then
        log_warning "Downloaded taxonomy did not pass validation."
        log_warning "Removing incomplete or invalid taxonomy files."

        rm -rf "$taxonomy_dir"
        rm -f "$complete_marker"
        return 1
    fi

    date -u '+%Y-%m-%dT%H:%M:%SZ' > "$complete_marker"

    log_ok "Shared Kraken2 taxonomy downloaded and validated."
    log_info "Taxonomy cache:"
    log_info "  $KRAKEN_TAXONOMY_CACHE"

    return 0
}

prepare_shared_kraken2_taxonomy() {
    local taxonomy_dir="$KRAKEN_TAXONOMY_CACHE/taxonomy"
    local complete_marker="$KRAKEN_TAXONOMY_CACHE/.mtd_taxonomy_complete"

    if [[ -f "$complete_marker" ]] &&
       validate_kraken2_taxonomy_dir "$taxonomy_dir"; then

        log_ok "Using the existing validated Kraken2 taxonomy cache."
        log_info "Taxonomy cache:"
        log_info "  $KRAKEN_TAXONOMY_CACHE"
        return 0
    fi

    if [[ -e "$taxonomy_dir" || -e "$complete_marker" ]]; then
        log_warning "Existing taxonomy cache is incomplete or invalid."
        log_warning "It will be removed and downloaded again."
        rm -rf "$taxonomy_dir"
        rm -f "$complete_marker"
    else
        log_info "No valid shared Kraken2 taxonomy cache was found."
    fi

    if ! retry_until_success \
        "Shared NCBI taxonomy for all Kraken2 databases" \
        download_shared_kraken2_taxonomy_once; then

        log_error "The shared Kraken2 taxonomy could not be prepared."
        exit 1
    fi
}

install_shared_kraken2_taxonomy() {
    local database="$1"
    local source_taxonomy="$KRAKEN_TAXONOMY_CACHE/taxonomy"
    local destination_taxonomy="$database/taxonomy"

    prepare_shared_kraken2_taxonomy

    log_info "Installing cached NCBI taxonomy into:"
    log_info "  $database"

    if ! mkdir -p "$database"; then
        log_error "Could not create Kraken2 database directory:"
        log_error "  $database"
        exit 1
    fi

    rm -rf "$destination_taxonomy"

    if ! cp -a "$source_taxonomy" "$database/"; then
        log_error "Could not copy the shared taxonomy into:"
        log_error "  $database"
        exit 1
    fi

    if ! validate_kraken2_taxonomy_dir "$destination_taxonomy"; then
        log_error "Copied taxonomy failed validation:"
        log_error "  $destination_taxonomy"
        exit 1
    fi

    log_ok "Cached taxonomy installed and validated:"
    log_ok "  $destination_taxonomy"
}

build_kraken2_database() {
    local database="$1"

    kraken2-build \
        --build \
        --threads "$threads" \
        --db "$database" \
        "${kraken_build_opts[@]}"
}

install_kraken_helper() {
    local source_file="$1"
    local target_name="$2"

    if [[ ! -d "$KRAKEN_ENV_LIBEXEC" ]]; then
        log_error "Kraken2 environment libexec directory not found:"
        log_error "  $KRAKEN_ENV_LIBEXEC"
        exit 1
    fi

    copy_required_file "$source_file" "$KRAKEN_ENV_LIBEXEC/$target_name"
    chmod +x "$KRAKEN_ENV_LIBEXEC/$target_name"

    # Preserve the original package-cache patch when that exact Kraken2 package
    # directory exists. The active environment copy above is the required one.
    if [[ -d "$KRAKEN_PKG_LIBEXEC" ]]; then
        cp -f "$source_file" "$KRAKEN_PKG_LIBEXEC/$target_name"
        chmod +x "$KRAKEN_PKG_LIBEXEC/$target_name"
    else
        log_warning "Kraken2 package-cache libexec not found; skipped optional copy:"
        log_warning "  $KRAKEN_PKG_LIBEXEC"
    fi
}

restore_default_rsync_helper() {
    install_kraken_helper "$dir/Installation/rsync_from_ncbi.pl" "rsync_from_ncbi.pl"
}

restore_default_genomic_library_helper() {
    install_kraken_helper \
        "$dir/Installation/download_genomic_library.sh" \
        "download_genomic_library.sh"
}

patch_perl_local_download_dir() {
    local perl_script="$1"
    local local_directory="$2"

    sed -i \
        's|^[[:space:]]*my \$local_download_dir = .*;|my $local_download_dir = "'"$local_directory"'";|' \
        "$perl_script"

    if ! grep -Fq 'my $local_download_dir = "'"$local_directory"'";' "$perl_script"; then
        log_error "Could not patch Perl local_download_dir in: $perl_script"
        exit 1
    fi
}

copy_manifest_with_offline_folder() {
    local source_script="$1"
    local destination_script="$2"

    copy_required_file "$source_script" "$destination_script"
    sed -i \
        "s|^offline_files_folder=.*|offline_files_folder=$offline_files_folder|" \
        "$destination_script"
    run_required_script "$destination_script"
}


# ------------------------------------------------------------------------------
# Installation stages
# ------------------------------------------------------------------------------

initialize_installation() {
    cd "$dir" || exit 1
    printf '%s\n' "$condapath" > "$dir/condaPath"
    printf '%s\n' "$offline_files_folder" > "$dir/offlineCachePath"
    # shellcheck disable=SC1090
    source "$condapath/etc/profile.d/conda.sh"

    ORIGINAL_CHANNEL_PRIORITY="$(conda config --show channel_priority 2>/dev/null | awk '{print $2}')"
    : "${ORIGINAL_CHANNEL_PRIORITY:=flexible}"

    log_info "Original Conda channel priority: $ORIGINAL_CHANNEL_PRIORITY"
    run_required_command \
        "Setting Conda channel priority to flexible for legacy MTD environments" \
        conda config --set channel_priority flexible
}

install_system_dependencies() {
    local -a system_packages=(
        libgeos-dev
        libharfbuzz-dev
        libfribidi-dev
        libfreetype6-dev
        libpng-dev
        libtiff5-dev
        libjpeg-dev
        rsync
        pigz
        curl
        wget
        aria2
        parallel
        ca-certificates
        build-essential
        pkg-config
        libssl-dev
        sra-toolkit
    )

    log_info "Installing system dependencies..."

    run_required_command \
        "Updating APT package indexes" \
        run_as_root apt-get update

    run_required_command \
        "Installing required system packages" \
        run_as_root apt-get install -y "${system_packages[@]}"

    log_ok "System dependencies installed successfully."
}

create_conda_environments() {
    safe_conda_deactivate

    run_required_command \
        "Creating MTD-fastp environment" \
        conda env create -f "$dir/Installation/MTD_fastp.yml"

    run_required_command \
        "Creating MTD environment" \
        conda env create -f "$dir/Installation/MTD.yml"

    run_required_command \
        "Updating MTD environment with R additions" \
        conda env update -n MTD -f "$dir/Installation/MTD_R_additions.yml"

    run_required_command \
        "Installing R packages in the MTD environment" \
        conda run -n MTD bash "$dir/update_fix/Install.R.packages.MTD.sh"

    run_required_command \
        "Checking R packages in the MTD environment" \
        conda run -n MTD bash "$dir/update_fix/check_R_pkg.MTD.sh"

    log_info "Creating Python 2 and HAllA environments..."
    sed -i 's/^rpy2[>=<]/# &/' "$dir/Installation/pip.requirements"

    run_required_command \
        "Creating py2 environment" \
        conda env create -f "$dir/Installation/py2.yml"

    run_required_command \
        "Creating halla0820 environment" \
        conda env create -f "$dir/Installation/halla0820.yml"

    safe_conda_deactivate

        run_required_command \
        "Creating initial R412 environment" \
        conda env create -f "$dir/Installation/R412.yml"

    # ----------------------------------------------------------
    # Dedicated environment for custom OrgDb construction
    # ----------------------------------------------------------

    run_required_command \
        "Creating dedicated MTD_orgdb environment" \
        conda env create -f "$dir/Installation/MTD_orgdb.yml"

    require_env_command MTD_orgdb Rscript
    require_env_command MTD_orgdb jq
    require_env_command MTD_orgdb yq

    run_required_command \
        "Validating MTD_orgdb R packages" \
        conda run -n MTD_orgdb \
        Rscript "$dir/Installation/check_MTD_orgdb.R"

    sed -i '/^# *rpy2/s/^# *//' "$dir/Installation/pip.requirements"

    chmod +x "$dir/aux_scripts/ssGSEA/resolve_ssgsea_go_terms.py"
    run_required_command \
        "Checking ssGSEA GO resolver syntax" \
        python3 -m py_compile "$dir/aux_scripts/ssGSEA/resolve_ssgsea_go_terms.py"

    require_env_command MTD kraken2-build
    require_env_command MTD bracken-build
    require_env_command MTD humann
    require_env_command MTD hisat2-build
    require_env_command MTD makeblastdb
}

install_halla_dependencies() {
    activate_required_env halla0820

    conda install -n halla0820 -y -c conda-forge pkg-config
    conda install -n halla0820 -y -c conda-forge ca-certificates openssl libcurl curl
    conda install -n halla0820 -y -c conda-forge libuv

    R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/lattice/lattice_0.22-7.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
    R -e "install.packages('$dir/update_fix/pvr_pkg/Matrix_1.6-5.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
    R -e "install.packages('$dir/update_fix/pvr_pkg/mnormt_2.1.0.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
    R -e "install.packages('$dir/update_fix/pvr_pkg/nlme_3.1-167.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
    R -e "install.packages('$dir/update_fix/pvr_pkg/GPArotation_2024.3-1.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
    R -e "install.packages('$dir/update_fix/pvr_pkg/psych_2.5.3.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
    R -e "install.packages('$dir/update_fix/pvr_pkg/foreign_0.8-89.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
    R -e "install.packages('$dir/update_fix/pvr_pkg/R.methodsS3_1.8.2.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
    R -e "install.packages('$dir/update_fix/pvr_pkg/R.oo_1.27.0.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
    R -e "install.packages('$dir/update_fix/pvr_pkg/rtf_0.4-14.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
    R -e "install.packages('$dir/update_fix/pvr_pkg/psychTools_2.4.3.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
    R -e "install.packages('https://cran.r-project.org/src/contrib/XICOR_0.4.1.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
    R -e "install.packages('https://cran.r-project.org/src/contrib/mclust_6.1.2.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
    R -e 'install.packages("BiocManager", repos = "https://cloud.r-project.org")'
    R -e "install.packages('$dir/update_fix/pvr_pkg/MASS_7.3-60.tar.gz', repos=NULL, type='source')"
    R -e "install.packages('$dir/update_fix/pvr_pkg/preprocessCore_1.72.0.tar.gz', repos=NULL, type='source')"
    R -e 'install.packages("remotes", repos="https://cloud.r-project.org")'
    R -e 'remotes::install_url("https://cran.r-project.org/src/contrib/EnvStats_3.1.0.tar.gz", dependencies=TRUE)'
    R -e 'remotes::install_version("Hmisc", version = "4.8-0", repos = "https://cloud.r-project.org")'
#    R -e "install.packages('https://cran.r-project.org/src/contrib/00Archive/eva/eva_0.2.6.tar.gz', repos=NULL, type='source')"
    R -e "install.packages('$dir/update_fix/pvr_pkg/eva_0.2.6.tar.gz', repos=NULL, type='source', Ncpus=$threads)"

    conda run -n halla0820 "$dir/update_fix/check_R_pkg.halla0820.sh"
    safe_conda_deactivate
}

install_mtd_extra_tools() {
    activate_required_env MTD

    run_required_command \
        "Installing pkg-config in MTD" \
        conda install -n MTD -y -c conda-forge pkg-config
    run_required_command \
        "Installing NCBI Datasets CLI in MTD" \
        conda install -n MTD -y -c conda-forge ncbi-datasets-cli
    run_required_command \
        "Installing eggNOG-mapper and DIAMOND in MTD" \
        conda install -n MTD -y -c conda-forge -c bioconda eggnog-mapper diamond

    require_env_command MTD kraken2-build
}

prepare_virome_files() {
    copy_required_file \
        "$offline_files_folder/Ref_genomes/MTD_virus/virushostdb.genomic.fna.gz" \
        "$dir/virushostdb.genomic.fna.gz"

    unpigz -f "$dir/virushostdb.genomic.fna.gz"
    cat \
        "$dir/Installation/M33262_SIVMM239.fa" \
        "$dir/virushostdb.genomic.fna" \
        > "$dir/viruses4kraken.fa"

    log_info "Preparing additional NCBI RefSeq viral files..."
    copy_manifest_with_offline_folder \
        "$dir/manifest.virus.sh" \
        "$offline_files_folder/Kraken2DB_micro/library/manifest.virus.sh"

    sed -i -E \
        's/^>([^ ]+) (.+)/>\1 [\1] \2./' \
        "$offline_files_folder/Kraken2DB_micro/library/viral/all_viral_genomes.fna"

    # The original installer intentionally leaves this append operation disabled:
    # cat "$offline_files_folder/Kraken2DB_micro/library/viral/all_viral_genomes.fna" >> "$dir/viruses4kraken.fa"
}

install_default_kraken_helpers() {
    restore_default_rsync_helper
    restore_default_genomic_library_helper
}

prepare_microbiome_manifests() {
    copy_manifest_with_offline_folder \
        "$dir/manifest.bacteria.sh" \
        "$offline_files_folder/Kraken2DB_micro/library/manifest.bacteria.sh"

    copy_required_file \
        "$dir/manifest.sh" \
        "$offline_files_folder/Kraken2DB_micro/library/manifest.sh"

    sed -i \
        "s|^offline_files_folder=.*|offline_files_folder=$offline_files_folder|" \
        "$offline_files_folder/Kraken2DB_micro/library/manifest.sh"
}

add_local_archaea_library() {
    local database="$1"
    local helper="$KRAKEN_ENV_LIBEXEC/rsync_from_ncbi.pl"

    copy_manifest_with_offline_folder \
        "$dir/manifest.archea.sh" \
        "$offline_files_folder/Kraken2DB_micro/library/manifest.archea.sh"

    install_kraken_helper \
        "$dir/Installation/rsync_from_ncbi_archaea.pl" \
        "rsync_from_ncbi.pl"

    patch_perl_local_download_dir \
        "$helper" \
        "$offline_files_folder/Kraken2DB_micro/library/archaea/all/"

    copy_required_file \
        "$offline_files_folder/Kraken2DB_micro/library/archaea/assembly_summary_archaea.txt" \
        "$offline_files_folder/Kraken2DB_micro/library/archaea/assembly_summary.txt"

    chmod +x "$helper"
    log_info "Adding local archaeal sequences to Kraken2 database..."
    download_kraken2_library_until_success "$database" "archaea" --use-ftp

    restore_default_rsync_helper
}

add_local_bacteria_library() {
    local database="$1"
    local helper="$KRAKEN_ENV_LIBEXEC/rsync_from_ncbi.pl"

    install_kraken_helper \
        "$dir/Installation/rsync_from_ncbi_bacteria.pl" \
        "rsync_from_ncbi.pl"

    patch_perl_local_download_dir \
        "$helper" \
        "$offline_files_folder/Kraken2DB_micro/library/bacteria/all/"

    copy_required_file \
        "$offline_files_folder/Kraken2DB_micro/library/bacteria/assembly_summary_bacteria.txt" \
        "$offline_files_folder/Kraken2DB_micro/library/bacteria/assembly_summary.txt"

    chmod +x "$helper"
    log_info "Adding local bacterial sequences to Kraken2 database..."
    download_kraken2_library_until_success "$database" "bacteria" --use-ftp

    restore_default_rsync_helper
}

add_local_plasmid_library() {
    local database="$1"
    local manifest_destination
    local custom_helper
    local plasmid_cache
    local plasmid_cache_count
    local download_status=0

    manifest_destination="$offline_files_folder/Kraken2DB_micro/library/manifest.plasmid.sh"
    custom_helper="$dir/Installation/download_genomic_library_plasmid.sh"
    plasmid_cache="$offline_files_folder/Kraken2DB_micro/library/plasmid"

    log_info "Preparing cached plasmid genomic sequences..."

    copy_required_file \
        "$dir/manifest.plasmid.sh" \
        "$manifest_destination"

    sed -i \
        "s|^LOCAL_DIR=.*|LOCAL_DIR=$plasmid_cache/|" \
        "$manifest_destination"

    run_required_script "$manifest_destination"

    if [[ ! -d "$plasmid_cache" ]]; then
        log_error "Plasmid cache directory not found:"
        log_error "  $plasmid_cache"
        exit 1
    fi

    plasmid_cache_count="$(
        find "$plasmid_cache" \
            -maxdepth 1 \
            -type f \
            -name '*.genomic.fna.gz' \
            -size +0c |
        wc -l
    )"

    if (( plasmid_cache_count == 0 )); then
        log_error "No usable genomic plasmid FASTA archives were found:"
        log_error "  $plasmid_cache"
        exit 1
    fi

    log_ok "Usable genomic plasmid FASTA archives: $plasmid_cache_count"
    log_info "Plasmid cache:"
    log_info "  $plasmid_cache"

    export MTD_KRAKEN2_PLASMID_CACHE="$plasmid_cache"

    install_kraken_helper \
        "$custom_helper" \
        "download_genomic_library.sh"

    log_info "Adding local plasmid sequences to Kraken2 database..."

    retry_until_success \
        "Kraken2 plasmid library for database '$database'" \
        kraken2-build \
        --download-library plasmid \
        --threads "$threads" \
        --db "$database" ||
        download_status=$?

    log_info "Restoring the standard Kraken2 genomic-library helper..."

    restore_default_genomic_library_helper

    unset MTD_KRAKEN2_PLASMID_CACHE

    if (( download_status != 0 )); then
        log_error "The Kraken2 plasmid library could not be prepared."
        exit "$download_status"
    fi

    local plasmid_library_dir="$database/library/plasmid"
    local required_file

    for required_file in \
        library.fna \
        manifest.txt \
        prelim_map.txt
    do
        if [[ ! -s "$plasmid_library_dir/$required_file" ]]; then
            log_error "Required plasmid-library output is missing or empty:"
            log_error "  $plasmid_library_dir/$required_file"
            exit 1
        fi
    done

    log_ok "Plasmid library prepared and validated:"
    log_ok "  $plasmid_library_dir"
}

build_microbiome_kraken_database() {
    local database="$dir/kraken2DB_micro"

    prepare_microbiome_manifests

    chmod +x "$dir/kraken2-build-download-taxonomy"
    log_info "Preparing shared NCBI taxonomy for the microbiome database..."
    install_shared_kraken2_taxonomy "$database"

    add_local_archaea_library "$database"
    add_local_bacteria_library "$database"

    log_info "Downloading RefSeq Protozoa library..."
    download_kraken2_library_until_success "$database" "protozoa" --use-ftp

    log_info "Downloading RefSeq Fungi library..."
    download_kraken2_library_until_success "$database" "fungi" --use-ftp

    add_local_plasmid_library "$database"

    log_info "Downloading UniVec_Core library..."
    download_kraken2_library_until_success "$database" "UniVec_Core" --use-ftp

    log_info "Adding custom viral sequences to Kraken2 database..."
    kraken2-build \
        --add-to-library "$dir/viruses4kraken.fa" \
        --threads "$threads" \
        --db "$database"

    log_info "Building final Kraken2 microbiome database..."
    build_kraken2_database "$database"
}

build_bracken_database() {
    local database="$dir/kraken2DB_micro"

    if [[ -z "$kmer" ]]; then
        bracken-build \
            -d "$database" \
            -t "$threads" \
            -l "$read_len"
    else
        bracken-build \
            -d "$database" \
            -t "$threads" \
            -l "$read_len" \
            -k "$kmer"
    fi
}

install_humann_databases() {
    local humann_dir="$dir/HUMAnN/ref_database"

    mkdir -p "$humann_dir/chocophlan"
    mkdir -p "$humann_dir/uniref"
    mkdir -p "$humann_dir/utility_mapping"

    copy_required_file \
        "$offline_files_folder/HUMAnN/full_chocophlan.v201901_v31.tar.gz" \
        "$humann_dir/full_chocophlan.v201901_v31.tar.gz"

    copy_required_file \
        "$offline_files_folder/HUMAnN/uniref90_annotated_v201901b_full.tar.gz" \
        "$humann_dir/uniref90_annotated_v201901b_full.tar.gz"

    copy_required_file \
        "$offline_files_folder/HUMAnN/full_mapping_v201901b.tar.gz" \
        "$humann_dir/full_mapping_v201901b.tar.gz"

    tar xzvf \
        "$humann_dir/full_chocophlan.v201901_v31.tar.gz" \
        -C "$humann_dir/chocophlan/"

    tar xzvf \
        "$humann_dir/uniref90_annotated_v201901b_full.tar.gz" \
        -C "$humann_dir/uniref/"

    tar xzvf \
        "$humann_dir/full_mapping_v201901b.tar.gz" \
        -C "$humann_dir/utility_mapping/"

    humann_config --update database_folders nucleotide "$humann_dir/chocophlan"
    humann_config --update database_folders protein "$humann_dir/uniref"
    humann_config --update database_folders utility_mapping "$humann_dir/utility_mapping"
}

install_r412_and_annotation_packages() {
    safe_conda_deactivate

    conda install -n py2 -y -c conda-forge pkg-config
    activate_required_env R412
    run_required_command \
        "Installing pkg-config in R412" \
        conda install -n R412 -y -c conda-forge pkg-config

    # Preserved libcurl troubleshooting step from the original installer.
    wget \
        -T 300 \
        -t 5 \
        -N \
        --no-if-modified-since \
        https://cran.r-project.org/src/contrib/Archive/curl/curl_4.3.2.tar.gz

    if [[ -f /usr/lib/x86_64-linux-gnu/pkgconfig/libcurl.pc ]]; then
        locate_lib=/usr/lib/x86_64-linux-gnu/pkgconfig
    else
        locate_lib="$(dirname "$(locate libcurl 2>/dev/null | grep '\.pc' | head -n 1)")"
    fi

    # `locate_lib` is intentionally retained for compatibility with the old
    # troubleshooting block, although it is not consumed later in this script.
    : "${locate_lib:=}"

    cd "$dir" || exit 1
    safe_conda_deactivate

    conda env remove -n R412 -y
    rm -rf "$condapath/envs/R412"
    rm -rf "$condapath/pkgs/r-base-4.1.2-hde4fec0_0"
    rm -f "$condapath"/pkgs/r-base-4.1.2-hde4fec0_0*.tar.bz2
    rm -f "$condapath"/pkgs/r-base-4.1.2-hde4fec0_0*.conda

    conda clean --packages --tarballs -y
    run_required_command \
        "Setting Conda channel priority to strict for R412 recreation" \
        conda config --set channel_priority strict
    run_required_command \
        "Recreating R412 environment" \
        conda env create -f "$dir/Installation/R412.yml"
    activate_required_env R412

    run_required_script "$dir/update_fix/Install.R.packages.R412_optimized.sh"

    Rscript -e \
        'install.packages("https://bioconductor.org/packages/3.19/bioc/src/contrib/UCSC.utils_1.0.0.tar.gz", repos=NULL, type="source", dependencies=FALSE); library(UCSC.utils); packageVersion("UCSC.utils")'

    run_required_script "$dir/update_fix/check_R_pkg.R412.sh"
    safe_conda_deactivate

    log_info "OrgDb annotation packages are managed by the dedicated MTD_orgdb environment."

    run_required_command \
        "Rechecking MTD_orgdb environment" \
        conda run -n MTD_orgdb \
        Rscript "$dir/Installation/check_MTD_orgdb.R"
}

show_r_package_versions() {
    echo "${g}"
    echo "*********************************"
    echo "R packages version for conda envs"
    echo "*********************************"
    echo "${w}"

    conda run -n MTD "$dir/update_fix/check_R_pkg.MTD.sh"
    conda run -n R412 "$dir/update_fix/check_R_pkg.R412.sh"
    conda run -n halla0820 "$dir/update_fix/check_R_pkg.halla0820.sh"

    echo
    echo "${g}MTD_orgdb environment:${w}"

    run_required_command \
        "Reporting MTD_orgdb package versions" \
        conda run -n MTD_orgdb \
        Rscript "$dir/Installation/check_MTD_orgdb.R"
}

# ------------------------------------------------------------------------------
# Main installation sequence
# ------------------------------------------------------------------------------

main() {
    init_colors

    parse_arguments "$@"
    validate_arguments

    install_system_dependencies
    install_miniconda

    configure_paths_and_options
    initialize_installation

    show_progress ">" "5%" "Preparing persistent installation cache..."
    prepare_installation_cache

    create_conda_environments

    show_progress ">>" "10%" "Installing HAllA dependencies..."
    install_halla_dependencies

    show_progress ">>>" "15%" "Preparing virome database and MTD tools..."
    install_mtd_extra_tools
    prepare_virome_files
    install_default_kraken_helpers

    show_progress ">>>>" "20%" \
        "Preparing microbiome (virus, bacteria, archaea, protozoa, fungi, plasmid, UniVec_Core) database..."
    build_microbiome_kraken_database

    show_progress ">>>>>>>>>" "45%" "Building Bracken database..."
    build_bracken_database

    show_progress ">>>>>>>>>>>" "55%" "Installing HUMAnN3 databases..."
    install_humann_databases

    show_progress ">>>>>>>>>>>>>>>>>>" "90%" "Installing R packages..."
    install_r412_and_annotation_packages
    show_r_package_versions

echo
print_rule
log_ok "MTD Explorer installation completed successfully."
log_info "No default host database was installed."
log_info "A host database must be created or installed before analyzing real data."
log_info "Custom host setup:"
log_info "  $dir/Create_custom_host.sh"
print_rule
}

main "$@"
