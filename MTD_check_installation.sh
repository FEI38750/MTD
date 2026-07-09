#!/usr/bin/env bash

# ==============================================================================
# MTD final installation checker
# Version 2026.07.08-r7 - tolerant residual bacterial-download validation
# ==============================================================================
# Read-only validation of:
#   - system dependencies
#   - Conda environments and YAML package inventories
#   - critical executables and Python imports
#   - R/Bioconductor packages and load health
#   - Kraken2/Bracken databases and restored helper scripts
#   - HUMAnN databases/configuration
#   - host GTF/FASTA references, HISAT2 indexes, and BLAST databases
#   - installer/helper script syntax and known portability hazards
#
# Exit codes:
#   0 = no FAIL results (warnings may exist)
#   1 = one or more FAIL results
#   2 = warnings only and --strict was requested
#
# This script does not install, update, repair, download, remove, or overwrite
# any MTD component. Only reports and temporary files are written.
# ============================================================================== 

set -uo pipefail

CHECKER_VERSION="2026.07.08-r7"
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
MTD_DIR="$SCRIPT_DIR"
CONDA_PATH=""
OFFLINE_DIR=""
READ_LEN=75
MODE="full"
STRICT=0
REPORT_DIR=""
KEEP_TEMP=0

PASS_COUNT=0
WARN_COUNT=0
FAIL_COUNT=0
SKIP_COUNT=0
CHECK_COUNT=0

RESULTS_TSV=""
FULL_LOG=""
SUMMARY_TXT=""
DEPENDENCY_TSV=""
TMP_WORK=""

COLOR_RESET=""
COLOR_BOLD=""
COLOR_PASS=""
COLOR_WARN=""
COLOR_FAIL=""
COLOR_SKIP=""

usage() {
    cat <<'USAGE'
Usage:
  MTD_check_installation.sh [options]

Options:
  -m, --mtd-dir PATH       MTD installation directory
                           [default: directory containing this checker]
  -p, --conda-path PATH    Conda installation directory
                           [default: value in MTD/condaPath, then ~/miniconda3]
  -o, --offline-dir PATH   Persistent installation cache
                           [default: value in MTD/offlineCachePath]
  -r, --read-length INT    Bracken read length [default: 75]
  --mode MODE              quick, full, or deep [default: full]
  --report-dir PATH        Output report directory
  --strict                 Return status 2 when warnings exist but no failures
  --keep-temp              Keep temporary test files inside the report directory
  -h, --help               Show this help
  --version                Show checker version

Modes:
  quick  Structural checks, package inventories, command presence, and syntax.
  full   quick + R package loading, HISAT2/BLAST/HUMAnN/Kraken inspection.
  deep   full + tiny Kraken2 classifications and full offline gzip validation.

Examples:
  ./MTD_check_installation.sh
  ./MTD_check_installation.sh -p ~/miniconda3 -m ~/MTD
  ./MTD_check_installation.sh -p ~/miniconda3 -m ~/MTD --mode deep
  ./MTD_check_installation.sh -p ~/miniconda3 -m ~/MTD -o /path/MTD_Offline_Install_files
USAGE
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -m|--mtd-dir)
            MTD_DIR="${2:-}"
            shift 2
            ;;
        -p|--conda-path)
            CONDA_PATH="${2:-}"
            shift 2
            ;;
        -o|--offline-dir)
            OFFLINE_DIR="${2:-}"
            shift 2
            ;;
        -r|--read-length)
            READ_LEN="${2:-}"
            shift 2
            ;;
        --mode)
            MODE="${2:-}"
            shift 2
            ;;
        --report-dir)
            REPORT_DIR="${2:-}"
            shift 2
            ;;
        --strict)
            STRICT=1
            shift
            ;;
        --keep-temp)
            KEEP_TEMP=1
            shift
            ;;
        --version)
            printf '%s\n' "$CHECKER_VERSION"
            exit 0
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            printf 'ERROR: Unknown option: %s\n' "$1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

case "$MODE" in
    quick|full|deep) ;;
    *)
        printf 'ERROR: --mode must be quick, full, or deep. Received: %s\n' "$MODE" >&2
        exit 1
        ;;
esac

if ! [[ "$READ_LEN" =~ ^[1-9][0-9]*$ ]]; then
    printf 'ERROR: --read-length must be a positive integer.\n' >&2
    exit 1
fi

if [[ -z "$MTD_DIR" ]]; then
    printf 'ERROR: MTD directory is empty.\n' >&2
    exit 1
fi
MTD_DIR="$(readlink -f "$MTD_DIR" 2>/dev/null || printf '%s' "$MTD_DIR")"

if [[ -z "$CONDA_PATH" && -s "$MTD_DIR/condaPath" ]]; then
    CONDA_PATH="$(head -n 1 "$MTD_DIR/condaPath" | tr -d '\r\n')"
fi
if [[ -z "$CONDA_PATH" ]]; then
    CONDA_PATH="$HOME/miniconda3"
fi
CONDA_PATH="$(readlink -f "$CONDA_PATH" 2>/dev/null || printf '%s' "$CONDA_PATH")"

if [[ -z "$OFFLINE_DIR" && -s "$MTD_DIR/offlineCachePath" ]]; then
    OFFLINE_DIR="$(head -n 1 "$MTD_DIR/offlineCachePath" | tr -d '\r\n')"
fi
if [[ -n "$OFFLINE_DIR" ]]; then
    OFFLINE_DIR="$(readlink -m "$OFFLINE_DIR" 2>/dev/null || printf '%s' "$OFFLINE_DIR")"
fi

if [[ -z "$REPORT_DIR" ]]; then
    REPORT_DIR="$MTD_DIR/installation_check_$(date +%Y%m%d_%H%M%S)"
fi
mkdir -p "$REPORT_DIR" || {
    printf 'ERROR: Could not create report directory: %s\n' "$REPORT_DIR" >&2
    exit 1
}
REPORT_DIR="$(readlink -f "$REPORT_DIR")"
RESULTS_TSV="$REPORT_DIR/MTD_installation_check.tsv"
FULL_LOG="$REPORT_DIR/MTD_installation_check.log"
SUMMARY_TXT="$REPORT_DIR/MTD_installation_summary.txt"
DEPENDENCY_TSV="$REPORT_DIR/MTD_installer_dependencies.tsv"
TMP_WORK="$REPORT_DIR/.tmp"
mkdir -p "$TMP_WORK"

printf 'Status\tSection\tCheck\tDetails\n' > "$RESULTS_TSV"
printf 'Status\tRelative_path\tReferenced_by\n' > "$DEPENDENCY_TSV"
: > "$FULL_LOG"

cleanup() {
    if [[ "$KEEP_TEMP" -eq 0 ]]; then
        rm -rf "$TMP_WORK"
    fi
}

on_signal() {
    cleanup
    exit 130
}

trap cleanup EXIT
trap on_signal INT TERM

init_colors() {
    if [[ -t 1 && -z "${NO_COLOR:-}" ]] && command -v tput >/dev/null 2>&1; then
        COLOR_RESET="$(tput sgr0 2>/dev/null || true)"
        COLOR_BOLD="$(tput bold 2>/dev/null || true)"
        COLOR_PASS="$(tput setaf 2 2>/dev/null || true)"
        COLOR_WARN="$(tput setaf 3 2>/dev/null || true)"
        COLOR_FAIL="$(tput setaf 1 2>/dev/null || true)"
        COLOR_SKIP="$(tput setaf 6 2>/dev/null || true)"
    fi
}

status_color() {
    case "$1" in
        PASS) printf '%s' "$COLOR_PASS" ;;
        WARN) printf '%s' "$COLOR_WARN" ;;
        FAIL) printf '%s' "$COLOR_FAIL" ;;
        SKIP) printf '%s' "$COLOR_SKIP" ;;
        *) printf '%s' "$COLOR_RESET" ;;
    esac
}

init_colors

sanitize_tsv() {
    printf '%s' "$1" | tr '\t\r\n' '   ' | sed -E 's/[[:space:]]+/ /g; s/^ //; s/ $//'
}

record() {
    local status="$1"
    local section="$2"
    local check="$3"
    local details="${4:-}"

    CHECK_COUNT=$((CHECK_COUNT + 1))
    case "$status" in
        PASS) PASS_COUNT=$((PASS_COUNT + 1)) ;;
        WARN) WARN_COUNT=$((WARN_COUNT + 1)) ;;
        FAIL) FAIL_COUNT=$((FAIL_COUNT + 1)) ;;
        SKIP) SKIP_COUNT=$((SKIP_COUNT + 1)) ;;
    esac

    local status_colour
    status_colour="$(status_color "$status")"
    printf '%s[%-4s]%s %-18s | %-42s | %s\n' \
        "$status_colour" "$status" "$COLOR_RESET" "$section" "$check" "$details"
    printf '%s\t%s\t%s\t%s\n' \
        "$(sanitize_tsv "$status")" \
        "$(sanitize_tsv "$section")" \
        "$(sanitize_tsv "$check")" \
        "$(sanitize_tsv "$details")" \
        >> "$RESULTS_TSV"
}

log_raw() {
    printf '%s\n' "$*" >> "$FULL_LOG"
}

short_output() {
    local file="$1"
    if [[ ! -s "$file" ]]; then
        printf 'no output'
        return
    fi
    head -n 4 "$file" | tr '\t\r\n' '   ' | sed -E 's/[[:space:]]+/ /g; s/^ //; s/ $//' | cut -c1-260
}

run_test() {
    local severity="$1"
    local section="$2"
    local label="$3"
    local timeout_seconds="$4"
    shift 4

    local out="$TMP_WORK/run_${CHECK_COUNT}_$RANDOM.log"
    log_raw "============================================================"
    log_raw "CHECK: $section :: $label"
    printf 'COMMAND:' >> "$FULL_LOG"
    printf ' %q' "$@" >> "$FULL_LOG"
    printf '\n' >> "$FULL_LOG"

    if timeout "$timeout_seconds" "$@" >"$out" 2>&1; then
        cat "$out" >> "$FULL_LOG"
        record PASS "$section" "$label" "$(short_output "$out")"
        return 0
    else
        local rc=$?
        cat "$out" >> "$FULL_LOG"
        if [[ "$severity" == "WARN" ]]; then
            record WARN "$section" "$label" "exit=$rc; $(short_output "$out")"
        else
            record FAIL "$section" "$label" "exit=$rc; $(short_output "$out")"
        fi
        return "$rc"
    fi
}

check_path() {
    local severity="$1"
    local section="$2"
    local label="$3"
    local path="$4"
    local kind="${5:-file}"

    local ok=0
    case "$kind" in
        file) [[ -f "$path" ]] && ok=1 ;;
        nonempty) [[ -s "$path" ]] && ok=1 ;;
        dir) [[ -d "$path" ]] && ok=1 ;;
        nonempty_dir) [[ -d "$path" ]] && find "$path" -mindepth 1 -maxdepth 1 -print -quit 2>/dev/null | grep -q . && ok=1 ;;
        executable) [[ -x "$path" ]] && ok=1 ;;
    esac

    if [[ "$ok" -eq 1 ]]; then
        local detail="$path"
        if [[ -f "$path" ]]; then
            detail="$path ($(du -h "$path" 2>/dev/null | awk '{print $1}' || true))"
        fi
        record PASS "$section" "$label" "$detail"
        return 0
    fi

    if [[ "$severity" == "WARN" ]]; then
        record WARN "$section" "$label" "missing/invalid: $path"
    else
        record FAIL "$section" "$label" "missing/invalid: $path"
    fi
    return 1
}

check_global_command() {
    local severity="$1"
    local cmd="$2"
    local section="System"
    local found=""
    found="$(command -v "$cmd" 2>/dev/null || true)"
    if [[ -n "$found" ]]; then
        record PASS "$section" "command: $cmd" "$found"
    elif [[ "$severity" == "WARN" ]]; then
        record WARN "$section" "command: $cmd" "not found in PATH"
    else
        record FAIL "$section" "command: $cmd" "not found in PATH"
    fi
}

check_apt_package() {
    local severity="$1"
    local package="$2"
    if ! command -v dpkg-query >/dev/null 2>&1; then
        record SKIP "System packages" "$package" "dpkg-query unavailable"
        return
    fi
    local status=""
    status="$(dpkg-query -W -f='${Status}' "$package" 2>/dev/null || true)"
    if [[ "$status" == "install ok installed" ]]; then
        local version
        version="$(dpkg-query -W -f='${Version}' "$package" 2>/dev/null || true)"
        record PASS "System packages" "$package" "$version"
    elif [[ "$severity" == "WARN" ]]; then
        record WARN "System packages" "$package" "not installed according to dpkg"
    else
        record FAIL "System packages" "$package" "not installed according to dpkg"
    fi
}

prefix_for_env() {
    local env="$1"
    if [[ "$env" == "base" ]]; then
        printf '%s' "$CONDA_PATH"
    else
        printf '%s/envs/%s' "$CONDA_PATH" "$env"
    fi
}

check_env_exists() {
    local env="$1"
    local prefix
    prefix="$(prefix_for_env "$env")"
    if [[ -d "$prefix/conda-meta" && ( -x "$prefix/bin/python" || -x "$prefix/bin/R" ) ]]; then
        record PASS "Conda" "environment: $env" "$prefix"
        return 0
    fi
    record FAIL "Conda" "environment: $env" "missing or incomplete: $prefix"
    return 1
}

check_env_command() {
    local severity="$1"
    local env="$2"
    local cmd="$3"
    local prefix
    prefix="$(prefix_for_env "$env")"
    local path="$prefix/bin/$cmd"
    if [[ -x "$path" ]]; then
        record PASS "Tools/$env" "$cmd" "$path"
    elif [[ "$severity" == "WARN" ]]; then
        record WARN "Tools/$env" "$cmd" "not executable: $path"
    else
        record FAIL "Tools/$env" "$cmd" "not executable: $path"
    fi
}

check_env_command_any() {
    local severity="$1"
    local env="$2"
    local label="$3"
    shift 3
    local prefix
    prefix="$(prefix_for_env "$env")"
    local cmd
    for cmd in "$@"; do
        if [[ -x "$prefix/bin/$cmd" ]]; then
            record PASS "Tools/$env" "$label" "$prefix/bin/$cmd"
            return 0
        fi
    done
    if [[ "$severity" == "WARN" ]]; then
        record WARN "Tools/$env" "$label" "none found: $*"
    else
        record FAIL "Tools/$env" "$label" "none found: $*"
    fi
    return 1
}

conda_pkg_version() {
    local env="$1"
    local pkg="$2"
    local prefix
    prefix="$(prefix_for_env "$env")"
    "$CONDA_PATH/bin/python" - "$prefix" "$pkg" <<'PY' 2>/dev/null
import glob, json, os, sys
prefix, wanted = sys.argv[1], sys.argv[2].lower()
for f in glob.glob(os.path.join(prefix, "conda-meta", "*.json")):
    try:
        with open(f) as h:
            x = json.load(h)
        if str(x.get("name", "")).lower() == wanted:
            print(x.get("version", ""))
            raise SystemExit(0)
    except Exception:
        pass
raise SystemExit(1)
PY
}

check_conda_version() {
    local severity="$1"
    local env="$2"
    local pkg="$3"
    local expected_prefix="$4"
    local version=""
    version="$(conda_pkg_version "$env" "$pkg" || true)"
    if [[ -z "$version" ]]; then
        if [[ "$severity" == "WARN" ]]; then
            record WARN "Versions/$env" "$pkg" "package not found"
        else
            record FAIL "Versions/$env" "$pkg" "package not found"
        fi
        return 1
    fi
    if [[ "$version" == "$expected_prefix"* ]]; then
        record PASS "Versions/$env" "$pkg" "$version (expected $expected_prefix*)"
    else
        if [[ "$severity" == "WARN" ]]; then
            record WARN "Versions/$env" "$pkg" "$version; expected $expected_prefix*"
        else
            record FAIL "Versions/$env" "$pkg" "$version; expected $expected_prefix*"
        fi
        return 1
    fi
}

check_yaml_inventory() {
    local env="$1"
    shift
    local prefix
    prefix="$(prefix_for_env "$env")"
    local out="$TMP_WORK/inventory_${env}.tsv"

    if [[ ! -x "$CONDA_PATH/bin/python" ]]; then
        record SKIP "Inventory/$env" "YAML package inventory" "base Conda Python unavailable"
        return
    fi

    "$CONDA_PATH/bin/python" - "$prefix" "$@" > "$out" 2>>"$FULL_LOG" <<'PY'
import glob, json, os, re, subprocess, sys
prefix = sys.argv[1]
yamls = sys.argv[2:]


def canon(name):
    """PEP 503-like canonical name: dots, underscores, and dashes are equal."""
    return re.sub(r"[-_.]+", "-", str(name).strip().lower())


def clean_item(item):
    item = item.split("#", 1)[0].strip()
    # Be defensive around nested/accidental YAML quotes and trailing commas.
    item = item.strip().strip("'\"").strip().rstrip(",").strip().strip("'\"")
    return item


expected_conda = {}
expected_pip = {}
for path in yamls:
    if not os.path.isfile(path):
        continue
    in_deps = False
    in_pip = False
    with open(path, errors="replace") as h:
        for raw in h:
            line = raw.rstrip("\n")
            if re.match(r"^dependencies:\s*$", line):
                in_deps = True
                in_pip = False
                continue
            if not in_deps:
                continue
            m = re.match(r"^  -\s+(.+?)\s*$", line)
            if m:
                item = clean_item(m.group(1))
                if not item:
                    continue
                if item == "pip:":
                    in_pip = True
                    expected_conda.setdefault("pip", path)
                    continue
                in_pip = False
                item = item.split("::")[-1]
                name = clean_item(re.split(r"[<>=!~ ]", item, 1)[0])
                name = canon(name)
                if name:
                    expected_conda[name] = path
                continue
            m = re.match(r"^    -\s+(.+?)\s*$", line)
            if in_pip and m:
                item = clean_item(m.group(1))
                name = clean_item(re.split(r"[<>=!~ ;\[]", item, 1)[0])
                name = canon(name)
                if name:
                    expected_pip[name] = path

installed_conda = set()
for f in glob.glob(os.path.join(prefix, "conda-meta", "*.json")):
    try:
        with open(f) as h:
            x = json.load(h)
        name = canon(x.get("name", ""))
        if name:
            installed_conda.add(name)
    except Exception:
        pass

installed_pip = set()
py = os.path.join(prefix, "bin", "python")
if os.path.isfile(py):
    try:
        p = subprocess.run([py, "-m", "pip", "list", "--format=json"],
                           capture_output=True, text=True, timeout=120)
        if p.returncode == 0:
            installed_pip = {canon(x.get("name", "")) for x in json.loads(p.stdout)}
    except Exception as e:
        print("PIP_CHECK_ERROR\tpip\t%s" % str(e).replace("\t", " "))

# Conda package names that are legitimately supplied by equivalent pip
# distributions in this historical environment. Only accept the equivalent
# when the pip distribution is both declared by the YAML and actually installed.
conda_to_pip = {
    "matplotlib-base": "matplotlib",
    "seaborn-base": "seaborn",
}

for name in sorted(expected_conda):
    if name in installed_conda:
        continue
    equivalent = conda_to_pip.get(name)
    if equivalent and equivalent in expected_pip and equivalent in installed_pip:
        print("CONDA_EQUIVALENT\t%s\t%s via pip\t%s" %
              (name, equivalent, expected_conda[name]))
    else:
        print("CONDA_MISSING\t%s\t%s" % (name, expected_conda[name]))

for name in sorted(expected_pip):
    if name not in installed_pip:
        print("PIP_MISSING\t%s\t%s" % (name, expected_pip[name]))

print("SUMMARY\t%d\t%d" % (len(expected_conda), len(expected_pip)))
PY

    local missing=0
    local conda_expected=0
    local pip_expected=0
    while IFS=$'\t' read -r kind name source extra; do
        case "$kind" in
            CONDA_EQUIVALENT)
                record PASS "Inventory/$env" "package provider equivalent: $name" "$source; expected by $(basename "$extra")"
                ;;
            CONDA_MISSING)
                missing=$((missing + 1))
                record WARN "Inventory/$env" "missing Conda package: $name" "expected by $(basename "$source")"
                ;;
            PIP_MISSING)
                missing=$((missing + 1))
                record WARN "Inventory/$env" "missing pip package: $name" "expected by $(basename "$source")"
                ;;
            PIP_CHECK_ERROR)
                record WARN "Inventory/$env" "pip inventory" "$source"
                ;;
            SUMMARY)
                conda_expected="$name"
                pip_expected="$source"
                ;;
        esac
    done < "$out"

    if [[ "$missing" -eq 0 ]]; then
        record PASS "Inventory/$env" "YAML package inventory" "$conda_expected Conda + $pip_expected pip requirements satisfied"
    else
        record WARN "Inventory/$env" "YAML package inventory summary" "$missing package(s) missing; expected $conda_expected Conda + $pip_expected pip"
    fi
}

check_python_distribution_version() {
    local severity="$1"
    local env="$2"
    local distribution="$3"
    local expected="$4"
    local prefix version
    prefix="$(prefix_for_env "$env")"
    version="$("$prefix/bin/python" -c 'import pkg_resources,sys; print(pkg_resources.get_distribution(sys.argv[1]).version)' "$distribution" 2>/dev/null || true)"
    if [[ -z "$version" ]]; then
        if [[ "$severity" == "WARN" ]]; then
            record WARN "Versions/$env" "$distribution" "Python distribution not found"
        else
            record FAIL "Versions/$env" "$distribution" "Python distribution not found"
        fi
    elif [[ "$version" == "$expected" ]]; then
        record PASS "Versions/$env" "$distribution" "$version"
    elif [[ "$severity" == "WARN" ]]; then
        record WARN "Versions/$env" "$distribution" "$version; expected $expected"
    else
        record FAIL "Versions/$env" "$distribution" "$version; expected $expected"
    fi
}

check_python_imports() {
    local severity="$1"
    local env="$2"
    shift 2
    local prefix
    prefix="$(prefix_for_env "$env")"
    local imports="$(IFS=,; printf '%s' "$*")"
    local script="$TMP_WORK/import_${env}_$RANDOM.py"
    cat > "$script" <<'PY'
from __future__ import print_function
import importlib
import os
mods = [x for x in os.environ.get("MTD_IMPORTS", "").split(",") if x]
failed = []
for name in mods:
    try:
        importlib.import_module(name)
    except Exception as e:
        failed.append("%s: %s" % (name, e))
if failed:
    print(" | ".join(failed))
    raise SystemExit(1)
print("Imported: " + ", ".join(mods))
PY
    run_test "$severity" "Python/$env" "module imports" 180 \
        env MTD_IMPORTS="$imports" PATH="$prefix/bin:$PATH" \
        LD_LIBRARY_PATH="$prefix/lib:${LD_LIBRARY_PATH:-}" \
        "$prefix/bin/python" "$script"
}

check_pip_health() {
    local env="$1"
    local prefix
    prefix="$(prefix_for_env "$env")"
    if [[ ! -x "$prefix/bin/python" ]]; then
        record SKIP "Python/$env" "pip dependency health" "Python unavailable"
        return
    fi
    # A metadata conflict is important but not automatically a runtime failure.
    # Critical tools are tested separately below.
    run_test WARN "Python/$env" "pip dependency health" 180 \
        "$prefix/bin/python" -m pip check
}

resolve_eggnog_db_for_check() {
    local candidate=""
    if [[ -n "${EGGNOG_DATA_DIR:-}" ]]; then
        candidate="$EGGNOG_DATA_DIR"
    elif [[ -n "$OFFLINE_DIR" ]]; then
        candidate="$OFFLINE_DIR/eggNOG/emapperdb-5.0.2"
    elif [[ -d "$MTD_DIR/eggnog_db" ]]; then
        candidate="$MTD_DIR/eggnog_db"
    fi
    if [[ -n "$candidate" ]]; then
        readlink -f "$candidate" 2>/dev/null || printf '%s' "$candidate"
    fi
}

check_eggnog_database() {
    local dbdir
    dbdir="$(resolve_eggnog_db_for_check)"
    if [[ -z "$dbdir" || ! -d "$dbdir" ]]; then
        record FAIL "eggNOG DB" "operational database directory" \
            "not found through EGGNOG_DATA_DIR, offlineCachePath, or legacy MTD/eggnog_db"
        return
    fi

    record PASS "eggNOG DB" "operational database directory" "$dbdir"
    local required
    local missing=0
    for required in \
        eggnog.db \
        eggnog_proteins.dmnd \
        eggnog.taxa.db \
        eggnog.taxa.db.traverse.pkl; do
        if [[ -s "$dbdir/$required" ]]; then
            record PASS "eggNOG DB" "$required" "$dbdir/$required ($(du -h "$dbdir/$required" 2>/dev/null | awk '{print $1}'))"
        else
            missing=$((missing + 1))
            record FAIL "eggNOG DB" "$required" "missing or empty: $dbdir/$required"
        fi
    done
    [[ "$missing" -eq 0 ]] || return

    local pipeline="$MTD_DIR/MTD_SE.sh"
    if [[ -s "$pipeline" ]] && grep -Fq 'resolve_eggnog_db_dir' "$pipeline" \
       && grep -Fq -- '--eggnog-db "$EGGNOG_DB_DIR"' "$pipeline"; then
        record PASS "eggNOG DB" "MTD_SE cache-aware resolution" \
            "uses resolved EGGNOG_DB_DIR instead of requiring a symlink"
    elif [[ -s "$pipeline" ]] && grep -Fq -- '--eggnog-db "$MTDIR/eggnog_db"' "$pipeline"; then
        record WARN "eggNOG DB" "MTD_SE cache-aware resolution" \
            "still uses legacy MTD/eggnog_db; cache database itself is valid"
    else
        record WARN "eggNOG DB" "MTD_SE cache-aware resolution" \
            "could not confirm the expected resolved EGGNOG_DB_DIR call"
    fi

    if [[ "$MODE" != "quick" ]]; then
        run_test FAIL "eggNOG DB" "emapper database discovery" 180 \
            "$CONDA_PATH/envs/MTD/bin/emapper.py" --data_dir "$dbdir" --version
        run_test FAIL "eggNOG DB" "DIAMOND database information" 600 \
            "$CONDA_PATH/envs/MTD/bin/diamond" dbinfo --db "$dbdir/eggnog_proteins.dmnd"

        local sqlite_test="$TMP_WORK/check_eggnog_sqlite.py"
        cat > "$sqlite_test" <<'PY'
import sqlite3, sys
path = sys.argv[1]
con = sqlite3.connect("file:%s?mode=ro" % path, uri=True)
count = con.execute("SELECT COUNT(*) FROM sqlite_master WHERE type='table'").fetchone()[0]
con.close()
print("SQLite tables:", count)
if count < 1:
    raise SystemExit("database contains no tables")
PY
        run_test FAIL "eggNOG DB" "eggnog.db read-only SQLite open" 600 \
            "$CONDA_PATH/envs/MTD/bin/python" "$sqlite_test" "$dbdir/eggnog.db"
    fi
}

extract_r_packages() {
    local out="$1"
    shift
    local py="$CONDA_PATH/bin/python"
    [[ -x "$py" ]] || py="$(command -v python3 2>/dev/null || true)"
    if [[ -z "$py" ]]; then
        : > "$out"
        return
    fi

    "$py" - "$@" > "$out" <<'PY_RPKG'
import os, re, sys

files = [x for x in sys.argv[1:] if os.path.isfile(x)]
found = set()
valid = re.compile(r"^[A-Za-z][A-Za-z0-9.]*$")

def add(x):
    x = x.strip().strip('"\'')
    if x.endswith(('.tar.gz', '.tgz')):
        x = os.path.basename(x)
        x = re.sub(r"_[-0-9].*$", "", x)
    if valid.match(x):
        found.add(x)

for path in files:
    text = open(path, errors="replace").read()
    for m in re.finditer(r"(?:pkgs|packages|cran_packages|bioc_packages|core_pkgs|cran_foundation|bioc_pkgs(?:_phase[0-9]+)?|[A-Z][A-Z0-9_]*(?:PKG|PACKAGES|STACK|CRAN|BIOC)[A-Z0-9_]*)\s*(?:<-|=)\s*c?\((.*?)\)", text, re.S):
        for q in re.findall(r'["\']([^"\']+)["\']', m.group(1)):
            add(q)
    for m in re.finditer(r"(?m)^\s*([A-Za-z_][A-Za-z0-9_]*)=\((.*?)^\s*\)", text, re.S | re.M):
        name, body = m.group(1), m.group(2)
        if not re.search(r"PKG|PACKAGE|STACK|CRAN|BIOC|PATCH|DEPEND", name, re.I):
            continue
        body = re.sub(r"(?m)#.*$", "", body)
        for q in re.findall(r'["\']([^"\']+)["\']', body):
            add(q)
        body = re.sub(r'["\'][^"\']+["\']', ' ', body)
        for token in re.findall(r"\b[A-Za-z][A-Za-z0-9.]*\b", body):
            add(token)
    call_patterns = [
        r'(?:install\.packages|BiocManager::install|remotes::install_version|remotes::install_bioc|library|requireNamespace|require)\s*\(\s*["\']([^"\']+)["\']',
        r'pkg\s*=\s*["\']([^"\']+)["\']',
    ]
    for pat in call_patterns:
        for x in re.findall(pat, text):
            add(x)
    for x in re.findall(r"([A-Za-z][A-Za-z0-9.]+)_[0-9][A-Za-z0-9.+-]*\.tar\.gz", text):
        add(x)

for x in sorted(found, key=str.lower):
    print(x)
PY_RPKG
}

rscript_for_env() {
    local env="$1"
    local prefix
    prefix="$(prefix_for_env "$env")"
    if [[ -x "$prefix/bin/Rscript" ]]; then
        printf '%s' "$prefix/bin/Rscript"
    elif [[ "$env" == "base" ]]; then
        command -v Rscript 2>/dev/null || true
    fi
}

check_r_package_set() {
    local severity="$1"
    local env="$2"
    local package_file="$3"
    local label="$4"
    local prefix rscript_bin
    prefix="$(prefix_for_env "$env")"
    rscript_bin="$(rscript_for_env "$env")"

    if [[ -z "$rscript_bin" || ! -x "$rscript_bin" ]]; then
        record FAIL "R/$env" "$label" "Rscript unavailable"
        return
    fi
    if [[ ! -s "$package_file" ]]; then
        record WARN "R/$env" "$label" "package list is empty"
        return
    fi

    local rscript="$TMP_WORK/check_r_packages.R"
    cat > "$rscript" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
pkg_file <- args[[1]]
out_file <- args[[2]]
mode <- args[[3]]
pkgs <- unique(trimws(readLines(pkg_file, warn = FALSE)))
pkgs <- pkgs[nzchar(pkgs)]
ip <- installed.packages()
clip <- function(x, n=240) {
  x <- paste(x, collapse=" | ")
  x <- gsub("[\\r\\n\\t]+", " ", x)
  if (nchar(x) > n) paste0(substr(x, 1, n-3), "...") else x
}
rows <- lapply(pkgs, function(pkg) {
  if (!pkg %in% rownames(ip)) {
    return(data.frame(Package=pkg, Version="-", Status="NOT_INSTALLED", Message="", stringsAsFactors=FALSE))
  }
  ver <- ip[pkg, "Version"]
  if (identical(mode, "quick")) {
    return(data.frame(Package=pkg, Version=ver, Status="OK", Message="", stringsAsFactors=FALSE))
  }
  expr <- sprintf("suppressPackageStartupMessages(library(%s)); cat('LOAD_OK')", deparse(pkg))
  outfile <- tempfile("r_load_")
  status <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", "-e", shQuote(expr)), stdout=outfile, stderr=outfile))
  msg <- if (file.exists(outfile)) paste(readLines(outfile, warn=FALSE), collapse=" | ") else ""
  unlink(outfile)
  if (!identical(status, 0L)) {
    data.frame(Package=pkg, Version=ver, Status="LOAD_FAIL", Message=clip(msg), stringsAsFactors=FALSE)
  } else {
    data.frame(Package=pkg, Version=ver, Status="OK", Message="", stringsAsFactors=FALSE)
  }
})
res <- do.call(rbind, rows)
write.table(res, out_file, sep="\t", quote=FALSE, row.names=FALSE)
RS

    local out="$TMP_WORK/r_${env}_$(echo "$label" | tr -cs 'A-Za-z0-9' '_').tsv"
    local load_mode="$MODE"
    local checker_rc=0
    if [[ "$env" == "base" ]]; then
        env -u LD_LIBRARY_PATH PATH="$prefix/bin:$PATH" \
            "$rscript_bin" --vanilla "$rscript" "$package_file" "$out" "$load_mode" \
            >>"$FULL_LOG" 2>&1 || checker_rc=$?
    else
        env PATH="$prefix/bin:$PATH" \
            LD_LIBRARY_PATH="$prefix/lib:${LD_LIBRARY_PATH:-}" \
            "$rscript_bin" --vanilla "$rscript" "$package_file" "$out" "$load_mode" \
            >>"$FULL_LOG" 2>&1 || checker_rc=$?
    fi
    if [[ "$checker_rc" -ne 0 ]]; then
        record FAIL "R/$env" "$label" "R package checker itself failed (exit=$checker_rc)"
        return
    fi

    local total=0
    local bad=0
    while IFS=$'\t' read -r pkg version status message; do
        [[ "$pkg" == "Package" ]] && continue
        total=$((total + 1))
        case "$status" in
            OK) ;;
            NOT_INSTALLED|LOAD_FAIL)
                bad=$((bad + 1))
                if [[ "$severity" == "WARN" ]]; then
                    record WARN "R/$env" "$pkg" "$status; $message"
                else
                    record FAIL "R/$env" "$pkg" "$status; $message"
                fi
                ;;
        esac
    done < "$out"

    if [[ "$bad" -eq 0 ]]; then
        local health_detail="presence checked"
        [[ "$MODE" != "quick" ]] && health_detail="isolated library() loading checked"
        record PASS "R/$env" "$label" "$total package(s); $health_detail"
    else
        if [[ "$severity" == "WARN" ]]; then
            record WARN "R/$env" "$label summary" "$bad of $total package(s) problematic"
        else
            record FAIL "R/$env" "$label summary" "$bad of $total package(s) problematic"
        fi
    fi
}

check_r_version_and_bioc() {
    local env="$1"
    local expected_r="$2"
    local expected_bioc="$3"
    local prefix
    prefix="$(prefix_for_env "$env")"
    if [[ ! -x "$prefix/bin/Rscript" ]]; then
        record FAIL "R/$env" "R version" "Rscript unavailable"
        return
    fi
    local rver
    rver="$(env PATH="$prefix/bin:$PATH" LD_LIBRARY_PATH="$prefix/lib:${LD_LIBRARY_PATH:-}" "$prefix/bin/Rscript" --vanilla -e 'cat(as.character(getRversion()))' 2>/dev/null || true)"
    if [[ "$rver" == "$expected_r"* ]]; then
        record PASS "R/$env" "R version" "$rver (expected $expected_r*)"
    else
        record FAIL "R/$env" "R version" "$rver; expected $expected_r*"
    fi

    local bver
    bver="$(env PATH="$prefix/bin:$PATH" LD_LIBRARY_PATH="$prefix/lib:${LD_LIBRARY_PATH:-}" "$prefix/bin/Rscript" --vanilla -e 'if(requireNamespace("BiocManager",quietly=TRUE)) cat(as.character(BiocManager::version()))' 2>/dev/null || true)"
    if [[ -z "$bver" ]]; then
        record FAIL "R/$env" "Bioconductor version" "BiocManager unavailable"
    elif [[ "$bver" == "$expected_bioc"* ]]; then
        record PASS "R/$env" "Bioconductor version" "$bver (expected $expected_bioc*)"
    else
        record WARN "R/$env" "Bioconductor version" "$bver; expected $expected_bioc*"
    fi
}

check_r_runtime() {
    local env="$1"
    local prefix
    prefix="$(prefix_for_env "$env")"
    if [[ ! -x "$prefix/bin/Rscript" ]]; then
        record FAIL "R/$env" "R runtime" "Rscript unavailable"
        return
    fi
    local rver
    rver="$(env PATH="$prefix/bin:$PATH" LD_LIBRARY_PATH="$prefix/lib:${LD_LIBRARY_PATH:-}" "$prefix/bin/Rscript" --vanilla -e 'cat(as.character(getRversion()))' 2>/dev/null || true)"
    if [[ -n "$rver" ]]; then
        record PASS "R/$env" "R runtime" "$rver"
    else
        record FAIL "R/$env" "R runtime" "could not query R version"
    fi
}

check_script_syntax() {
    local path="$1"
    local type="$2"
    local severity="${3:-FAIL}"
    if [[ ! -f "$path" ]]; then
        if [[ "$severity" == "WARN" ]]; then
            record WARN "Installer source" "$(basename "$path")" "missing: $path"
        else
            record FAIL "Installer source" "$(basename "$path")" "missing: $path"
        fi
        return
    fi
    case "$type" in
        bash) run_test "$severity" "Installer source" "syntax: $(basename "$path")" 30 bash -n "$path" ;;
        perl) run_test "$severity" "Installer source" "syntax: $(basename "$path")" 30 perl -c "$path" ;;
        python) run_test "$severity" "Installer source" "compile: $(basename "$path")" 30 python3 -m py_compile "$path" ;;
    esac
}

check_installer_dependency_graph() {
    local installer="$MTD_DIR/Install.sh"
    local py="$CONDA_PATH/bin/python"
    [[ -x "$py" ]] || py="$(command -v python3 2>/dev/null || true)"

    if [[ ! -f "$installer" ]]; then
        record FAIL "Dependency graph" "Install.sh" "missing: $installer"
        return
    fi
    if [[ -z "$py" ]]; then
        record SKIP "Dependency graph" "recursive source discovery" "Python 3 unavailable"
        return
    fi

    printf 'Status\tRelative_path\tReferenced_by\n' > "$DEPENDENCY_TSV"
    local raw="$TMP_WORK/dependency_graph.raw.tsv"

    "$py" - "$MTD_DIR" > "$raw" <<'PY_DEPGRAPH'
from pathlib import Path
import re, sys

root = Path(sys.argv[1]).resolve()
queue = [(Path("Install.sh"), "entrypoint")]
seen = set()
rows = []
prefixes = ("Installation/", "update_fix/", "aux_scripts/")
top_names = re.compile(r"^(?:manifest(?:\.[A-Za-z0-9_-]+)?\.sh|kraken2-build-download-taxonomy)$")
text_ext = {".sh", ".pl", ".py", ".R", ".r", ".yml", ".yaml", ".txt", ".requirements"}
patterns = [
    re.compile(r"\$(?:\{)?(?:dir|MTD_DIR|SCRIPT_DIR)(?:\})?/([A-Za-z0-9_.+@/-]+)"),
    re.compile(r"(?<![A-Za-z0-9_.-])((?:Installation|update_fix|aux_scripts)/[A-Za-z0-9_.+@/-]+)"),
    re.compile(r"(?<![A-Za-z0-9_.-])(manifest(?:\.[A-Za-z0-9_-]+)?\.sh|kraken2-build-download-taxonomy)\b"),
]

def clean(candidate):
    candidate = candidate.strip().rstrip("'\";,):]")
    candidate = re.sub(r"\$\{?[A-Za-z_][A-Za-z0-9_]*\}?", "", candidate)
    candidate = re.sub(r"/+", "/", candidate).lstrip("./")
    if not candidate or ".." in Path(candidate).parts:
        return None
    if re.search(r"(^|/)[^/]*logs?$", candidate, re.I):
        return None
    if not (candidate.startswith(prefixes) or top_names.match(candidate)):
        return None
    return Path(candidate)

while queue:
    rel, caller = queue.pop(0)
    key = rel.as_posix()
    if key in seen:
        continue
    seen.add(key)
    path = root / rel
    if path.is_file():
        status = "EXISTS" if path.stat().st_size > 0 else "MISSING"
    elif path.is_dir():
        try:
            status = "EXISTS" if any(path.iterdir()) else "MISSING"
        except OSError:
            status = "MISSING"
    else:
        status = "MISSING"
    rows.append((status, key, caller))
    if status != "EXISTS" or path.is_dir():
        continue
    if path.stat().st_size > 8_000_000:
        continue
    if path.suffix not in text_ext and path.name != "kraken2-build-download-taxonomy":
        continue
    try:
        text = path.read_text(errors="replace")
    except Exception:
        continue
    text = re.sub(r"(?m)^[ \t]*#.*$", "", text)
    for pat in patterns:
        for match in pat.findall(text):
            child = clean(match)
            if child is not None and child.as_posix() not in seen:
                queue.append((child, key))

for row in rows:
    print("\t".join(row))
PY_DEPGRAPH

    local total=0 missing=0 existing=0
    while IFS=$'\t' read -r status rel caller; do
        [[ -z "$status" ]] && continue
        total=$((total + 1))
        printf '%s\t%s\t%s\n' "$status" "$rel" "$caller" >> "$DEPENDENCY_TSV"
        if [[ "$status" == "MISSING" ]]; then
            missing=$((missing + 1))
            record FAIL "Dependency graph" "$rel" "referenced by $caller but missing/empty"
        else
            existing=$((existing + 1))
        fi
    done < "$raw"

    if [[ "$missing" -eq 0 ]]; then
        record PASS "Dependency graph" "Install.sh recursive dependencies" "$existing source/configuration file(s) found"
    else
        record FAIL "Dependency graph" "Install.sh recursive dependencies" "$missing of $total referenced file(s) missing"
    fi
    record PASS "Dependency graph" "dependency manifest" "$DEPENDENCY_TSV"

    while IFS=$'\t' read -r status rel caller; do
        [[ "$status" == "EXISTS" ]] || continue
        case "$rel" in
            *.sh) check_script_syntax "$MTD_DIR/$rel" bash FAIL ;;
            *.pl) check_script_syntax "$MTD_DIR/$rel" perl FAIL ;;
            *.py) check_script_syntax "$MTD_DIR/$rel" python FAIL ;;
            kraken2-build-download-taxonomy)
                if head -n 1 "$MTD_DIR/$rel" 2>/dev/null | grep -qi perl; then
                    check_script_syntax "$MTD_DIR/$rel" perl FAIL
                else
                    check_script_syntax "$MTD_DIR/$rel" bash FAIL
                fi
                ;;
        esac
    done < <(tail -n +2 "$DEPENDENCY_TSV")
}

check_installer_contract() {
    local installer="$MTD_DIR/Install.sh"
    check_script_syntax "$installer" bash FAIL

    for required_function in run_as_root ensure_sudo_credentials prepare_installation_cache prepare_shared_kraken2_taxonomy; do
        if grep -Eq "^[[:space:]]*(function[[:space:]]+)?${required_function}[[:space:]]*(\\(\\))?[[:space:]]*\\{" "$installer"; then
            record PASS "Installer contract" "function: $required_function" "defined in Install.sh"
        else
            record FAIL "Installer contract" "function: $required_function" "required current installer helper is not defined"
        fi
    done

    if grep -Eq '(^|[[:space:]])sudo_with_pass([[:space:]]|$)' "$installer" && \
       ! grep -Eq '^[[:space:]]*(function[[:space:]]+)?sudo_with_pass[[:space:]]*(\(\))?[[:space:]]*\{' "$installer"; then
        record WARN "Installer contract" "undefined sudo_with_pass call" "Install.sh calls sudo_with_pass but defines run_as_root instead; source issue for future installs"
    else
        record PASS "Installer contract" "sudo helper calls" "no undefined sudo_with_pass invocation detected"
    fi

    for marker in condaPath offlineCachePath; do
        if grep -Fq "$marker" "$installer"; then
            record PASS "Installer contract" "writes $marker" "referenced by Install.sh"
        else
            record WARN "Installer contract" "writes $marker" "not found in Install.sh"
        fi
    done

    if grep -Fq 'src/contrib/00Archive/' "$installer"; then
        record WARN "Installer contract" "invalid CRAN archive URL"             "Install.sh contains src/contrib/00Archive/; source issue for future installs"
    else
        record PASS "Installer contract" "CRAN archive URL layout" "no src/contrib/00Archive/ path detected"
    fi

    if grep -Fq 'Installation/pip.requirements' "$installer"; then
        if grep -Eq '(^|[[:space:]])(python[0-9.]*[[:space:]]+-m[[:space:]]+pip|pip[0-9.]*)[[:space:]]+install[^#\n]*-[rR][[:space:]]+[^#\n]*pip\.requirements' "$installer"; then
            record PASS "Installer contract" "pip.requirements consumed" "pip install -r invocation detected"
        else
            record WARN "Installer contract" "pip.requirements consumed"                 "Install.sh edits/references Installation/pip.requirements but no pip install -r invocation was detected"
        fi
    fi
}

check_required_local_tarballs() {
    local base="$MTD_DIR/update_fix/pvr_pkg"
    local files=(
        Matrix_1.6-5.tar.gz
        mnormt_2.1.0.tar.gz
        nlme_3.1-167.tar.gz
        GPArotation_2024.3-1.tar.gz
        psych_2.5.3.tar.gz
        foreign_0.8-89.tar.gz
        R.methodsS3_1.8.2.tar.gz
        R.oo_1.27.0.tar.gz
        rtf_0.4-14.tar.gz
        psychTools_2.4.3.tar.gz
        MASS_7.3-60.tar.gz
        preprocessCore_1.72.0.tar.gz
    )
    local file
    for file in "${files[@]}"; do
        check_path FAIL "Installer source" "required HAllA tarball: $file" "$base/$file" nonempty
    done
}

check_known_source_hazards() {
    local bacteria="$MTD_DIR/manifest.bacteria.sh"
    local virus="$MTD_DIR/manifest.virus.sh"
    local plasmid_helper="$MTD_DIR/Installation/download_genomic_library_plasmid.sh"
    local tax_helper="$MTD_DIR/kraken2-build-download-taxonomy"

    for f in "$bacteria" "$virus"; do
        if [[ -f "$f" ]] && grep -q 'ftp\.ncbi\.nih\.gov' "$f"; then
            record WARN "Installer source" "invalid NCBI hostname in $(basename "$f")" "ftp.ncbi.nih.gov must be ftp.ncbi.nlm.nih.gov; source issue for future cache refreshes"
        elif [[ -f "$f" ]] && grep -q 'ftp\.ncbi\.nlm\.nih\.gov' "$f"; then
            record PASS "Installer source" "NCBI hostname in $(basename "$f")" "ftp.ncbi.nlm.nih.gov"
        fi
    done

    if [[ -f "$plasmid_helper" ]] && grep -Eq '^dir=~/?MTD' "$plasmid_helper"; then
        if [[ "$MTD_DIR" != "$HOME/MTD" ]]; then
            record WARN "Installer source" "plasmid helper portability" "hardcodes ~/MTD, but installation is $MTD_DIR"
        else
            record WARN "Installer source" "plasmid helper portability" "hardcodes ~/MTD; works only at $HOME/MTD"
        fi
    fi

    if [[ -f "$plasmid_helper" ]] && grep -Fq 'LIBRARY_DIR="$dir/$KRAKEN2_DB_NAME/library"' "$plasmid_helper"; then
        record FAIL "Installer source" "plasmid helper database path" \
            'prefixes KRAKEN2_DB_NAME with $dir; an absolute --db path becomes duplicated'
    elif [[ -f "$plasmid_helper" ]] && grep -Fq 'LIBRARY_DIR="$KRAKEN2_DB_NAME/library"' "$plasmid_helper"; then
        record PASS "Installer source" "plasmid helper database path" \
            'uses $KRAKEN2_DB_NAME/library directly'
    fi

    for f in "$MTD_DIR/Installation/rsync_from_ncbi_archaea.pl" "$MTD_DIR/Installation/rsync_from_ncbi_bacteria.pl"; do
        if [[ -f "$f" ]] && grep -q '\$ENV{HOME}/MTD/kraken2DB_micro' "$f"; then
            if [[ "$MTD_DIR" != "$HOME/MTD" ]]; then
                record WARN "Installer source" "$(basename "$f") portability" "hardcodes HOME/MTD while installation is $MTD_DIR"
            else
                record WARN "Installer source" "$(basename "$f") portability" "hardcodes HOME/MTD"
            fi
        fi
    done

    if [[ -f "$tax_helper" ]]; then
        local hardcoded
        hardcoded="$(awk -F"'" '/^my [$]KRAKEN2_DIR = / {print $2; exit}' "$tax_helper")"
        local expected="$CONDA_PATH/envs/MTD/libexec"
        if [[ -n "$hardcoded" && "$hardcoded" != "$expected" ]]; then
            record WARN "Installer source" "taxonomy helper libexec path" "hardcoded $hardcoded; active path is $expected"
        elif [[ -n "$hardcoded" ]]; then
            record PASS "Installer source" "taxonomy helper libexec path" "$hardcoded"
        fi
    fi

    if [[ -f "$MTD_DIR/manifest.sh" ]] && grep -q '^URLS_FILE="/Kraken2DB_micro' "$MTD_DIR/manifest.sh"; then
        record WARN "Installer source" "legacy manifest.sh path" "URLS_FILE starts at filesystem root; script is copied but not executed by current installer"
    fi
}

check_helper_restoration() {
    local env_libexec="$CONDA_PATH/envs/MTD/libexec"
    local src_rsync="$MTD_DIR/Installation/rsync_from_ncbi.pl"
    local dst_rsync="$env_libexec/rsync_from_ncbi.pl"
    local src_genomic="$MTD_DIR/Installation/download_genomic_library.sh"
    local dst_genomic="$env_libexec/download_genomic_library.sh"

    for tuple in \
        "$src_rsync|$dst_rsync|rsync_from_ncbi.pl" \
        "$src_genomic|$dst_genomic|download_genomic_library.sh"; do
        IFS='|' read -r src dst name <<< "$tuple"
        if [[ ! -s "$src" || ! -s "$dst" ]]; then
            record FAIL "Kraken helpers" "$name restored" "source or active helper missing"
            continue
        fi
        local a b
        a="$(sha256sum "$src" | awk '{print $1}')"
        b="$(sha256sum "$dst" | awk '{print $1}')"
        if [[ "$a" == "$b" ]]; then
            record PASS "Kraken helpers" "$name restored" "active helper matches default source"
        else
            record FAIL "Kraken helpers" "$name restored" "active helper differs from default; a temporary local helper may still be installed"
        fi
        if [[ -x "$dst" ]]; then
            record PASS "Kraken helpers" "$name executable" "$dst"
        else
            record FAIL "Kraken helpers" "$name executable" "not executable: $dst"
        fi
    done
}

check_taxonomy_directory() {
    local label="$1"
    local dir="$2"
    local severity="${3:-FAIL}"
    local missing=0
    local required=(names.dmp nodes.dmp nucl_gb.accession2taxid nucl_wgs.accession2taxid)
    local f
    for f in "${required[@]}"; do
        if [[ -s "$dir/$f" ]]; then
            record PASS "Kraken taxonomy" "$label/$f" "$dir/$f"
        else
            missing=$((missing + 1))
            if [[ "$severity" == "WARN" ]]; then
                record WARN "Kraken taxonomy" "$label/$f" "missing/empty: $dir/$f"
            else
                record FAIL "Kraken taxonomy" "$label/$f" "missing/empty: $dir/$f"
            fi
        fi
    done
    [[ "$missing" -eq 0 ]] || return 1

    if [[ "$MODE" != "quick" ]]; then
        if head -n 1 "$dir/nucl_gb.accession2taxid" 2>/dev/null | grep -q $'^accession\taccession.version\ttaxid'; then
            record PASS "Kraken taxonomy" "$label nucl_gb header" "valid"
        else
            record FAIL "Kraken taxonomy" "$label nucl_gb header" "unexpected header"
        fi
        if head -n 1 "$dir/nucl_wgs.accession2taxid" 2>/dev/null | grep -q $'^accession\taccession.version\ttaxid'; then
            record PASS "Kraken taxonomy" "$label nucl_wgs header" "valid"
        else
            record FAIL "Kraken taxonomy" "$label nucl_wgs header" "unexpected header"
        fi
        for f in names.dmp nodes.dmp; do
            if head -n 10 "$dir/$f" 2>/dev/null | grep -q $'\t|\t'; then
                record PASS "Kraken taxonomy" "$label $f format" "valid NCBI dump structure"
            else
                record FAIL "Kraken taxonomy" "$label $f format" "unexpected NCBI dump structure"
            fi
        done
    fi
}

check_shared_taxonomy_cache() {
    [[ -n "$OFFLINE_DIR" ]] || {
        record SKIP "Kraken taxonomy" "shared taxonomy cache" "offlineCachePath unavailable"
        return
    }
    local root="$OFFLINE_DIR/Kraken2_taxonomy_cache"
    local taxonomy="$root/taxonomy"
    check_path FAIL "Kraken taxonomy" "shared cache directory" "$root" dir
    check_path FAIL "Kraken taxonomy" "shared cache completion marker" "$root/.mtd_taxonomy_complete" nonempty
    check_taxonomy_directory "shared-cache" "$taxonomy" FAIL

    local db name f cache_size db_size
    for name in microbiome human mouse rhesus; do
        case "$name" in
            microbiome) db="$MTD_DIR/kraken2DB_micro" ;;
            human) db="$MTD_DIR/kraken2DB_human" ;;
            mouse) db="$MTD_DIR/kraken2DB_mice" ;;
            rhesus) db="$MTD_DIR/kraken2DB_rhesus" ;;
        esac
        for f in names.dmp nodes.dmp nucl_gb.accession2taxid nucl_wgs.accession2taxid; do
            [[ -s "$taxonomy/$f" && -s "$db/taxonomy/$f" ]] || continue
            cache_size="$(stat -c%s "$taxonomy/$f" 2>/dev/null || echo 0)"
            db_size="$(stat -c%s "$db/taxonomy/$f" 2>/dev/null || echo 1)"
            if [[ "$cache_size" == "$db_size" ]]; then
                record PASS "Kraken taxonomy" "$name copy: $f" "size matches shared cache ($cache_size bytes)"
            else
                record WARN "Kraken taxonomy" "$name copy: $f" "database size $db_size differs from cache $cache_size"
            fi
        done
    done
}

check_kraken_standard_library() {
    local db_label="$1"
    local db="$2"
    local library="$3"
    local require_manifest="${4:-yes}"
    local libdir="$db/library/$library"

    if [[ ! -d "$libdir" ]]; then
        record FAIL "Kraken library" "$db_label/$library directory" "missing: $libdir"
        return
    fi
    record PASS "Kraken library" "$db_label/$library directory" "$libdir"

    if [[ ! -s "$libdir/library.fna" ]]; then
        record FAIL "Kraken library" "$db_label/$library library.fna" "missing or empty"
    elif LC_ALL=C grep -a -m1 -q '^>kraken:taxid|' "$libdir/library.fna" 2>/dev/null; then
        record PASS "Kraken library" "$db_label/$library library.fna" \
            "non-empty FASTA with Kraken taxid headers ($(du -h "$libdir/library.fna" | awk '{print $1}'))"
    elif LC_ALL=C grep -a -m1 -q '^>' "$libdir/library.fna" 2>/dev/null; then
        record WARN "Kraken library" "$db_label/$library library.fna" \
            "FASTA detected, but headers do not show the expected kraken:taxid prefix"
    else
        record FAIL "Kraken library" "$db_label/$library library.fna" "no FASTA header detected"
    fi

    if [[ -s "$libdir/prelim_map.txt" ]]; then
        record PASS "Kraken library" "$db_label/$library prelim_map.txt" \
            "present ($(du -h "$libdir/prelim_map.txt" | awk '{print $1}'))"
    else
        record FAIL "Kraken library" "$db_label/$library prelim_map.txt" "missing or empty"
    fi

    if [[ "$require_manifest" == "yes" ]]; then
        check_path FAIL "Kraken library" "$db_label/$library manifest.txt" "$libdir/manifest.txt" nonempty
    fi
}

check_kraken_added_library() {
    local db_label="$1"
    local db="$2"
    local description="$3"
    local added="$db/library/added"
    if [[ ! -d "$added" ]]; then
        record FAIL "Kraken library" "$db_label added library" "missing: $added ($description)"
        return
    fi
    local count size
    count="$(find "$added" -type f -size +0c 2>/dev/null | wc -l)"
    size="$(du -sh "$added" 2>/dev/null | awk '{print $1}')"
    if [[ "$count" -gt 0 ]]; then
        record PASS "Kraken library" "$db_label added library" "$count non-empty file(s), $size; $description"
    else
        record FAIL "Kraken library" "$db_label added library" "directory exists but contains no non-empty files; $description"
    fi
}

# A tiny residual deficit in the very large bacterial RefSeq download set can
# persist after retries because individual upstream files may be temporarily
# unavailable or change while the assembly manifest is being refreshed.
# Treat only a very small bacterial deficit as WARN: at most 50 genomes and
# no more than 0.05% of the expected set. All other libraries remain strict.
is_tolerable_bacterial_shortfall() {
    local label="$1"
    local missing="$2"
    local expected="$3"

    [[ "$label" == *bacteria* ]] || return 1
    [[ "$missing" =~ ^[0-9]+$ ]] || return 1
    [[ "$expected" =~ ^[1-9][0-9]*$ ]] || return 1

    (( missing > 0 )) || return 1
    (( missing <= 50 )) || return 1
    (( missing * 2000 <= expected )) || return 1
}

check_manifest_cache_set() {
    local label="$1"
    local manifest="$2"
    local cache_dir="$3"
    local report="$TMP_WORK/cache_manifest_${label//[^A-Za-z0-9]/_}.txt"

    if [[ ! -s "$manifest" ]]; then
        record FAIL "Kraken cache" "$label manifest" "missing or empty: $manifest"
        return
    fi
    if [[ ! -d "$cache_dir" ]]; then
        record FAIL "Kraken cache" "$label genome directory" "missing: $cache_dir"
        return
    fi

    "$CONDA_PATH/bin/python" - "$manifest" "$cache_dir" > "$report" 2>>"$FULL_LOG" <<'PY'
import os, sys
manifest, cache = sys.argv[1:]
expected = set()
with open(manifest, errors="replace") as h:
    for raw in h:
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        expected.add(os.path.basename(line))
local = {name for name in os.listdir(cache)
         if name.endswith(".gz") and os.path.isfile(os.path.join(cache, name))
         and os.path.getsize(os.path.join(cache, name)) > 0}
missing = sorted(expected - local)
extra = sorted(local - expected)
print("expected=%d" % len(expected))
print("local=%d" % len(local))
print("missing=%d" % len(missing))
print("extra=%d" % len(extra))
for name in missing[:10]:
    print("MISSING\t" + name)
for name in extra[:10]:
    print("EXTRA\t" + name)
raise SystemExit(0 if expected and not missing else 1)
PY
    local rc=$?
    local expected local_count missing extra
    expected="$(awk -F= '$1=="expected"{print $2}' "$report")"
    local_count="$(awk -F= '$1=="local"{print $2}' "$report")"
    missing="$(awk -F= '$1=="missing"{print $2}' "$report")"
    extra="$(awk -F= '$1=="extra"{print $2}' "$report")"
    cat "$report" >> "$FULL_LOG"
    if [[ "$rc" -eq 0 ]]; then
        record PASS "Kraken cache" "$label manifest completeness" \
            "$local_count local genome(s); $expected expected; missing=0; extra=$extra"
    else
        if is_tolerable_bacterial_shortfall "$label" "${missing:-}" "${expected:-}"; then
            record WARN "Kraken cache" "$label manifest completeness" \
                "$local_count local; $expected expected; missing=$missing; residual upstream shortfall within tolerance (<=0.05%, max 50); see full log"
        else
            record FAIL "Kraken cache" "$label manifest completeness" \
                "$local_count local; $expected expected; missing=${missing:-unknown}; see full log"
        fi
    fi
}

check_installed_manifest_alignment() {
    local label="$1"
    local installed_manifest="$2"
    local expected_manifest="$3"
    local cache_dir="$4"
    local mode="${5:-manifest}"
    local report="$TMP_WORK/installed_alignment_${label//[^A-Za-z0-9]/_}.txt"

    if [[ ! -s "$installed_manifest" ]]; then
        record FAIL "Kraken library" "$label installed manifest" "missing or empty: $installed_manifest"
        return
    fi

    "$CONDA_PATH/bin/python" - "$installed_manifest" "$expected_manifest" "$cache_dir" "$mode" > "$report" 2>>"$FULL_LOG" <<'PY'
import os, sys
installed_path, expected_path, cache_dir, mode = sys.argv[1:]

def names_from_file(path):
    out = set()
    if os.path.isfile(path):
        with open(path, errors="replace") as h:
            for raw in h:
                line = raw.strip()
                if line and not line.startswith("#"):
                    out.add(os.path.basename(line))
    return out

installed = names_from_file(installed_path)
if mode == "plasmid_genomic_cache":
    # Only genomic nucleotide FASTA archives are inputs to the Kraken2
    # plasmid library. Protein, GenBank, RNA, and WGS-master archives in
    # the same NCBI cache directory are intentionally excluded.
    expected = {n for n in os.listdir(cache_dir)
                if n.endswith(".genomic.fna.gz")
                and os.path.isfile(os.path.join(cache_dir, n))
                and os.path.getsize(os.path.join(cache_dir, n)) > 0}
elif mode == "cache":
    expected = {n for n in os.listdir(cache_dir)
                if n.endswith(".gz") and os.path.isfile(os.path.join(cache_dir, n))
                and os.path.getsize(os.path.join(cache_dir, n)) > 0}
else:
    expected = names_from_file(expected_path)
missing = sorted(expected - installed)
extra = sorted(installed - expected)
print("expected=%d" % len(expected))
print("installed=%d" % len(installed))
print("missing=%d" % len(missing))
print("extra=%d" % len(extra))
for n in missing[:10]: print("MISSING\t" + n)
raise SystemExit(0 if expected and not missing else 1)
PY
    local rc=$?
    cat "$report" >> "$FULL_LOG"
    local exp inst miss extra
    exp="$(awk -F= '$1=="expected"{print $2}' "$report")"
    inst="$(awk -F= '$1=="installed"{print $2}' "$report")"
    miss="$(awk -F= '$1=="missing"{print $2}' "$report")"
    extra="$(awk -F= '$1=="extra"{print $2}' "$report")"
    if [[ "$rc" -eq 0 ]]; then
        record PASS "Kraken library" "$label installed/cache alignment" \
            "$inst installed manifest entries; $exp expected; missing=0; extra=$extra"
    else
        record FAIL "Kraken library" "$label installed/cache alignment" \
            "$inst installed; $exp expected; missing=${miss:-unknown}; see full log"
    fi
}

check_kraken_database() {
    local name="$1"
    local db="$2"
    local severity="${3:-FAIL}"
    local missing=0
    for f in hash.k2d opts.k2d taxo.k2d; do
        if [[ ! -s "$db/$f" ]]; then
            missing=$((missing + 1))
        fi
    done
    if [[ "$missing" -eq 0 ]]; then
        record PASS "Kraken DB" "$name core files" "$db"
    elif [[ "$severity" == "WARN" ]]; then
        record WARN "Kraken DB" "$name core files" "$missing Kraken2 core file(s) missing in $db"
        return
    else
        record FAIL "Kraken DB" "$name core files" "$missing Kraken2 core file(s) missing in $db"
        return
    fi

    check_taxonomy_directory "$name" "$db/taxonomy" "$severity"

    if [[ "$MODE" != "quick" ]]; then
        local inspect="$CONDA_PATH/envs/MTD/bin/kraken2-inspect"
        local out="$TMP_WORK/inspect_${name}.log"
        log_raw "============================================================"
        log_raw "CHECK: Kraken DB :: inspect $name"
        local -a inspect_args=(--db "$db")
        if "$inspect" --help 2>&1 | grep -q -- '--skip-counts'; then
            inspect_args+=(--skip-counts)
        fi
        set +o pipefail
        timeout 300 "$inspect" "${inspect_args[@]}" 2>&1 | head -n 8 > "$out"
        local rc=${PIPESTATUS[0]}
        set -o pipefail
        cat "$out" >> "$FULL_LOG"
        if [[ "$rc" -eq 0 || "$rc" -eq 141 ]]; then
            record PASS "Kraken DB" "inspect $name" "$(short_output "$out")"
        elif [[ "$rc" -eq 124 ]] && grep -q '^# Database options:' "$out" && grep -q '^# Total taxonomy nodes:' "$out"; then
            record PASS "Kraken DB" "inspect $name" "database opened and header validated; full inspection exceeded timeout"
        else
            record FAIL "Kraken DB" "inspect $name" "exit=$rc; $(short_output "$out")"
        fi
    fi

    if [[ "$MODE" == "deep" ]]; then
        local fasta="$TMP_WORK/kraken_smoke.fa"
        local output="$TMP_WORK/kraken_${name}.out"
        local report="$TMP_WORK/kraken_${name}.report"
        printf '>mtd_smoke\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n' > "$fasta"
        run_test FAIL "Kraken DB" "tiny classification: $name" 600 \
            "$CONDA_PATH/envs/MTD/bin/kraken2" --db "$db" --threads 1 \
            --output "$output" --report "$report" "$fasta"
    fi
}

check_bracken_database() {
    local db="$MTD_DIR/kraken2DB_micro"
    local distrib="$db/database${READ_LEN}mers.kmer_distrib"
    if [[ -s "$distrib" ]]; then
        record PASS "Bracken DB" "read length $READ_LEN distribution" "$distrib ($(du -h "$distrib" | awk '{print $1}'))"
    else
        record FAIL "Bracken DB" "read length $READ_LEN distribution" "missing: $distrib"
    fi
    if [[ -s "$db/database.kraken" ]]; then
        record PASS "Bracken DB" "database.kraken" "$db/database.kraken"
    else
        record WARN "Bracken DB" "database.kraken" "not found; some Bracken builds may not retain this file"
    fi
}

check_humann() {
    local root="$MTD_DIR/HUMAnN/ref_database"
    local nuc="$root/chocophlan"
    local prot="$root/uniref"
    local util="$root/utility_mapping"
    check_path FAIL "HUMAnN" "ChocoPhlAn directory" "$nuc" nonempty_dir
    check_path FAIL "HUMAnN" "UniRef directory" "$prot" nonempty_dir
    check_path FAIL "HUMAnN" "utility mapping directory" "$util" nonempty_dir

    local count
    count="$(find "$nuc" -type f 2>/dev/null | wc -l)"
    [[ "$count" -gt 0 ]] && record PASS "HUMAnN" "ChocoPhlAn files" "$count file(s)" || record FAIL "HUMAnN" "ChocoPhlAn files" "none found"
    count="$(find "$prot" -type f \( -name '*.dmnd' -o -name '*.faa*' \) 2>/dev/null | wc -l)"
    [[ "$count" -gt 0 ]] && record PASS "HUMAnN" "UniRef database files" "$count candidate file(s)" || record FAIL "HUMAnN" "UniRef database files" "no .dmnd/.faa files found"
    count="$(find "$util" -type f 2>/dev/null | wc -l)"
    [[ "$count" -gt 0 ]] && record PASS "HUMAnN" "utility mapping files" "$count file(s)" || record FAIL "HUMAnN" "utility mapping files" "none found"

    if [[ "$MODE" != "quick" ]]; then
        local cfg="$TMP_WORK/humann_config.txt"
        if timeout 60 "$CONDA_PATH/envs/MTD/bin/humann_config" --print > "$cfg" 2>&1; then
            cat "$cfg" >> "$FULL_LOG"
            local miss=0
            grep -Fq "$nuc" "$cfg" || miss=$((miss + 1))
            grep -Fq "$prot" "$cfg" || miss=$((miss + 1))
            grep -Fq "$util" "$cfg" || miss=$((miss + 1))
            if [[ "$miss" -eq 0 ]]; then
                record PASS "HUMAnN" "configured database paths" "all three paths point inside $root"
            else
                record FAIL "HUMAnN" "configured database paths" "$miss expected path(s) absent from humann_config"
            fi
        else
            record FAIL "HUMAnN" "humann_config --print" "command failed: $(short_output "$cfg")"
        fi
    fi
}

check_hisat_index() {
    local species="$1"
    local dir="$2"
    local prefix="$dir/genome_tran"
    local shards
    shards="$(find "$dir" -maxdepth 1 -type f \( -name 'genome_tran.*.ht2' -o -name 'genome_tran.*.ht2l' \) 2>/dev/null | wc -l)"
    if [[ "$shards" -ge 8 ]]; then
        record PASS "HISAT2" "$species index shards" "$shards files in $dir"
    else
        record FAIL "HISAT2" "$species index shards" "$shards found; expected at least 8"
        return
    fi
    for f in genome.fa genome.gtf genome.ss genome.exon; do
        check_path FAIL "HISAT2" "$species $f" "$dir/$f" nonempty
    done

    if [[ "$MODE" != "quick" ]]; then
        run_test FAIL "HISAT2" "inspect $species index" 300 \
            "$CONDA_PATH/envs/MTD/bin/hisat2-inspect" -s "$prefix"
    fi
}

check_blast_database() {
    local species="$1"
    local base="$2"
    local nfiles
    nfiles="$(find "$(dirname "$base")" -maxdepth 1 -type f -name "$(basename "$base").n*" 2>/dev/null | wc -l)"
    if [[ "$nfiles" -gt 0 ]]; then
        record PASS "BLAST DB" "$species database files" "$nfiles file(s) for $base"
    else
        record FAIL "BLAST DB" "$species database files" "no nucleotide BLAST database files for $base"
        return
    fi
    if [[ "$MODE" != "quick" ]]; then
        run_test FAIL "BLAST DB" "blastdbcmd info: $species" 180 \
            "$CONDA_PATH/envs/MTD/bin/blastdbcmd" -db "$base" -info
    fi
}

check_reference_gzip() {
    local label="$1"
    local path="$2"
    if [[ ! -s "$path" ]]; then
        record FAIL "Host references" "$label" "missing: $path"
        return
    fi
    if [[ "$MODE" == "quick" ]]; then
        record PASS "Host references" "$label" "$path"
    else
        run_test FAIL "Host references" "gzip integrity: $label" 300 gzip -t "$path"
    fi
}

check_small_script_smokes() {
    local resolver="$MTD_DIR/aux_scripts/ssGSEA/resolve_ssgsea_go_terms.py"
    local exons="$MTD_DIR/Installation/hisat2_extract_exons.py"
    local splice="$MTD_DIR/Installation/hisat2_extract_splice_sites.py"

    if [[ -f "$resolver" ]]; then
        local gct="$TMP_WORK/smoke.gct"
        local outg="$TMP_WORK/smoke_named.gct"
        local outm="$TMP_WORK/smoke_map.tsv"
        printf '#1.2\n1\t1\nName\tDescription\tS1\nGO:0008150\ttest\t1\n' > "$gct"
        run_test FAIL "Script smoke" "ssGSEA GO resolver" 60 \
            "$CONDA_PATH/envs/MTD/bin/python" "$resolver" \
            --scores "$gct" --out-gct "$outg" --out-map "$outm" \
            --go-db no --quickgo no --quickgo-search no
        [[ -s "$outg" && -s "$outm" ]] \
            && record PASS "Script smoke" "ssGSEA resolver outputs" "GCT and map generated" \
            || record FAIL "Script smoke" "ssGSEA resolver outputs" "expected outputs missing"
    fi

    local gtf="$TMP_WORK/smoke.gtf"
    printf 'chr1\ttest\texon\t1\t100\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\nchr1\ttest\texon\t201\t300\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n' > "$gtf"
    if [[ -f "$exons" ]]; then
        run_test FAIL "Script smoke" "HISAT2 exon extractor" 30 python3 "$exons" "$gtf"
    fi
    if [[ -f "$splice" ]]; then
        run_test FAIL "Script smoke" "HISAT2 splice-site extractor" 30 python3 "$splice" "$gtf"
    fi
}

validate_cache_archive() {
    local label="$1"
    local path="$2"
    [[ -s "$path" ]] || return
    [[ "$MODE" == "quick" ]] && return
    case "$path" in
        *.tar.gz) run_test FAIL "Installation cache" "archive integrity: $label" 900 tar -tzf "$path" ;;
        *.gz) run_test FAIL "Installation cache" "gzip integrity: $label" 900 gzip -t "$path" ;;
    esac
}

check_installation_cache() {
    if [[ -z "$OFFLINE_DIR" ]]; then
        record FAIL "Installation cache" "persistent cache" \
            "no -o value and MTD/offlineCachePath is absent; current installer requires a persistent cache"
        return
    fi
    if [[ ! -d "$OFFLINE_DIR" ]]; then
        record FAIL "Installation cache" "cache directory" "missing: $OFFLINE_DIR"
        return
    fi
    record PASS "Installation cache" "cache directory" "$OFFLINE_DIR"

    if [[ -s "$MTD_DIR/offlineCachePath" ]]; then
        local saved
        saved="$(head -n 1 "$MTD_DIR/offlineCachePath" | tr -d '\r\n')"
        saved="$(readlink -m "$saved" 2>/dev/null || printf '%s' "$saved")"
        if [[ "$saved" == "$OFFLINE_DIR" ]]; then
            record PASS "Paths" "MTD/offlineCachePath" "$saved"
        else
            record WARN "Paths" "MTD/offlineCachePath" "contains $saved; checker uses $OFFLINE_DIR"
        fi
    else
        record FAIL "Paths" "MTD/offlineCachePath" "missing or empty"
    fi

    local files=(
        "Ref_genomes/MTD_virus/virushostdb.genomic.fna.gz"
        "Ref_genomes/Mus_musculus/GCF_000001635.27_GRCm39_genomic.fna.gz"
        "Ref_genomes/Mus_musculus/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
        "Ref_genomes/Mus_musculus/Mus_musculus.GRCm39.104.gtf.gz"
        "Ref_genomes/Macaca_mulatta/Macaca_mulatta.Mmul_10.104.gtf.gz"
        "Ref_genomes/Macaca_mulatta/Macaca_mulatta.Mmul_10.dna.toplevel.fa.gz"
        "Ref_genomes/Macaca_mulatta/GCF_003339765.1_Mmul_10_genomic.fna.gz"
        "Ref_genomes/Homo_sapiens/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
        "Ref_genomes/Homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz"
        "HUMAnN/full_mapping_v201901b.tar.gz"
        "HUMAnN/uniref90_annotated_v201901b_full.tar.gz"
        "HUMAnN/full_chocophlan.v201901_v31.tar.gz"
    )
    local rel
    for rel in "${files[@]}"; do
        check_path FAIL "Installation cache" "$rel" "$OFFLINE_DIR/$rel" nonempty
        validate_cache_archive "$rel" "$OFFLINE_DIR/$rel"
    done

    # Exact offline completeness for the three URL-manifest libraries.
    check_manifest_cache_set viral \
        "$OFFLINE_DIR/Kraken2DB_micro/library/viral/manifest_viral.list.txt" \
        "$OFFLINE_DIR/Kraken2DB_micro/library/viral/all"
    check_manifest_cache_set bacteria \
        "$OFFLINE_DIR/Kraken2DB_micro/library/bacteria/manifest_bacteria.list.txt" \
        "$OFFLINE_DIR/Kraken2DB_micro/library/bacteria/all"
    check_manifest_cache_set archaea \
        "$OFFLINE_DIR/Kraken2DB_micro/library/archaea/manifest_archaea.list.txt" \
        "$OFFLINE_DIR/Kraken2DB_micro/library/archaea/all"

    check_path FAIL "Installation cache" "viral combined FASTA" \
        "$OFFLINE_DIR/Kraken2DB_micro/library/viral/all_viral_genomes.fna" nonempty

    for rel in \
        "Kraken2DB_micro/library/viral/assembly_summary_viral.txt" \
        "Kraken2DB_micro/library/viral/manifest_viral.list.txt" \
        "Kraken2DB_micro/library/bacteria/assembly_summary_bacteria.txt" \
        "Kraken2DB_micro/library/bacteria/manifest_bacteria.list.txt" \
        "Kraken2DB_micro/library/archaea/assembly_summary_archaea.txt" \
        "Kraken2DB_micro/library/archaea/manifest_archaea.list.txt"; do
        check_path FAIL "Installation cache" "$rel" "$OFFLINE_DIR/$rel" nonempty
    done

    local plasmid_cache="$OFFLINE_DIR/Kraken2DB_micro/library/plasmid"
    check_path FAIL "Installation cache" "Kraken2DB_micro/library/plasmid" "$plasmid_cache" nonempty_dir
    local plasmid_count
    plasmid_count="$(find "$plasmid_cache" -maxdepth 1 -type f -name '*.genomic.fna.gz' -size +0c 2>/dev/null | wc -l)"
    if [[ "$plasmid_count" -gt 0 ]]; then
        record PASS "Kraken cache" "plasmid genomic FASTA archives"             "$plasmid_count non-empty *.genomic.fna.gz file(s)"
    else
        record FAIL "Kraken cache" "plasmid genomic FASTA archives"             "no non-empty *.genomic.fna.gz files found in $plasmid_cache"
    fi

    local eggnog="$OFFLINE_DIR/eggNOG/emapperdb-5.0.2"
    check_path FAIL "Installation cache" "eggNOG cache directory" "$eggnog" dir
    for rel in eggnog.db eggnog_proteins.dmnd eggnog.taxa.db eggnog.taxa.db.traverse.pkl; do
        check_path FAIL "Installation cache" "eggNOG/$rel" "$eggnog/$rel" nonempty
    done
    check_path FAIL "Installation cache" "custom host cache directory" "$OFFLINE_DIR/Customized_hosts" dir

    local failure_file found_failure=0
    while IFS= read -r -d '' failure_file; do
        found_failure=1
        if [[ -s "$failure_file" ]]; then
            record FAIL "Installation cache" "failed-download record" \
                "$failure_file is non-empty; database cache cannot be considered complete"
        else
            record PASS "Installation cache" "failed-download record" "$failure_file is empty"
        fi
    done < <(find "$OFFLINE_DIR" -maxdepth 2 -type f -name 'failed_downloads*.txt' -print0 2>/dev/null)
    [[ "$found_failure" -eq 1 ]] || record SKIP "Installation cache" "failed-download records" "none present"

    check_shared_taxonomy_cache

    if [[ "$MODE" == "deep" ]]; then
        local bad=0 total=0 gz
        while IFS= read -r -d '' gz; do
            total=$((total + 1))
            if ! gzip -t "$gz" >/dev/null 2>&1; then
                bad=$((bad + 1))
                record FAIL "Installation cache" "corrupted gzip" "$gz"
            fi
        done < <(find "$OFFLINE_DIR/Kraken2DB_micro/library" -type f -name '*.gz' -print0 2>/dev/null)
        if [[ "$bad" -eq 0 && "$total" -gt 0 ]]; then
            record PASS "Installation cache" "all Kraken library gzip files" "$total file(s) passed gzip -t"
        elif [[ "$total" -eq 0 ]]; then
            record FAIL "Installation cache" "Kraken library gzip summary" "no gzip files found"
        else
            record FAIL "Installation cache" "Kraken library gzip summary" "$bad of $total corrupted"
        fi
    fi
}

print_header() {
    cat <<EOF
============================================================
MTD final installation check v$CHECKER_VERSION
============================================================
MTD directory : $MTD_DIR
Conda path    : $CONDA_PATH
Cache path    : ${OFFLINE_DIR:-not detected}
Mode          : $MODE
Bracken length: $READ_LEN
Report folder : $REPORT_DIR
Dependency map: $DEPENDENCY_TSV
============================================================
EOF
}

print_header | tee -a "$FULL_LOG"

# ------------------------------------------------------------------------------
# 1. Root paths and system dependencies
# ------------------------------------------------------------------------------
check_path FAIL "Paths" "MTD directory" "$MTD_DIR" dir
check_path FAIL "Paths" "current Install.sh" "$MTD_DIR/Install.sh" nonempty
check_path FAIL "Paths" "Conda executable" "$CONDA_PATH/bin/conda" executable
check_path FAIL "Paths" "Conda initialization" "$CONDA_PATH/etc/profile.d/conda.sh" nonempty

if [[ -s "$MTD_DIR/condaPath" ]]; then
    saved_conda="$(head -n1 "$MTD_DIR/condaPath" | tr -d '\r\n')"
    saved_conda="$(readlink -f "$saved_conda" 2>/dev/null || printf '%s' "$saved_conda")"
    if [[ "$saved_conda" == "$CONDA_PATH" ]]; then
        record PASS "Paths" "MTD/condaPath" "$saved_conda"
    else
        record WARN "Paths" "MTD/condaPath" "contains $saved_conda; checker uses $CONDA_PATH"
    fi
else
    record WARN "Paths" "MTD/condaPath" "file missing or empty"
fi

for cmd in bash awk sed grep find xargs gzip gunzip tar wget curl rsync pigz unpigz perl python3 sha256sum timeout pkg-config; do
    check_global_command FAIL "$cmd"
done
check_global_command WARN aria2c
check_global_command WARN parallel

for package in \
    libgeos-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev \
    libpng-dev libtiff5-dev libjpeg-dev rsync pigz curl wget \
    ca-certificates build-essential pkg-config libssl-dev; do
    check_apt_package FAIL "$package"
done

for pc in geos harfbuzz fribidi freetype2 libpng libtiff-4; do
    if pkg-config --exists "$pc" >/dev/null 2>&1; then
        record PASS "System libraries" "$pc" "$(pkg-config --modversion "$pc" 2>/dev/null || true)"
    else
        record WARN "System libraries" "$pc" "not visible through pkg-config"
    fi
done

# ------------------------------------------------------------------------------
# 2. Conda environments, package inventories, versions, and tools
# ------------------------------------------------------------------------------
for f in \
    "$MTD_DIR/Installation/MTD_fastp.yml" \
    "$MTD_DIR/Installation/MTD.yml" \
    "$MTD_DIR/Installation/MTD_R_additions.yml" \
    "$MTD_DIR/Installation/py2.yml" \
    "$MTD_DIR/Installation/halla0820.yml" \
    "$MTD_DIR/Installation/R412.yml" \
    "$MTD_DIR/Installation/pip.requirements" \
    "$MTD_DIR/Installation/M33262_SIVMM239.fa"; do
    check_path FAIL "Installer source" "$(basename "$f")" "$f" nonempty
done
check_path WARN "Installer source" "local R tarball directory" "$MTD_DIR/update_fix/pvr_pkg" nonempty_dir
check_required_local_tarballs

for env in base MTD_fastp MTD py2 halla0820 R412; do
    check_env_exists "$env"
done

check_yaml_inventory MTD_fastp "$MTD_DIR/Installation/MTD_fastp.yml"
check_yaml_inventory MTD "$MTD_DIR/Installation/MTD.yml" "$MTD_DIR/Installation/MTD_R_additions.yml"
check_yaml_inventory py2 "$MTD_DIR/Installation/py2.yml"
check_yaml_inventory halla0820 "$MTD_DIR/Installation/halla0820.yml"
check_yaml_inventory R412 "$MTD_DIR/Installation/R412.yml"

check_conda_version WARN MTD_fastp fastp 1.3.3
check_conda_version WARN MTD r-base 4.0.3
check_conda_version WARN MTD kraken2 2.1.2
check_conda_version WARN MTD bracken 2.6.0
check_conda_version WARN MTD hisat2 2.2.1
check_conda_version WARN MTD humann 3.1.1
check_conda_version FAIL py2 python 2.7
check_conda_version FAIL halla0820 python 3.10
check_conda_version FAIL halla0820 r-base 4.1.2
check_conda_version FAIL R412 r-base 4.1.2
check_python_distribution_version FAIL halla0820 halla 0.8.20
check_python_distribution_version WARN halla0820 rpy2 3.4.5

check_env_command FAIL MTD_fastp fastp
check_env_command FAIL MTD_fastp R
check_env_command FAIL MTD_fastp Rscript

for cmd in \
    python R Rscript fastp kraken2 kraken2-build kraken2-inspect \
    bracken bracken-build hisat2 hisat2-build hisat2-inspect bowtie2 \
    samtools featureCounts makeblastdb blastdbcmd blastn magicblast \
    humann humann_config metaphlan diamond emapper.py datasets \
    STAR rsem-calculate-expression nextflow ktImportTaxonomy kraken-biom \
    mafft mash trimal prefetch fasterq-dump parallel; do
    check_env_command FAIL MTD "$cmd"
done
check_env_command_any FAIL MTD "IQ-TREE" iqtree2 iqtree
check_env_command_any FAIL MTD "FastTree" FastTree fasttree

check_env_command FAIL py2 python
check_env_command FAIL py2 hclust2.py
check_path FAIL "Tools/py2" "export2graphlan.py" "$MTD_DIR/Tools/export2graphlan/export2graphlan.py" nonempty
check_path FAIL "Tools/MTD" "graphlan.py" "$MTD_DIR/Tools/graphlan/graphlan.py" nonempty
check_path FAIL "Tools/MTD" "graphlan_annotate.py" "$MTD_DIR/Tools/graphlan/graphlan_annotate.py" nonempty
check_path FAIL "Tools/MTD" "verify_and_correct_annotations.py" "$MTD_DIR/Tools/graphlan/verify_and_correct_annotations.py" nonempty
if [[ "$MODE" != "quick" ]]; then
    run_test FAIL "Tools/py2" "export2graphlan runtime" 120 \
        "$CONDA_PATH/envs/py2/bin/python" "$MTD_DIR/Tools/export2graphlan/export2graphlan.py" --help
    run_test FAIL "Tools/MTD" "GraPhlAn renderer runtime" 120 \
        "$CONDA_PATH/envs/MTD/bin/python" "$MTD_DIR/Tools/graphlan/graphlan.py" --help
    run_test FAIL "Tools/MTD" "GraPhlAn annotator runtime" 120 \
        "$CONDA_PATH/envs/MTD/bin/python" "$MTD_DIR/Tools/graphlan/graphlan_annotate.py" --help
fi
check_env_command FAIL halla0820 python
check_env_command FAIL halla0820 Rscript
check_env_command FAIL halla0820 halla
check_env_command FAIL R412 Rscript
check_env_command FAIL R412 kraken2

check_python_imports FAIL MTD Bio biom numpy pandas scipy sklearn yaml rpy2
check_python_imports FAIL py2 numpy pandas scipy matplotlib hclust2 biom
check_python_imports FAIL halla0820 halla rpy2 numpy pandas scipy sklearn jinja2
check_pip_health MTD

# ------------------------------------------------------------------------------
# 3. R package health
# ------------------------------------------------------------------------------
MTD_R_LIST="$TMP_WORK/MTD_R_required.txt"
R412_R_LIST="$TMP_WORK/R412_R_required.txt"
HALLA_R_LIST="$TMP_WORK/HALLA_R_required.txt"
BASE_R_LIST="$TMP_WORK/BASE_R_required.txt"
FASTP_R_LIST="$TMP_WORK/MTD_fastp_R_required.txt"
MTD_R_OPTIONAL="$TMP_WORK/MTD_R_optional.txt"

extract_r_packages "$MTD_R_LIST" \
    "$MTD_DIR/update_fix/check_R_pkg.MTD.sh" \
    "$MTD_DIR/update_fix/Install.R.packages.MTD.sh"
grep -Ev '^(UCSC\.utils|rgeos)$' "$MTD_R_LIST" > "$MTD_R_LIST.tmp" 2>/dev/null || true
mv "$MTD_R_LIST.tmp" "$MTD_R_LIST"
printf '%s\n' RCurl >> "$MTD_R_LIST"
sort -fu "$MTD_R_LIST" -o "$MTD_R_LIST"
printf '%s\n' UCSC.utils rgeos > "$MTD_R_OPTIONAL"

extract_r_packages "$R412_R_LIST" \
    "$MTD_DIR/update_fix/check_R_pkg.R412.sh" \
    "$MTD_DIR/update_fix/Install.R.packages.R412_optimized.sh"
printf '%s\n' mgcv nlme survival treeio ggtree enrichplot mia cplm Maaslin2 Seurat SeuratObject sctransform flowCore cytolib cmapR tximeta hdf5r Rgraphviz pathview pbkrtest >> "$R412_R_LIST"
sort -fu "$R412_R_LIST" -o "$R412_R_LIST"

printf "%s\n" ggVennDiagram eulerr > "$FASTP_R_LIST"

cat > "$HALLA_R_LIST" <<'EOF'
lattice
MASS
mnormt
nlme
GPArotation
psych
foreign
R.methodsS3
R.oo
rtf
psychTools
XICOR
mclust
BiocManager
preprocessCore
remotes
EnvStats
Hmisc
eva
Matrix
EOF

cat > "$BASE_R_LIST" <<'EOF'
BiocManager
GenomeInfoDb
optparse
readr
httr
AnnotationForge
biomaRt
EOF

check_r_runtime MTD_fastp
check_r_version_and_bioc MTD 4.0 3.12
check_r_version_and_bioc R412 4.1 3.14
check_r_version_and_bioc halla0820 4.1 3.14

check_r_package_set FAIL MTD_fastp "$FASTP_R_LIST" "Venn/Euler packages"
check_r_package_set FAIL MTD "$MTD_R_LIST" "required packages"
check_r_package_set WARN MTD "$MTD_R_OPTIONAL" "optional compatibility packages"
check_r_package_set FAIL R412 "$R412_R_LIST" "required packages"
check_r_package_set FAIL halla0820 "$HALLA_R_LIST" "required packages"
check_r_package_set FAIL base "$BASE_R_LIST" "annotation packages in base"

# ------------------------------------------------------------------------------
# 4. Installer/helper source health
# ------------------------------------------------------------------------------
check_installer_contract
check_installer_dependency_graph

for f in \
    "$MTD_DIR/manifest.virus.sh" \
    "$MTD_DIR/manifest.bacteria.sh" \
    "$MTD_DIR/manifest.archea.sh" \
    "$MTD_DIR/manifest.plasmid.sh" \
    "$MTD_DIR/manifest.sh" \
    "$MTD_DIR/Installation/download_genomic_library.sh" \
    "$MTD_DIR/Installation/download_genomic_library_plasmid.sh" \
    "$MTD_DIR/update_fix/Install.R.packages.MTD.sh" \
    "$MTD_DIR/update_fix/check_R_pkg.MTD.sh" \
    "$MTD_DIR/update_fix/Install.R.packages.R412_optimized.sh" \
    "$MTD_DIR/update_fix/check_R_pkg.R412.sh" \
    "$MTD_DIR/update_fix/check_R_pkg.halla0820.sh" \
    "$MTD_DIR/update_fix/Install.R.AnnotPackages.base.sh"; do
    check_script_syntax "$f" bash FAIL
done

for f in \
    "$MTD_DIR/Installation/rsync_from_ncbi.pl" \
    "$MTD_DIR/Installation/rsync_from_ncbi_archaea.pl" \
    "$MTD_DIR/Installation/rsync_from_ncbi_bacteria.pl" \
    "$MTD_DIR/kraken2-build-download-taxonomy"; do
    check_script_syntax "$f" perl FAIL
done

for f in \
    "$MTD_DIR/aux_scripts/ssGSEA/resolve_ssgsea_go_terms.py" \
    "$MTD_DIR/Installation/hisat2_extract_exons.py" \
    "$MTD_DIR/Installation/hisat2_extract_splice_sites.py"; do
    check_script_syntax "$f" python FAIL
done

for f in \
    "$MTD_DIR/update_fix/check_R_pkg.MTD.sh" \
    "$MTD_DIR/update_fix/check_R_pkg.halla0820.sh" \
    "$MTD_DIR/aux_scripts/ssGSEA/resolve_ssgsea_go_terms.py" \
    "$MTD_DIR/kraken2-build-download-taxonomy"; do
    check_path FAIL "Installer source" "executable: $(basename "$f")" "$f" executable
done

check_known_source_hazards
check_helper_restoration
check_small_script_smokes

# ------------------------------------------------------------------------------
# 5. Kraken2, Bracken, and eggNOG databases
# ------------------------------------------------------------------------------
check_path FAIL "Kraken DB" "custom viral FASTA" "$MTD_DIR/viruses4kraken.fa" nonempty

# Every library explicitly requested by Install.sh is mandatory.
for lib in bacteria archaea protozoa fungi plasmid; do
    check_kraken_standard_library microbiome "$MTD_DIR/kraken2DB_micro" "$lib" yes
done
check_kraken_standard_library microbiome "$MTD_DIR/kraken2DB_micro" UniVec_Core no
check_kraken_added_library microbiome "$MTD_DIR/kraken2DB_micro" "custom viral FASTA added to the database"

check_kraken_standard_library human "$MTD_DIR/kraken2DB_human" human yes
check_path FAIL "Kraken DB" "mouse source FASTA" \
    "$MTD_DIR/kraken2DB_mice/GCF_000001635.27_GRCm39_genomic.fa" nonempty
check_path FAIL "Kraken DB" "rhesus source FASTA" \
    "$MTD_DIR/kraken2DB_rhesus/GCF_003339765.1_Mmul_10_genomic.fa" nonempty
check_kraken_added_library mouse "$MTD_DIR/kraken2DB_mice" "mouse reference genome"
check_kraken_added_library rhesus "$MTD_DIR/kraken2DB_rhesus" "rhesus reference genome"

# Compare cached manifests with files actually incorporated by the local helpers.
if [[ -n "$OFFLINE_DIR" ]]; then
    check_installed_manifest_alignment "microbiome/bacteria" \
        "$MTD_DIR/kraken2DB_micro/library/bacteria/manifest.txt" \
        "$OFFLINE_DIR/Kraken2DB_micro/library/bacteria/manifest_bacteria.list.txt" \
        "$OFFLINE_DIR/Kraken2DB_micro/library/bacteria/all" cache
    check_installed_manifest_alignment "microbiome/archaea" \
        "$MTD_DIR/kraken2DB_micro/library/archaea/manifest.txt" \
        "$OFFLINE_DIR/Kraken2DB_micro/library/archaea/manifest_archaea.list.txt" \
        "$OFFLINE_DIR/Kraken2DB_micro/library/archaea/all" cache
    check_installed_manifest_alignment "microbiome/plasmid" \
        "$MTD_DIR/kraken2DB_micro/library/plasmid/manifest.txt" \
        /dev/null \
        "$OFFLINE_DIR/Kraken2DB_micro/library/plasmid" plasmid_genomic_cache
fi

check_kraken_database microbiome "$MTD_DIR/kraken2DB_micro"
check_kraken_database human "$MTD_DIR/kraken2DB_human"
check_kraken_database mouse "$MTD_DIR/kraken2DB_mice"
check_kraken_database rhesus "$MTD_DIR/kraken2DB_rhesus"
check_bracken_database
check_eggnog_database

# ------------------------------------------------------------------------------
# 6. HUMAnN databases
# ------------------------------------------------------------------------------
check_humann

# ------------------------------------------------------------------------------
# 7. Host references, HISAT2, and BLAST
# ------------------------------------------------------------------------------
check_reference_gzip "rhesus GTF" "$MTD_DIR/ref_rhesus/Macaca_mulatta.Mmul_10.104.gtf.gz"
check_reference_gzip "human GTF" "$MTD_DIR/ref_human/Homo_sapiens.GRCh38.104.gtf.gz"
check_reference_gzip "mouse GTF" "$MTD_DIR/ref_mouse/Mus_musculus.GRCm39.104.gtf.gz"

check_hisat_index rhesus "$MTD_DIR/hisat2_index_rhesus"
check_hisat_index mouse "$MTD_DIR/hisat2_index_mouse"
check_hisat_index human "$MTD_DIR/hisat2_index_human"

check_blast_database rhesus "$MTD_DIR/rhesus_blastdb/rhesus_blastdb"
check_blast_database mouse "$MTD_DIR/mouse_blastdb/mouse_blastdb"
check_blast_database human "$MTD_DIR/human_blastdb/human_blastdb"

# ------------------------------------------------------------------------------
# 8. Optional offline source folder
# ------------------------------------------------------------------------------
check_installation_cache

# ------------------------------------------------------------------------------
# Final report
# ------------------------------------------------------------------------------
{
    echo "============================================================"
    echo "MTD installation check summary v$CHECKER_VERSION"
    echo "============================================================"
    echo "Date:        $(date '+%F %T %Z')"
    echo "Host:        $(hostname)"
    echo "MTD:         $MTD_DIR"
    echo "Conda:       $CONDA_PATH"
    echo "Mode:        $MODE"
    echo "Cache:       ${OFFLINE_DIR:-not detected}"
    echo "Dependencies: $DEPENDENCY_TSV"
    echo "Checks:      $CHECK_COUNT"
    echo "PASS:        $PASS_COUNT"
    echo "WARN:        $WARN_COUNT"
    echo "FAIL:        $FAIL_COUNT"
    echo "SKIP:        $SKIP_COUNT"
    echo "Results TSV: $RESULTS_TSV"
    echo "Full log:    $FULL_LOG"
    echo "============================================================"
} | tee "$SUMMARY_TXT"

if [[ "$FAIL_COUNT" -gt 0 ]]; then
    printf 'FINAL STATUS: FAIL\n' >> "$SUMMARY_TXT"
    printf '%s%sFINAL STATUS: FAIL%s\n' "$COLOR_BOLD" "$COLOR_FAIL" "$COLOR_RESET"
elif [[ "$WARN_COUNT" -gt 0 ]]; then
    printf 'FINAL STATUS: PASS WITH WARNINGS\n' >> "$SUMMARY_TXT"
    printf '%s%sFINAL STATUS: PASS WITH WARNINGS%s\n' "$COLOR_BOLD" "$COLOR_WARN" "$COLOR_RESET"
else
    printf 'FINAL STATUS: PASS\n' >> "$SUMMARY_TXT"
    printf '%s%sFINAL STATUS: PASS%s\n' "$COLOR_BOLD" "$COLOR_PASS" "$COLOR_RESET"
fi

if [[ "$FAIL_COUNT" -gt 0 ]]; then
    exit 1
fi
if [[ "$STRICT" -eq 1 && "$WARN_COUNT" -gt 0 ]]; then
    exit 2
fi
exit 0
