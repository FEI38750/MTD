# Verify the installation

MTD Explorer includes a dedicated verification program that checks whether
the software environments, commands, packages, reference files, indexes,
and databases required by the pipeline are available and usable.

The verification program is:

```text
MTD_check_installation.sh
```

!!! important

```
Finishing `Install.sh` does not by itself guarantee that every
environment, package, index, and reference database is ready for use.

Run the installation checker before analyzing a real dataset.
```

## Display the checker help

From the MTD Explorer directory:

```bash
cd ~/MTD
bash MTD_check_installation.sh --help
```

The general syntax is:

```text
MTD_check_installation.sh [options]
```

## Basic verification

The checker can usually determine the MTD Explorer directory, Conda
installation, and persistent cache paths automatically.

Run the default full verification with:

```bash
cd ~/MTD
bash MTD_check_installation.sh
```

The default verification mode is:

```text
full
```

## Verification modes

MTD Explorer provides three verification levels.

| Mode    | Description                                                                                                                 |
| ------- | --------------------------------------------------------------------------------------------------------------------------- |
| `quick` | Checks directory structure, package inventories, command availability, and script syntax                                    |
| `full`  | Performs all quick checks and additionally tests R package loading and inspects HISAT2, BLAST, HUMAnN, and Kraken resources |
| `deep`  | Performs all full checks and additionally runs small Kraken2 classifications and validates all cached gzip files            |

### Quick mode

Use quick mode for a fast structural check:

```bash
bash MTD_check_installation.sh --mode quick
```

This mode is useful after:

* updating scripts;
* changing file permissions;
* modifying software environments;
* transferring an existing installation;
* checking whether the expected commands are available.

### Full mode

Full mode is the default and is recommended after a normal installation:

```bash
bash MTD_check_installation.sh --mode full
```

The following two commands are therefore equivalent:

```bash
bash MTD_check_installation.sh
```

```bash
bash MTD_check_installation.sh --mode full
```

### Deep mode

Deep mode performs the most comprehensive validation:

```bash
bash MTD_check_installation.sh --mode deep
```

In addition to the full checks, deep mode performs small Kraken2
classification tests and validates all gzip files in the persistent
installation cache.

!!! note "Deep verification"

```
Deep mode may take considerably longer than quick or full mode,
particularly when the installation cache contains many large
compressed files.

It is especially useful after:

- a clean installation;
- an interrupted database download;
- moving the cache to another disk;
- installing MTD Explorer on another computer;
- suspected cache corruption.
```

## Checker options

| Option                | Argument | Description                                                        |
| --------------------- | -------- | ------------------------------------------------------------------ |
| `-m`, `--mtd-dir`     | `PATH`   | MTD Explorer installation directory                                |
| `-p`, `--conda-path`  | `PATH`   | Conda installation directory                                       |
| `-o`, `--offline-dir` | `PATH`   | Persistent installation cache                                      |
| `-r`, `--read-length` | `INT`    | Bracken read length; default: `75`                                 |
| `--mode`              | `MODE`   | Verification mode: `quick`, `full`, or `deep`                      |
| `--report-dir`        | `PATH`   | Directory where verification reports are written                   |
| `--strict`            | —        | Return status `2` when warnings exist but no failures are detected |
| `--keep-temp`         | —        | Preserve temporary test files inside the report directory          |
| `-h`, `--help`        | —        | Display the checker help and exit                                  |
| `--version`           | —        | Display the checker version and exit                               |

## Automatic path detection

### MTD Explorer directory

When `--mtd-dir` is not provided, the checker uses the directory containing
`MTD_check_installation.sh`.

For a standard installation, this is normally sufficient:

```bash
cd ~/MTD
bash MTD_check_installation.sh
```

A different installation directory can be specified explicitly:

```bash
bash MTD_check_installation.sh \
  --mtd-dir /path/to/MTD
```

### Conda installation

When `--conda-path` is not provided, the checker searches for the Conda
installation in the following order:

1. the path recorded in `MTD/condaPath`;
2. `$HOME/miniconda3`.

A path can also be supplied explicitly:

```bash
bash MTD_check_installation.sh \
  --conda-path /home/user/miniconda3
```

### Persistent installation cache

When `--offline-dir` is not provided, the checker uses the path recorded in:

```text
MTD/offlineCachePath
```

The cache can be specified manually:

```bash
bash MTD_check_installation.sh \
  --offline-dir /path/to/MTD_install_cache
```

An explicit cache path is particularly useful when:

* the cache was moved;
* the installation directory was copied from another computer;
* more than one cache is available;
* `offlineCachePath` is missing or outdated.

## Complete explicit command

All primary paths can be supplied explicitly:

```bash
bash MTD_check_installation.sh \
  --mtd-dir /home/user/MTD \
  --conda-path /home/user/miniconda3 \
  --offline-dir /path/to/MTD_install_cache \
  --read-length 75 \
  --mode full
```

Replace the example paths with the paths used by the current installation.

## Bracken read length

The default Bracken read length checked by the program is:

```text
75
```

A different value can be selected with:

```bash
bash MTD_check_installation.sh \
  --read-length 100
```

The value should match the Bracken read length used during installation and
database preparation.

For example, when the installer was run with:

```bash
bash Install.sh \
  -o /path/to/MTD_install_cache \
  -r 100
```

the checker should also use:

```bash
bash MTD_check_installation.sh \
  --read-length 100
```

## Report directory

A specific output directory can be selected with `--report-dir`:

```bash
bash MTD_check_installation.sh \
  --mode full \
  --report-dir ./MTD_installation_check
```

Using a dedicated report directory is recommended when:

* testing a clean installation;
* comparing different computers;
* preparing benchmark records;
* reporting an installation problem;
* preserving results from multiple checker runs.

A timestamped report directory can be created with:

```bash
bash MTD_check_installation.sh \
  --mode full \
  --report-dir "MTD_check_$(date +%Y%m%d_%H%M%S)"
```

## Preserve temporary test files

Temporary files created during verification are normally removed.

Use `--keep-temp` to preserve them inside the report directory:

```bash
bash MTD_check_installation.sh \
  --mode deep \
  --report-dir ./MTD_deep_check \
  --keep-temp
```

This option is mainly useful for debugging failed tests.

## Strict mode

By default, warnings do not necessarily produce a failing process status.

With `--strict`, the checker returns exit status `2` when warnings are
present but no failures are detected:

```bash
bash MTD_check_installation.sh \
  --mode full \
  --strict
```

Strict mode is useful for:

* automated validation;
* continuous integration;
* installation benchmarking;
* detecting warnings in shell scripts;
* requiring a completely clean verification report.

The exit status can be inspected immediately after the checker finishes:

```bash
echo $?
```

## Recommended validation workflow

After a clean installation, the following sequence is recommended.

### Step 1: Quick structural check

```bash
bash MTD_check_installation.sh \
  --mode quick \
  --report-dir ./MTD_check_quick
```

### Step 2: Full functional check

```bash
bash MTD_check_installation.sh \
  --mode full \
  --report-dir ./MTD_check_full
```

### Step 3: Deep cache and Kraken2 check

```bash
bash MTD_check_installation.sh \
  --mode deep \
  --report-dir ./MTD_check_deep
```

Running the three modes separately makes it easier to identify whether a
problem involves:

* directory structure or syntax;
* software environments or packages;
* reference databases or compressed cache files;
* Kraken2 database functionality.

## Save the terminal output

The checker output can also be saved with `tee`:

```bash
bash MTD_check_installation.sh \
  --mode full \
  --report-dir ./MTD_check_full \
  2>&1 | tee MTD_check_full.log
```

For a timestamped log:

```bash
bash MTD_check_installation.sh \
  --mode full \
  --report-dir "MTD_check_$(date +%Y%m%d_%H%M%S)" \
  2>&1 | tee "MTD_check_$(date +%Y%m%d_%H%M%S).log"
```

## Understanding the results

The checker reports conditions using statuses such as:

| Status | Meaning                                                           |
| ------ | ----------------------------------------------------------------- |
| `PASS` | The component was found and passed the corresponding check        |
| `WARN` | The component may be usable, but the condition requires attention |
| `FAIL` | A required component is missing, invalid, incomplete, or unusable |

A warning should be reviewed before running a real analysis.

A failure should be resolved before using MTD Explorer.

!!! warning "Do not evaluate only the final line"

```
Review the complete checker output and generated reports.

An installation can contain several valid components while still
having an incomplete environment, database, index, or reference file.
```

## Information to preserve

When reporting an installation problem, preserve:

* the complete installer log;
* the checker terminal output;
* the checker report directory;
* temporary test files when relevant;
* the MTD Explorer Git commit;
* the current Git working-tree state;
* operating-system information;
* hardware information;
* the Conda installation path;
* the persistent cache path;
* the Bracken read length;
* the verification mode used.

Record the current commit with:

```bash
git log -1 --oneline
```

Record the repository state with:

```bash
git status
```

Record basic system information with:

```bash
cat /etc/os-release
uname -a
nproc
free -h
df -h
```

## Next step

After the required checks pass, continue to the
[Quick start guide](quick-start.md).

