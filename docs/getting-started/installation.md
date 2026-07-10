## Install Git

Git is required to download the MTD Explorer source code and to keep the
installation updated.

Check whether Git is already installed:

```bash
git --version
```

When Git is available, the command displays its installed version.

If Git is not installed, install it using the package manager of your Linux
distribution.

### Ubuntu and Linux Mint

```bash
sudo apt update
sudo apt install -y git
```

Confirm the installation:

```bash
git --version
```

## Clone the MTD Explorer repository

After installing Git, download the MTD Explorer repository.

### Using HTTPS

```bash
cd ~
git clone https://github.com/patrick-douglas/MTD.git
cd MTD
```

HTTPS is the recommended method for users who have not configured an SSH
key for GitHub.

Confirm that the repository was downloaded correctly:

```bash
git status
git log -1 --oneline
```

The installation commands in the following sections should be executed
from inside the MTD Explorer repository directory.

## Installer command

The MTD Explorer installer is executed through `Install.sh`.

Display the current command-line help with:

```bash
bash Install.sh -h
```

The general syntax is:

```text
Install.sh -o <installation-cache> [options]
```

## Required option

### Installation cache

The `-o` option specifies the persistent installation cache:

```text
-o PATH
```

This option is required.

The directory is created and populated automatically when it does not already exist. Files already present in a valid cache may be reused during subsequent or interrupted installations.

Example:

```bash
bash Install.sh \
  -o /home/user/MTD_install_cache
```

The cache may also be located on a separate mounted disk:

```bash
bash Install.sh \
  -o /path/to/large/storage/MTD_install_cache
```

!!! important "Use an absolute path"

````
An absolute cache path is recommended.

Correct:

```text
/home/user/MTD_install_cache
```

Correct:

```text
/media/user/storage/MTD_install_cache
```

Incorrect:

```text
/home/user/~/MTD_install_cache
```

The `~` character is expanded only when it appears at the beginning of
a shell path.
````

## Installer options

| Option | Argument | Description                                                          |
| ------ | -------- | -------------------------------------------------------------------- |
| `-o`   | `PATH`   | Persistent installation cache. This option is required.              |
| `-p`   | `PATH`   | Miniconda installation directory. The default is `$HOME/miniconda3`. |
| `-k`   | `INT`    | Kraken2 k-mer length used when building databases.                   |
| `-m`   | `INT`    | Kraken2 minimizer length used when building databases.               |
| `-s`   | `INT`    | Kraken2 minimizer-spaces value used when building databases.         |
| `-r`   | `INT`    | Bracken read length. The default is `75`.                            |
| `-h`   | —        | Display the installer help message and exit.                         |

!!! note "Kraken2 database parameters"

```
Most users should leave `-k`, `-m`, and `-s` unset.

These parameters alter how Kraken2 databases are built and should be
changed only when there is a specific technical reason to use
non-default database settings.
```

## Standard installation

For a standard installation using the default Miniconda location and a
Bracken read length of 75:

```bash
cd ~/MTD

bash Install.sh \
  -o /home/user/MTD_install_cache
```

Replace `/home/user/MTD_install_cache` with the desired cache location.

## Custom Miniconda location

By default, Miniconda is installed at:

```text
$HOME/miniconda3
```

A different location can be selected with `-p`:

```bash
bash Install.sh \
  -o /path/to/MTD_install_cache \
  -p /path/to/miniconda3
```

The selected parent directory must be writable by the current user.

## Existing Miniconda installation

The installer automatically downloads and installs Miniconda.

If the selected Miniconda directory already exists, the installer displays
a warning and requests explicit confirmation before permanently deleting
that directory.

!!! danger "Existing Miniconda directory"

```
Confirming the removal of an existing Miniconda directory permanently
deletes the environments and packages stored inside it.

Before confirming:

1. verify that the displayed path is correct;
2. confirm that no important Conda environments are stored there;
3. back up anything that must be preserved.

Do not confirm the deletion when the installer displays an unexpected
directory.
```

## Bracken read length

The default Bracken read length is:

```text
75
```

A different value can be selected using `-r`:

```bash
bash Install.sh \
  -o /path/to/MTD_install_cache \
  -r 100
```

The selected value should correspond to the read length intended for
Bracken abundance estimation.

Changing this option affects the Bracken files generated during database
preparation.

## Save the installation log

The complete installer output should be saved, especially during clean
installations or tests on a new computer.

```bash
bash Install.sh \
  -o /path/to/MTD_install_cache \
  2>&1 | tee MTD_installation.log
```

This displays the installation output in the terminal while also saving it
to:

```text
MTD_installation.log
```

For a timestamped log:

```bash
bash Install.sh \
  -o /path/to/MTD_install_cache \
  2>&1 | tee "MTD_installation_$(date +%Y%m%d_%H%M%S).log"
```

## During installation

The installer may display status messages such as:

| Status    | Meaning                                              |
| --------- | ---------------------------------------------------- |
| `INFO`    | Describes the current installation step              |
| `PASS`    | A required step or validation completed successfully |
| `WARNING` | A non-fatal condition requires attention             |
| `RETRY`   | A failed operation is being attempted again          |
| `ERROR`   | A required operation failed                          |

Large reference databases may require considerable download and processing
time.

The installation cache should not be deleted after a successful
installation. It can be reused when:

* reinstalling MTD Explorer;
* recovering from an interrupted installation;
* installing the pipeline on another computer;
* validating downloaded reference files;
* avoiding repeated large downloads.

## After installation

Do not begin a real analysis solely because `Install.sh` reached the end.

First, run the dedicated installation verification procedure described in
[Verify the installation](verify-installation.md).

