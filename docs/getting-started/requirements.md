# System requirements

This page describes the operating system, hardware, storage, and network requirements for installing and running MTD Explorer.

!!! warning "Requirements under validation"

```
The minimum and recommended hardware requirements are currently being validated using clean installations and benchmark runs on different computer configurations.
```

## Supported operating systems

MTD Explorer is developed and tested on 64-bit GNU/Linux systems.

### Currently tested

* Linux Mint
* Ubuntu-based Linux distributions

### Not currently supported

* Native Microsoft Windows
* Native macOS

Windows users may be able to use a Linux server, virtual machine, or WSL2, but these configurations have not yet been formally validated.

## Hardware requirements

| Resource            |          Minimum |      Recommended | Notes                                                              |
| ------------------- | ---------------: | ---------------: | ------------------------------------------------------------------ |
| CPU architecture    |           x86-64 |           x86-64 | ARM64 is not currently validated                                   |
| CPU threads         | Under validation |       16 or more | More threads reduce installation and analysis time                 |
| RAM                 | Under validation | Under validation | Database construction is one of the most memory-intensive steps    |
| Free storage        |  1TB             | 2TB              | Depends on databases, cache, host references, and analysis outputs |
| Internet connection |         Required | Stable broadband | Large reference databases are downloaded during installation       |

!!! note "Minimum versus recommended"

```
A computer may be able to run MTD Explorer with less memory or fewer CPU threads, but installation and analysis may take considerably longer.
```

## Storage planning

Storage should be available for:

1. MTD Explorer source code;
2. Conda environments;
3. Kraken2 databases;
4. host reference genomes and indexes;
5. eggNOG and HUMAnN databases;
6. the optional reusable installation cache;
7. input FASTQ files;
8. analysis outputs and temporary files.

The reusable installation cache may be stored on a separate disk.

## Required permissions

The user should have:

* permission to write inside the MTD Explorer directory;
* permission to write to the selected cache directory;
* permission to create files in the selected output directory;
* `sudo` access when system packages need to be installed.

MTD Explorer itself should not normally be executed as the root user.

## Check your computer

Use the following commands before installation.

### Operating system

```bash
cat /etc/os-release
```

### CPU architecture

```bash
uname -m
```

### Available CPU threads

```bash
nproc
```

### Available memory

```bash
free -h
```

### Available disk space

```bash
df -h
```

### Git availability

```bash
git --version
```

## Next step

Continue to the [installation guide](installation.md).

