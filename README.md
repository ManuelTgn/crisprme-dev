# CRISPRme2

High-performance and scalable tool for genome-wide off-target assessment in CRISPR-Cas systems. It supports variant-aware and haplotype-aware predictions, integrating SNVs, indels, and population-scale haplotypes with orthogonal genomic annotations to prioritize off-targets across personal and population genomes

## Table of Contents

0 [System Requirements](#0-system-requirements)
<br>1 [Installation](#1-installation)
<br>&nbsp;&nbsp;1.1 [Install CRISPRme2 from Mamba/Conda](#11-install-crisprme2-from-mambaconda)
<br>&nbsp;&nbsp;&nbsp;&nbsp;1.1.1 [Install Conda or Mamba](#111-install-conda-or-mamba)
<br>&nbsp;&nbsp;&nbsp;&nbsp;1.1.2 [Install CRISPRme2](#112-install-crisprme2)
<br>&nbsp;&nbsp;1.2 [Install CRISPRme2 from Docker](#12-install-crisprme2-from-docker)
<br>&nbsp;&nbsp;1.3 [Install CRISPRme2 from PyPI](#13-install-crisprme2-from-pypi)
<br>&nbsp;&nbsp;1.4 [Install CRISPRme2 from Source Code](#14-install-crisprme2-from-source-code)
<br>2 [Usage](#2-usage)
<br>&nbsp;&nbsp;2.1 [General Syntax](#21-general-syntax)
<br>&nbsp;&nbsp;2.2 [Complete Search](#22-complete-search)
<br>3 [Test](#3-test)
<br>4 [Citation](#4-citation)
<br>5 [Contacts](#5-contacts)
<br>6 [License](#6-license)

## 0 System Requirements

To ensure optimal performance, CRISPRme2 requires the following system specifications:

- **Operating System**:
<br>macOS or any modern Linux distribution (e.g., Ubuntu, CentOS)

- **Required Disk Space**:
<br> 100 MB

- **Minimum RAM**:
<br>32 GB — sufficient for standard use cases and small to medium-sized datasets

- **Recommended RAM for Large-Scale Analyses**:
<br>128 GB or more — recommended for memory-intensive tasks such as:

    - Processing large-scale variant datasets (e.g., gnomAD data, All of US data, etc.)

> ⚠️ **Note**: For optimal performance and stability, especially when dealing with large-scale variant datasets, ensure that your system meets or exceeds the recommended specifications.

## 1 Installation

This section provides step-by-step instructions to install CRISPR-HAWK and external dependencies. Choose the method that best suits your environment and preferences:

- **[Install CRISPRme2 from Mamba/conda](#11-install-crisprme2-from-mambaconda)** (recommended)
<br>Best for users seeking an isolated and reproducible environment with minimal manual dependency handling.

- **[Install CRISPRme2 from Docker](#12-install-crisprme2-from-docker)**
<br>Ideal for users who prefer containerized deployments or want to avoid configuring the environment manually.

- **[Install CRISPRme2 from PyPI](#13-install-crisprme2-from-pypi)**
<br>Quick option for Python users already working within a virtual environment. May require manual handling of some dependencies.

- **[Install CRISPRme2 from source code](#14-install-crisprme2-from-source-code)**
<br>Suitable for developers or contributors who want full control over the codebase or plan to customize CRISPRme2.

> ⚠️ **Note:** We recommend using the Mamba/Conda or Docker installation methods for most users, as they ensure the highest compatibility and stability across systems.

### 1.1 Install CRISPRme2 from Mamba/Conda

#### 1.1.1 Install Conda or Mamba

Before installing CRISPRme2, ensure that either **Conda** or **Mamba** is installed on your system. Based on recommendations from the Bioconda community and performance testing during CRISPR-HAWK development, we **recommend using [Mamba](https://mamba.readthedocs.io/en/latest/index.html)** over Conda. Mamba is a fast, efficient drop-in replacement for Conda, built with a high-performance dependency solver in C++.

**Installation Steps**

**1. Install Conda or Mamba**

* To install **Conda**, follow the official instructions:
<br>[Conda Installation Guide](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)

* To install **Mamba**, follow the official instructions:
<br>[Mamba Installation Guide](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html)

**2. Configure Bioconda Channels**
Once Mamba (or Conda) is installed, configure your environment with the appropriate channels used by CRISPR-HAWK:
```bash
mamba config --add channels bioconda
mamba config --add channels defaults
mamba config --add channels conda-forge
mamba config --set channel_priority strict
```

> 💡 **Tip**: If you are using Conda instead of Mamba, simply replace `mamba` with `conda` in the commands above.

By completing these steps your system will be correctly configured to install CRISPR-HAWK and all required dependencies via Bioconda.

**Apple Silicon (M1/M2/M3) Support**

If you're using a Mac with Apple Silicon, follow these additional steps to ensure compatibility with Bioconda packages (which are primarily built for Intel)

> 💡 **Tip**: Not sure if your Mac uses Apple Silicon (M1, M2, or M3)? You can check by visiting Apple’s official support page: [Identify your Mac model and chip](https://support.apple.com/en-us/116943)

**System-wide (Recommended)**

Make sure [Rosetta](https://support.apple.com/en-us/102527) is installed:
```zsh
softwareupdate --install-rosetta
```

Configure Mamba (or Conda) to prefer Intel (x86_64) builds:
```zsh
mamba config --add subdirs osx-64
```

This will allow Bioconda to fetch compatible packages globally across all environments.

**Environment-specific (Alternative)**

You can also enable Intel compatibility in a specific environment only:
```zsh
CONDA_SUBDIR=osx-64 mamba create -n crisprme2-dev -c bioconda crisprme2
```

> ⚠️ **Note**: If you use this method, remember to prepend CONDA_SUBDIR=osx-64 to every future conda install command within this environment — or set the variable globally in your shell profile.


#### 1.1.2 Install CRISPRme2

TBA

### 1.2 Install CRISPRme2 from Docker

TBA

### 1.3 Install CRISPRme2 from PyPI

TBA

### 1.4 Install CRISPRme2 from Source Code

Installing CRISPRme2 from source is ideal for developers, contributors, or users who wish to inspect or customize the codebase.

This method assumes you already have **Python 3.8** installed and accessible from your system’s environment.

**Prerequisites**

- Python **3.8** (strictly required)

- `git`

- A virtual environment (optional but recommended)

**Installation Steps**

**1. Clone the Repository**
```bash
git clone https://github.com/ManuelTgn/CRISPRme2.git
cd CRISPRme2
```

**2. (Optional) Create and Activate a Virtual Environment**
```bash
mamba create -n crisprme2-env python=3.8 -y
mamba activate crisprme2-env
```

**3. Install CRISPR-HAWK and Its Dependencies**
```bash
pip install .  # regular installation
pip install -e .  # development-mode installation
```

The `.` tells `pip` to install the current directory as a package, including all dependencies specified in `setup.py` or `pyproject.toml`.

**Quick Test**

Once installation is complete, verify that the command-line interface is working:
```bash
crisprme2 -h
```

If the help message is displayed correctly, CRISPR-HAWK is successfully installed and callable from any directory in your system.

## 2 Usage

CRISPRme2 provides multiple functionalities designed to support variant- and haplotype-aware CRISPR off-targets nomination, off-targets risk evaluation, and integration with downstream analysis pipelines. Each command serves a distinct role in the workflow.

### 2.1 General Syntax
```bash
crisprme2 <command> [options]
```

To view available commands:
```bash
crisprme2 --help
```

To check version:
```bash
crisprme2 --version
```

### 2.2 Complete Search

The `crisprme2 complete-search` command is the core functionality of CRISPRme2, designed to identify and annotate candidate off-targets in both reference and variant genomes.
It integrates variant-aware search, functional annotation, and predictive scoring to help you prioritize the most robust and context-aware guides for CRISPR editing.

The search includes:

* Support for any Cas system (Cas9, Cpf1, SaCas9, etc.)
* Compatibility with custom PAM sequences and guide lengths
* Variant-aware design from individual or population-level VCF files (SNVs and indels)
* Scoring using **CFD**, **CRISTA**, **CRISPR-bulge**, and **Elevation**
* Functional and gene annotation using user-specified BED files 
* Output in detailed and structured reports (off-target tables, graphical reports)

Usage:
```bash
crisprme2 search --genome <genome-dir> --vcf <vcf-dir> --outdir <output-dir>
```

> ⚠️ **Note**: All FASTA files in `<genome-dir>` must be one per chromosome (e.g., chr1.fa, chr2.fa, etc.).

---

#### Required Arguments

#### Optional Arguments

## 3 Test

## 4 Citation

If you use CRISPRme2 in your research, please cite:

```bibtex
@software{crisprme22025,
  title = {CRISPRme2: High-performance and scalable tool for variant- and haplotype-aware genome-wide off-target assessment in CRISPR-Cas systems},
  author = {Manuel Tognon},
  year = {2025},
  url = {https://github.com/ManuelTgn/CRISPRme2}
}
```

## 5 Contacts

* Manuel Tognon
  <br>manuel.tognon@univr.it

* Rosalba Giugno
  <br>rosalba.giugno@univr.it

* Luca Pinello
  <br>lpinello@mgh.harvard.edu

## 6 License

CRISPRme2 is licensed under the AGPL-3.0 license, which permits its use for academic research purposes only.

For any commercial or for-profit use, please contact the authors.

