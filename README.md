# MIDORI BLAST Parser

Python tools and notebooks for downloading and parsing MIDORI reference databases, with support for selecting a specific or the latest GenBank release and exporting sequence tables for a selected taxonomic class.

## Overview

This repository contains a small workflow for working with the [MIDORI reference database](https://www.reference-midori.info/download.php).

The main script:

- parses the MIDORI download page
- finds `.fasta.zip` archives in `longest/fasta` and `uniq/fasta`
- selects either a user-defined GenBank release or the latest available release
- extracts FASTA records directly from ZIP archives
- filters sequences by taxonomic class
- saves final tabular datasets with species name, gene, and sequence

It is especially useful for preparing class-specific datasets, for example for **Aves**.

---

## Repository structure

```text
.
├── data/
│   ├── longest_fasta/          # downloaded MIDORI longest FASTA ZIP archives
│   ├── uniq_fasta/             # downloaded MIDORI uniq FASTA ZIP archives
│   ├── longest.csv        # example output table
│   └── uniq.csv           # example output table
│
├── script/
│   └── midori_blast_parse.py   # main MIDORI parsing script
│
├── notebook/
│   └── *.ipynb                 # notebooks for testing, exploration, and documentation
│
└── README.md
