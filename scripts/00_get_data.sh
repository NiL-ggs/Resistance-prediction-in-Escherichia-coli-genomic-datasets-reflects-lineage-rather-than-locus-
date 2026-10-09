#!/usr/bin/env bash
# Download the two external inputs. Neither is redistributed in this repository.
set -euo pipefail
mkdir -p giessen reference

# 1. Giessen SNP matrix and phenotypes (Ren et al. 2022)
curl -L -o giessen/Giessen_dataset.zip \
  https://raw.githubusercontent.com/YunxiaoRen/ML-iAMR/main/Giessen_dataset.zip
unzip -o -q giessen/Giessen_dataset.zip -d giessen/
#   -> giessen/Giessen_dataset/cip_ctx_ctz_gen_multi_data.csv   (99 MB)
#   -> giessen/Giessen_dataset/cip_ctx_ctz_gen_pheno.csv

# 2. E. coli K-12 MG1655 reference annotation (GenBank, ~12 MB)
curl -L -o reference/NC_000913.3.gb \
  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000913.3&rettype=gbwithparts&retmode=text"

echo "Checksums of downloaded inputs:"
md5sum giessen/Giessen_dataset/*.csv reference/NC_000913.3.gb
