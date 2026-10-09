#!/usr/bin/env bash
# Maximum-likelihood phylogeny from the SNP alignment written by 02_prepare_matrix.py.
# Runtime is roughly 25 minutes single-threaded for 809 taxa x 16,435 sites.
set -euo pipefail
FastTree -nt -gtr -gamma -nosupport -log ftlog.txt snps.fa > tree.nwk
