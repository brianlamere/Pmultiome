#!/bin/bash
# links.sh - create input directory structure and symlinks for remaining samples

# LG26 (ATPSY 13, S8/S1)
mkdir -p input/LG26/GEX input/LG26/ATAC
ln -s ../../../originfiles/13-3GEX_S8_L002_R1_001.fastq.gz  input/LG26/GEX/LG26_GEX_S8_L002_R1_001.fastq.gz
ln -s ../../../originfiles/13-3GEX_S8_L002_R2_001.fastq.gz  input/LG26/GEX/LG26_GEX_S8_L002_R2_001.fastq.gz
ln -s ../../../originfiles/13-ATAC_S1_L004_R1_001.fastq.gz  input/LG26/ATAC/LG26_ATAC_S1_L004_R1_001.fastq.gz
ln -s ../../../originfiles/13-ATAC_S1_L004_R2_001.fastq.gz  input/LG26/ATAC/LG26_ATAC_S1_L004_R2_001.fastq.gz
ln -s ../../../originfiles/13-ATAC_S1_L004_R3_001.fastq.gz  input/LG26/ATAC/LG26_ATAC_S1_L004_R3_001.fastq.gz

# LG31 (ATPSY 14, S9/S2)
mkdir -p input/LG31/GEX input/LG31/ATAC
ln -s ../../../originfiles/14-3GEX_S9_L002_R1_001.fastq.gz  input/LG31/GEX/LG31_GEX_S9_L002_R1_001.fastq.gz
ln -s ../../../originfiles/14-3GEX_S9_L002_R2_001.fastq.gz  input/LG31/GEX/LG31_GEX_S9_L002_R2_001.fastq.gz
ln -s ../../../originfiles/14-ATAC_S2_L004_R1_001.fastq.gz  input/LG31/ATAC/LG31_ATAC_S2_L004_R1_001.fastq.gz
ln -s ../../../originfiles/14-ATAC_S2_L004_R2_001.fastq.gz  input/LG31/ATAC/LG31_ATAC_S2_L004_R2_001.fastq.gz
ln -s ../../../originfiles/14-ATAC_S2_L004_R3_001.fastq.gz  input/LG31/ATAC/LG31_ATAC_S2_L004_R3_001.fastq.gz

# LG22 (ATPSY 17, S10/S3)
mkdir -p input/LG22/GEX input/LG22/ATAC
ln -s ../../../originfiles/17-3GEX_S10_L002_R1_001.fastq.gz input/LG22/GEX/LG22_GEX_S10_L002_R1_001.fastq.gz
ln -s ../../../originfiles/17-3GEX_S10_L002_R2_001.fastq.gz input/LG22/GEX/LG22_GEX_S10_L002_R2_001.fastq.gz
ln -s ../../../originfiles/17-ATAC_S3_L004_R1_001.fastq.gz  input/LG22/ATAC/LG22_ATAC_S3_L004_R1_001.fastq.gz
ln -s ../../../originfiles/17-ATAC_S3_L004_R2_001.fastq.gz  input/LG22/ATAC/LG22_ATAC_S3_L004_R2_001.fastq.gz
ln -s ../../../originfiles/17-ATAC_S3_L004_R3_001.fastq.gz  input/LG22/ATAC/LG22_ATAC_S3_L004_R3_001.fastq.gz

# LG38 (ATPSY 20, S11/S4)
mkdir -p input/LG38/GEX input/LG38/ATAC
ln -s ../../../originfiles/20-3GEX_S11_L002_R1_001.fastq.gz input/LG38/GEX/LG38_GEX_S11_L002_R1_001.fastq.gz
ln -s ../../../originfiles/20-3GEX_S11_L002_R2_001.fastq.gz input/LG38/GEX/LG38_GEX_S11_L002_R2_001.fastq.gz
ln -s ../../../originfiles/20-ATAC_S4_L004_R1_001.fastq.gz  input/LG38/ATAC/LG38_ATAC_S4_L004_R1_001.fastq.gz
ln -s ../../../originfiles/20-ATAC_S4_L004_R2_001.fastq.gz  input/LG38/ATAC/LG38_ATAC_S4_L004_R2_001.fastq.gz
ln -s ../../../originfiles/20-ATAC_S4_L004_R3_001.fastq.gz  input/LG38/ATAC/LG38_ATAC_S4_L004_R3_001.fastq.gz

# LG300 (ATPSY 23, S12/S5)
mkdir -p input/LG300/GEX input/LG300/ATAC
ln -s ../../../originfiles/23-3GEX_S12_L002_R1_001.fastq.gz input/LG300/GEX/LG300_GEX_S12_L002_R1_001.fastq.gz
ln -s ../../../originfiles/23-3GEX_S12_L002_R2_001.fastq.gz input/LG300/GEX/LG300_GEX_S12_L002_R2_001.fastq.gz
ln -s ../../../originfiles/23-ATAC_S5_L004_R1_001.fastq.gz  input/LG300/ATAC/LG300_ATAC_S5_L004_R1_001.fastq.gz
ln -s ../../../originfiles/23-ATAC_S5_L004_R2_001.fastq.gz  input/LG300/ATAC/LG300_ATAC_S5_L004_R2_001.fastq.gz
ln -s ../../../originfiles/23-ATAC_S5_L004_R3_001.fastq.gz  input/LG300/ATAC/LG300_ATAC_S5_L004_R3_001.fastq.gz

# LG301 (ATPSY 25, S13/S6)
mkdir -p input/LG301/GEX input/LG301/ATAC
ln -s ../../../originfiles/25-3GEX_S13_L002_R1_001.fastq.gz input/LG301/GEX/LG301_GEX_S13_L002_R1_001.fastq.gz
ln -s ../../../originfiles/25-3GEX_S13_L002_R2_001.fastq.gz input/LG301/GEX/LG301_GEX_S13_L002_R2_001.fastq.gz
ln -s ../../../originfiles/25-ATAC_S6_L004_R1_001.fastq.gz  input/LG301/ATAC/LG301_ATAC_S6_L004_R1_001.fastq.gz
ln -s ../../../originfiles/25-ATAC_S6_L004_R2_001.fastq.gz  input/LG301/ATAC/LG301_ATAC_S6_L004_R2_001.fastq.gz
ln -s ../../../originfiles/25-ATAC_S6_L004_R3_001.fastq.gz  input/LG301/ATAC/LG301_ATAC_S6_L004_R3_001.fastq.gz
