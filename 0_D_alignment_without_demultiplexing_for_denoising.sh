# NOTE (2026-11-06): Subject labels H and 594 were found to be
# swapped. All labels in this repository, the config CSVs, the FASTQ file names, and the data folders have
# been corrected, so the labels below refer to the correct samples. The original sequencing library
# IDs (CKDL...) stay attached to their data. Corrected files are marked with an r.
#
# Previously, I ran cellranger multi (to combine gex and plx) for each *sample* individually
# (executed in previous notebooks, 0_A_demultiplex_cellranger_multi_pbmc.sh, 0_B_bamtofastq.sh, 
# and 0_C_align_cellranger_multi_pbmc.sh).
# 
# However, for the purpose of denoising with SoupX (RNA) Mand DSB (ADT), I want raw and filtered GEX and PLX counts
# matrices that have not yet been demultiplexed. So, in this notebook, we perform alignment without demultiplexing.

# Full path to cellranger, because the tmux session may not load the PATH from your .bashrc.
cellranger=/home/jelrod/software/cellranger-8.0.1/bin/cellranger

# For each library: ask for a server, then ssh there and start cellranger in a detached (-d)
# tmux session named after the sample. (ssh -n stops ssh from reading what you type at the
# next prompt.)

# Sample 591

#############
## SAMPLE 515
#############

read -p "Server for 515: " server
ssh -n jelrod@leviathan$server.stat.cmu.edu tmux new-session -d -s 515 \
"'$cellranger multi --id=p515 \
--csv=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_csvs/config_csvs_denoise/515_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_output/cellranger_denoise_2026/p515'"

#############
## SAMPLE 591
#############

read -p "Server for 591: " server
ssh -n jelrod@leviathan$server.stat.cmu.edu tmux new-session -d -s 591 \
"'$cellranger multi --id=p591 \
--csv=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_csvs/config_csvs_denoise/591_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_output/cellranger_denoise_2026/p591'"

###############
## SAMPLE 594_r
###############

read -p "Server for 594_r: " server
ssh -n jelrod@leviathan$server.stat.cmu.edu tmux new-session -d -s 594_r \
"'$cellranger multi --id=p594_r \
--csv=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_csvs/config_csvs_denoise/594_r_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_output/cellranger_denoise_2026/p594_r'"

##############
## SAMPLE H_r 
##############

read -p "Server for H_r: " server
ssh -n jelrod@leviathan$server.stat.cmu.edu tmux new-session -d -s H_r \
"'$cellranger multi --id=H_r \
--csv=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_csvs/config_csvs_denoise/H_r_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_output/cellranger_denoise_2026/H_r'"
