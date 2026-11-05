# NOTE (2026-11-05): Subject labels H and 594 were found to be
# swapped. All labels in this repository, the config CSVs, the FASTQ file names, and the data folders have
# been corrected, so the labels below refer to the correct samples. The original sequencing library
# IDs (CKDL...) stay attached to their data. Corrected files are marked with an r.
#
# In this script, we demultiplex each PBMC library using cellranger multi.
#
# We run cellranger on a distinct remote server for each library. This script asks which
# leviathan to use for each library (type 5 for leviathan5), then starts cellranger there
# in a detached tmux session, so the run keeps going after you log off. Run it with:
#   bash 0_A_demultiplex_cellranger_multi_pbmc.sh
# To check on a run, ssh to its server and run `tmux attach -t 591` (Ctrl+b then d to leave).
# The tmux session closes when cellranger finishes; cellranger's own log is in the _log file
# in that sample's output folder.
#
# Note that the paths and folders referenced here refer to remote paths with a different
# file structure than appears in this repository. Config CSVs referenced here are found in
# cellranger_csvs/demultiplexing. The HTO barcodes used for demultiplexing are defined in
# 5p_hashing_dmux_cmo-set.csv, which is found in the same subfolder. The hashtags used to
# label subject 594 at six and twelve months were swapped, and this is reflected in the
# config csv.
#
# This script is followed up by 0_B_align_cellranger_multi_pbmc.sh, where we will merge our
# protein (ADT), T cell receptor (TCR), and B cell receptor (BCR) data with our demultiplexed
# gene expression (GEX data).
#
# A third script, 0_C_alignment_without_demultiplexing_for_denoising.sh, uses cellranger multi
# to merge RNA, ADT, TCR, and BCR data without demultiplexing. This is necessary for ambient RNA
# correction with SoupX and ADT denoising with the DSB algorithm.
#
# On 2026-11-05 this script was changed with help from Claude (Anthropic, Claude Opus 5.5):
# previously each library was started by hand (log in to the server, start tmux, run
# cellranger); now the script asks for the servers and starts every run itself. The
# cellranger commands are otherwise unchanged.

# Full path to cellranger, because the tmux session may not load the PATH from your .bashrc.
cellranger=/home/jelrod/software/cellranger-8.0.1/bin/cellranger

# For each library: ask for a server, then ssh there and start cellranger in a detached (-d)
# tmux session named after the sample. (ssh -n stops ssh from reading what you type at the
# next prompt.)

# Sample 591

read -p "Server for 591: " server
ssh -n jelrod@leviathan$server.stat.cmu.edu tmux new-session -d -s 591 \
"'$cellranger multi --id=demultiplexed_samples \
--csv=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_csvs/591_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_output/cellranger_final_2026/591'"

# Sample 515

read -p "Server for 515: " server
ssh -n jelrod@leviathan$server.stat.cmu.edu tmux new-session -d -s 515 \
"'$cellranger multi --id=demultiplexed_samples_2 \
--csv=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_csvs/515_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_output/cellranger_final_2026/515'"

# Sample 594
# (r marks corrected sample IDs)

read -p "Server for 594_r: " server
ssh -n jelrod@leviathan$server.stat.cmu.edu tmux new-session -d -s 594_r \
"'$cellranger multi --id=demultiplexed_samples_3 \
--csv=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_csvs/594_r_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_output/cellranger_final_2026/594_r'"

# Healthy controls
# (r marks corrected sample IDs)

read -p "Server for H_r: " server
ssh -n jelrod@leviathan$server.stat.cmu.edu tmux new-session -d -s H_r \
"'$cellranger multi --id=demultiplexed_samples_4 \
--csv=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_csvs/H_r_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_output/cellranger_final_2026/H_r'"
