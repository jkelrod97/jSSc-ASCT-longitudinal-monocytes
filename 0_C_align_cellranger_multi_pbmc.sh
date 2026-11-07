# NOTE (2026-11-06): Subject labels H and 594 were found to be
# swapped. All labels in this repository, the config CSVs, the FASTQ file names, and the data folders have
# been corrected, so the labels below refer to the correct samples. The original sequencing library
# IDs (CKDL...) stay attached to their data. Corrected files are marked with an r.
#
# This script is a follow-up to 0_A_demultiplex_cellranger_multi_pbmc.sh.
#
# We run cellranger on a distinct remote server for each library. This script asks which
# leviathan to use for each library (type 5 for leviathan5), then starts cellranger there
# in a detached tmux session, so the run keeps going after you log off. Run it with:
#   bash 0_B_align_cellranger_multi_pbmc.sh
# To check on a run, ssh to its server and run `tmux attach -t 515bas` (Ctrl+b then d to leave).
# The tmux session closes when cellranger finishes; cellranger's own log is in the _log file
# in that sample's output folder.
# 
# In the alignment step, we match up the protein (ADT), T cell receptor (TCR), and B cell receptor (BCR)
# data to the per-sample gene expression (GEX) data. 
#
# This requires 15 separate runs of cellranger multi, one per sample.
#
# Note that the paths and folders referenced here refer to remote paths with a different
# file structure than appears in this repository. Config CSVs referenced here are found in
# cellranger_csvs/alignment.
#
# A third script, 0_C_alignment_without_demultiplexing_for_denoising.sh, uses cellranger multi
# to merge RNA, ADT, TCR, and BCR data without demultiplexing. This is necessary for ambient RNA
# correction with SoupX and ADT denoising with the DSB algorithm. 


# Full path to cellranger, because the tmux session may not load the PATH from your .bashrc.
cellranger=/home/jelrod/software/cellranger-8.0.1/bin/cellranger

# For each library: ask for a server, then ssh there and start cellranger in a detached (-d)
# tmux session named after the sample. (ssh -n stops ssh from reading what you type at the
# next prompt.)

##############
## SUBJECT 515
##############

# Baseline
read -p "Server for 515_baseline: " server
ssh -n jelrod@leviathan$server.stat.cmu.edu tmux new-session -d -s 515bas \
"'$cellranger multi --id=515_baseline \
--csv=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_csvs/final_configs/515_baseline_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_output/cellranger_final_2026_step2/515_baseline'"

# Six mo.
read -p "Server for 515_six: " server
ssh -n jelrod@leviathan$server.stat.cmu.edu tmux new-session -d -s 515six \
"'$cellranger multi --id=515_six \
--csv=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_csvs/final_configs/515_six_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_output/cellranger_final_2026_step2/515_six'"

# Twenty-four mo.
read -p "Server for 515_twen4: " server
ssh -n jelrod@leviathan$server.stat.cmu.edu tmux new-session -d -s 515twen4 \
"'$cellranger multi --id=515_twen4 \
--csv=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_csvs/final_configs/515_twen4_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_output/cellranger_final_2026_step2/515_twen4'"

##############
## SUBJECT 591
##############

# Baseline
read -p "Server for 591_baseline: " server
ssh -n jelrod@leviathan$server.stat.cmu.edu tmux new-session -d -s 591bas \
"'$cellranger multi --id=591_baseline \
--csv=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_csvs/final_configs/591_baseline_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_output/cellranger_final_2026_step2/591_baseline'"

# Six mo.
read -p "Server for 591_six: " server
ssh -n jelrod@leviathan$server.stat.cmu.edu tmux new-session -d -s 591six \
"'$cellranger multi --id=591_six \
--csv=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_csvs/final_configs/591_six_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_output/cellranger_final_2026_step2/591_six'"

# Twelve mo.
read -p "Server for 591_twelve: " server
ssh -n jelrod@leviathan$server.stat.cmu.edu tmux new-session -d -s 591twelve \
"'$cellranger multi --id=591_twelve \
--csv=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_csvs/final_configs/591_twelve_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_output/cellranger_final_2026_step2/591_twelve'"

# Twenty-four mo.
read -p "Server for 591_twen4: " server
ssh -n jelrod@leviathan$server.stat.cmu.edu tmux new-session -d -s 591twen4 \
"'$cellranger multi --id=591_twen4 \
--csv=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_csvs/final_configs/591_twen4_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_output/cellranger_final_2026_step2/591_twen4'"


##############
## SUBJECT 594
##############

# Baseline
read -p "Server for 594_r_baseline: " server
ssh -n jelrod@leviathan$server.stat.cmu.edu tmux new-session -d -s 594_r_bas \
"'$cellranger multi --id=594_r_baseline \
--csv=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_csvs/final_configs/594_r_baseline_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_output/cellranger_final_2026_step2/594_r_baseline'"

# Six mo.
read -p "Server for 594_r_six: " server
ssh -n jelrod@leviathan$server.stat.cmu.edu tmux new-session -d -s 594_r_six \
"'$cellranger multi --id=594_r_six \
--csv=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_csvs/final_configs/594_r_six_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_output/cellranger_final_2026_step2/594_r_six'"

# Twelve mo.
read -p "Server for 594_r_twelve: " server
ssh -n jelrod@leviathan$server.stat.cmu.edu tmux new-session -d -s 594_r_twelve \
"'$cellranger multi --id=594_r_twelve \
--csv=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_csvs/final_configs/594_r_twelve_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_output/cellranger_final_2026_step2/594_r_twelve'"

# Twenty-four mo.
read -p "Server for 594_r_twen4: " server
ssh -n jelrod@leviathan$server.stat.cmu.edu tmux new-session -d -s 594_r_twen4 \
"'$cellranger multi --id=594_r_twen4 \
--csv=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_csvs/final_configs/594_r_twen4_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_output/cellranger_final_2026_step2/594_r_twen4'"


#############
## HEALTHY ##
#############

# 1
read -p "Server for H1_r: " server
ssh -n jelrod@leviathan$server.stat.cmu.edu tmux new-session -d -s H1r \
"'$cellranger multi --id=H1_r \
--csv=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_csvs/final_configs/H1_r_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_output/cellranger_final_2026_step2/H1_r'"

# 2
read -p "Server for H2_r: " server
ssh -n jelrod@leviathan$server.stat.cmu.edu tmux new-session -d -s H2r \
"'$cellranger multi --id=H2_r \
--csv=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_csvs/final_configs/H2_r_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_output/cellranger_final_2026_step2/H2_r'"

# 3
read -p "Server for H3_r: " server
ssh -n jelrod@leviathan$server.stat.cmu.edu tmux new-session -d -s H3r \
"'$cellranger multi --id=H3_r \
--csv=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_csvs/final_configs/H3_r_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_output/cellranger_final_2026_step2/H3_r'"

# 4
read -p "Server for H4_r: " server
ssh -n jelrod@leviathan$server.stat.cmu.edu tmux new-session -d -s H4r \
"'$cellranger multi --id=H4_r \
--csv=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_csvs/final_configs/H4_r_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_output/cellranger_final_2026_step2/H4_r'"
