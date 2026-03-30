# This script is a follow-up to 0_A_demultiplex_cellranger_multi_pbmc.sh.
# 
# In the alignment step, we match up the protein (ADT), T cell receptor (TCR), and B cell receptor (BCR)
# data to the per-sample gene expression (GEX) data. 
#
# This requires 15 separate runs of cellranger multi, one per sample.
#
# We are starting with sample 515 (recall we are missing the twelve-month benchmark for this sample):
#
# Once again, we use a distinct remote server and tmux session for each cellranger run.
#
# Note that the paths and folders referenced here refer to remote paths with a different
# file structure than appears in this repository. Config CSVs referenced here are found in
# cellranger_csvs/alignment.
#
# A third script, 0_C_alignment_without_demultiplexing_for_denoising.sh, uses cellranger multi
# to merge RNA, ADT, TCR, and BCR data without demultiplexing. This is necessary for ambient RNA
# correction with SoupX and ADT denoising with the DSB algorithm. 

##############
## SUBJECT 515
##############

# Baseline
# Server alias
leviathan1
# Start tmux session
tmux
cellranger multi --id=515_baseline --csv=/home/jelrod/scl-CITE-SEQ/cellranger_csvs/final_configs/515_baseline_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/cellranger_output/final_analysis/output/515_baseline

# Six mo.
# Server alias
leviathan2
# Start tmux session
tmux
cellranger multi --id=515_six --csv=/home/jelrod/scl-CITE-SEQ/cellranger_csvs/final_configs/515_six_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/cellranger_output/final_analysis/output/515_six

# Twenty-four mo.
# Server alias
leviathan3
# Start tmux session
tmux
cellranger multi --id=515_twen4 --csv=/home/jelrod/scl-CITE-SEQ/cellranger_csvs/final_configs/515_twen4_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/cellranger_output/final_analysis/output/515_twen4

##############
## SUBJECT 591
##############

# Baseline
# Server alias
leviathan4
# Start tmux session
tmux
cellranger multi --id=591_baseline --csv=/home/jelrod/scl-CITE-SEQ/cellranger_csvs/final_configs/591_baseline_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/cellranger_output/final_analysis/output/591_baseline

# Six mo.
# Server alias
leviathan5
# Start tmux session
tmux
cellranger multi --id=591_six --csv=/home/jelrod/scl-CITE-SEQ/cellranger_csvs/final_configs/591_six_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/cellranger_output/final_analysis/output/591_six

# Twelve mo.
# Server alias
leviathan6
# Start tmux session
tmux
cellranger multi --id=591_twelve --csv=/home/jelrod/scl-CITE-SEQ/cellranger_csvs/final_configs/591_twelve_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/cellranger_output/final_analysis/output/591_twelve


# Twenty-four mo.
# Server alias
leviathan7
# Start tmux session
tmux
# Delete & recreate output directory (to get rid of failed pipestance)
rm -rf /home/jelrod/scl-CITE-SEQ/cellranger_output/final_analysis/output/591_twen4
mkdir /home/jelrod/scl-CITE-SEQ/cellranger_output/final_analysis/output/591_twen4
cellranger multi --id=591_twen4 --csv=/home/jelrod/scl-CITE-SEQ/cellranger_csvs/final_configs/591_twen4_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/cellranger_output/final_analysis/output/591_twen4


##############
## SUBJECT 594
##############

# Baseline
# Server alias
leviathan8
# Start tmux session
tmux
cellranger multi --id=594_baseline --csv=/home/jelrod/scl-CITE-SEQ/cellranger_csvs/final_configs/594_baseline_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/cellranger_output/final_analysis/output/594_baseline

# Six mo.
# Server alias
leviathan9
# Start tmux session
tmux
cellranger multi --id=594_six --csv=/home/jelrod/scl-CITE-SEQ/cellranger_csvs/final_configs/594_six_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/cellranger_output/final_analysis/output/594_six

# Twelve mo.
# Server alias
leviathan10
# Start tmux session
tmux
cellranger multi --id=594_twelve --csv=/home/jelrod/scl-CITE-SEQ/cellranger_csvs/final_configs/594_twelve_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/cellranger_output/final_analysis/output/594_twelve

# Twenty-four mo.
# Server alias
leviathan11
# Start tmux session
tmux
cellranger multi --id=594_twen4 --csv=/home/jelrod/scl-CITE-SEQ/cellranger_csvs/final_configs/594_twen4_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/cellranger_output/final_analysis/output/594_twen4


#############
## HEALTHY ##
#############

# Baseline
# Server alias
leviathan1
# Start tmux session
tmux
cellranger multi --id=H1 --csv=/home/jelrod/scl-CITE-SEQ/cellranger_csvs/final_configs/H1_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/cellranger_output/final_analysis/output/H1

# Six mo.
# Server alias
leviathan2
# Start tmux session
tmux
cellranger multi --id=H2 --csv=/home/jelrod/scl-CITE-SEQ/cellranger_csvs/final_configs/H2_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/cellranger_output/final_analysis/output/H2

# Twelve mo.
# Server alias
leviathan3
# Start tmux session
tmux
cellranger multi --id=H3 --csv=/home/jelrod/scl-CITE-SEQ/cellranger_csvs/final_configs/H3_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/cellranger_output/final_analysis/output/H3

# Twenty-four mo.
# Server alias
leviathan4
# Start tmux session
tmux
cellranger multi --id=H4 --csv=/home/jelrod/scl-CITE-SEQ/cellranger_csvs/final_configs/H4_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/cellranger_output/final_analysis/output/H4
