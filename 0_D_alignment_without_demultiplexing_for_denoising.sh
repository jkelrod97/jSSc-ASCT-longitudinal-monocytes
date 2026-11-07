# NOTE (2026-11-03): Sample labels H1-H4 and 594 baseline/six/twelve/twen4 were found to be
# swapped (whole pool: the library originally labelled "594" is the healthy pool and vice versa).
# All labels in this repository, the config CSVs, the FASTQ file names, and the data folders have
# been corrected, so the labels below refer to the correct samples. The original sequencing library
# IDs (CKDL...) stay attached to their data. See sample_swap.csv and sample_swap.py.
#
# Previously, I ran cellranger multi (to combine gex and plx) for each *sample* individually 
# However, for the purpose of denoising with SoupX (RNA) Mand DSB (ADT), I want raw and filtered GEX and PLX counts
# matrices that have not yet been demultiplexed.

#############
## SAMPLE 515
#############

# Server alias
leviathan1
# Start tmux session
tmux
cellranger multi --id=p515 --csv=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_csvs/config_csvs_denoise/515_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_output/redo_for_denoise/p515

#############
## SAMPLE 591
#############

# Server alias
leviathan2
# Start tmux session
tmux
cellranger multi --id=p591 --csv=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_csvs/config_csvs_denoise/591_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_output/redo_for_denoise/p591

#############
## SAMPLE 594
#############

# Server alias
leviathan3
# Start tmux session
tmux
cellranger multi --id=p594 --csv=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_csvs/config_csvs_denoise/594_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_output/redo_for_denoise/p594

#############
## SAMPLE H #
#############

# Server alias
leviathan4
# Start tmux session
tmux
cellranger multi --id=H --csv=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_csvs/config_csvs_denoise/H_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_output/redo_for_denoise/H
