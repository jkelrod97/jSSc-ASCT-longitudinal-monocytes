# NOTE (2026-11-05): Subject labels H and 594 were found to be
# swapped. All labels in this repository, the config CSVs, the FASTQ file names, and the data folders have
# been corrected, so the labels below refer to the correct samples. The original sequencing library
# IDs (CKDL...) stay attached to their data. Corrected files are marked with an r.
#
# In this script, we convert the demultiplexed BAM files made by
# 0_A_demultiplex_cellranger_multi_pbmc.sh back to FASTQ files with bamtofastq, so that each
# sample's gene expression reads can be aligned on their own in the next step.
#
# We have 15 valid BAM files: 11 patient/time point combos and 4 healthy controls. There is
# no BAM for 515 at twelve months, because no cells were mapped to that sample.
# --reads-per-fastq=2200000000 keeps the output FASTQs from being split up.
#
# We run bamtofastq on a distinct remote server for each sample. This script asks which
# leviathan to use for each sample (type 5 for leviathan5), then starts bamtofastq there
# in a detached tmux session, so the run keeps going after you log off. Run it with:
#   bash 0_B_bamtofastq_copy.sh
# To check on a run, ssh to its server and run `tmux attach -t 515_baseline` (Ctrl+b then d
# to leave). The tmux session closes when bamtofastq finishes; its output is also saved to
# <sample>_bamtofastq.log in the FASTQs folder.
#
# Input BAMs are read from:
#   cellranger_final_2026/<library>/outs/per_sample_outs/<time point>/count/sample_alignments.bam
# Output FASTQs are written to:
#   bamtofastq_2026/FASTQs/<sample>
# Both folders are under /home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_output.
# Note that these are remote paths with a different file structure than appears in this
# repository.
#
# This script was changed on 2026-11-07 with help from Claude
# (Anthropic, Claude Opus 5.5): previously each sample was started by hand (log in to a hydra,
# start tmux, run bamtofastq); now the script asks for the servers and starts every run itself,
# in the same way as 0_A_demultiplex_cellranger_multi_pbmc.sh. The file paths were also updated
# to the new folder layout. The bamtofastq options are otherwise unchanged.

# Full path to cellranger, because the tmux session may not load the PATH from your .bashrc.
# bamtofastq is run as a cellranger subcommand (cellranger bamtofastq).
cellranger=/home/jelrod/software/cellranger-8.0.1/bin/cellranger

bam_dir=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_output/cellranger_final_2026
fastq_dir=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_output/bamtofastq_2026/FASTQs

# For one sample: ask for a server, then ssh there and start bamtofastq in a detached (-d)
# tmux session named after the sample. (ssh -n stops ssh from reading what you type at the
# next prompt.) The FASTQs folder is created first if it does not exist yet; the sample's own
# output folder must not exist yet, or bamtofastq will refuse to run.
#   $1 = sample name (used for the output folder and the tmux session)
#   $2 = library folder in cellranger_final_2026
#   $3 = time point folder in per_sample_outs
run_bamtofastq() {
  sample=$1
  library=$2
  timepoint=$3
  read -p "Server for $sample: " server
  ssh -n jelrod@leviathan$server.stat.cmu.edu tmux new-session -d -s $sample \
  "'mkdir -p $fastq_dir && $cellranger bamtofastq --reads-per-fastq=2200000000 \
$bam_dir/$library/outs/per_sample_outs/$timepoint/count/sample_alignments.bam \
$fastq_dir/$sample 2>&1 | tee $fastq_dir/${sample}_bamtofastq.log'"
}

# To skip a sample, comment out its line.

# Sample 515
# (515 / twelve months has no BAM: no cells were mapped to this sample)
run_bamtofastq 515_baseline 515 baseline
run_bamtofastq 515_six      515 six_months
run_bamtofastq 515_twen4    515 twen4_months

# # Sample 591
# run_bamtofastq 591_baseline 591 baseline
# run_bamtofastq 591_six      591 six_months
# run_bamtofastq 591_twelve   591 twelve_months
# run_bamtofastq 591_twen4    591 twen4_months

# # Sample 594
# # (r marks corrected sample IDs)
# run_bamtofastq 594_r_baseline 594_r baseline
# run_bamtofastq 594_r_six      594_r six_months
# run_bamtofastq 594_r_twelve   594_r twelve_months
# run_bamtofastq 594_r_twen4    594_r twen4_months

# # Healthy controls
# # (r marks corrected sample IDs; in H_r_config.csv the four controls are named 1-4)
# run_bamtofastq H_r_1 H_r 1
# run_bamtofastq H_r_2 H_r 2
# run_bamtofastq H_r_3 H_r 3
# run_bamtofastq H_r_4 H_r 4
