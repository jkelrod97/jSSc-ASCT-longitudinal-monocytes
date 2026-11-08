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
# bamtofastq is fast, so all 15 samples run one after another on a single remote server.
# This script asks which leviathan to use (type 5 for leviathan5). It copies the commands
# below to run_bamtofastq.sh in the FASTQs folder on the server, then runs that file there
# in one detached tmux session named bamtofastq, so the runs keep going after you log off.
# (tmux limits how long a command can be, so the commands are sent as a file.)
# Run it with:
#   bash 0_B_bamtofastq_copy.sh
# To check on the runs, ssh to that server and run `tmux attach -t bamtofastq` (Ctrl+b then d
# to leave). The tmux session closes when the last sample finishes.
#
# If one sample fails, the script moves on to the next one. Each sample's bamtofastq output
# is saved to <sample>_bamtofastq.log in the FASTQs folder, and bamtofastq_progress.log in the
# same folder lists when each sample started and whether it finished or failed.
# bamtofastq refuses to write into an output folder that already exists, so to re-run a failed
# sample, delete its partial output folder and comment out the samples you don't need.
#
# Input BAMs are read from:
#   cellranger_final_2026/<library>/outs/per_sample_outs/<time point>/count/sample_alignments.bam
# Output FASTQs are written to:
#   bamtofastq_2026/FASTQs/<sample>
# Both folders are under /home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_output.
# Note that these are remote paths with a different file structure than appears in this
# repository.
#
# This script is a copy of 0_B_bamtofastq.sh, changed on 2026-11-07 with help from Claude
# (Anthropic, Claude Opus 5.5): previously each sample was started by hand (log in to a hydra,
# start tmux, run bamtofastq); now the script asks for a server and runs every sample there
# itself. The file paths were also updated to the new folder layout. The bamtofastq options
# are otherwise unchanged.

fastq_dir=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_output/bamtofastq_2026/FASTQs

read -p "Server for bamtofastq: " server

# Step 1: create the FASTQs folder on the server and copy the commands between the two
# EOF lines into run_bamtofastq.sh there. ('EOF' is quoted, so nothing in it is changed
# on the way; it all runs on the server.)
ssh jelrod@leviathan$server.stat.cmu.edu "mkdir -p $fastq_dir && cat > $fastq_dir/run_bamtofastq.sh" <<'EOF'
# Full path to cellranger, because the tmux session may not load the PATH from your .bashrc.
# bamtofastq is run as a cellranger subcommand (cellranger bamtofastq).
cellranger=/home/jelrod/software/cellranger-8.0.1/bin/cellranger

bam_dir=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_output/cellranger_final_2026
fastq_dir=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_output/bamtofastq_2026/FASTQs
progress_log=$fastq_dir/bamtofastq_progress.log

# Samples to convert, one per line: sample name (used for the output folder), library folder
# in cellranger_final_2026, and time point folder in per_sample_outs.
# To skip a sample, comment out its line.
samples=(
  # Sample 515
  # (515 / twelve months has no BAM: no cells were mapped to this sample)
  "515_baseline 515 baseline"
  "515_six      515 six_months"
  "515_twen4    515 twen4_months"

  # Sample 591
  "591_baseline 591 baseline"
  "591_six      591 six_months"
  "591_twelve   591 twelve_months"
  "591_twen4    591 twen4_months"

  # Sample 594
  # (r marks corrected sample IDs)
  "594_r_baseline 594_r baseline"
  "594_r_six      594_r six_months"
  "594_r_twelve   594_r twelve_months"
  "594_r_twen4    594_r twen4_months"

  # Healthy controls
  # (r marks corrected sample IDs; in H_r_config.csv the four controls are named 1-4)
  "H_r_1 H_r 1"
  "H_r_2 H_r 2"
  "H_r_3 H_r 3"
  "H_r_4 H_r 4"
)

# Run every sample in turn. A failed sample is logged and does not stop the ones after it.
for entry in "${samples[@]}"; do
  read -r sample library timepoint <<< "$entry"
  echo "$(date): starting $sample" | tee -a $progress_log
  if $cellranger bamtofastq --reads-per-fastq=2200000000 \
    $bam_dir/$library/outs/per_sample_outs/$timepoint/count/sample_alignments.bam \
    $fastq_dir/$sample > $fastq_dir/${sample}_bamtofastq.log 2>&1; then
    echo "$(date): $sample finished" | tee -a $progress_log
  else
    echo "$(date): $sample FAILED, see ${sample}_bamtofastq.log" | tee -a $progress_log
  fi
done
echo "$(date): all samples done" | tee -a $progress_log
EOF

# Step 2: start run_bamtofastq.sh on the server in a detached (-d) tmux session.
# (ssh -n stops ssh from reading from the keyboard.)
ssh -n jelrod@leviathan$server.stat.cmu.edu tmux new-session -d -s bamtofastq "'bash $fastq_dir/run_bamtofastq.sh'"
