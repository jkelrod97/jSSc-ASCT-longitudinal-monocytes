# NOTE (2026-11-05): Subject labels H and 594 were found to be
# swapped. All labels in this repository, the config CSVs, the FASTQ file names, and the data folders have
# been corrected, so the labels below refer to the correct samples. The original sequencing library
# IDs (CKDL...) stay attached to their data. Corrected files are marked with an r.
#
# In this script, we convert the demultiplexed BAM files made by
# 0_A_demultiplex_cellranger_multi_pbmc.sh back to FASTQ files with bamtofastq, so that each
# sample's gene expression reads can be aligned on their own in the next step,
# 0_C_align_cellranger_multi_pbmc.sh.
#
# We have 15 valid BAM files: 11 patient/time point combos and 4 healthy controls. There is
# no BAM for 515 at twelve months, because no cells were mapped to that sample.
# --reads-per-fastq=2200000000 keeps the output FASTQs from being split up.
#
# bamtofastq writes two subfolders per sample: demultiplexed_samples_[N_]0_1_<flowcell> holds
# the gene expression FASTQs and demultiplexed_samples_[N_]1_1_<flowcell> holds the multiplexing
# capture FASTQs. N is 2, 3 or 4 depending on the library, and is absent for 591. This script
# moves the gene expression FASTQs directly into the sample's folder and deletes both subfolders.
#
# bamtofastq is fast, so all 15 samples run one after another on a single remote server.
# This script asks which leviathan to use (type 5 for leviathan5). It copies the commands
# below to run_bamtofastq.sh in the bamtofastq_2026 folder on the server, then runs that file
# there in one detached tmux session named bamtofastq, so the runs keep going after you log off.
# (tmux limits how long a command can be, so the commands are sent as a file.)
# Run it with:
#   bash 0_B_bamtofastq.sh
# To check on the runs, ssh to that server and run `tmux attach -t bamtofastq` (Ctrl+b then d
# to leave). The tmux session closes when the last sample finishes.
#
# If one sample fails, the script moves on to the next one. No log files are saved, so to confirm
# that every sample worked, check that each sample folder directly contains its FASTQ files.
# A sample that failed will not have them.
# bamtofastq refuses to write into an output folder that already exists, so to re-run a failed
# sample, delete its partial output folder and comment out the samples you don't need.
#
# Input BAMs are read from:
#   cellranger_final_2026/<library>/outs/per_sample_outs/<time point>/count/sample_alignments.bam
# Output FASTQs are written to:
#   bamtofastq_2026/<sample>
# Both folders are under /home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_output.
# Note that these are remote paths with a different file structure than appears in this
# repository.
#
# On 2026-11-08 this script was changed with help from Claude
# (Anthropic, Claude Opus 5.5): sample folders now sit directly in bamtofastq_2026 instead of a
# FASTQs subfolder, only the gene expression FASTQs are kept, and log files are no longer saved.
# The bamtofastq options are unchanged.

out_dir=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_output/bamtofastq_2026

read -p "Server for bamtofastq: " server

# Step 1: create the output folder on the server and copy the commands between the two
# EOF lines into run_bamtofastq.sh there. ('EOF' is quoted, so nothing in it is changed
# on the way; it all runs on the server.)
ssh jelrod@leviathan$server.stat.cmu.edu "mkdir -p $out_dir && cat > $out_dir/run_bamtofastq.sh" <<'EOF'
# Full path to cellranger, because the tmux session may not load the PATH from your .bashrc.
# bamtofastq is run as a cellranger subcommand (cellranger bamtofastq).
cellranger=/home/jelrod/software/cellranger-8.0.1/bin/cellranger

bam_dir=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_output/cellranger_final_2026
out_dir=/home/jelrod/scl-CITE-SEQ/preprocessing/alignment/cellranger/cellranger_output/bamtofastq_2026

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

# Run every sample in turn. A failed sample does not stop the ones after it.
for entry in "${samples[@]}"; do
  read -r sample library timepoint <<< "$entry"
  echo "$(date): starting $sample"
  if ! $cellranger bamtofastq --reads-per-fastq=2200000000 \
    $bam_dir/$library/outs/per_sample_outs/$timepoint/count/sample_alignments.bam \
    $out_dir/$sample; then
    echo "$(date): $sample FAILED"
    continue
  fi

  # Keep only the gene expression FASTQs. The subfolders are deleted only if exactly one
  # gene expression subfolder was found and its FASTQs moved successfully.
  gex_dirs=($out_dir/$sample/demultiplexed_samples_*0_1_*)
  if [ ${#gex_dirs[@]} -eq 1 ] && [ -d "${gex_dirs[0]}" ] \
    && mv "${gex_dirs[0]}"/*.fastq.gz $out_dir/$sample/; then
    rm -r $out_dir/$sample/demultiplexed_samples_*
    echo "$(date): $sample finished"
  else
    echo "$(date): $sample FAILED: could not move the gene expression FASTQs"
  fi
done
echo "$(date): all samples done"
EOF

# Step 2: start run_bamtofastq.sh on the server in a detached (-d) tmux session.
# (ssh -n stops ssh from reading from the keyboard.)
ssh -n jelrod@leviathan$server.stat.cmu.edu tmux new-session -d -s bamtofastq "'bash $out_dir/run_bamtofastq.sh'"
