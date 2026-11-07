# We have 16 BAM files--12 patient/time point combos and 4 healthy controls.
# We use bamtofastq --reads-per-fastq=2200000000 to convert these to FASTQ files
# First run--stored in folder cellranger_multi_with_bam
#
# Updates: We only have 15 valid BAM files (none for sample 15, twelve months)
# Added bamtofastq to path, so I was able to remove cellranger prefix
# Added in --reads-per-fastq=2200000000 so that output FASTQs will not be split up
# Updating available hydras and tmux sessions as I go
# Second run--stored in folder cellranger_multi_with_bam_2

# HYDRA 1
ssh jelrod@hydra1.stat.cmu.edu
# TMUX SESSION 14
tmux
# 515 / BASELINE
bamtofastq --reads-per-fastq=2200000000 ~/scl-CITE-SEQ/cellranger_output/cellranger_multi_with_bam/515/outs/per_sample_outs/baseline/count/sample_alignments.bam \
~/scl-CITE-SEQ/cellranger_output/cellranger_multi_with_bam_2/FASTQs/515_baseline
# TMUX SESSION 15
tmux
# 515 / SIX MONTHS
bamtofastq --reads-per-fastq=2200000000 ~/scl-CITE-SEQ/cellranger_output/cellranger_multi_with_bam/515/outs/per_sample_outs/six_months/count/sample_alignments.bam \
~/scl-CITE-SEQ/cellranger_output/cellranger_multi_with_bam_2/FASTQs/515_six


# 515 / TWELVE MONTHS 
# NO CELLS WERE MAPPED TO THIS SAMPLE

# HYDRA 3
ssh jelrod@hydra3.stat.cmu.edu
# TMUX SESSION 19
tmux
# 515 / TWENTY-FOUR MONTHS
bamtofastq --reads-per-fastq=2200000000 ~/scl-CITE-SEQ/cellranger_output/cellranger_multi_with_bam/515/outs/per_sample_outs/twen4_months/count/sample_alignments.bam \
~/scl-CITE-SEQ/cellranger_output/cellranger_multi_with_bam_2/FASTQs/515_twen4
# TMUX SESSION 20
tmux
# 591 / BASELINE
bamtofastq --reads-per-fastq=2200000000 ~/scl-CITE-SEQ/cellranger_output/cellranger_multi_with_bam/591/outs/per_sample_outs/baseline/count/sample_alignments.bam \
~/scl-CITE-SEQ/cellranger_output/cellranger_multi_with_bam_2/FASTQs/591_baseline

# HYDRA 4
ssh jelrod@hydra4.stat.cmu.edu
# TMUX SESSION 2
tmux
# 591 / SIX MONTHS
bamtofastq --reads-per-fastq=2200000000 ~/scl-CITE-SEQ/cellranger_output/cellranger_multi_with_bam/591/outs/per_sample_outs/six_months/count/sample_alignments.bam \
~/scl-CITE-SEQ/cellranger_output/cellranger_multi_with_bam_2/FASTQs/591_six
# TMUX SESSION 3
tmux
# 591 / TWELVE MONTHS
bamtofastq --reads-per-fastq=2200000000 ~/scl-CITE-SEQ/cellranger_output/cellranger_multi_with_bam/591/outs/per_sample_outs/twelve_months/count/sample_alignments.bam \
~/scl-CITE-SEQ/cellranger_output/cellranger_multi_with_bam_2/FASTQs/591_twelve

# HYDRA 6
ssh jelrod@hydra6.stat.cmu.edu
# TMUX SESSION 1
tmux
# 591 / TWENTY-FOUR MONTHS
bamtofastq --reads-per-fastq=2200000000 ~/scl-CITE-SEQ/cellranger_output/cellranger_multi_with_bam/591/outs/per_sample_outs/twen4_months/count/sample_alignments.bam \
~/scl-CITE-SEQ/cellranger_output/cellranger_multi_with_bam_2/FASTQs/591_twen4
# TMUX SESSION 2
tmux
# 594 / BASELINE
bamtofastq --reads-per-fastq=2200000000 ~/scl-CITE-SEQ/cellranger_output/cellranger_multi_with_bam/594/outs/per_sample_outs/baseline/count/sample_alignments.bam \
~/scl-CITE-SEQ/cellranger_output/cellranger_multi_with_bam_2/FASTQs/594_baseline
# TMUX SESSION 3
tmux
# 594 / SIX MONTHS
bamtofastq --reads-per-fastq=2200000000 ~/scl-CITE-SEQ/cellranger_output/cellranger_multi_with_bam/594/outs/per_sample_outs/six_months/count/sample_alignments.bam \
~/scl-CITE-SEQ/cellranger_output/cellranger_multi_with_bam_2/FASTQs/594_six

# HYDRA 1
ssh jelrod@hydra1.stat.cmu.edu
# TMUX SESSION 16
tmux
# 594 / TWELVE MONTHS
bamtofastq --reads-per-fastq=2200000000 ~/scl-CITE-SEQ/cellranger_output/cellranger_multi_with_bam/594/outs/per_sample_outs/twelve_months/count/sample_alignments.bam \
~/scl-CITE-SEQ/cellranger_output/cellranger_multi_with_bam_2/FASTQs/594_twelve

# HYDRA 10
ssh jelrod@hydra10.stat.cmu.edu
# TMUX SESSION 2
tmux
# 594 / TWENTY-FOUR MONTHS
bamtofastq --reads-per-fastq=2200000000 ~/scl-CITE-SEQ/cellranger_output/cellranger_multi_with_bam/594/outs/per_sample_outs/twen4_months/count/sample_alignments.bam \
~/scl-CITE-SEQ/cellranger_output/cellranger_multi_with_bam_2/FASTQs/594_twen4

# HYDRA 11
ssh jelrod@hydra11.stat.cmu.edu
# TMUX SESSION 1
tmux
# H / BASELINE (HEALTHY CONTROL 1)
bamtofastq --reads-per-fastq=2200000000 ~/scl-CITE-SEQ/cellranger_output/cellranger_multi_with_bam/H/outs/per_sample_outs/baseline/count/sample_alignments.bam \
~/scl-CITE-SEQ/cellranger_output/cellranger_multi_with_bam_2/FASTQs/H_1

# HYDRA 12
ssh jelrod@hydra12.stat.cmu.edu
# TMUX SESSION 1
tmux
# H / SIX MONTHS (HEALTHY CONTROL 2)
bamtofastq --reads-per-fastq=2200000000 ~/scl-CITE-SEQ/cellranger_output/cellranger_multi_with_bam/H/outs/per_sample_outs/six_months/count/sample_alignments.bam \
~/scl-CITE-SEQ/cellranger_output/cellranger_multi_with_bam_2/FASTQs/H_2

# HYDRA 13
ssh jelrod@hydra13.stat.cmu.edu
# TMUX SESSION 2
tmux
# H / TWELVE MONTHS (HEALTHY CONTROL 3)
bamtofastq --reads-per-fastq=2200000000 ~/scl-CITE-SEQ/cellranger_output/cellranger_multi_with_bam/H/outs/per_sample_outs/twelve_months/count/sample_alignments.bam \
~/scl-CITE-SEQ/cellranger_output/cellranger_multi_with_bam_2/FASTQs/H_3
# TMUX SESSION 3
tmux
# H / TWENTY-FOUR MONTHS (HEALTHY CONTROL 4)
bamtofastq --reads-per-fastq=2200000000 ~/scl-CITE-SEQ/cellranger_output/cellranger_multi_with_bam/H/outs/per_sample_outs/twen4_months/count/sample_alignments.bam \
~/scl-CITE-SEQ/cellranger_output/cellranger_multi_with_bam_2/FASTQs/H_4
