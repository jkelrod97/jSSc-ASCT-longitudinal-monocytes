Code accompanying "Longitudinal single-cell modeling reveals monocyte reprogramming in juvenile systemic sclerosis following autologous stem cell transplantation" (Elrod et al. 2026)

Preprint: https://www.biorxiv.org/content/10.64898/2026.08.21.738279v1

Torok lab, UPMC Children's hospital.

Townes lab, Carnegie Mellon Statistics & Data Science.

Files are executed in this order:
1. **0_A_demultiplex_cellranger_multi_pbmc.sh**: Demultiplex and align PBMC RNA samples with cellranger multi. 
2. **0_B_align_cellranger_multi_pbmc.sh**: Merge demultiplexed PBMC RNA data with other modalities: protein (ADT), T cell receptor (TCR), and B cell receptor (BCR).
3. **0_C_alignment_without_demultiplexing_for_denoising.sh**: Align PBMC RNA data with other modalities without demultiplexing first. Output used for RNA and protein denoising.
4. **1_pbmc_preprocessing.Rmd**: Preprocess multimodal PBMC data using Seurat and standard single-cell pipeline.
5. **1a_composite_pbmc.py**: Run scomposite algorithm for multiplet detection. https://github.com/CHPGenetics/COMPOSITE
6. **2_pbmc_cell_type_annotation.Rmd**: Perform multimodal clustering and cell type annotation on PBMCs.
7. **3_monocyte_subtype_annotation.Rmd**: Sub-cluster monocytes and annotate subtypes.
8. **4_monocyte_longitudinal_analysis.Rmd**: Implement longitudinal analysis to identify differentially express genes in jSSc monocytes over time since ASCT. 
9. **5_generate_supplemental_spreadsheet.Rmd**: Generate Excel workbook with supplemental tables.
10. **6_generate_supporting_data_values**: Generate supporting data values XLSX workbook.

## Sample label correction (2026-11-03)

Healthy controls H1–H4 and patient 594 (baseline, 6, 12, and 24 months) were found to have their labels swapped. The whole pool was affected: the library originally labelled "594" is the healthy pool, and the other way round. Every label in this repository, the cellranger config CSVs, the FASTQ file names, and the data folders has been corrected, so the labels used throughout refer to the correct samples. The sequencing library IDs (`CKDL…`) stay attached to their data. Cellranger was **not** rerun; its existing outputs were renamed.

- `sample_swap.csv`: every rename rule (sample, pool, library/FASTQ prefix, and bamtofastq folder level).
- `sample_swap.py`: applies the rules.
  - `python3 sample_swap.py rename <dir> --dry-run` previews renaming files and folders under `<dir>`. Drop `--dry-run` to apply.
  - Renames go through temporary names in two passes. md5 checksums are recorded before and after, and a log and manifest are written to `<dir>/sample_swap_logs/`.
  - `python3 sample_swap.py text <files>` swaps labels inside config files.

Downstream objects saved before the correction (h5seurat, SoupX, COMPOSITE, DEG/GSEA outputs, supplement) are out of date. Regenerate them by rerunning steps 4–10 above.
