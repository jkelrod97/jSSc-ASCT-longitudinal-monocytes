Code accompanying "Longitudinal single-cell modeling reveals monocyte reprogramming in juvenile systemic sclerosis following autologous stem cell transplantation" (Elrod et al. 2026)

Preprint: https://www.biorxiv.org/content/10.64898/2026.08.21.738279v1

Torok lab, UPMC Children's hospital.

Townes lab, Carnegie Mellon Statistics & Data Science.

Files are executed in this order:
1. **0_A_demultiplex_cellranger_multi_pbmc.sh**: Demultiplex RNA samples with cellranger multi.
2. **0_B_bamtofastq.sh**: Convert bam files from previous step into FASTQs. 
3. **0_C_align_cellranger_multi_pbmc.sh**: Merge demultiplexed RNA data with other modalities: protein (ADT), T cell receptor (TCR), and B cell receptor (BCR) and align to reference.
4. **0_D_alignment_without_demultiplexing_for_denoising.sh**: Align RNA data with other modalities without demultiplexing first. Output used for RNA and protein denoising.
5. **1_pbmc_preprocessing.Rmd**: Preprocess multimodal PBMC data using Seurat and standard single-cell pipeline.
6. **1a_composite_pbmc.py**: Run scomposite algorithm for multiplet detection. https://github.com/CHPGenetics/COMPOSITE
7. **2_pbmc_cell_type_annotation.Rmd**: Perform multimodal clustering and cell type annotation on PBMCs.
8. **3_monocyte_subtype_annotation.Rmd**: Sub-cluster monocytes and annotate subtypes.
9. **4_monocyte_longitudinal_analysis.Rmd**: Implement longitudinal analysis to identify differentially express genes in jSSc monocytes over time since ASCT. 
10. **5_generate_supplemental_spreadsheet.Rmd**: Generate Excel workbook with supplemental tables.
11. **6_generate_supporting_data_values**: Generate supporting data values XLSX workbook.
