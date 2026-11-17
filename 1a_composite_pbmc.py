# In this Python script, we run the multimodal sccomposite multiplet detection algorithm for our PBMCs,
# using MTX files exported from R (1_pbmc_preprocessing.Rmd).
#
# COMPOSITE is run separately for each library (10x lane: 515, 591, 594_r, H_r), since multiplets
# form within a lane and each lane has its own loading density and sequencing depth.
# (Running on all libraries pooled made multiplet calls track library-level ADT depth.)
#
# Inputs (per library, from 1_pbmc_preprocessing.Rmd; copy data/composite/by_library/ to the server):
#   RNA_<lib>.mtx, ADT_<lib>.mtx (features x cells), cells_<lib>.csv (cell names, in mtx column order)
# Outputs (per library; copy back to data/composite/by_library/):
#   multiplet_prediction_<lib>.csv with columns cell_id, multiplet_classification
#
# COMPOSITE prints an RNA and ADT goodness-of-fit score for each library
# (<3: poor fit, 3~5: moderate fit, >5: good fit). Record them here after running:
#   515:   RNA = , ADT =
#   591:   RNA = , ADT =
#   594_r: RNA = , ADT =
#   H_r:   RNA = , ADT =

# Load packages
# See https://github.com/CHPGenetics/COMPOSITE
import sccomposite
import pandas as pd

# Set project path
project_path = "/home/jelrod/scl-CITE-SEQ/preprocessing/multiplet_detection/"
# Folder with per-library input and output
lib_dir = project_path + "PBMC/by_library/"

# Libraries (one 10x lane each)
libraries = ["515", "591", "594_r", "H_r"]

for lib in libraries:
    print(f"\n===== Running COMPOSITE for library {lib} =====")

    # Paths to RNA and ADT matrices (exported from R)
    rna_path = lib_dir + f"RNA_{lib}.mtx"
    adt_path = lib_dir + f"ADT_{lib}.mtx"

    # Cell names, in the same order as the mtx columns
    cell_ids = pd.read_csv(lib_dir + f"cells_{lib}.csv")["cell_id"]

    # Run COMPOSITE algorithm for multiplet detection
    multiplet_classification, multiplet_probability = sccomposite.Multiomics.composite_multiomics(
        RNA=rna_path, ADT=adt_path
    )

    # Make sure there is one prediction per cell
    if len(multiplet_classification) != len(cell_ids):
        raise ValueError(
            f"{lib}: {len(multiplet_classification)} predictions but {len(cell_ids)} cells"
        )

    # Save predictions keyed by cell name
    out = pd.DataFrame({
        "cell_id": cell_ids,
        "multiplet_classification": multiplet_classification,
    })
    out.to_csv(lib_dir + f"multiplet_prediction_{lib}.csv", index=False)

    print(f"{lib}: {int(out['multiplet_classification'].sum())} of {len(out)} cells called multiplets "
          f"({100 * out['multiplet_classification'].mean():.1f}%)")
