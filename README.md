# Metabolic control of macrophage crystal surveillance

Analysis code for the paper **"Metabolic control of macrophage crystal surveillance"**.

The repository contains the single-cell RNA-seq and 10x Genomics Visium spatial
transcriptomics analyses of mouse kidneys with glyoxylate-induced calcium
oxalate crystal nephropathy (**GLY**) and of untreated controls (**NC**).

## Data availability

| Dataset | Accession | Use in this repository |
| --- | --- | --- |
| scRNA-seq and Visium spatial transcriptomics of NC / GLY mouse kidney | [GSE269465](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE269465) | all main analyses (scripts 01–08) |
| Independent mouse kidney Visium section (GSM6250307) | [GSE206306](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206306) | external validation of the metabolic scoring (script 09) |

Raw and processed data are **not** included in the repository. Download them
from GEO and arrange them as described below.

## Repository layout

```
.
├── R/                     shared configuration and helper functions
│   ├── config.R           paths, sample groups, gene sets, analysis parameters
│   ├── utils.R            object loading, PC selection, plot saving helpers
│   ├── plots.R            reusable plot builders (violin, bubble, UMAP, spatial)
│   ├── plot_volcano.R     volcano plot
│   └── setup.R            single entry point sourced by every script
├── scripts/               numbered analysis steps, run in order
│   ├── 01_scrna_kidney_overview.R
│   ├── 02_scrna_myeloid_subsets.R
│   ├── 03_scrna_crystal_adhesion_genes.R
│   ├── 04_scrna_de_and_enrichment.R
│   ├── 05_spatial_preprocess.R
│   ├── 06_spatial_feature_maps.R
│   ├── 07_spatial_deconvolution.R
│   ├── 08_spatial_cellchat.R
│   └── 09_spatial_metabolism.R
├── data/                  input data (not tracked)
└── results/               figures and tables (not tracked)
```

## Expected data layout

```
data/
├── scRNA/
│   ├── IntegrateData_3.Rdata          annotated whole-kidney object (scRNA1, labels in cell_type)
│   └── IntegrateData_3_Macro.Rdata    re-clustered myeloid object   (scRNA1, labels in celltype)
├── spatial/
│   ├── NC/
│   │   ├── filtered_feature_bc_matrix.h5
│   │   └── spatial/                   tissue images and scalefactors_json.json
│   └── GLY/
│       ├── filtered_feature_bc_matrix.h5
│       └── spatial/
├── spatial_validation/
│   └── GSM6250307/                    same structure as above (GSE206306)
└── processed/                         intermediate objects written by the pipeline
```

The two `.Rdata` files hold the integrated Seurat objects obtained from the
GSE269465 count matrices with the standard Seurat workflow (QC filtering,
`SCTransform`, integration, PCA, clustering, manual annotation).

## Analysis steps

| Script | Content |
| --- | --- |
| `01_scrna_kidney_overview.R` | marker violin plots of the whole-kidney annotation, cell composition per group, `Triobp` expression across compartments and in endothelium |
| `02_scrna_myeloid_subsets.R` | PC selection, UMAP and marker bubble plot of the ten myeloid / lymphoid subsets, subset composition per group |
| `03_scrna_crystal_adhesion_genes.R` | UMAP feature plots of the crystal adhesion genes (`Cd81`, `Itgb1`, `Itgb5`, `Rac1`, `Actr2`, `Actr3`), `Itgb5` in Cx3cr1+ macrophages (NC vs GLY) |
| `04_scrna_de_and_enrichment.R` | differential expression (Cx3cr1+ vs other macrophages; NC vs GLY within Cx3cr1+), volcano plots, GO biological process enrichment |
| `05_spatial_preprocess.R` | `SCTransform`, PCA, clustering and UMAP of the NC and GLY Visium sections |
| `06_spatial_feature_maps.R` | spatial expression of the fatty-acid transport, crystal adhesion, SPP1 / MIF / CCL and metabolic gene sets, crystal adhesion module score |
| `07_spatial_deconvolution.R` | RCTD deconvolution against the scRNA-seq reference: dominant cell type per spot (doublet mode) and cell type weights (full mode, scatterpies) |
| `08_spatial_cellchat.R` | spatially informed CellChat analysis of cell-cell contact signalling, aggregated networks and spatial maps of selected ligand-receptor pairs |
| `09_spatial_metabolism.R` | scMetabolism (VISION / KEGG) pathway scoring of the independent GSE206306 section |

Scripts 05 → 07 → 08 must be run in that order, because 07 consumes the
clustered objects written by 05 and 08 consumes the annotated objects written
by 07. The scRNA-seq scripts (01–04) are independent of each other.

## Running the analysis

All scripts assume the repository root as working directory and read their
configuration from `R/config.R`.

```bash
Rscript scripts/01_scrna_kidney_overview.R
Rscript scripts/02_scrna_myeloid_subsets.R
Rscript scripts/03_scrna_crystal_adhesion_genes.R
Rscript scripts/04_scrna_de_and_enrichment.R
Rscript scripts/05_spatial_preprocess.R
Rscript scripts/06_spatial_feature_maps.R
Rscript scripts/07_spatial_deconvolution.R
Rscript scripts/08_spatial_cellchat.R
Rscript scripts/09_spatial_metabolism.R
```

Each script writes its figures and tables into its own sub-directory of
`results/`, so nothing is overwritten between steps.

## Dependencies

R (>= 4.1) and the following packages.

CRAN:

```r
install.packages(c("Seurat", "ggplot2", "dplyr", "tidyr", "reshape2",
                   "RColorBrewer", "ggsignif", "ggrepel", "ggsci",
                   "cols4all", "stringr", "jsonlite", "Matrix"))
```

Bioconductor:

```r
BiocManager::install(c("clusterProfiler", "org.Mm.eg.db", "GOplot"))
```

GitHub:

```r
remotes::install_github("dmcable/spacexr")            # RCTD deconvolution
remotes::install_github("JEFworks-Lab/STdeconvolve")  # scatterpie visualisation
remotes::install_github("jinworks/CellChat")          # cell-cell communication
remotes::install_github("YosefLab/VISION")            # required by scMetabolism
remotes::install_github("wu-yc/scMetabolism")         # metabolic pathway scoring
```

`SCTransform` requires `glmGamPoi` for reasonable runtimes
(`BiocManager::install("glmGamPoi")`).

## Notes

* All analyses are performed on mouse data, so the mouse annotation resources
  are used throughout (`org.Mm.eg.db`, `CellChatDB.mouse`, `PPI.mouse`).
* The number of principal components used downstream (`n_pcs_scrna`,
  `n_pcs_spatial` in `R/config.R`) was chosen once with `choose_pcs()` and is
  kept fixed; the scripts still report the suggested value for transparency.
* RCTD is run with `max_cores = 1` for reproducibility; increase it in
  `scripts/07_spatial_deconvolution.R` to speed up the deconvolution.

## Citation

If you use this code, please cite:

> Metabolic control of macrophage crystal surveillance.
