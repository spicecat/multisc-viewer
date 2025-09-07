# MultiSC-Data

Datasets and publications data for [MultiSC-Viewer](https://git.jasonxu.dev/JasonXu/plot-viewer).

## Datasets

Add datasets to the `data/datasets` directory and `data/datasets/meta.json` file. Each dataset is a directory that includes:
- `data.rds` — a Seurat object
- `genes.json` — an array of genes
- `cluster.colors.rds` — a vector mapping cluster labels to colors
- `genotype.colors.rds` — a vector mapping genotypes/conditions to colors

Generate `genes.json` by running the `genes.R` script from each dataset directory (where `data.rds` is located):
```bash
cd path/to/My_Dataset
Rscript ../../../multisc-daemon/genes.R
```

Example `cluster.colors.rds`
```r
cluster_colors <- c(
	"K_RORB+" = "#00B4F0",
	"K_FAT2+" = "royalblue3",
	"K_TH_SOX6hi" = "#D19300",
	"K_TH_SOX6lo" = "#F17D50"
	# ... add all clusters present in your dataset
)
saveRDS(cluster_colors, "cluster.colors.rds")
```

Example `genotype.colors.rds`
```r
genotype_colors <- c(
	CTRL = "green",
	DLBD = "orange1",
	PD   = "blue",
	PDD  = "purple2",
	HD   = "midnightblue",
	AD   = "red",
	NHD  = "sienna4"
)
saveRDS(genotype_colors, "genotype.colors.rds")
```

- Use `data/datasets/example_meta.json` as a template for `data/datasets/meta.json`.
- Ensure each entry’s `id` matches the dataset directory name (e.g., `data/datasets/<id>`).
- Required fields: `id`, `title`, `year`, `authors`, `PMID`, `region`, `disease`, `cellType`, `species`.
- Optional field: `defaultGene`.

## Publications
Add publications to the `data/publications/meta.json` file. 
- Use `data/publications/example_meta.json` as a template for `data/publications/meta.json`.
- Required fields: `id`, `title`, `year`, `authors`, `PMID`, `journal`, `abstract`, `datasets`.
  - `abstract` is raw HTML.
  - `datasets` is an array of dataset IDs included in the publication.


Example data layout:
```
plot-viewer
├── data
│   ├── datasets
│   │   ├── My_Dataset
│   │   │   ├── cluster.colors.rds
│   │   │   ├── data.rds
│   │   │   ├── genes.json
│   │   │   └── genotype.colors.rds
│   │   ├── Another_Dataset
│   │   │   ├── cluster.colors.rds
│   │   │   ├── data.rds
│   │   │   ├── genes.json
│   │   │   └── genotype.colors.rds
│   │   ├── meta.json
│   └── publications
│       └── meta.json
├── multisc-daemon
│   ├── daemon.R
│   ├── genes.R
│   ├── plumber.R
│   └── ...
├── multisc-viewer
│   └── ...
└── ...
```