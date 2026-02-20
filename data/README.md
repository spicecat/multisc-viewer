# MultiSC-Data

Datasets and publications data for [MultiSC-Viewer](https://github.com/spicecat/multisc-viewer/tree/main).

## Data Directory Layout

```
data
├── datasets
│   ├── index.json
│   └── <Dataset>
│       ├── metadata.json
│       ├── data.rds
│       ├── genes.json
│       ├── degs.json
│       ├── cluster.colors.rds
│       └── genotype.colors.rds
└── publications
    ├── index.json
    └── <Publication>
        └── metadata.json
```

---

## 1. Datasets

Single-cell RNA-seq datasets are stored as folders containing a serialized Seurat object, dataset metadata, and gene lists.

### 1.1 Datasets Index (`data/datasets/index.json`)

`data/datasets/index.json` is a map of dataset `id` to a path of the dataset folder (relative to `data/datasets`).
All `id` should be unique.

Example `index.json`:

```jsonc
{
  "example-ds-1": "Dataset-1",
  "example-ds-2": "./Dataset-2",
  "example-ds-3": "/absolute/path/to/EXAMPLE-datasets/Dataset-3"
}
```

---

### 1.2 Dataset Folder

Each dataset is stored as its own folder.
Dataset folders names may be arbitrary, but must match the path specified in the datasets index ([1.1](#11-datasets-index-datadatasetsindexjson)).

#### 1.2.1 Dataset Folder Layout

```
<Dataset>
├── metadata.json
├── data.rds
├── genes.json
├── degs.json
├── cluster.colors.rds
└── genotype.colors.rds
```

Required Files

| File                  | Description                                                                                    |
| --------------------- | ---------------------------------------------------------------------------------------------- |
| `metadata.json`       | JSON object of dataset metadata ([1.2.2](#122-dataset-metadata-metadatajson))                  |
| `data.rds`            | Serialized R Seurat object                                                                     |
| `cluster.colors.rds`  | Serialized R named vector mapping cluster label to colors ([1.2.5](#125-r-color-vectors))      |
| `genotype.colors.rds` | Serialized R named vector mapping genotype/condition to colors ([1.2.5](#125-r-color-vectors)) |

Optional Files

| File         | Description                                                                                                    |
| ------------ | -------------------------------------------------------------------------------------------------------------- |
| `genes.json` | Array of dataset genes ([1.2.3](#123-dataset-genes-genesjson))                                                 |
| `degs.json`  | Object mapping DEG set id to an object of DEGs ([1.2.4](#124-dataset-differentially-expressed-genes-degsjson)) |

#### 1.2.2 Dataset Metadata (`metadata.json`)

Dataset metadata schema is sourced from the [Dataset component](https://github.com/spicecat/multisc-viewer/tree/main/MultiSCViewerR/openapi.json#/components/schemas/Dataset) of the [MultiSC-Daemon API](https://github.com/spicecat/multisc-viewer/tree/main/MultiSCViewerR/openapi.json).

`required`: `id`

All other properties are optional, but recommended when applicable.

| Field             | Type                | Description                                                                                    |
| ----------------- | ------------------- | ---------------------------------------------------------------------------------------------- |
| `id`             | string              | Unique dataset id (key used in datasets index [1.1](#11-datasets-index-datadatasetsindexjson)) |
| `author`          | Author              | Array of author objects ([1.2.2.1](#1221-author))                                              |
| `citation`        | Citation            | Array of citation objects ([1.2.2.2](#1222-citation))                                          |
| `publicationDate`            | string (YYYY-MM-DD) | Release publicationDate                                                                                   |
| `defaultGene`     | string              | Default gene to plot                                                                           |
| `description`     | string (HTML)       | Abstract (supports raw HTML)                                                                   |
| `displayName`     | string              | Short display label                                                                            |
| `healthCondition` | OntologyTerm        | Array of diseases or health conditions ([1.2.2.3](#1223-ontologyterm))                         |
| `identifier`      | string              | External accession (e.g., GEO / SRA / PMID)                                                    |
| `name`            | string (HTML)       | Dataset title (supports raw HTML)                                                              |
| `species`         | OntologyTerm        | Array of species objects ([1.2.2.3](#1223-ontologyterm))                                       |
| `cellType`        | OntologyTerm        | Array of cell type objects ([1.2.2.3](#1223-ontologyterm))                                     |
| `tissue`          | OntologyTerm        | Array of tissue or anatomical source objects ([1.2.2.3](#1223-ontologyterm))                   |
| `url`             | string (URI)        | Link to dataset source (if `identifier` is provided)                                           |

Full dataset `metadata.json` examples: [`Example-datasets/Dataset-1/metadata.json`](EXAMPLE-datasets/Dataset-1/metadata.json) and [`EXAMPLE-datasets/Dataset-2/metadata.json`](EXAMPLE-datasets/Dataset-2/metadata.json).

Example dataset folder `metadata.json`:

```jsonc
{
  "id": "example-ds-1",
  "name": "Example Dataset 1 with <i>HTML</i>",
  "displayName": "Example Dataset 1",
  "description": "This is a summary of this dataset.",
  "author": [{ "name": "Doe J" }, { "name": "Smith A" }],
  "citation": [
    {
      "identifier": "PMID:00000000",
      "name": "Primary Study Title",
      "url": "https://pubmed.ncbi.nlm.nih.gov/00000000/"
    }
  ],
  "publicationDate": "2024-01-15",
  "identifier": "GSE000000",
  "defaultGene": "APOE",
  "healthCondition": [
    {
      "name": "Alzheimer Disease",
      "displayName": "AD",
      "identifier": "mesh:D000544"
    }
  ],
  "cellType": [{ "name": "Astrocyte" }],
  "tissue": [{ "name": "Middle Temporal Gyrus" }],
  "species": [{ "name": "Homo sapiens", "displayName": "Human" }],
  "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE000000"
}
```

##### 1.2.2.1 Author

Array of objects with at least `name`.

| Field  | Type   | Description |
| ------ | ------ | ----------- |
| `name` | string | Author name |

Example `Author`:

```jsonc
[{ "name": "Doe J" }, { "name": "Smith A" }]
```

##### 1.2.2.2 Citation

Array of objects with at least one of `name`, `identifier`, or `url`.

| Field        | Type         | Description                      |
| ------------ | ------------ | -------------------------------- |
| `name`       | string       | Publication title                |
| `identifier` | string       | External id (e.g., PMID:#######) |
| `url`        | string (URI) | Link to citation source          |

Example `citation`:

```jsonc
[
  {
    "identifier": "PMID:00000000",
    "name": "Primary Study Title",
    "url": "https://pubmed.ncbi.nlm.nih.gov/00000000/"
  }
]
```

##### 1.2.2.3 OntologyTerm

Array of objects with at least `name`.

The ontology term schema is used for the `healthCondition`, `species`, `cellType`, and `tissue` fields of a dataset `metadata.json` file ([1.2.2](#122-dataset-metadata-metadatajson)).

| Field         | Type         | Description                                          |
| ------------- | ------------ | ---------------------------------------------------- |
| `name`        | string       | Term name                                            |
| `displayName` | string       | Short label (e.g., abbreviation)                     |
| `identifier`  | string       | External id (e.g., `mesh:D000544`, `UBERON:0002435`) |
| `url`         | string (URI) | Link to term source                                  |

Example `healthCondition`:

```jsonc
[
  {
    "name": "Alzheimer Disease",
    "displayName": "AD",
    "identifier": "mesh:D000544",
    "url": "http://id.nlm.nih.gov/mesh/D000544"
  }
]
```

Example `species`:

```jsonc
[
  {
    "displayName": "Human | Homo sapiens",
    "identifier": "9606",
    "name": "Homo sapiens",
    "url": "https://www.uniprot.org/taxonomy/9606"
  }
]
```

Example `cellType`:

```jsonc
[
  {
    "url": "http://id.nlm.nih.gov/mesh/D001253",
    "name": "Astrocytes",
    "identifier": "mesh:D001253"
  }
]
```

Example `tissue`:

```jsonc
[
  {
    "url": "http://id.nlm.nih.gov/mesh/D013378",
    "name": "Substantia Nigra",
    "identifier": "mesh:D013378"
  }
]
```

#### 1.2.3 Dataset Genes (`genes.json`)

Array of strings representing genes in the dataset.

If `genes.json` is not provided, [MultiSC-Daemon](https://github.com/spicecat/multisc-viewer/tree/main/MultiSCViewerR) will generate a `genes.json` file from the Seurat object in `data.rds`.

`genes.json` can be generated manually in a Dataset folder with the R script:

```r
library(SeuratObject)
rds <- readRDS(file = "data.rds")
path <- "genes.json"
genes <- rownames(rds@assays$RNA@counts)
jsonlite::write_json(genes, path)
```

Example `genes.json`:

```jsonc
["SNCA", "MAPT", "LRRK2", "GBA", "APP", "PSEN1", "PSEN2"]
```

#### 1.2.4 Dataset Differentially Expressed Genes (`degs.json`)

Object whose keys are DEG set ids for the dataset. Each value is an object with:

| Field  | Type     | Description                             |
| ------ | -------- | --------------------------------------- |
| `id`  | string   | Unique DEG set id                       |
| `gene` | string[] | Array of differentially expressed genes |
| `name` | string   | Display name                            |

Example `degs.json`:

```jsonc
{
  "example-ds-condition-1-deg": {
    "id": "example-ds-deg",
    "name": "Dataset DEGs",
    "gene": ["SNCA", "MAPT", "LRRK2"]
  },
  "example-ds-condition-2-deg": {
    "id": "example-ds-deg",
    "name": "Dataset DEGs",
    "gene": ["LRRK2", "GBA", "APP"]
  }
}
```

#### 1.2.5 R Color Vectors

`cluster.colors.rds` is a serialized R named vector mapping cluster labels to colors.

`cluster.colors.rds` can be generated manually in a Dataset folder with the R script:

```r
cluster_colors <- c(
	"K_RORB+" = "#00B4F0",
	"K_FAT2+" = "royalblue3",
	"K_TH_SOX6hi" = "#D19300",
	"K_TH_SOX6lo" = "#F17D50"
	# ...
)
saveRDS(cluster_colors, "cluster.colors.rds")
```

`genotype.colors.rds` is a serialized R named vector mapping genotypes/conditions to colors.

`genotype.colors.rds` can be generated manually in a Dataset folder with the R script:

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

---

## 3. Publications

Publications represent papers that reference one or more datasets.
Publications are stored as folders containing publication metadata.

### 3.1 Publications Index (`data/publications/index.json`)

`data/publications/index.json` is a map of publication `id` to a path of the publication folder (relative to `data/publications`).
All `id` should be unique.

Example `index.json`:

```jsonc
{
  "example-pub-1": "Publication-1",
  "example-pub-2": "./Publication-2",
  "example-pub-3": "./relative/path/to/EXAMPLE-publications/Publication-3",
  "example-pub-4": "/absolute/path/to/EXAMPLE-publications/Publication-4"
}
```

### 3.2 Publication Folder

Each publication is stored as its own folder.
Publication folder names may be arbitrary, but must match the path specified in the publications index ([3.1](#31-publications-index-datapublicationsindexjson)).

#### 3.2.1 Publication Folder Layout

```
<Publication>
└── metadata.json
```

#### 3.2.2 Publication Metadata (`metadata.json`)

Publication metadata schema is sourced from the [Publication component](https://github.com/spicecat/multisc-viewer/tree/main/MultiSCViewerR/openapi.json#/components/schemas/Publication) of the [MultiSC-Daemon API](https://github.com/spicecat/multisc-viewer/tree/main/MultiSCViewerR/openapi.json).

`required`: `id`, `datasets`

| Field         | Type                | Description                                                        |
| ------------- | ------------------- | ------------------------------------------------------------------ |
| `id`         | string              | Unique publication id                                              |
| `author`      | Author              | Array of author objects (same as Dataset, [1.2.2.1](#1221-author)) |
| `publicationDate`        | string (YYYY-MM-DD) | Publication publicationDate                                                   |
| `description` | string (HTML)       | Abstract (supports raw HTML)                                       |
| `journalName` | string              | Journal title                                                      |
| `identifier`  | string              | External id (e.g., PMID:#######)                                   |
| `name`        | string (HTML)       | Publication title (supports raw HTML)                              |
| `url`         | string (URI)        | Link to publication source (if `identifier` is provided)           |
| `datasets`    | string[]            | Array of dataset `id`s referenced                                 |

Example publication folder `metadata.json`:

```jsonc
{
  "id": "example-pub-1",
  "name": "Multi-omic integration of astrocyte states",
  "journalName": "Nature",
  "identifier": "PMID:00000000",
  "publicationDate": "2025-02-10",
  "author": [{ "name": "Doe J" }, { "name": "Smith A" }],
  "description": "This is an abstract of this publication.",
  "url": "https://pubmed.ncbi.nlm.nih.gov/00000000/",
  "datasets": ["example-ds-1", "example-ds-2"]
}
```

---

## Related

- MultiSC: [https://github.com/spicecat/multisc-viewer/tree/main/multisc-viewer](https://github.com/spicecat/multisc-viewer/tree/main/README.md)
- Add MultiSC-Data: [https://github.com/spicecat/multisc-viewer/tree/main/data](https://github.com/spicecat/multisc-viewer/tree/main/data/README.md)
- Run MultiSC-Daemon: [https://github.com/spicecat/multisc-viewer/tree/main/MultiSCViewerR](https://github.com/spicecat/multisc-viewer/tree/main/MultiSCViewerR/README.md)
- Run MultiSC-Viewer: [https://github.com/spicecat/multisc-viewer/tree/main/multisc-viewer](https://github.com/spicecat/multisc-viewer/tree/main/multisc-viewer/README.md)
