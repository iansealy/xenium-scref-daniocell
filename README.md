# Daniocell Single Cell Reference for Xenium Custom Panel Design

Download Daniocell and annotation data to `/data/scratch/bty107/daniocell`

From https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE223922:

- GSE223922_Sur2023_counts.mtx.gz
- GSE223922_Sur2023_counts_cols_cells.txt.gz
- GSE223922_Sur2023_counts_rows_genes.txt.gz
- GSE223922_Sur2023_metadata.tsv.gz

From https://daniocell.nichd.nih.gov/:

- cluster_annotations.csv
- Daniocell2023_SeuratV4.rds

From https://www.umassmed.edu/lawson-lab/reagents/zebrafish-transcriptome/:

- v4_3_2geneinfo.txt

Run RStudio 2024 (2024.09.0-375 & R 4.4.1 (Rocky9), 1 core, 128GB, high memory node) on Open OnDemand::

```
install.packages("tidyverse")
install.packages("Seurat")
install.packages("R.utils")

library(Matrix)
library(Seurat)
library(tidyverse)

setwd("/data/scratch/bty107/daniocell")

# Load the data from GEO and make a Seurat object
counts <- as(readMM("GSE223922_Sur2023_counts.mtx.gz"), "CsparseMatrix")
colnames(counts) <- read.delim("GSE223922_Sur2023_counts_cols_cells.txt.gz", header=FALSE)$V1
rownames(counts) <- read.delim("GSE223922_Sur2023_counts_rows_genes.txt.gz", header=FALSE)$V1
meta <- read.delim("GSE223922_Sur2023_metadata.tsv.gz")
daniocell <- CreateSeuratObject(counts, project="Daniocell", assay="RNA", meta.data=meta)

# Update to the latest metadata
daniocell.annot <- read.csv("cluster_annotations.csv")
rownames(daniocell.annot) <- daniocell.annot$clust
daniocell.annot <- daniocell.annot[,c("tissue", "identity.super", "identity.sub", "identity.super.short", "identity.sub.short", "zfin")]
daniocell@meta.data <- cbind(daniocell@meta.data, daniocell.annot[daniocell@meta.data$clust,])

# Save full object
SaveSeuratRds(daniocell, "Daniocell-Seurat-counts.rds")

# Subset to required tissue and stage
cells.to.keep <- rownames(daniocell@meta.data)[daniocell@meta.data[[26]] %in% c("glial", "neural")]
glial_neural <- subset(daniocell, cells = cells.to.keep)
glial_neural_5dpf <- subset(glial_neural, subset = stage.group == "120")
SaveSeuratRds(glial_neural, "Daniocell-Seurat-glial-neural-counts.rds")
SaveSeuratRds(glial_neural_5dpf, "Daniocell-Seurat-glial-neural-5dpf-counts.rds")

# Output gene name metadata
write_tsv(as.data.frame(rownames(daniocell)), "gene_names.tsv", col_names=FALSE)
```

The Daniocell paper says:

> We mapped reads to the GRCz11 genome, annotated using the Lawson Lab Zebrafish Transcriptome Annotation (v4.3.2) that harmonizes Ensembl and Refseq annotations, includes improved 3' UTR models, and proposes additional gene models.

But the Daniocell data only contains gene names and these don't map uniquely to the v4.3.2 annotation, so perhaps another version was used

Additionally, some names aren't unique (but have had a suffix added in Daniocell) and can't be mapped back to a single Ensembl ID:

join -t $'\t' -j1 -a 1 <(sort -k 1b,1 gene_names.tsv) <(cut -f6,8 v4_3_2geneinfo.txt | sort -k 1b,1 | grep ENSDARG) \
  | sort -u | cut -f1 | uniq -c | sed -e 's/^ *//' | sed -e 's/ /\t/' | awk -F"\t" '$1 > 1 { print $2 }' \
  | join -t $'\t' -j1 - <(cut -f6,8 v4_3_2geneinfo.txt | sort -k 1b,1 | grep ENSDARG) \
  > problematic-genes.txt

Ignore these genes when mapping back to Ensembl IDs:

join -t $'\t' -j1 -a 1 <(sort -k 1b,1 gene_names.tsv) <(cut -f6,8 v4_3_2geneinfo.txt | sort -k 1b,1 | grep ENSDARG | grep -vf <(cut -f2 problematic-genes.txt | sort -u)) \
  | sort -u | awk -F"\t" -v OFS="\t" '{ if ($2 !~ /ENS/) $2 = "FALSE"; print $0 }' | sed -E 's/\.[0-9]+$/\tTRUE/' | sed -e 's/FALSE$/\tFALSE/' \
  > gene_names_ens_unsorted.tsv

Sort into original order:

awk -F"\t" 'FNR == NR { lineno[$1] = NR; next} {print lineno[$1] "\t" $0;}' gene_names.tsv gene_names_ens_unsorted.tsv | sort -k 1,1n | cut -f2- \
  > gene_names_ens_sorted.tsv

28769 of the 36250 gene names have been mapped uniquely to an Ensembl ID

Check the list of genes supplied by Amanda:

sed -e 's/.*,//' gene_list.csv | grep ENSDARG | sort -u | wc -l
334

Two genes are repeated in the list:

grep ENSDARG00000077840 gene_list.csv
95,meis2a,ENSDARG00000077840
96,meis2b,ENSDARG00000077840

grep ENSDARG00000056218 gene_list.csv
13,tnika,ENSDARG00000056218
267,tnika,ENSDARG00000056218

sed -e 's/.*,//' gene_list.csv | grep ENSDARG | grep -c -f - gene_names_ens_sorted.tsv
329

So 5 genes missing:

for gene in `sed -e 's/.*,//' gene_list.csv | grep ENSDARG`; do echo -ne "$gene\t"; grep -c $gene gene_names_ens_sorted.tsv; done | grep 0$ | cut -f1 | grep -f - gene_list.csv
6,crhr2,ENSDARG00000092918
49,wnt8a,ENSDARG00000052910
206,atp2b4,ENSDARG00000044902
241,her4.2,ENSDARG00000056729
323,vav3b,ENSDARG00000073713

3 are mapped to a different Ensembl ID, one is a problematic gene, and another is missing entirely from the Lawson Lab annotation:

for gene in `sed -e 's/.*,//' gene_list.csv | grep ENSDARG`; do echo -ne "$gene\t"; grep -c $gene gene_names_ens_sorted.tsv; done | grep 0$ | cut -f1 \
  | grep -f - gene_list.csv | sed -e 's/,ENSDARG.*//' | sed -e 's/.*,//' | grep -if - gene_names_ens_sorted.tsv
her4.2	ENSDARG00000094426	TRUE
vav3b	ENSDARG00000075962	TRUE
wnt8a	ENSDARG00000078507	TRUE
crhr2		FALSE
crhr2.1		FALSE
crhr2.2		FALSE

for gene in `sed -e 's/.*,//' gene_list.csv | grep ENSDARG`; do echo -ne "$gene\t"; grep -c $gene gene_names_ens_sorted.tsv; done | grep 0$ | cut -f1 | grep -f - gene_list.csv | sed -e 's/,ENSDARG.*//' | sed -e 's/.*,//' | grep -if - v4_3_2geneinfo.txt | cut -f1-10
chr23	21455152	21471022	-	LL0000012050	her4.2	"hairy-related 4, tandem duplicate 2"	ENSDARG00000094426.4	30301	ZDB-GENE-060815-1
chr2	15776156	16159491	-	LL0000031813	vav3b	vav 3 guanine nucleotide exchange factor b	ENSDARG00000075962.6	559145	ZDB-GENE-070912-251
chr14	34490445	34497724	+	LL0000032103	wnt8a	"wingless-type MMTV integration site family, member 8a"	ENSDARG00000078507.6	553976	ZDB-GENE-980526-332
chr14	32783973	32788046	+	LL0000039783	crhr2	corticotropin releasing hormone receptor 2 [Source:NCBI gene;Acc:100002312]	ENSDARG00000092918.3	100002312
chr2	2503397	2511944	+	LL0000039784	crhr2	corticotropin releasing hormone receptor 2 [Source:NCBI gene;Acc:100002312]	ENSDARG00000103704.2	100002312
chr2	2563193	2608112	+	LL0000039785	crhr2	corticotropin releasing hormone receptor 2 [Source:NCBI gene;Acc:100002312]	ENSDARG00000062377.6	100335005

Back to RStudio:

```
# Subset to only genes with Ensembl IDs
ens.meta <- read_tsv("gene_names_ens_sorted.tsv", col_names=c("gene_name", "ensembl_id", "ensembl_mapped"))
genes.to.keep <- ens.meta$gene_name[ens.meta$ensembl_mapped]
glial_neural_5dpf_ens <- subset(glial_neural_5dpf, features = genes.to.keep)

# Function for creating MEX files
writeCounts <- function(out_dir, data, barcodes = colnames(data), gene.id = rownames(data), gene.symbol = rownames(data), feature.type = "Gene Expression", subset = 1:length(barcodes)) {
  require("R.utils")
  require("Matrix")

  if (file.exists(out_dir) || (dir.exists(out_dir) && length(list.files(out_dir)) > 0)) {
    stop("The specified output directory already exists! Not overwriting")
  }
  dir.create(out_dir, recursive = TRUE)

  write.table(
    data.frame(barcode = barcodes[subset]),
    gzfile(file.path(out_dir, "barcodes.tsv.gz")),
    sep = "\t", quote = FALSE,
    col.names = FALSE, row.names = FALSE
  )

  write.table(
    data.frame(
      gene_id = gene.id,
      gene_symbol = gene.symbol,
      feature_type = feature.type
    ),
    gzfile(file.path(out_dir, "features.tsv.gz")),
    sep = "\t", quote = FALSE,
    col.names = FALSE, row.names = FALSE
  )

  Matrix::writeMM(data[, subset], file.path(out_dir, "matrix.mtx"))
  R.utils::gzip(file.path(out_dir, "matrix.mtx"), remove = TRUE)
}

# Create MEX files
writeCounts(
  "reference_data",
  GetAssayData(glial_neural_5dpf_ens, assay="RNA", slot="counts"),
  gene.id = ens.meta$ensembl_id[ens.meta$ensembl_mapped],
  gene.symbol = rownames(glial_neural_5dpf_ens)
)

# Function for saving cell type annotation and bundling reference data
bundleOutputs <- function(out_dir, data, barcodes = colnames(data), cell_type = "cell_type", subset = 1:length(barcodes)) {
  data.table::fwrite(
    data.table::data.table(
      barcode = barcodes,
      annotation = unlist(data[[cell_type]])
    )[subset, ],
    file.path(out_dir, "annotations.csv")
  )

  bundle <- file.path(out_dir, paste0(basename(out_dir), ".zip"))

  utils::zip(
    bundle,
    list.files(out_dir, full.names = TRUE),
    zip = "zip"
  )

  if (file.info(bundle)$size / 1e6 > 2000) {
    warning("The output file is more than 2G and will need to be subset further.")
  }
}

# Create cell type annotation and bundle reference data
bundleOutputs(out_dir = "reference_data", data = glial_neural_5dpf_ens, cell_type="identity.sub")
```
