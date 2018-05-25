## geneConvert

**Note**: This is a WIP.

geneConvert is a simple approach to managing and converting gene annotations. geneConvert integrates SQL and parses large databases (such as NCBI) with XML to fetch gene annotation and store them locally in an SQL database. Unlike other tools, the philosophy of geneConvert is to only fetch annotation once and store them instead of constantly querying databases for the same genes over and over again. When queries are run again, geneConvert checks the SQL table for the existence of the gene before scraping a webpage for annotations.

## Installation
Using the R package, ``devtools`` run ``devtools::install_github("Dudemanguy911/RNASeqSuite")``.

## Quick Usage
At the moment, geneConvert only supports fetching human annotations from NCBI, but this is planned to be expanded upon in the future. The main function, ``convert``, can either return a character vector of matches or a full data frame.
```
library(geneConvert)
genes <- readLines("gene_list_names.txt")
output <- convert(genes, "Homo sapiens", "hgnc", "geneid")
```

## License
GPLv2 or later.
