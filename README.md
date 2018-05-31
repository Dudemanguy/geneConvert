## geneConvert

**Note**: This is a WIP.

geneConvert is a simple approach to managing and converting gene annotations. geneConvert integrates SQL and parses NCBI's gene database with XML to fetch gene annotation and store them locally in an SQL database. Unlike other tools, the philosophy of geneConvert is to only fetch annotation once and store them instead of constantly querying databases for the same genes over and over again. When queries are run again, geneConvert checks the SQL table for the existence of the gene before scraping a webpage for annotations.

## Installation
Using the R package, ``devtools``, run ``devtools::install_github("Dudemanguy911/RNASeqSuite")``.

## Quick Usage
By default, geneConvert comes with tables for human, mouse and rat. However, you can add your own organism with the ``addNewOrganism`` function. Note that in order for the scraper to work properly, you must use the correct, scientific name. The main function, ``convert``, can return a dataframe with just the desired output annotation or it can return full data frame showing all annotations for  the  queried genes.
```
library(geneConvert)
genes <- readLines("gene_list_names.txt")
output <- convert(genes, "human", "symbol", "geneid")
```

## License
GPLv2 or later.
