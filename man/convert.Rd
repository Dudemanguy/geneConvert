\name{convert}
\alias{convert}

\title{Convert gene annotation}
\description{
Main function used to convert the inputed genes from one set of annotations to another.
}

\usage{
convert(genes, "Homo sapiens", "symbol", "description", full=FALSE)
}

\arguments{
  \item{genes}{A character vector of genes}
  \item{organism}{The name of the organism}
  \item{input}{The annotation type of the inputted genes}
  \item{output}{The annotation type of the output}
  \item{full}{A boolean controlling whether to return all annotations}
}

\details{
\code{organism} must be the proper, scientific name. However, there are three shortcuts hardcoded for \code{homo_sapiens}, \code{mus_musculus}, and \code{rattus_norvegicus}. These would be \code{human}, \code{mouse}, and \code{rat} respectively. If you want to convert annotation for a different organism, please add it first using the \code{addNewOrganism} function.

\code{input} and \code{output} must match the field names in the SQL table. Valid names are \code{symbol}, \code{geneid}, \code{description}, \code{geneloc}, \code{refseq}, \code{protein}, \code{ensembl}, and \code{date}. Note that \code{output} can accept multilength character vectors.

If enabled, \code{full} will return all annotations for the queried genes instead of just the specified output.
}

\examples{
output <- convert(genes, "human", "symbol", "ensembl")
}
