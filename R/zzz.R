.onLoad <- function(libname, pkgname) {
	from <- system.file("extdata", "annotations.sqlite", package="geneConvert")
	directory <- file.path(path.expand("~"), ".config/geneConvert")
	dir.create(directory)
	to <- file.path(directory, "annotations.sqlite")
	file.copy(from, to, copy.mode=FALSE)
}
