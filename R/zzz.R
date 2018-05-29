.onLoad <- function(libname, pkgname) {
	from <- system.file("extdata", "annotations.sqlite", package="geneConvert")
	directory <- file.path(path.expand("~"), ".config/geneConvert")
	if (!file.exists(directory)) {
		dir.create(directory)
	}
	to <- file.path(directory, "annotations.sqlite")
	if (!file.exists(to)) {
		file.copy(from, to, copy.mode=FALSE)
	}
}
