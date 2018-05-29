argumentHandling <- function(argument, values) {
	valid <- colnames(values)
	while (!(argument %in% valid)) {
		argument_grep <- grep(argument, valid, ignore.case=TRUE)
		valid <- valid[argument_grep]
		if (length(valid) == 0) {
			stop("No results found for the entered annotation. Please try again")
		} else {
			print(valid)
			argument <- readline("Type in the name of the input id from the list. \n")
		}
	}
	argument
}

convert <- function(genes, organism, input, output, full_table=FALSE) {
	path <- file.path(path.expand("~"), ".config/geneConvert/annotations.sqlite")
	con <- dbConnect(RSQLite::SQLite(), path)
	values <- dbReadTable(con, organism)
	if (!(input %in% colnames(values))) {
		input <- argumentHandling(input, values)
	}
	if (!(output %in% colnames(values))) {
		output <- argumentHandling(output, values)
	}
	new_genes <- unique(genes[!(genes %in% values[[input]])])
	if (length(new_genes) > 0) {
		scraped_genes <- scraper(new_genes, organism)
		values <- rbind(values, scraped_genes)
	}

	dbDisconnect(con)

	if (full_table) {
		return (values)
	} else {
		m <- match(genes, values[[input]])
		output_values <- values[[output]][m]
		return (output_values)
	}
}

scraper <- function(genes, organism) {
	total_list <- list()
	path <- file.path(path.expand("~"), ".config/geneConvert/annotations.sqlite")
	con <- dbConnect(RSQLite::SQLite(), path)
	for (i in seq_along(genes)) {
		fileURL <- paste0("https://www.ncbi.nlm.nih.gov/gene?term=(", genes[[i]], "[gene])%20AND%20(", organism, "[orgn])")
		fileURL <- gsub(" ", "%20", fileURL)
		xData <- getURL(fileURL)
		doc <- htmlParse(xData, encoding="UTF-8")
		print(paste0("Scraping ", genes[[i]]))
		valid_URL <- grepl("Full Report", xmlValue(getNodeSet(doc, "/*")[[1]]))
		search_results <- grepl("Search results", xmlValue(getNodeSet(doc, "/*")[[1]]))
		if (!valid_URL & !search_results) {
			message("Invalid gene or organism name inputted; skipping")
			next
		}
		if (search_results) {
			results <- getNodeSet(doc, "//a[contains(@href, '/gene/')]")
			results_list <- list()
			for (n in seq_along(results)) {
				results_list[[n]] <- as(results[[n]], "character")
			}
			results_grep <- as.character(results_list[grep(genes[[i]], results_list)])
			results_split <- strsplit(results_grep, "ref")[[1]][2]
			id <- as.character(gsub("\\D", "", results_split))
			newURL <- paste0("https://www.ncbi.nlm.nih.gov/gene/", id)
			xData <- getURL(newURL)
			doc <- htmlParse(xData, encoding="UTF-8")
		}
		description <- xmlValue(getNodeSet(doc, "//title")[[1]][1]$text)
		description <- trimws(gsub("\\[.*", "", description))
		date <- xmlValue(getNodeSet(doc, "//*[@class='geneid']")[[1]][1]$text)
		date <- trimws(gsub(".*\n", "", date))
		hgnc <- xmlValue(getNodeSet(doc, "//span[@class='gn']")[[1]][1]$text)
		geneid <- xmlValue(getNodeSet(doc, "//*[@class='geneid']")[[1]][1]$text)
		geneid <- trimws(gsub(".*\\:", "", geneid))
		geneid <- gsub(",.*", "", geneid)
		description <- xmlValue(getNodeSet(doc, "//title")[[1]][1]$text)
		description <- trimws(gsub("\\[.*", "", description))
		geneloc <- xmlValue(getNodeSet(doc, "//p[@class='withnote margin_t2em']/strong")[[1]][1]$text)
		geneloc <- trimws(gsub(".*-", "", geneloc))
		refseq_length <- getNodeSet(doc, "//p/a[contains(@href, 'NM')]")
		if (length(refseq_length) > 0) {
			refseq <- list()
			for (j in seq_along(refseq_length)) {
				refseq[[j]] <- xmlValue(refseq_length[[j]][1]$text)
			}
			refseq <- unlist(refseq)
		} else {
			refseq <- NA
			refseq_length <- NA
		}
		protein_length <- getNodeSet(doc, "//p/a[contains(@href, 'protein/NP_')]")
		if (length(protein_length) > 0) {
			protein <- list()
			for (k in seq_along(protein_length)) {
				protein[[k]] <- xmlValue(protein_length[[k]][1]$text)
			}
			protein <- unlist(protein)
		} else {
			protein <- NA
			protein_length <- NA
		}
		ensembl_check <- getNodeSet(doc, "//dd/a[@class='genome-browser-link']")
		if (length(ensembl_check) > 0) {
			ensembl <- xmlValue(ensembl_check[[1]][1]$text)
			ensembl <- gsub(".*\\:", "", ensembl)
		} else {
			ensembl <- NA
		}
		total_list[[i]] <- data.frame(hgnc, geneid, description, geneloc, refseq, protein, ensembl, date)
	}
	total_frame <- do.call(rbind, total_list)
	if (length(total_frame) > 0) {
		dbWriteTable(con, organism, total_frame, overwrite=FALSE, append=TRUE)
	}
	total_frame
}
