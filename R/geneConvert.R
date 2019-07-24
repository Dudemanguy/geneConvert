addNewOrganism <- function(organism) {
	if (class(organism) != "character" || length(organism) != 1) {
		stop("The organism name must be a character vector of length 1.") 
	}
	path <- file.path(path.expand("~"), ".config/geneConvert/annotations.sqlite")
	con <- dbConnect(RSQLite::SQLite(), path)
	tables <- dbListTables(con)
	if (!(organism %in% tables)) {
		templateFrame <- data.frame(geneid=numeric(), symbol=character(), description=character(),
									geneloc=character(), transcript=character(), protein=character(),
									ensembl=character(), goterm=character(), date=character())
		dbWriteTable(con, organism, templateFrame)
	} else {
		stop(paste0("Table named ", organism, " already exists. Quitting."))
	}
	dbDisconnect(con)
}

argumentHandling <- function(argument, values) {
	valid_names <- colnames(values)
	argument_grep <- list()
	valid_return <- list()
	for (i in seq_along(argument)) {
		while (!(argument[[i]] %in% valid_names)) {
			argument_grep[[i]] <- grep(argument[[i]], valid_names, ignore.case=TRUE)
			valid_return[[i]] <- valid_names[argument_grep[[i]]]
			if (identical(valid_return[[i]], character(0))) {
				stop("No results found for the entered annotation. Please try again")
			} else {
				print(valid_return[[i]])
				argument[[i]] <- readline("Type in the name of the input id from the list. \n")
			}
		}
	}
	argument
}

convert <- function(genes, organism, input, output, scrape=TRUE, force=FALSE, full=FALSE, no_version=TRUE, query=3000) {
	path <- file.path(path.expand("~"), ".config/geneConvert/annotations.sqlite")
	con <- dbConnect(RSQLite::SQLite(), path)
	organism <- organismSelect(organism)
	values <- dbReadTable(con, organism)
	if (!(input %in% colnames(values))) {
		input <- argumentHandling(input, values)
	}
	if (no_version && (identical(input, "geneloc") || identical(input, "protein")
			|| identical(input, "transcript"))) {
		geneloc <- gsub("\\.\\d+", "", values[["geneloc"]])
		protein <- gsub("\\.\\d+", "", values[["protein"]])
		transcript <- gsub("\\.\\d+", "", values[["transcript"]])
		values[["geneloc"]] <- geneloc
		values[["protein"]] <- protein
		values[["transcript"]] <- transcript
	}
		
	if (any(!(output %in% colnames(values)))) {
		output <- argumentHandling(output, values)
	}
	if (identical(input, "symbol")) {
		genes <- matchCase(genes, organism)
	}
	if (identical(input, "transcript")) {
		genes <- genes[!(grepl("dup", genes))]
		genes <- genes[!(grepl("chr", genes))]
	}
	if (identical(input, "transcript") || identical(input, "protein")) {
		new_genes <- character()
		split_genes <- list()
		for (i in seq_along(genes)) {
			if (grepl(",", genes[i])) {
				gene_split <- unlist(strsplit(genes[i], ","))
				split_genes[[i]] <- gene_split
			} else if (any(grepl(genes[i], values[[input]]))) {
				next
			} else {
				new_genes[i] <- genes[i]
			}
		}
		split_genes <- unique(unlist(split_genes))
		for (i in seq_along(split_genes)) {
			if (any(grepl(split_genes[i], values[[input]]))) {
				next
			} else {
				new_genes[i+length(new_genes)] <- split_genes[i]
			}
		}
		new_genes <- new_genes[!is.na(new_genes)]
	} else {
		new_genes <- genes[!(genes %in% values[[input]])]
	}

	if (identical(input, "goterm")) {
		warning("Please be sure you have run gotermGrab to store the needed goterms before using convert.")
	}

	if (length(new_genes) > 0 && identical(scrape, TRUE) && identical(force, FALSE)) {
		scraped_genes <- scraper(new_genes, input, organism, query)
		values <- dbReadTable(con, organism)
	}
	if (identical(force, TRUE)) {
		scraped_genes <- scraper(genes, input, organism, query)
		values <- dbReadTable(con, organism)
	}

	if (identical(input, "transcript") || identical(input, "protein") || identical(input, "goterm")) {
		value_list <- list()
		split_list <- list()
		for (i in seq_along(genes)) {
			if (grepl(",", genes[i])) {
				gene_split <- unlist(strsplit(genes[i], ","))
				split_list[[i]] <- gene_split
			} else {
				value_list[[i]] <- values[grepl(genes[i], values[[input]]),]
			}
		}
		gene_split <- unique(unlist(split_list))
		for (i in seq_along(gene_split)) {
				 value_list[[i+length(value_list)]] <- values[grepl(gene_split[i], values[[input]]),]
		}
		values <- do.call(rbind, unique(value_list))
	} else {
		values <- values[values[[input]] %in% genes,]
	}
	dbDisconnect(con)

	if (no_version) {
		geneloc <- gsub("\\.\\d+", "", values[["geneloc"]])
		protein <- gsub("\\.\\d+", "", values[["protein"]])
		transcript <- gsub("\\.\\d+", "", values[["transcript"]])
		values[["geneloc"]] <- geneloc
		values[["protein"]] <- protein
		values[["transcript"]] <- transcript
	}

	if (full) {
		return (values)
	} else {
		return (values[c(input, output)])
	}
}

deleteOrganism <- function(organism) {
	if (class(organism) != "character" || length(organism) != 1) {
		stop("The organism name must be a character vector of length 1.") 
	}
	path <- file.path(path.expand("~"), ".config/geneConvert/annotations.sqlite")
	con <- dbConnect(RSQLite::SQLite(), path)
	tables <- dbListTables(con)
	if (!(organism %in% tables)) {
		stop(paste0("Table named ", organism, " does not exist."))
	} else {
		confirmation <- readline(paste0("Are you sure you want to delete ", organism, " ? Type 'y' to confirm.\n"))
		if (identical(confirmation, "y")) {
			dbRemoveTable(con, organism)
		}
	}
	dbDisconnect(con)
}

forceUpdate <- function(organism, query=3000) {
	organism <- organismSelect(organism)
	confirmation <- readline(paste0("This will delete all records in the inputted table and then rescrape annotations. Type 'y' to confirm.\n"))
	if (identical(confirmation, "y")) {
		path <- file.path(path.expand("~"), ".config/geneConvert/annotations.sqlite")
		con <- dbConnect(RSQLite::SQLite(), path)
		values <- dbReadTable(con, organism)
		dbSendQuery(con, paste0("DELETE FROM ", organism))
		genes <- unique(values[["symbol"]])
		invisible(scraper(genes, "symbol", organism, query))
	}
}

gotermGrab <- function(genes, organism) {
	print("Note that this function requires biomaRt organism names.")
	if ((!"biomaRt" %in% installed.packages())) {
		warning("R package 'biomaRt' not found in the system. Please install the missing packages to use geneConvert functions with goterms.")
	}
	ensembl <- biomaRt::useMart("ensembl")
	ensembl <- biomaRt::useDataset(organism, mart=ensembl)
	output <- biomaRt::getBM(attributes=c("entrezgene", "go_id"), filters="entrezgene", values=genes, mart=ensembl)
	path <- file.path(path.expand("~"), ".config/geneConvert/annotations.sqlite")
	con <- dbConnect(RSQLite::SQLite(), path)
	sql_organism <- readline("Input regular organism name from database. \n")
	values <- dbReadTable(con, sql_organism)
	index <- unique(output[,1])
	for (i in seq_along(index)) {
		sub_output <- output[output[,1] == index[i],]
		go <- sub_output[,2]
		go <- go[go != ""]
		go <- paste(go, collapse=",")
		rs <- dbSendStatement(con, paste0("UPDATE ", sql_organism, " SET goterm = '", go, "' WHERE geneid = ", index[i], ";"))
		dbClearResult(rs)
	}
	dbDisconnect(con)
}

matchCase <- function(genes, organism) {
	if (identical(organism, "homo_sapiens")) {
		toupper(genes)
	}
	if (identical(organism, "mus_musculus") || identical(organism, "rattus_norvegicus")) {
		split <- strsplit(genes, " ")
		genes <- paste0(toupper(substring(split, 1, 1)), tolower(substring(split, 2)))
		if (any(grepl("-", genes))) {
			hyphen_genes <- genes[grep("-", genes)]
			no_hyphen_genes <- genes[genes != hyphen_genes]
			hyphen_genes <- toupper(hyphen_genes)
			genes <- c(no_hyphen_genes, hyphen_genes)
		}
	}
	genes
}

organismSelect <- function(organism) {
	if (identical(organism, "human")) {
		organism <- "homo_sapiens"
	}
	if (identical(organism, "mouse")) {
		organism <- "mus_musculus"
	}
	if (identical(organism, "rat")) {
		organism <- "rattus_norvegicus"
	}
	organism
}

retrieveURLs <- function(genes, input, organism) {
	searchURLs <- c()
	for (i in seq_along(genes)) {
		if (identical(input, "geneloc")) {
			warning("Unable to scrape missing genes based on only geneloc.")
			break
		}
		if (identical(input, "symbol") || identical(input, "description")) {
			searchURLs[[i]] <- paste0("https://www.ncbi.nlm.nih.gov/gene?term=(", genes[[i]], "[gene])%20AND%20(", organism, "[orgn])")
			searchURLs[[i]] <- gsub(" ", "%20", searchURLs[[i]])
		}
		if (identical(input, "geneid")) {
			searchURLs[[i]] <- paste0("https://www.ncbi.nlm.nih.gov/gene/", genes[[i]])
		}
		if (identical(input, "transcript") || identical(input, "protein") || identical(input, "ensembl")) {
			searchURLs[[i]] <- paste0("https://www.ncbi.nlm.nih.gov/gene?term=", genes[[i]])
		}
	}
	searchURLs
}

scraper <- function(genes, input, organism, query=3000) {
	path <- file.path(path.expand("~"), ".config/geneConvert/annotations.sqlite")
	con <- dbConnect(RSQLite::SQLite(), path)
	searchURLs <- retrieveURLs(genes, input, organism)
	data <- list()
	success <- function(res) {
		data <<- c(data, list(res))
	}
	failure <- function(msg) {
		cat("Request failed.", msg, "\n")
	}
	if (identical(query, "max")) {
		query <- length(searchURLs)
	}
	for (i in seq_along(searchURLs)) {
		curl_fetch_multi(searchURLs[[i]], success, failure)
		if (length(searchURLs) > query) {
			if (identical((i %% query), 0)) {
				print(paste("Queries", (i+1) - query, "-", i, "sent."))
				multi_run()
			}
			if (i > (floor(length(searchURLs)/query) * query)) {
				if (identical(i, length(searchURLs))) {
					print(paste("Queries", ((floor(length(searchURLs)/query)*query) + 1), "-", length(genes), "sent."))
					multi_run()
				}
			}
		} else {
			if (identical(i, length(searchURLs))) {
				print(paste(length(searchURLs), "queries sent."))
				multi_run()
			}
		}
	}
	for (i in seq_along(data)) {
		xdata <- rawToChar(data[[i]]$content)
		doc <- htmlParse(xdata, encoding="UTF-8")
		values <- dbReadTable(con, organism)
		print(paste0("Scraping ", genes[[i]]))
		exact_match <- grepl("Full Report", xpathApply(doc, "/*", xmlValue))
		search_results <- grepl("Search results", xpathApply(doc, "/*", xmlValue))
		if (!exact_match && !search_results) {
			message("Invalid gene or organism name inputted; skipping")
			next
		}
		if (search_results) {
			result_values <- xpathApply(doc, "//a[contains(@href, '/gene/')]", xmlValue)
			index <- match(genes[[i]], result_values)
			if (is.na(index)) {
				message("Gene not found. Skipping")
				next
			}
			result_xml <- as(xpathApply(doc, "//a[contains(@href, '/gene/')]")[[index]], "character")
			result_split <- strsplit(result_xml, "ref")[[1]][2]
			id <- as.character(gsub("\\D", "", result_split))
			newURL <- paste0("https://www.ncbi.nlm.nih.gov/gene/", id)
			xdata <- getURL(newURL)
			doc <- htmlParse(xdata, encoding="UTF-8")
		}
		geneid <- unlist(xpathApply(doc, "//*[@class='geneid']", xmlValue))
		geneid <- trimws(gsub(".*\\:", "", geneid))
		geneid <- gsub(",.*", "", geneid)
		if (geneid %in% values[["geneid"]]) {
			next
		}
		symbol <- unlist(xpathApply(doc, "//span[@class='gn']", xmlValue))
		description <- unlist(xpathApply(doc, "//title", xmlValue))
		description <- trimws(gsub("\\[.*", "", description))
		geneloc <- unlist(xpathApply(doc, "//p[@class='withnote margin_t2em']/strong", xmlValue))
		geneloc <- trimws(gsub(".*-", "", geneloc))
		if (identical(geneloc, character(0))) {
			geneloc <- NA
		}
		transcript_nm <- unlist(xpathApply(doc, "//p/a[contains(@href, 'NM')]", xmlValue))
		transcript_nr <- unlist(xpathApply(doc, "//p/a[contains(@href, 'NR')]", xmlValue))
		transcript <- paste(c(transcript_nm, transcript_nr), collapse=",")
		if (identical(transcript, "")) {
			transcript <- NA
		}
		protein <- paste(xpathApply(doc, "//p/a[contains(@href, 'protein/NP_')]", xmlValue), collapse=",")
		if (identical(protein, "")) {
			protein <- NA
		}
		ensembl <- unlist(xpathApply(doc, "//dd/a[@class='genome-browser-link']", xmlValue))
		ensembl <- gsub(".*\\:", "", ensembl)
		if (identical(ensembl, character(0))) {
			ensembl <- NA
		}
		date <- unlist(xpathApply(doc, "//*[@class='geneid']", xmlValue))
		date <- trimws(gsub(".*\n", "", date))
		values <- data.frame(geneid, symbol, description, geneloc, transcript,  protein, ensembl, date)
		dbWriteTable(con, organism, values, overwrite=FALSE, append=TRUE)
	}
	dbDisconnect(con)
}

updateFields <- function(organism) {
	organism <- organismSelect(organism)
	localpath <- file.path(path.expand("~"), ".config/geneConvert/annotations.sqlite")
	old <- dbConnect(RSQLite::SQLite(), localpath)
	sourcepath <- system.file("extdata/annotations.sqlite", package="geneConvert")
	new <- dbConnect(RSQLite::SQLite(), sourcepath)
	oldFields <- dbListFields(old, organism)
	newFields <- dbListFields(new, organism)
	diff <- newFields[!(newFields %in% oldFields)]
	if (!(identical(diff, character(0)))) {
		rs <- dbSendStatement(old, paste("ALTER TABLE", organism, "ADD", diff, ";"))
		dbClearResult(rs)
	}
	dbDisconnect(old)
	dbDisconnect(new)
}

updateTables <- function() {
	localpath <- file.path(path.expand("~"), ".config/geneConvert/annotations.sqlite")
	old <- dbConnect(RSQLite::SQLite(), localpath)
	oldTables <- dbListTables(old)
	sourcepath <- system.file("extdata/annotations.sqlite", package="geneConvert")
	new <- dbConnect(RSQLite::SQLite(), sourcepath)
	newTables <- dbListTables(new)
	updatedTables <- newTables[!(newTables %in% oldTables)]
	if (identical(length(updatedTables), 0)) {
		message("No new tables found")
	}
	if (length(updatedTables) > 0) {
		templateFrame <- data.frame(geneid=character(), symbol=character(), description=character(),
									geneloc=character(), transcript=character(), protein=character(),
									ensembl=character(), goterm=character(), date=character())
		for (i in seq_along(updatedTables)) {
			dbWriteTable(old, updatedTables[[i]], templateFrame)
		}
	}
	dbDisconnect(old)
	dbDisconnect(new)
}
