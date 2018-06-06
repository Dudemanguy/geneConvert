addNewOrganism <- function(organism) {
	if (class(organism) != "character" || length(organism) != 1) {
		stop("The organism name must be a character vector of length 1.") 
	}
	path <- file.path(path.expand("~"), ".config/geneConvert/annotations.sqlite")
	con <- dbConnect(RSQLite::SQLite(), path)
	tables <- dbListTables(con)
	if (!(organism %in% tables)) {
		templateFrame <- data.frame(symbol=character(), geneid=character(), description=character(),
									geneloc=character(), transcript=character(), protein=character(),
									ensembl=character(), date=character())
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

convert <- function(genes, organism, input, output, scrape=TRUE, full=FALSE, no_version=TRUE, query=3000) {
	path <- file.path(path.expand("~"), ".config/geneConvert/annotations.sqlite")
	con <- dbConnect(RSQLite::SQLite(), path)
	organism <- organismSelect(organism)
	values <- dbReadTable(con, organism)
	if (!(input %in% colnames(values))) {
		input <- argumentHandling(input, values)
	}
	if (no_version && (identical(input, "geneloc") || identical(input, "protein")
			|| identical(input, "transcript"))) {
		geneloc <- gsub("\\..*", "", values[["geneloc"]])
		protein <- gsub("\\..*", "", values[["protein"]])
		transcript <- gsub("\\..*", "", values[["transcript"]])
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
	new_genes <- unique(genes[!(genes %in% values[[input]])])
	if (length(new_genes) > 0 && identical(scrape, TRUE)) {
		scraped_genes <- scraper(new_genes, input, organism, query)
		values <- rbind(values, scraped_genes)
	}
	dbDisconnect(con)

	if (no_version) {
		geneloc <- gsub("\\..*", "", values[["geneloc"]])
		protein <- gsub("\\..*", "", values[["protein"]])
		transcript <- gsub("\\..*", "", values[["transcript"]])
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
	total_list <- list()
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
		if (identical((i %% query), 0)) {
			print(paste("Queries", (i+1) - query, "-", i, "sent."))
			multi_run()
		}
		if (identical(i, length(searchURLs))) {
			print(paste(length(searchURLs), "queries sent."))
		}
	}
	multi_run()
	for (i in seq_along(data)) {
		xdata <- rawToChar(data[[i]]$content)
		doc <- htmlParse(xdata, encoding="UTF-8")
		print(paste0("Scraping ", genes[[i]]))
		exact_match <- grepl("Full Report", xpathApply(doc, "/*", xmlValue))
		search_results <- grepl("Search results", xpathApply(doc, "/*", xmlValue))
		if (!exact_match & !search_results) {
			message("Invalid gene or organism name inputted; skipping")
			next
		}
		if (search_results) {
			result_values <- xpathApply(doc, "//a[contains(@href, '/gene/')]", xmlValue)
			index <- match(genes[[i]], result_values)
			result_xml <- as(xpathApply(doc, "//a[contains(@href, '/gene/')]")[[index]], "character")
			result_split <- strsplit(result_xml, "ref")[[1]][2]
			id <- as.character(gsub("\\D", "", result_split))
			newURL <- paste0("https://www.ncbi.nlm.nih.gov/gene/", id)
			xdata <- getURL(newURL)
			doc <- htmlParse(xdata, encoding="UTF-8")
		}
		symbol <- unlist(xpathApply(doc, "//span[@class='gn']", xmlValue))
		geneid <- unlist(xpathApply(doc, "//*[@class='geneid']", xmlValue))
		geneid <- trimws(gsub(".*\\:", "", geneid))
		geneid <- gsub(",.*", "", geneid)
		description <- unlist(xpathApply(doc, "//title", xmlValue))
		description <- trimws(gsub("\\[.*", "", description))
		geneloc <- unlist(xpathApply(doc, "//p[@class='withnote margin_t2em']/strong", xmlValue))
		geneloc <- trimws(gsub(".*-", "", geneloc))
		if (identical(geneloc, character(0))) {
			geneloc <-  NA
		}
		transcript <- unlist(xpathApply(doc, "//p/a[contains(@href, 'NM')]", xmlValue))
		if (identical(transcript, NULL)) {
			transcript <- unlist(xpathApply(doc, "//p/a[contains(@href, 'NR')]", xmlValue))
		}
		if (identical(transcript, NULL)) {
			transcript <- NA
		}
		protein <- unlist(xpathApply(doc, "//p/a[contains(@href, 'protein/NP_')]", xmlValue))
		if (identical(protein, NULL)) {
			protein <- NA
		}
		ensembl <- unlist(xpathApply(doc, "//dd/a[@class='genome-browser-link']", xmlValue))
		ensembl <- gsub(".*\\:", "", ensembl)
		if (identical(ensembl, character(0))) {
			ensembl <- NA
		}
		if (length(ensembl) != length(transcript) &&
			length(ensembl) > 1) {
			for (i in ((length(ensembl) + 1):length(transcript))) {
				ensembl[[i]] <- NA
			}
		}
		date <- unlist(xpathApply(doc, "//*[@class='geneid']", xmlValue))
		date <- trimws(gsub(".*\n", "", date))
		total_list[[i]] <- data.frame(symbol, geneid, description, geneloc, transcript, protein, ensembl, date)
	}
	total_frame <- do.call(rbind, total_list)
	if (length(total_frame) > 0) {
		dbWriteTable(con, organism, total_frame, overwrite=FALSE, append=TRUE)
	}
	dbDisconnect(con)
	total_frame
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
		templateFrame <- data.frame(symbol=character(), geneid=character(), description=character(),
									geneloc=character(), transcript=character(), protein=character(),
									ensembl=character(), date=character())
		for (i in seq_along(updatedTables)) {
			dbWriteTable(old, updatedTables[[i]], templateFrame)
		}
	}
	dbDisconnect(old)
	dbDisconnect(new)
}
