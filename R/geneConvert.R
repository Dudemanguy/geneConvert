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

convert <- function(genes, organism, input, output, full=FALSE, no_version=TRUE) {
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
	genes <- matchCase(genes, organism)
	new_genes <- unique(genes[!(genes %in% values[[input]])])
	if (length(new_genes) > 0) {
		scraped_genes <- scraper(new_genes, input, organism)
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
		if (confirmation == 'y') {
			dbRemoveTable(con, organism)
		}
	}
	dbDisconnect(con)
}

forceUpdate <- function(organism) {
	organism <- organismSelect(organism)
	confirmation <- readline(paste0("This will delete all records in the inputted table and then rescrape annotations. Type 'y' to confirm.\n"))
	if (confirmation == "y") {
		path <- file.path(path.expand("~"), ".config/geneConvert/annotations.sqlite")
		con <- dbConnect(RSQLite::SQLite(), path)
		values <- dbReadTable(con, organism)
		dbSendQuery(con, paste0("DELETE FROM ", organism))
		genes <- unique(values[["symbol"]])
		invisible(scraper(genes, "symbol", organism))
	}
}

matchCase <- function(genes, organism) {
	if (organism == "homo_sapiens") {
		toupper(genes)
	}
	if (identical(organism, "mus_musculus") || identical(organism, "rattus_norvegicus")) {
		split <- strsplit(genes, " ")
		genes <- paste0(toupper(substring(split, 1, 1)), tolower(substring(split, 2)))
	}
	genes
}

organismSelect <- function(organism) {
	if (organism == "human") {
		organism <- "homo_sapiens"
	}
	if (organism == "mouse") {
		organism <- "mus_musculus"
	}
	if (organism == "rat") {
		organism <- "rattus_norvegicus"
	}
	organism
}

scraper <- function(genes, input, organism) {
	total_list <- list()
	path <- file.path(path.expand("~"), ".config/geneConvert/annotations.sqlite")
	con <- dbConnect(RSQLite::SQLite(), path)
	for (i in seq_along(genes)) {
		if (input == "geneloc") {
			warning("Unable to scrape missing genes based on only geneloc.")
			break
		}
		if (input == "symbol" || input == "description") {
			searchURL <- paste0("https://www.ncbi.nlm.nih.gov/gene?term=(", genes[[i]], "[gene])%20AND%20(", organism, "[orgn])")
			searchURL <- gsub(" ", "%20", searchURL)
		}
		if (input == "geneid") {
			searchURL <- paste0("https://www.ncbi.nlm.nih.gov/gene/", genes[[i]])
		}
		if (input == "transcript" || input == "protein" || input == "ensembl") {
			searchURL <- paste0("https://www.ncbi.nlm.nih.gov/gene?term=", genes[[i]])
		}
		xdata <- getURL(searchURL)
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
		transcript <- unlist(xpathApply(doc, "//p/a[contains(@href, 'NM')]", xmlValue))
		if (identical(transcript, NULL)) {
			transcript <- unlist(xpathApply(doc, "//p/a[contains(@href, 'NR')]", xmlValue))
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
		date <- unlist(xpathApply(doc, "//*[@class='geneid']", xmlValue))
		date <- trimws(gsub(".*\n", "", date))
		total_list[[i]] <- data.frame(symbol, geneid, description, geneloc, transcript, protein, ensembl, date)
	}
	total_frame <- do.call(rbind, total_list)
	if (length(total_frame) > 0) {
		dbWriteTable(con, organism, total_frame, overwrite=FALSE, append=TRUE)
	}
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
	if (length(updatedTables) == 0) {
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
