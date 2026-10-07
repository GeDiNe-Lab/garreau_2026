# =============================================================================
# Fetch PubMed abstracts for the mechanism-oriented queries
# -----------------------------------------------------------------------------
# Input : PubMed_mechanism_queries.tsv  (Gene, Tier, Level, Query, ...)
# Output: PubMed_abstracts_<date>/
#           Direct/<GENE>_Direct.txt     PubMed "Abstract" text format,
#           Related/<GENE>_Related.txt   same layout as an export from the
#           General/<GENE>_General.txt   PubMed website (Save > Abstract)
#           search_log.tsv               gene, level, hits, date, sort, PMIDs
#           PubMed_mechanism_queries_filled.tsv  tracking grid with date + hits
# Records are ordered by relevance ("Best match"), like the website default.
# =============================================================================

library(rentrez)
library(data.table)

setwd(paste(here::here(), "projects", sep = "/"))

# Optional NCBI API key: 10 requests/s instead of 3
# Sys.setenv(ENTREZ_KEY = "your_key")
pause <- if (nzchar(Sys.getenv("ENTREZ_KEY"))) 0.12 else 0.35
sort_order <- "relevance" # "relevance" = Best match ; "pub_date" = Most recent
batch_size <- 200 # records per efetch call
max_hits <- 9999 # E-utilities limit for a single esearch

run_date <- format(Sys.Date())
out_dir <- paste0("article/transit_data/PubMed_abstracts_", run_date)
for (lv in c("Direct", "Related", "General")) {
    dir.create(file.path(out_dir, lv), recursive = TRUE, showWarnings = FALSE)
}

queries <- fread("article/transit_data/PubMed_mechanism_queries.tsv")

# efetch numbers records from 1 in each batch: shift the numbers of later batches
renumber <- function(txt, offset) {
    if (offset == 0) {
        return(txt)
    }
    m <- gregexpr("(^|\n\n\n)\\d+\\. ", txt)
    regmatches(txt, m) <- lapply(regmatches(txt, m), function(x) {
        n <- as.integer(sub("^\n*(\\d+)\\. $", "\\1", x))
        paste0(sub("\\d+\\. $", "", x), n + offset, ". ")
    })
    txt
}

fetch_abstracts <- function(ids) {
    chunks <- split(ids, ceiling(seq_along(ids) / batch_size))
    out <- character(length(chunks))
    offset <- 0
    for (k in seq_along(chunks)) {
        txt <- NULL
        for (attempt in 1:3) { # retry on transient NCBI errors
            txt <- tryCatch(
                entrez_fetch(
                    db = "pubmed", id = chunks[[k]],
                    rettype = "abstract", retmode = "text"
                ),
                error = function(e) NULL
            )
            Sys.sleep(pause)
            if (!is.null(txt)) break
            Sys.sleep(2)
        }
        if (is.null(txt)) stop("efetch failed for batch ", k)
        out[k] <- renumber(sub("\n+$", "", sub("^\n+", "", txt)), offset)
        offset <- offset + length(chunks[[k]])
    }
    paste(out, collapse = "\n\n\n")
}

log <- vector("list", nrow(queries))
for (i in seq_len(nrow(queries))) {
    g <- queries$Gene[i]
    lv <- queries$Level[i]
    message(sprintf("[%d/%d] %s - %s", i, nrow(queries), g, lv))

    s <- entrez_search(
        db = "pubmed", term = queries$Query[i],
        retmax = max_hits, sort = sort_order
    )
    Sys.sleep(pause)
    if (s$count > max_hits) {
        warning(
            g, " ", lv, ": ", s$count,
            " hits, only the first ", max_hits, " are saved"
        )
    }

    f <- file.path(out_dir, lv, paste0(g, "_", lv, ".txt"))
    if (length(s$ids) == 0) {
        writeLines(sprintf("No results (%s, PubMed query run on %s).", g, run_date), f)
    } else {
        writeLines(fetch_abstracts(s$ids), f)
    }

    log[[i]] <- data.table(
        Gene = g, Tier = queries$Tier[i], Level = lv,
        n_hits = s$count, n_saved = length(s$ids),
        date = run_date, sort = sort_order, file = f,
        pmids = paste(s$ids, collapse = ";"),
        query = queries$Query[i]
    )
}

log <- rbindlist(log)
fwrite(log, file.path(out_dir, "search_log.tsv"), sep = "\t")

# fill date and hit count in the tracking grid
# (empty columns are read as logical by fread: drop them and recreate them)
queries[, c("Date searched", "Hits") := NULL]
queries[, `Date searched` := run_date]
queries[log, Hits := i.n_hits, on = .(Gene, Level)]
setcolorder(queries, c(
    "Gene", "Tier", "Level", "Query", "PubMed link",
    "Date searched", "Hits"
))
fwrite(queries, file.path(out_dir, "PubMed_mechanism_queries_filled.tsv"), sep = "\t")
print(dcast(log, Gene ~ Level, value.var = "n_hits")[, .(Gene, Direct, Related, General)])
