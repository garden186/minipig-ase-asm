custom_lib <- Sys.getenv("ASM_GO_RLIBRARY", "")
if (nzchar(custom_lib)) .libPaths(c(custom_lib, .libPaths()))
options(stringsAsFactors=FALSE, warn=1, enrichment_force_universe=FALSE)
suppressPackageStartupMessages({
    library(jsonlite)
})
need <- function(ok, message) if (!isTRUE(ok)) stop(message, call.=FALSE)
read_tab <- function(path) {
    con <- if (grepl("[.]gz$", path)) gzfile(path, "rt") else file(path, "rt")
    on.exit(close(con)); read.delim(con, check.names=FALSE, na.strings="", quote="", comment.char="")
}
write_tab <- function(x, path) {
    con <- if (grepl("[.]gz$", path)) gzfile(path, "wt") else file(path, "wt")
    on.exit(close(con)); write.table(x, con, sep="\t", row.names=FALSE, quote=FALSE, na="NA")
}
write_json <- function(x, path) jsonlite::write_json(x, path, pretty=TRUE, auto_unbox=TRUE, na="null", digits=16)
sorted_unique <- function(x) sort(unique(as.character(x[!is.na(x) & nzchar(x)])))
empty_result <- function() data.frame(ID=character(), Description=character(), GeneRatio=character(),
    BgRatio=character(), pvalue=numeric(), p.adjust=numeric(), qvalue=numeric(), geneID=character(), Count=integer())
go_sets <- function(t2g, background) {
    t2g <- unique(t2g[t2g$gene %in% background, c("term", "gene")])
    lapply(split(t2g$gene, t2g$term), unique)
}

# This oracle uses R's separate phyper/p.adjust implementations and independently
# reconstructs every membership and denominator; no replacement enrichment code.
validate_result <- function(object, t2g, background, candidate, min_size=10, max_size=500) {
    bg <- intersect(sorted_unique(background), sorted_unique(t2g$gene))
    query <- intersect(sorted_unique(candidate), bg)
    sets <- go_sets(t2g, bg)
    sets <- sets[lengths(sets)>=min_size & lengths(sets)<=max_size]
    hits <- lapply(sets, intersect, y=query)
    expected_ids <- names(sets)[lengths(hits)>0]
    result <- if (is.null(object)) empty_result() else object@result
    need(!anyDuplicated(result$ID) && setequal(result$ID, expected_ids), "GO hypothesis family differs from expected positive-overlap package family")
    if (!is.null(object)) need(setequal(object@universe, bg), "Effective GO universe mismatch")
    max_p_error <- 0; max_q_error <- 0
    if (nrow(result)) {
        K <- lengths(sets[result$ID]); k <- lengths(hits[result$ID]); M <- length(bg); n <- length(query)
        need(all(result$Count==k), "Term overlap count mismatch")
        need(all(result$GeneRatio==paste0(k,"/",n)) && all(result$BgRatio==paste0(K,"/",M)), "GO ratio denominator mismatch")
        for (j in seq_len(nrow(result))) need(setequal(strsplit(result$geneID[j], "/", fixed=TRUE)[[1]], hits[[result$ID[j]]]), "Contributing gene mismatch")
        oracle_p <- phyper(k-1, K, M-K, n, lower.tail=FALSE)
        need(all(is.finite(result$pvalue) & result$pvalue>0 & result$pvalue<=1), "Invalid or underflowed GO P")
        max_p_error <- max(abs(result$pvalue-oracle_p))
        need(all(abs(result$pvalue-oracle_p)<=1e-12+1e-9*oracle_p), "Package P differs from independent hypergeometric calculation")
        oracle_q <- p.adjust(oracle_p, method="BH")
        max_q_error <- max(abs(result$p.adjust-oracle_q))
        need(all(abs(result$p.adjust-oracle_q)<=1e-12+1e-9*oracle_q), "Package BH differs from complete within-job BH replay")
    }
    list(status="PASS_ALL_TERM_COUNTS_P_AND_WITHIN_JOB_BH", tested_terms=nrow(result),
         size_eligible_terms=length(sets), zero_overlap_terms_not_in_package_BH=length(sets)-nrow(result),
         background_GO_genes=length(bg), selected_GO_genes=length(query),
         significant_BH05=sum(result$p.adjust<=.05), max_P_error=max_p_error,max_BH_error=max_q_error)
}
