# 10_prepare_gse79396.R
# Parse GSE79396 series matrix + GPL13158 annotation into a clean workspace dataset.
# Outputs: data/GSE79396/vaccine_data.rds, samples_all.csv, annot_clean.csv
suppressPackageStartupMessages({
  library(GEOquery)
  library(data.table)
})

base_dir <- Sys.getenv("TER_ROOT", getwd())
src_series <- Sys.getenv("GSE79396_SERIES_MATRIX",
                         file.path(base_dir, "data/GSE79396/GSE79396_series_matrix.txt.gz"))
if (!file.exists(src_series)) {
  stop("GSE79396 series matrix not found: ", src_series,
       "\nSet GSE79396_SERIES_MATRIX or place the file under data/GSE79396/.")
}
src_annot  <- file.path(base_dir, "data/GPL13158/GPL13158.annot.gz")
out_dir    <- file.path(base_dir, "data/GSE79396")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---- 1. expression + sample metadata ---------------------------------------
gse <- getGEO(filename = src_series)
if (is.list(gse) && !is(gse, "ExpressionSet")) gse <- gse[[1]]
expr <- exprs(gse)
pd   <- as.data.frame(pData(gse))
cat("expr dims:", dim(expr), "\n")
cat("pData columns:\n"); print(colnames(pd))

get_char <- function(x, pattern) {
  out <- rep(NA_character_, length(x))
  hit <- grepl(pattern, x)
  out[hit] <- sub(pattern, "", x[hit])
  trimws(out)
}
ch_cols <- grep("^characteristics_ch1", colnames(pd), value = TRUE)
ch <- do.call(cbind, lapply(ch_cols, function(cc) as.character(pd[[cc]])))
cat("\ncharacteristics columns:", length(ch_cols), "\n")
cat("example row:\n"); print(ch[1, ])

pick_field <- function(field) {
  vals <- rep(NA_character_, nrow(ch))
  for (j in seq_len(ncol(ch))) {
    pat <- paste0("^", field, ": ")
    hit <- grepl(pat, ch[, j])
    vals[hit] <- sub(pat, "", ch[hit, j])
  }
  vals
}

desc <- as.character(pd$description)
subject_from_desc <- sub("_.*$", "", desc)
day_from_desc <- as.numeric(sub("^.*_D", "", desc))

samples <- data.frame(
  gsm        = as.character(pd$geo_accession),
  title      = as.character(pd$title),
  source     = as.character(pd$source_name_ch1),
  subject_id = pick_field("subject id"),
  vaccine    = pick_field("vaccine"),
  cohort     = pick_field("cohort"),
  visit      = pick_field("visit"),
  sex        = pick_field("Sex"),
  age        = as.numeric(pick_field("age")),
  age_group  = pick_field("age group"),
  tissue     = pick_field("tissue/cell type"),
  day        = day_from_desc,
  subject_desc = subject_from_desc,
  stringsAsFactors = FALSE
)
# consistency checks between description and characteristics
stopifnot(all(samples$subject_id == samples$subject_desc))
stopifnot(!anyNA(samples$day))
cat("\n#samples:", nrow(samples),
    "| #subjects:", length(unique(samples$subject_id)), "\n")
print(table(samples$age_group, useNA = "ifany"))
print(table(samples$day))
print(table(samples$cohort))
print(table(samples$sex))

## ---- 2. annotation ----------------------------------------------------------
pre <- readLines(gzfile(src_annot), n = 40)
hdr <- which(startsWith(pre, "ID\t") | startsWith(pre, "\"ID\"\t"))
stopifnot(length(hdr) >= 1)
ann <- fread(src_annot, skip = hdr[1] - 1, sep = "\t", header = TRUE,
             quote = "", na.strings = c("", "NA", "---"))
cat("\nannot dims:", dim(ann), "\n")
cat("annot columns:", paste(head(colnames(ann), 8), collapse = " | "), "\n")
ann_clean <- data.frame(
  probe  = ann[["ID"]],
  symbol = ann[["Gene symbol"]],
  title  = ann[["Gene title"]],
  entrez = if ("Gene ID" %in% colnames(ann)) ann[["Gene ID"]] else NA,
  stringsAsFactors = FALSE
)
cat("annot rows:", nrow(ann_clean),
    "| with symbol:", sum(!is.na(ann_clean$symbol) & ann_clean$symbol != ""),
    "\n")
cat("probes in expr also in annot:",
    sum(rownames(expr) %in% ann_clean$probe), "/", nrow(expr), "\n")

## ---- 3. save ----------------------------------------------------------------
saveRDS(list(expr = expr, samples = samples, annot = ann_clean),
        file.path(out_dir, "vaccine_data.rds"))
fwrite(samples, file.path(out_dir, "samples_all.csv"))
fwrite(ann_clean, file.path(base_dir, "data/GPL13158/annot_clean.csv"))

## ---- 4. compare with v1 sample set ------------------------------------------
v1_path <- Sys.getenv("TER_V1_METADATA",
                      file.path(base_dir, "data/v1_tables/sample_metadata_complete.csv"))
if (file.exists(v1_path)) {
  v1 <- fread(v1_path)
  cat("\n--- v1 metadata comparison ---\n")
  cat("v1 samples:", nrow(v1), "| v1 subjects:", length(unique(v1$subject_id)), "\n")
  cat("same GSM set:", setequal(v1$geo_accession, samples$gsm), "\n")
  m <- merge(v1, samples, by.x = "geo_accession", by.y = "gsm",
             suffixes = c(".v1", ".new"))
  cat("merged rows:", nrow(m), "\n")
  cat("subject match:", all(m$subject_id.v1 == m$subject_id.new), "\n")
  cat("day match:", all(m$day.v1 == m$day.new), "\n")
  cat("age_group match:", all(m$age_group.v1 == m$age_group.new), "\n")
  cat("sex match:", all(tolower(m$sex.v1) == tolower(m$sex.new)), "\n")
  cat("cohort match:", all(m$cohort.v1 == m$cohort.new), "\n")
  cat("age match:", all(m$age.v1 == m$age.new), "\n")
}
cat("\nDONE\n")
