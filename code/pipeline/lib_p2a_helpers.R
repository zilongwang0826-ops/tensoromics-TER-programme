# p2a_helpers.R — Phase-4 P2A candidate helpers (audit-only, NOT part of the original pipeline).
# Used by the fixed scripts in audit/phase4-p2/code/ and by tests/run_helpers_tests.R.
# Purpose: idempotent, atomic table updates; chip-column write guard.

upsert_rows <- function(df, new_rows, keys) {
  stopifnot(is.data.frame(df), is.data.frame(new_rows))
  stopifnot(all(keys %in% names(df)), all(keys %in% names(new_rows)))
  key_df <- function(x) do.call(paste, c(lapply(keys, function(k) as.character(x[[k]])), sep = "\r"))
  keep <- !(key_df(df) %in% key_df(new_rows))
  out <- rbind(df[keep, , drop = FALSE], new_rows)
  rownames(out) <- NULL
  out
}

atomic_write_csv <- function(df, path) {
  dir <- dirname(path)
  tmp <- tempfile(tmpdir = dir, fileext = ".tmp.csv")
  data.table::fwrite(df, tmp)
  ok <- file.copy(tmp, path, overwrite = TRUE)
  unlink(tmp)
  if (!isTRUE(ok)) stop("atomic_write_csv: copy failed for ", path)
  invisible(path)
}

write_chip_idempotent <- function(samples, chip_vec, path) {
  stopifnot(length(chip_vec) == nrow(samples))
  if (file.exists(path)) {
    cur <- data.table::fread(path)
    if ("chip" %in% names(cur)) {
      if (identical(as.character(cur$chip), as.character(chip_vec))) {
        return(FALSE)  # unchanged: skip write (idempotent)
      }
      stop("write_chip_idempotent: existing chip column differs from derived values; ",
           "refusing to overwrite (inspect CEL-derived serials vs stored file)")
    }
  }
  samples$chip <- chip_vec
  atomic_write_csv(samples, path)
  TRUE
}
