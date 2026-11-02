.trace_input <- function(path) {
  if (is.character(path) && length(path) == 1L && file.exists(path)) {
    cat(normalizePath(path, winslash = "/"), "\n", file = Sys.getenv("P26_TRACE_FILE"), append = TRUE)
  }
}
if (nzchar(Sys.getenv("P26_TRACE_FILE"))) {
  invisible(trace("read.table", tracer = quote(.trace_input(file)), print = FALSE, where = asNamespace("utils")))
}
