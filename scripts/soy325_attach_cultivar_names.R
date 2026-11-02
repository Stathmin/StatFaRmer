# Attach proper cultivar names from this project's handmade onto `cultivar`.
# Keeps numeric id in `cultivar_id`. No-op if handmade / column missing.
.sf_soy325_attach_cultivar_names <- function(d, project_id = "project_2025-07-14_b51-91_soy") {
  hm_path <- file.path("data", project_id, paste0(project_id, "_handmade.csv"))
  if (!file.exists(hm_path) || is.null(d) || !is.data.frame(d) || !"cultivar" %in% names(d)) {
    return(d)
  }
  hm <- tryCatch(
    utils::read.csv(hm_path, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8-BOM"),
    error = function(e) {
      utils::read.csv(hm_path, stringsAsFactors = FALSE, check.names = FALSE, encoding = "UTF-8")
    }
  )
  names(hm) <- sub("^\ufeff", "", names(hm))
  if (!all(c("cultivar", "cultivar_name") %in% names(hm))) return(d)
  map <- unique(hm[, c("cultivar", "cultivar_name"), drop = FALSE])
  map$cultivar <- as.character(map$cultivar)
  map$cultivar_name <- as.character(map$cultivar_name)
  map <- map[!is.na(map$cultivar) & nzchar(map$cultivar) &
               !is.na(map$cultivar_name) & nzchar(map$cultivar_name), , drop = FALSE]
  if (nrow(map) < 1L) return(d)
  lookup <- setNames(map$cultivar_name, map$cultivar)
  d$cultivar_id <- as.character(d$cultivar)
  mapped <- unname(lookup[d$cultivar_id])
  ok <- !is.na(mapped) & nzchar(mapped)
  if (!any(ok)) return(d)
  d$cultivar <- as.character(d$cultivar)
  d$cultivar[ok] <- mapped[ok]
  d
}
