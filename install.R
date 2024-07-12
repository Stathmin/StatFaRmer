# renv -----
renv::activate()
renv::hydrate(prompt = FALSE)
if (.Platform$OS.type == "windows") {
  install.packages('flextable', type="binary")
}
