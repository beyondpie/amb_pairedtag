Sys.setenv("_R_USE_PIPEBIND_" = TRUE)
devtools::document()
devtools::install(pkg = ".", reload = TRUE, quick = FALSE,
  build = TRUE, upgrade = "never")
