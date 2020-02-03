local({
  lib <- file.path("lib", R.version$platform, substr(getRversion(), 1, 3))
  if (!file.exists(lib)) {
    dir.create(lib, recursive = TRUE)
  }
  .libPaths(c(lib, Sys.getenv("R_LIBS_USER")))
})
