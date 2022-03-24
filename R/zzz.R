.onLoad = function (libname, pkgname) {
  datafile = system.file("extdata", "orbis_sim_sample.csv", package = "revenueCES")
  assign('datafile', datafile, envir = topenv())
}
