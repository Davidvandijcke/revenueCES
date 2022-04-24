# set the "prod" object class #
setClass("prod", representation(Model = "list", Data = "list", Estimates = "list"))
# end of object classs setting #

# Show method #
setMethod('show', signature(object = 'prod'), function(object){

  Pars = object@Estimates$pars
  nPars = names(Pars)
  st.err = object@Estimates$std.errors
  Method = object@Model$method

  cat("\n-            Revenue CES Production Function Estimation           -")
  cat(paste("\n                   Method:   ", Method, "             "))
  cat(paste("\n Observations                         :  ", N, sep = "") )
  cat(paste("\n                      ", paste("  ", nPars, collapse = "    ", " ", sep = "")))
  cat(paste("\n ", paste(" ", round(Pars, digits = 3), collapse = "   ", " ", sep = "")))
  cat(paste("\n                      ", paste("(", round(st.err, digits = 3), collapse = "   ", ")", sep = "")))
  cat("\n-------------------------------------------------------")
})
# end of show method

# Show method #
setMethod('summary', signature(object = 'prod'), function(object){

  Pars = object@Estimates$pars
  namePars = names(Pars)
  st.err = object@Estimates$std.errors
  Method = object@Model$method
  N = nrow(object@Data$Y)
  Time = object@Model$elapsed.time
  nCores = object@Model$nCores
  opt = object@Model$opt
  R = object@Model$boot.repetitions
  theta0 = object@Model$theta0

  cat("\n-              Revenue CES Production Function Estimation              -")
  cat(paste("\n                   Method :   ", Method, "             "))
  cat(paste("\n Observations                         :  ", N, sep = "") )
  cat(paste("\n                            ", paste("  ", namePars, collapse = "    ", " ", sep = "")))
  cat(paste("\n     : ", paste(" ", round(Pars, digits = 3), collapse = "   ", " ", sep = "")))
  cat(paste("\n                            ", paste("(", round(st.err, digits = 3), collapse = "   ", ")", sep = "")))
  cat("\n-------------------------------------------------------------")
  cat(paste("\nElapsed Time              :  ", round(as.double(Time, units = 'mins'), digits = 2), " mins", sep = "") )
  if (!is.null(nCores)){
    cat(paste("\n# Cores                 : ", nCores, sep = "") )
  }
  cat("\n-------------------------------------------------------------")
})
# end of show method
