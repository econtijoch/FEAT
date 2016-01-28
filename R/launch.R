#' @export
launchFEAT <- function() {
  appDir <- system.file("shiny-FEAT", package = "FEAT")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `FEAT`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}
