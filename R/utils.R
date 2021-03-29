#' Check if required binaries exist in PATH
#'
#' Some of cragr's functionalities require certain third-party binaries to exist
#' in PATH. For example, selectively loading cfDNA fragments from a certain
#' region under the hood invokes tabix to complete the task. This function checks
#' whether certain binaries exist.
#' @param binaries Character vector indicating what binaries to check. Default
#'   is tabix.
#' @param verbose Logical value indicating whether to output detailed messages.
#'   Default is FALSE.
#' @param stop_on_error Stop the script completely if the result is `FALSE`.
#' @return Logical value. TRUE if all required binaries can be found in PATH.
#' @export
check_binaries <- function(binaries = "tabix", verbose = FALSE, stop_on_fail = FALSE) {
  check_results <- binaries %>% purrr::map_lgl(function(binary) {
    # check if binary is in path
    if ("" == Sys.which(binary)) {
      if (verbose) {
        message(str_interp("Checking ${binary}: FAIL"))
      }
      return(FALSE)
    }
    else {
      if (verbose) {
        message(str_interp("Checking ${binary}: PASS"))
      }
      return(TRUE)
    }
  })

  check_results <- all(check_results)
  if (stop_on_fail)
    assertthat::assert_that(check_results, msg = paste0("Failed in locating ", paste(binaries, collapse = ", ")))
  else
    check_results
}
