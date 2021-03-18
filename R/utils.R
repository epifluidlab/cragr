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
#' @return Logical value. TRUE if all required binaries can be found in PATH.
check_binaries <- function(binaries = "tabix", verbose = FALSE) {
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

  all(check_results)
}
