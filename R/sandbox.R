# #!/usr/bin/env Rscript
#
# suppressMessages(library(readr))
# suppressMessages(library(dplyr))
# suppressMessages(library(purrr))
# suppressMessages(library(stringr))
#
# # while(nrow(df <- data.table::fread('file:///dev/stdin')) > 0) {
# #
# # }
#
# lineno <- 0
#
# # f <- file("stdin")#, open = "rt")
# # open(f)
# # while (nrow(df <- tryCatch(read.delim(f, nrows = 150), error = function(e) NULL)) %||% 0 > 0) {
# #   # print(df)
# #   # lineno <- nrow(df)
# # }
# #
#
# while (TRUE) {
#   # lines <- readLines(con = f, n = -1)
#   # if (length(lines) == 0) {
#   #   break
#   # }
#
#   # message(str_interp("Chunk: ${length(lines)}"))
#   # lineno <<- lineno + length(lines)
#
#   # df <- data.table::fread(text = lines)
#   #
#   df <- suppressWarnings(data.table::fread("file:///dev/stdin"))
#   # # # df <- tryCatch(read.delim("/dev/stdin", nrows = 1500, comment.char = "#"), error = function(e) NULL)
#   # #
#   if (is.null(df) || nrow(df) == 0) {
#     print(df)
#     break
#   }
#   #
#   # colnames(df)[2] <- "start"
#   #
#   start_pos <- df$start[1]
#   end_pos <- tail(df$start, n = 1)
#   message(str_interp("Chunk: ${nrow(df)}"))
#   message(str_interp("${start_pos} -> ${end_pos}"))
#   # # print(df)
#   # # tryCatch(print(df), error = function(e) NULL)
#   lineno <<- lineno + nrow(df)
# }
#
# # close(f)
#
# message(str_interp("Total #: ${lineno}"))
#
# # while(nrow(df <- data.table::fread("file:///dev/stdin")) > 0) {
# #
# # }
# #
# # f <- file("stdin")
# # open(f)
# # while (length(lines <- readLines(f, n = -1)) > 0) {
# #   # conn <- textConnection(lines)
# #   # df <- read_tsv(file = conn, col_names = FALSE, col_types = cols())
# #   # df <- read.table(file = conn, sep = "\t") %>% as_tibble()
# #   # write(nrow(df), stderr())
# #   # print(df)
# #   # write(line, stderr())
# #   # process line
# # }
