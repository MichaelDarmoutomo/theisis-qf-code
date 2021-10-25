load_data <- function(filename, skip=0) {
  library(readxl)
  print(getwd())
  
  dir = "../data/"
  file_dir = paste0(dir, file)
  cat('Trying to read file', file_dir, '\n')
  cat(getwd())
  
  # Check if .csv or .xlsx
  if (grepl('csv', filename, fixed = TRUE)) {
    o = read.csv(file_dir, sep=';', header = FALSE)
  } else if (grepl('xlsx', filename, fixed = TRUE)) {
    o = read_excel(file_dir, skip=skip)
  } else {
    stop('File extension not recognized.')
  }
  o
}

