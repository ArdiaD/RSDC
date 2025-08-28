# Build the packaged dataset "greenbrown" from the Excel in inst/extdata

# Only needed for building the data:
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")

file <- system.file("extdata", "green-brown-ptf.xlsx", package = "RSDC")

greenbrown <- readxl::read_excel(file)

# Clean/types
greenbrown$DATE <- as.Date(greenbrown$DATE)

# Ensure it's a base data.frame
greenbrown <- as.data.frame(greenbrown)

# Order by DATE ascending
greenbrown <- greenbrown[order(greenbrown$DATE), ]

# Save as compressed .rda in data/
usethis::use_data(greenbrown, overwrite = TRUE, compress = "xz")
