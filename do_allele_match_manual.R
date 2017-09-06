# Load packages
library(allelematch)

# Import data (can be csv, even SQL statement)
# Implementing reading from a csv with decimal delimiter and tab as separator. 
xy <- read.table("./data/76311e7f71e2c1de0d8d2cd1c5917a64.csv", header = TRUE, sep = ";")
xy <- read.table("./data/testni_podatki.txt", header = TRUE, sep = ";")
xy <- read.table("./data/sample_data_with_id.csv", header = TRUE, sep = ";")
# podvojena imena vzorcev, vrže error
xy <- read.table("./data/518ac7b0a6ab2ead5323b6b65c16b480.csv", header = TRUE, sep = ";")
head(xy)
xy

xy <- reshape(xy, timevar = "marker", idvar = c("sample", "id"), direction = "wide")
# xy <- reshape(xy, timevar = "marker", idvar = "sample", direction = "wide")
names(xy) <- gsub("measured\\.", replacement = "", x = names(xy))

xy[1:5, 1:5]

image(t(as.matrix(xy[, c(-1, -2)])))

# imajo preveč NAjev?
unname(apply(xy, MARGIN = 1, FUN = function(x) {
  m <- table(is.na(x))
  m["TRUE"]/sum(m)
}))

# Convert the imported dataset into a proper structure for matching
d.xy <- amDataset(xy, missingCode = NA, indexColumn = "sample", metaDataColumn = "id")
# d.xy <- amDataset(xy, missingCode = NA, indexColumn = "sample")

# Do profiling. Figure is saved to working directory.
png(filename = "tmp.png",  width = 500, height = 500) # output figure size, adapt if needed
prof <- tryCatch(amUniqueProfile(amDatasetFocal = d.xy, verbose = TRUE),
                 error = function(e) e,
                 warning = function(w) w)

if (any(class(prof) %in% c("error", "warning"))) {
  plot.new()
  plot.window(xlim = c(-5, 5), ylim = c(-5, 5))
  text(x = 0, y = 1, labels = "Profiling failed with error message:")
  text(x = 0, y = 0, labels = prof$message)
}
dev.off()

cat("Done performing profiling, figure saved.", file = opt$verbose)

message("Figure saved.")
mp <- c("alleleMismatch", "matchThreshold", "cutHeight")
opt[mp] <- lapply(opt[mp], FUN = function(x) {
  if (is.na(x)) {
    x <- NULL
  }
  x
})

if (all(is.null(opt$alleleMismatch), is.null(opt$matchThreshold), is.null(opt$cutHeight))) {
  stop("Please specify at least one parameter (alleleMismatch, matchThreshold or cutHeight), see --help.")
}

cat("Performing matching.", file = opt$verbose)

message("Performing matching...")
result <- amUnique(amDatasetFocal = d.xy, 
                   alleleMismatch = opt$alleleMismatch, 
                   matchThreshold = opt$matchThreshold,
                   cutHeight = opt$cutHeight
)

cat("Matching done.", file = opt$verbose)

# Produce result as specified in output
# Catch csv, if not, assume html is requested
message("Writing output...")

if (grepl("\\.csv$", opt$output)) {
  suppressMessages(require(tidyr))
  
  cat("Writing result to CSV.", file = opt$verbose)
  
  amCSV.amUnique(x = result, csvFile = opt$output)
  
  # Reread the data and reflow it into a long format, subset only 
  # relevant columns and resave it
  rein <- read.table(opt$output, header = TRUE, sep = ",")
  rein <- gather(rein, key = markerName, value = value, -uniqueGroup, -rowType, -uniqueIndex,
                 -matchIndex, -nUniqueGroup, -alleleMismatch, -matchThreshold, -cutHeight,
                 -Psib, -score)[, c("uniqueGroup", "rowType", "uniqueIndex", "matchIndex", "Psib",
                                    "markerName", "value")]
  write.table(rein, file = opt$output, quote = FALSE, sep = ",",
              col.names = TRUE, row.names = FALSE)
  
  cat("Done writing CSV file.")
} else {
  cat("Writing data to HTML.", file = opt$verbose)
  amHTML.amUnique(x = result, htmlFile = opt$output)
  cat("Done writing data to HTML.", file = opt$verbose)
}
}
