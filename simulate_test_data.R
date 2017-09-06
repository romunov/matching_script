library(tidyr)

x <- read.table("test_genotipi.csv", header = TRUE, sep = ";")
x <- x[, c("sample", "marker", "Al1", "Al2")]

x <- gather(x, key = allel, value = measured, -sample, -marker)

x$marker <- paste(x$marker, "_", gsub(pattern = "^Al", replacement = "", x = x$allel), sep = "")
x$allel <- NULL

write.table(x, file = "testni_podatki.txt", sep = ";", quote = FALSE, row.names = FALSE,
            col.names = TRUE)
