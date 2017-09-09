xy <- read.table("./sandbox/gis_data_sim.csv", sep = ";", header = TRUE)
xy$id <- as.numeric(xy$sample)
write.table(xy, file = "./sandbox/test_genotipi_id.csv", sep = ";", col.names = TRUE,
            row.names = FALSE, quote = FALSE)

# Rscript do_allele_match.R --profile --input=gis_demo_gene.csv --fig=mojaslika.png
# Rscript do_allele_match.R --input=test_genotipi_id.csv --output=testout.csv --alleleMismatch=3
# Rscript do_allele_match.R --input=./data/sample_id.csv --output=sample_id_result.html --alleleMismatch=1
# Rscript do_allele_match.R --input=./data/sample_id.csv --output=sample_id_result.csv --alleleMismatch=1

# create sample data, matching and calibration table
library(readxl)
library(data.table)

xy <- as.data.table(read_excel("./data/DAB_Genetika_DemoData.RealnoMbase2017.xlsx"))

xy <- melt(xy, id.vars = c("Sample", "Mrkr", "Rel", "QualityIndex", "Animal"),
     measure.vars = c("Al1", "Al2"),
     value.name = "measured")[order(Sample, Mrkr), ]
xy[, marker := sprintf("%s_%s", Mrkr, gsub("^Al(\\d$)", "\\1", xy$variable))]
xy[1:10, ]

out <- xy[, .(Sample, marker, measured)]
names(out) <- c("sample", "marker", "measured")
fwrite(out, file = "./data/sample_real_dinalpbear_data.csv", sep = ";")

# Matching
library(allelematch)

mx <- xy[, .(Sample, marker, Animal, measured)] # subset only relevant variables
mx <- dcast(mx, Sample + Animal ~ marker, value.var = "measured") # cast from long to wide

d.mx <- amDataset(multilocusDataset = mx, missingCode = NA, indexColumn = "Sample",
                  metaDataColumn = "Animal")

prf <- amUniqueProfile(amDatasetFocal = d.mx, verbose = TRUE)
rs <- amUnique(amDatasetFocal = d.mx, alleleMismatch = 4)
amHTML.amUnique(x = rs, htmlFile = "result_dinalpbear.html")
amCSV.amUnique(x = rs, csvFile = "result_dinalpbear.csv")

# Calibration table
setkey(xy, Mrkr, measured)
out <- xy[, lapply(.SD, function(x) head(x, 1)), by = key(xy)]
out <- out[, .(Mrkr, measured, Sample)]
out <- out[!is.na(measured), ]
out[, sequence := ""]
out[, comment := ""]
fwrite(out, file = "./data/calibration_table_dab_testdata.csv", sep = "\t")
