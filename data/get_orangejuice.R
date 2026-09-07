# get_orangejuice.R -- fetch the orangeJuice data of the R package bayesm (Rossi; GPL >= 2)
# from CRAN and export its two data frames to CSV in the working directory.  Called by
# examples/00_get_data.do; run from the data/ folder.
repos <- "https://cloud.r-project.org"
if (!requireNamespace("bayesm", quietly = TRUE)) install.packages("bayesm", repos = repos)
library(bayesm)
cat("bayesm version:", as.character(packageVersion("bayesm")), "\n")
data(orangeJuice)
yx <- orangeJuice$yx
sd <- orangeJuice$storedemo
cat("yx: rows", nrow(yx), "cols", ncol(yx), "; stores", length(unique(yx$store)),
    "brands", length(unique(yx$brand)), "weeks", length(unique(yx$week)), "\n")
write.csv(yx, "orangeJuice_yx.csv", row.names = FALSE)
write.csv(sd, "orangeJuice_storedemo.csv", row.names = FALSE)
cat("written orangeJuice_yx.csv and orangeJuice_storedemo.csv\n")
