# Load required packages and set working directory.
rm(list=ls())
setwd("~/github/ecolog-plates/")
require(xlsx)


# Read in and parse .xlsx files from plate reader.
filename <- "0723_3C_EcoLog_07112015_RNF.xlsx"
sheet <- read.xlsx(paste(getwd(),filename, sep="/"), sheetIndex=1, header=FALSE)
plate <- as.matrix(sheet[(19:26),(3:14)])
rownames(plate) <- c("A", "B", "C", "D", "E", "F", "G", "H")
colnames(plate) <- c("1", "2", "3", "4", "5", "6", 
                     "7", "8", "9", "10", "11", "12")

# Now, we have a matrix of the plate, 
# indexible by row name ("A":"H") and column number (1:12).
plate

# Save a .csv of the plate output for reading in later.
write.csv(plate, file=sub(".xlsx", replacement=".csv", filename))

