library(stringr)
library(dplyr)

# expects hypsography to be in proportions of lake area
setwd("./data")
dir_name <- "Fourth_set_AVP"
Hypsos <- list.files(dir_name, pattern = "csv", full.names = T, recursive = F)

QCreport <- data.frame(File = NA, Flag = NA)


for (i in seq_along(Hypsos)) {
  QC_file <- data.frame(File = NA, Flag = NA)
  
  hypso_data <- read.csv(Hypsos[i])
  
  # Check file name
  filename <- str_split(Hypsos[i], "/", simplify = T)[,2] 
  
  # Check for missing underscore
  if (!grepl("_", filename, fixed = TRUE)) {
    QC_temp <- data.frame(File = Hypsos[i], Flag = "File name missing underscore")
    QC_file <- rbind(QC_file, QC_temp)
  }
  #Only relevant for MN hypsos
  # Check DOW is 8 digits
  if (grepl("_", filename, fixed = TRUE)) {
    DOW <- str_split(filename, "_", simplify = T)
    DOW <- DOW[,ncol(DOW)] %>% str_replace(".csv", "")
    if (!nchar(DOW) == 8) {
      QC_temp <- data.frame(File = Hypsos[i], Flag = "DOW incorrect")
      QC_file <- rbind(QC_file, QC_temp)
    }
  }
  # Check column names
  if (!"depth_feet" %in% colnames(hypso_data)) {
    QC_temp <- data.frame(File = Hypsos[i], Flag = "Depth column name")
    QC_file <- rbind(QC_file, QC_temp)
  }
  if (!"proportion_area" %in% colnames(hypso_data)) {
    QC_temp <- data.frame(File = Hypsos[i], Flag = "Area column name")
    QC_file <- rbind(QC_file, QC_temp)
  }
  
  # Check if there are any commas in the data (instead of decimal points)
  if (grepl(",", paste0(hypso_data[,1], hypso_data[,2], collapse = ""), fixed = TRUE)) {
    QC_temp <- data.frame(File = Hypsos[i], Flag = "Comma in data")
    QC_file <- rbind(QC_file, QC_temp)
  }
  
  # Check for extra columns
  if (ncol(hypso_data) > 2) {
    QC_temp <- data.frame(File = Hypsos[i], Flag = "Extra columns")
    QC_file <- rbind(QC_file, QC_temp)
  }  
  
  hypso_data2 <- hypso_data[,1:2]
  hypso_data2 <- hypso_data2[!is.na(hypso_data2[,1]),]
  
  # Check if columns are numeric
  if (!(is.numeric(hypso_data2[,1]) & is.numeric(hypso_data2[,2]))) {
    QC_temp <- data.frame(File = Hypsos[i], Flag = "Columns not numeric")
    QC_file <- rbind(QC_file, QC_temp)
  }
  
  if (sum(is.na(hypso_data2[,2])) > 0) {
    QC_temp <- data.frame(File = Hypsos[i], Flag = "Missing Areas")
    QC_file <- rbind(QC_file, QC_temp)
  }
  
  # Check for failures in monotonicity
  if (is.numeric(hypso_data2[,1])) {
    if (any(diff(hypso_data2[,1]) < 0)) {
      QC_temp <- data.frame(File = Hypsos[i], Flag = "Depth not increasing")
      QC_file <- rbind(QC_file, QC_temp)
    }
  }
  if (is.numeric(hypso_data2[,2]) & sum(is.na(hypso_data2[,2])) == 0) {
    if (any(diff(hypso_data2[,2]) > 0)) {
      QC_temp <- data.frame(File = Hypsos[i], Flag = "Area not decreasing")
      QC_file <- rbind(QC_file, QC_temp)
    }
  }
  
  # Check first row
  if (is.numeric(hypso_data2[,1]) & is.numeric(hypso_data2[,2])) {
    if (!(hypso_data[1,2] == 1 & hypso_data[1,1] == 0)) {
      QC_temp <- data.frame(File = Hypsos[i], Flag = "First row incorrect")
      QC_file <- rbind(QC_file, QC_temp)
    }
  }
  
  # If there are flags, add to file
  if (!is.na(QC_file$File[2])) QCreport <- rbind(QCreport, QC_file)
  
}


QCreport <- QCreport[!is.na(QCreport[,1]),]
QCreport

write.csv(QCreport, paste0(dir_name,"QC_report.csv"), row.names = FALSE)

