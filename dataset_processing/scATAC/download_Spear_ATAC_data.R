#Download SpearATAC Test Data
library(here)
data_dir = paste0(here::here(), "/scATAC/data/Spear_ATAC/")

aws <- "https://jeffgranja.s3.amazonaws.com/SpearATAC_2020"



files = c()

data_file_types <- c(".fragments.tsv.gz", ".fragments.tsv.gz.tbi", ".sgRNA.rds", ".singlecell.csv")

replicates <- 1:4
for (replicate in replicates) {
    for (data_file_type in data_file_types) {
        filename = paste0("GM-LargeScreen-R", replicate, data_file_type, collapse = ", ")
        files = c(files, filename)
    }
}

replicates <- 1:6
for (replicate in replicates) {
    for (data_file_type in data_file_types) {
        filename = paste0("K562-LargeScreen-R", replicate, data_file_type, collapse = ", ")
        files = c(files, filename) 
    }
}


replicates <- 1:4
for (replicate in replicates) {
    for (data_file_type in data_file_types) {
        filename = paste0("MCF7-LargeScreen-R", replicate, data_file_type, collapse = ", ")
        files = c(files, filename) 
    }
}


dir.create(data_dir, showWarnings = FALSE)

for(i in seq_along(files)){
    if(!file.exists(file.path(data_dir, files[i]))){
        download.file(
            url = file.path(aws, files[i]),
            destfile = file.path(data_dir, files[i])
        )
    }
}