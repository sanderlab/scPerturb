library(here)
data_dir = paste0(here::here(), "/scATAC/data/")

download.file(
    url = "https://github.com/GreenleafLab/SpearATAC_MS_2021/raw/main/AnalyzeSpearATAC/data/Vierstra-Human-Motifs.rds",
    destfile = file.path(data_dir, "Vierstra-Human-Motifs.rds")
)