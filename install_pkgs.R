user_lib <- Sys.getenv("R_LIBS_USER")
if (!dir.exists(user_lib)) dir.create(user_lib, recursive = TRUE)

local({r <- getOption('repos'); r['CRAN'] <- 'https://cloud.r-project.org'; options(repos=r)})

if (!requireNamespace('BiocManager', quietly = TRUE)) {
    install.packages('BiocManager', lib = user_lib)
}

library(BiocManager, lib.loc = user_lib)
BiocManager::install(c('ComplexHeatmap', 'circlize', 'tidyverse'), lib = user_lib, ask = FALSE)
