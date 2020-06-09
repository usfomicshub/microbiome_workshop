# 6/3/2020: source-script for setting up environment to run the USF Omics Hub Microbiome Data-Analysis Workshop pipeline(s).
# Author: J. Oberstaller

# install necessary general packages/dependencies if they aren't already installed (combining AS, JG, JO, and TEK package-requirements)
packages <- c("BiocManager",
              "knitr",
              "caTools",
              "ggpubr",
              "dplyr",
              "RColorBrewer",
              "reshape2",
              "ggplot2",
              "hrbrthemes",
              "gcookbook",
              "RCurl",
              "fs",
              "readr",
              "GUniFrac",
              "simEd",
              "tictoc",
              "vegan"
             )
if (length(setdiff(packages,
                   rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages,
                           rownames(installed.packages())),
                   repos='https://cloud.r-project.org') 
}

# install our primary analysis-packages from bioconductor if they aren't already installed
bioconductor.packages <- c("dada2",
                           "microbiome",
                           "phyloseq",
                           "ALDEx2",
                           "SIAMCAT",
                           "metagenomeSeq")
if (length(setdiff(bioconductor.packages, rownames(installed.packages()))) > 0){
  BiocManager::install(c("dada2",
                         "microbiome",
                         "phyloseq",
                         "ALDEx2",
                         "SIAMCAT",
                         "metagenomeSeq"),
                       ask = FALSE,
                       quiet = TRUE,
                       verbose = FALSE)
}

# clean up leftover variables from above
rm("packages", "bioconductor.packages")

# add these load-package parts to the markdown for the appropriate day of the course instead of this source-document
# load required packages (AS)
library(dada2)
library(RCurl)
library(fs)
library(dplyr)
library(GUniFrac)
library(simEd)
library(tictoc)

# load required packages (JG)
library(microbiome)
library(phyloseq)
library(vegan)
library(metagenomeSeq)
library(knitr)
library(ggpubr)
library(dplyr)
library(RColorBrewer)
library(reshape2)
library(ggplot2)
library(hrbrthemes)
library(gcookbook)

# load required packages (TEK)
library(ALDEx2)
library(readr)
library(dplyr)
library(SIAMCAT)
