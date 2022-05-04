#!/usr/bin/env R
## Installation script for MaGplotR

if (!("optparse" %in% installed.packages())) { 
  install.packages("optparse");
}

if (!("ggplot2" %in% installed.packages())) { 
  install.packages("ggplot2");
}

if (!("tidyr" %in% installed.packages())) { 
  install.packages("tidyr");
}

if (!("reshape2" %in% installed.packages())) { 
  install.packages("reshape2");
}

if (!("dplyr" %in% installed.packages())) { 
  install.packages("dplyr");
}

if (!("stringr" %in% installed.packages())) { 
  install.packages("stringr");
}

if (!("tidyverse" %in% installed.packages())) { 
  install.packages("tidyverse");
}

if (!("tidyverse" %in% installed.packages())) { 
  install.packages("tidyverse");
}

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")
BiocManager::install("ReactomePA")
BiocManager::install("clusterProfiler")
