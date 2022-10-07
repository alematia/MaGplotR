#!/usr/bin/env R
## Installation script for MaGplotR

if (!("optparse" %in% installed.packages())) { 
  install.packages("optparse");
}

if (!("tidyverse" %in% installed.packages())) { 
  install.packages("tidyverse");
}

if (!("reshape2" %in% installed.packages())) { 
  install.packages("reshape2");
}

if (!("BiocManager" %in% installed.packages())) { 
  install.packages("BiocManager");
}

if (!("org.Hs.eg.db" %in% installed.packages())) { 
  BiocManager::install("org.Hs.eg.db");
}

if (!("ReactomePA" %in% installed.packages())) { 
  BiocManager::install("ReactomePA");
}

if (!("clusterProfiler" %in% installed.packages())) { 
  BiocManager::install("clusterProfiler");
}

if (!("matrixStats" %in% installed.packages())) { 
  BiocManager::install("matrixStats");
}
