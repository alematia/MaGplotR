#!/usr/bin/env R

## Welcome to ScrispR

options(warn=-1)

## Load libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyverse))


# option parsing #

parser <- OptionParser()

parser <- add_option(parser, 
                     opt_str = c("-f", "--format"), 
                     type = "character",
                     dest = 'format',
                     help="Input file format (mageck -mg- or screenbeam -sb-). mageck by default."
)

parser <- add_option(parser, 
                     opt_str = c("-i", "--inputdirectory"), 
                     type = "character",
                     dest = 'input.directory',
                     help="Input directory (required). Path to directory where test files are saved."
)

parser <- add_option(parser, 
                     opt_str = c("-c", "--control"), 
                     type = "character",
                     default = NA,
                     dest = 'control.file',
                     help="A file path to the control experiment."
)


parser <- add_option(parser, 
                     opt_str = c("-o", "--outputdirectory"), 
                     type = "character",
                     dest = 'output.directory',
                     default = getwd(),
                     help="Output directory (required). Path to directory where images will be saved."
)

parser <- add_option(parser, 
                     opt_str = c("-p", "--plotformat"), 
                     type = "character",
                     default = 'png',
                     dest = 'plot.format',
                     help="A string: png (default), pdf, ps, jpeg, tiff, bmp are supported."
)




opt = parse_args(parser)



# Set input path
first_dir <- getwd()  # Get wd.
input_dir_path = opt$input.directory  # Input test files directory.
setwd(input_dir_path)  # Set wd to read files.


# Set output directory
if(is.null(opt$output.directory)){ 
  output.directory = file.path(dirname(input_dir_path))
} else {
  output.directory = file.path(opt$output.directory)
}


# Select plot format
if (opt$plot.format %in% c('png', 'pdf', 'ps', 'jpeg', 'tiff', 'bmp')){
  plot.format = opt$plot.format
} else {
  flog.error("The plot format (-p) must be 'png', 'pdf', 'ps', 'jpeg', 'tiff', or 'bmp'.")
  stop()
}


# Set control file path. If no control file is loaded, pass.
if(is.null(opt$control.file)){ 
} else {
  control_file_path <- file.path(opt$control.file)
}

# Read MaGeCK or ScreenBEAM files.
if (opt$format == "mg") {
  MaGeCK_files = list.files(pattern="*.txt")
  input_files_txt = lapply(MaGeCK_files, read.delim)
  print("Analyzing MaGeCK data...")
  if(is.na(opt$control.file)){ 
    control_file_mg <- NULL
    print("No control file detected.")
  } else {
    control_file_mg <- read.delim(file.path(control_file_path))
    print("Control file detected")
  }
} else if (opt$format == "sb") {
  ScreenBEAM_files = list.files(pattern="*.csv")
  input_files_csv <- lapply(ScreenBEAM_files, read.csv)
  print("Analyzing ScreenBEAM data...")
} else {print("Error: wrong input format")}


setwd(first_dir)


## MaGeCK:
#
# gene id and its lfc is substracted from each df.
id_lfc_maker <- function(input_files_txt){
  num <- 1
  sub_input_files_txt <- c()
  for (i in input_files_txt){
    sub_input_files_txt[[num]] <- data.frame(i$id, i$pos.lfc)
    colnames(sub_input_files_txt[[num]]) <- c(paste0("id", num), "LFC")
    num <- num + 1
  }
  return(sub_input_files_txt)
}


# Melts 
melter <- function(sub_input_files_txt){
  melted_sub_input_files <- c()
  num <- 1
  for (i in sub_input_files_txt){
    melted_sub_input_files[[num]] <- reshape2::melt(i, id.vars = 2)
    num <- num + 1
  }
  return(melted_sub_input_files)
}


# Creates a vector of basic x axis labels (T1, T2...)
#xlabs <- function(){
#  vector <- c()
#  num <- 1
#  for (i in 1:9999){
#    vector <- c(vector, paste0("T", num))
#    num <- num + 1
#  }
#  return(vector)
#}

## Creates sub dataset with gene id and rank
id_rank_maker_pos <- function(input_files_txt){
  num <- 1
  sub_mg_rank_files <- c()
  for (i in input_files_txt){
    sub_mg_rank_files[[num]] <- data.frame(i$id, i$pos.rank)
    colnames(sub_mg_rank_files[[num]]) <- c("id", paste0("Rank", num))
    num <- num + 1
  }
  return(sub_mg_rank_files)
}


id_rank_maker_neg <- function(input_files_txt){
  num <- 1
  sub_mg_rank_files <- c()
  for (i in input_files_txt){
    sub_mg_rank_files[[num]] <- data.frame(i$id, i$neg.rank)
    colnames(sub_mg_rank_files[[num]]) <- c("id", paste0("Rank", num))
    num <- num + 1
  }
  return(sub_mg_rank_files)
}


main_mg <- function(x = input_files_txt, y = control_file_mg){
  sub_input_files_txt <- id_lfc_maker(input_files_txt)
  melted_sub_input_files <- melter(sub_input_files_txt)
  bound <- bind_rows(as.vector(melted_sub_input_files)) # Binds all id_lfc dfs
  #experiments_labs <- xlabs()
  
  ## 1. Boxplot with control if supplied
  if (is.null(control_file_mg)){
    boxplot <- ggplot(bound, aes(variable, y= LFC))+
      geom_boxplot(aes(colour=variable), outlier.shape = NA)+
      geom_jitter(position = position_jitter(width = 0.25) , size=0.005, alpha = 0.05, aes(colour=variable))+
      xlab("Test experiments")+
      ylab("Log2 Fold Change (LFC)")+
      theme(legend.position = "None", axis.text.x = element_text(size = 10))#+
    #scale_x_discrete(labels=experiments_labs, )
  } else {
    sub_control_mg <- data.frame(Control = control_file_mg$id, LFC = control_file_mg$pos.lfc)
    melted_control_mg <- reshape2::melt(sub_control_mg, id.vars=2)
    all_data_boxplot_mg <- bind_rows(bound, melted_control_mg)
    boxplot <- ggplot(all_data_boxplot_mg, aes(variable, y= LFC))+
      geom_boxplot(aes(colour=variable), outlier.shape = NA)+
      geom_jitter(position = position_jitter(width = 0.25) , size=0.005, alpha = 0.05, aes(colour=variable))+
      xlab("Test experiments")+
      ylab("Log2 Fold Change (LFC)")+
      theme(legend.position = "None", axis.text.x = element_text(size = 10))
  }
  suppressMessages(ggsave(path = output.directory, filename = paste0("boxplot.", plot.format), plot = boxplot, device = plot.format))
  print("Boxplot saved in output directory.")
  
  ## 2. Heatmaps + record csv files (pos and neg).
  sub_mg_rank_files_pos <- id_rank_maker_pos(input_files_txt)
  ## Merge all data frames according to the gene id.
  merged_mg_pos <- sub_mg_rank_files_pos %>% reduce(inner_join, by = "id")
  ## Add the mean of all exp ranks to each gene as a new column.
  rank_data_mg_pos <- merged_mg_pos %>% 
    select(matches("Rank"))
  merged_mg_pos$RankMeans <- rowMeans(rank_data_mg_pos)
  ## Generate csv file with rank records.
  sorted_whole_mg_pos <- merged_mg_pos[order(merged_mg_pos$RankMeans),]
  write.csv(sorted_whole_mg_pos, paste0(output.directory, "/MaGeCK_ranked_genes_pos.csv"), row.names = FALSE)
  print("Positive selection ranked genes (.csv file) saved in output directory.")
  ## Head top 25
  top_25_mg_pos <- head(sorted_whole_mg_pos, 25)
  ## Melt df to tidy large format
  melted_mg_pos <- melt(top_25_mg_pos,id.vars = c("id", "RankMeans"))
  heatmap_mg_pos <- ggplot(melted_mg_pos, aes(x=variable, y=reorder(id, -RankMeans), fill=value))+
    geom_tile(colour="white", size=.2)+
    ggtitle("Gene position in rank")+
    theme(panel.background = element_blank(),
          legend.position = "left",
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())+
    scale_fill_distiller(palette = "RdYlGn", direction = 1, trans= "reverse", 
                         name="", labels = c("Highest\nposition\nin rank", "Lowest\nposition\nin rank"),
                         breaks=seq(min(melted_mg_pos$value), max(melted_mg_pos$value),
                                    by=(max(melted_mg_pos$value)-min(melted_mg_pos$value))))  # Legend scale from highest to lowest position in rank
  suppressMessages(ggsave(path = output.directory, filename = paste0("heatmap_mageck_pos.", plot.format), plot = heatmap_mg_pos, device = plot.format))
  print("Heatmap (positive selection) saved in output directory.")
  
  ## Same with negative ranks
  sub_mg_rank_files_neg <- id_rank_maker_neg(input_files_txt)
  merged_mg_neg <- sub_mg_rank_files_neg %>% reduce(inner_join, by = "id")
  rank_data_mg_neg <- merged_mg_neg %>% 
    select(matches("Rank"))
  merged_mg_neg$RankMeans <- rowMeans(rank_data_mg_neg)
  sorted_whole_mg_neg <- merged_mg_neg[order(merged_mg_neg$RankMeans),]
  write.csv(sorted_whole_mg_neg, paste0(output.directory, "/MaGeCK_ranked_genes_neg.csv"), row.names = FALSE)
  print("Negative selection ranked genes (.csv file) saved in output directory.")
  top_25_mg_neg <- head(sorted_whole_mg_neg, 25)
  melted_mg_neg <- melt(top_25_mg_neg,id.vars = c("id", "RankMeans"))
  heatmap_mg_neg <- ggplot(melted_mg_neg, aes(x=variable, y=reorder(id, -RankMeans), fill=value))+
    geom_tile(colour="white", size=.2)+
    ggtitle("Gene position in rank")+
    theme(panel.background = element_blank(),
          legend.position = "left",
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())+
    scale_fill_distiller(palette = "RdYlGn", direction = 1, trans= "reverse", 
                         name="", labels = c("Highest\nposition\nin rank", "Lowest\nposition\nin rank"),
                         breaks=seq(min(melted_mg_neg$value), max(melted_mg_neg$value),
                                    by=(max(melted_mg_neg$value)-min(melted_mg_neg$value))))  # Legend scale from highest to lowest position in rank
  suppressMessages(ggsave(path = output.directory, filename = paste0("heatmap_mageck_neg.", plot.format), plot = heatmap_mg_neg, device = plot.format))
  print("Heatmap (negative selection) saved in output directory.")
  
  ## Control plots for heatmaps if control file is given
  if(is.null(control_file_mg)){
    control_file_mg <- NULL
  } else {
    ## Control LFC plot for positive heatmap
    simple_top_pos <- data.frame(Control = top_25_mg_pos$id, rank = 1:length(top_25_mg_pos$id))
    control_merge_pos <- merge(x=simple_top_pos, y=sub_control_mg , by = "Control")
    control_merge_pos <- control_merge_pos[order(control_merge_pos$rank),]
    self_plot_pos <- ggplot(control_merge_pos, aes(x = LFC , y = reorder(Control, -rank)))+
      geom_point(size=1.5)+
      xlim(floor(min(control_merge_pos$LFC)), max(control_merge_pos$LFC))+
      theme(text = element_text(size=12), legend.position = "none", panel.grid.major.y = element_line(colour="black"), 
            panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
            axis.line.x = element_blank(), axis.line.y = element_blank(), axis.text.x= element_text(),
            axis.ticks.x = element_line(), axis.ticks.y = element_blank(), axis.title.y = element_blank())+
      geom_vline(xintercept=0, linetype="dashed", color = "red")
    suppressMessages(ggsave(path = output.directory, filename = paste0("self_enrichment_pos.", plot.format), plot = self_plot_pos, device = plot.format))
    print("Self enrichment plot (positive selection) saved in output directory.")
    
    
    ## Control LFC plot for negative heatmap
    simple_top_neg <- data.frame(Control = top_25_mg_neg$id, rank = 1:length(top_25_mg_neg$id))
    control_merge_neg <- merge(x=simple_top_neg, y=sub_control_mg , by = "Control")
    control_merge_neg <- control_merge_neg[order(control_merge_neg$rank),]
    self_plot_neg <- ggplot(control_merge_neg, aes(x = LFC , y = reorder(Control, -rank)))+
      geom_point(size=1.5)+
      xlim(floor(min(control_merge_neg$LFC)), max(control_merge_neg$LFC))+
      theme(text = element_text(size=12), legend.position = "none", panel.grid.major.y = element_line(colour="black"), 
            panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
            axis.line.x = element_blank(), axis.line.y = element_blank(), axis.text.x= element_text(),
            axis.ticks.x = element_line(), axis.ticks.y = element_blank(), axis.title.y = element_blank())+
      geom_vline(xintercept=0, linetype="dashed", color = "red")
    suppressMessages(ggsave(path = output.directory, filename = paste0("self_enrichment_neg.", plot.format), plot = self_plot_neg, device = plot.format))
    print("Self enrichment plot (negative selection) saved in output directory.")
  }
  
  
}

if (opt$format == "mg"){
  main_mg(input_files_txt, control_file_mg)
} #else if (opt$format == "sb"){
#main_sb()
#}
