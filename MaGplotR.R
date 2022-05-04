#!/usr/bin/env R

## Welcome to MaGplotR

options(warn=-1)

## Load libraries
print("Loading libraries...")
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ReactomePA))


# option parsing #

parser <- OptionParser()


parser <- add_option(parser, 
                     opt_str = c("-i", "--inputdirectory"), 
                     type = "character",
                     dest = 'input.directory',
                     help="Input directory (required). Path to directory where gene summary test files are saved."
)


parser <- add_option(parser, 
                     opt_str = c("-c", "--control"), 
                     type = "character",
                     default = NA,
                     dest = 'control.file',
                     help="A file path to the control experiment (optional)."
)


parser <- add_option(parser, 
                     opt_str = c("-t", "--top-cutoff"), 
                     type = "character",
                     dest = 'top.cutoff',
                     help="A number for heatmap visualization (top hits). 25 is set as default."
)


parser <- add_option(parser, 
                     opt_str = c("-o", "--outputdirectory"), 
                     type = "character",
                     dest = 'output.directory',
                     default = getwd(),
                     help="Output directory (optional). Path to directory where images will be saved. Working directory is default."
)


parser <- add_option(parser, 
                     opt_str = c("-g", "--sgrna-inputdirectory"), 
                     type = "character",
                     dest = 'sgrna.input.directory',
                     help="sgRNA input directory (optional). Path to directory where sgRNA summary test files are saved."
)


parser <- add_option(parser, 
                     opt_str = c("-p", "--plotformat"), 
                     type = "character",
                     default = 'png',
                     dest = 'plot.format',
                     help="A string: png (default), pdf, ps, jpeg, tiff, bmp are supported."
)




opt = parse_args(parser)


first_dir <- getwd()  # Get wd.

# Set input directory (gene summary files)
input_dir_path = opt$input.directory



# Set output directory
if(is.null(opt$output.directory)){ 
  output.directory = file.path(dirname(input_dir_path))
} else {
  output.directory = file.path(opt$output.directory)
}


# Set control file path. If no control file is loaded, pass.
if(!is.null(opt$control.file)){control_file_path <- file.path(opt$control.file)}


# Set number of hits shown in heatmap. 25 is default.
if(!is.null(opt$top.cutoff)){
  top_cutoff <- opt$top.cutoff
} else {
  top_cutoff <- 25
}


# Set sgRNA input directory
if(!is.null(opt$sgrna.input.directory)){sgrna_input_dir_path = opt$sgrna.input.directory}


# Select plot format
if (opt$plot.format %in% c('png', 'pdf', 'ps', 'jpeg', 'tiff', 'bmp')){
  plot.format = opt$plot.format
} else {
  flog.error("The plot format (-p) must be 'png', 'pdf', 'ps', 'jpeg', 'tiff', or 'bmp'.")
  stop()
}





# Read MaGeCK files
setwd(input_dir_path)  # Set wd to read gene summary files.
MaGeCK_files = list.files(pattern="*.txt")
input_files_txt = lapply(MaGeCK_files, read.delim)
print("MaGeCK gene summary data detected.")
if(is.na(opt$control.file)){ 
  control_file <- NULL
  print("No control file detected.")
} else {
  control_file <- read.delim(file.path(control_file_path))
  print("Control file detected")
}


# Read sgRNA summary files.
if(!is.null(opt$sgrna.input.directory)){
  print("MaGeCK sgRNA data detected.")
  setwd(sgrna_input_dir_path)  # Set wd to read sgRNA summary files.
  sgMaGeCK_files = list.files(pattern="*.txt")
  sginput_files_txt = lapply(sgMaGeCK_files, read.delim)
  #if(is.na(opt$control.file)){ 
  #  control_file <- NULL
  #  print("No control file detected.")
  #}
}



setwd(first_dir)


## MaGeCK gene summary analysis:
#
# gene id and LFC is obtained from each df
genes_shortener <- function(input_files){
  num <- 1
  sub_input_files <- c()
  for (i in input_files_txt){
    sub_input_files[[num]] <- data.frame(i$id, i$pos.lfc)
    colnames(sub_input_files[[num]]) <- c(gsub(".gene_summary.txt", "", MaGeCK_files[num]), "LFC")
    num <- num + 1
  }
  return(sub_input_files)
}


# Melts gene summary sub dfs
melter <- function(sub_input_files){
  melted_sub_input_files <- c()
  num <- 1
  for (i in sub_input_files){
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
    colnames(sub_mg_rank_files[[num]]) <- c("id", gsub(".gene_summary.txt", " ", MaGeCK_files[num]))
    num <- num + 1
  }
  return(sub_mg_rank_files)
}


id_rank_maker_neg <- function(input_files_txt){
  num <- 1
  sub_mg_rank_files <- c()
  for (i in input_files_txt){
    sub_mg_rank_files[[num]] <- data.frame(i$id, i$neg.rank)
    colnames(sub_mg_rank_files[[num]]) <- c("id", gsub(".gene_summary.txt", " ", MaGeCK_files[num]))
    num <- num + 1
  }
  return(sub_mg_rank_files)
}


id_LFC_maker <- function(input_files_txt){
  num <- 1
  sub_mg_LFC_files <- c()
  for (i in input_files_txt){
    sub_mg_LFC_files[[num]] <- data.frame(i$id, i$pos.lfc)
    colnames(sub_mg_LFC_files[[num]]) <- c("id", gsub(".gene_summary.txt", " ", MaGeCK_files[num]))
    num <- num + 1
  }
  return(sub_mg_LFC_files)
}

# gene analysis 
gene_analysis <- function(x = input_files_txt, y = control_file){
  print("Analyzing MaGeCK gene summary data...")
  sub_input_files <- genes_shortener(input_files_txt)
  melted_sub_input_files <- melter(sub_input_files)
  bound <- bind_rows(as.vector(melted_sub_input_files)) # Binds all id_lfc dfs
  #experiments_labs <- xlabs()
  
  ## 1. Boxplot with control if supplied
  if (is.null(control_file)){
    boxplot <- ggplot(bound, aes(variable, y = LFC))+
      geom_boxplot(aes(colour=variable), outlier.shape = NA)+
      geom_jitter(position = position_jitter(width = 0.25) , size = 0.005, alpha = 0.05, aes(colour=variable))+
      xlab("Test experiments")+
      ylab("Genes Log2 Fold Change (LFC)")+
      theme(legend.position = "None", axis.text.x = element_text(size = 10, angle = 45, hjust = 1))#+
    #scale_x_discrete(labels=experiments_labs, )
  } else {
    sub_control_mg <- data.frame(id = control_file$id, LFC = control_file$pos.lfc)
    sub_control_mg2 <- data.frame(Control = control_file$id, LFC = control_file$pos.lfc)
    melted_control_mg <- reshape2::melt(sub_control_mg2, id.vars = 2)
    all_data_boxplot_mg <- bind_rows(bound, melted_control_mg)
    boxplot <- ggplot(all_data_boxplot_mg, aes(variable, y= LFC))+
      geom_boxplot(aes(colour=variable), outlier.shape = NA)+
      geom_jitter(position = position_jitter(width = 0.25) , size=0.005, alpha = 0.05, aes(colour=variable))+
      xlab("Test experiments")+
      ylab("Genes Log2 Fold Change (LFC)")+
      theme(legend.position = "None", axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
  }
  suppressMessages(ggsave(path = output.directory, filename = paste0("genes_boxplot.", plot.format), plot = boxplot, device = plot.format))
  print("Genes boxplot saved in output directory.")
  
  ## 2. Heatmaps + record csv files (pos and neg).
  sub_mg_rank_files_pos <- id_rank_maker_pos(input_files_txt)
  ## Merge all data frames according to the gene id.
  merged_mg_pos <- sub_mg_rank_files_pos %>% reduce(inner_join, by = "id")
  ## Add the mean of all exp ranks to each gene as a new column.
  rank_data_mg_pos <- merged_mg_pos %>% 
    select(matches(" "))
  merged_mg_pos$RankMeans <- rowMeans(rank_data_mg_pos)
  ## Generate csv file with rank records.
  sorted_whole_mg_pos <- merged_mg_pos[order(merged_mg_pos$RankMeans),]
  write.csv(sorted_whole_mg_pos, paste0(output.directory, "/MaGeCK_ranked_genes_pos.csv"), row.names = FALSE)
  print("Positive selection ranked genes (.csv file) saved in output directory.")
  ## Head top 25
  top_mg_pos <- head(sorted_whole_mg_pos, as.numeric(top_cutoff))
  ## Melt df to tidy large format
  melted_mg_pos <- melt(top_mg_pos,id.vars = c("id", "RankMeans"))
  heatmap_mg_pos <- ggplot(melted_mg_pos, aes(x=variable, y=reorder(id, -RankMeans), fill=value))+
    geom_tile(colour="white", size=.2)+
    ggtitle("Gene position in positive rank")+
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
    select(matches(" "))
  merged_mg_neg$RankMeans <- rowMeans(rank_data_mg_neg)
  sorted_whole_mg_neg <- merged_mg_neg[order(merged_mg_neg$RankMeans),]
  write.csv(sorted_whole_mg_neg, paste0(output.directory, "/MaGeCK_ranked_genes_neg.csv"), row.names = FALSE)
  print("Negative selection ranked genes (.csv file) saved in output directory.")
  top_mg_neg <- head(sorted_whole_mg_neg, as.numeric(top_cutoff))
  melted_mg_neg <- melt(top_mg_neg,id.vars = c("id", "RankMeans"))
  heatmap_mg_neg <- ggplot(melted_mg_neg, aes(x=variable, y=reorder(id, -RankMeans), fill=value))+
    geom_tile(colour="white", size=.2)+
    ggtitle("Gene position in negative rank")+
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
  if(is.null(control_file)){
    control_file <- NULL
  } else {
    # Obtain LFC mean for every gene in all experiments:
    sub_mg_LFC_files <- id_LFC_maker(input_files_txt)
    ## Merge all data frames according to the gene id:
    merged_mg_LFC <- sub_mg_LFC_files %>% reduce(inner_join, by = "id")
    ## Add the mean of all exp LFCs to each gene as a new column.
    LFC_data_mg <- merged_mg_LFC %>% 
      select(matches(" "))
    merged_mg_LFC$LFCMeans <- rowMeans(LFC_data_mg)
    all_LFC <- data.frame(id = merged_mg_LFC$id, expLFCMeans = merged_mg_LFC$LFCMeans)
    
    ## Control LFC plot for positive heatmap
    simple_top_pos <- data.frame(id = top_mg_pos$id, rank = 1:length(top_mg_pos$id))
    control_merge_pos <- merge(x=simple_top_pos, y=sub_control_mg , by = "id")
    control_merge_pos <- merge(x=control_merge_pos, y=all_LFC , by = "id")
    control_merge_pos <- control_merge_pos[order(control_merge_pos$rank),]
    trial_plot_pos <- control_merge_pos %>% pivot_longer(cols = c(LFC, expLFCMeans), names_to = "lfcs")  #pivot plot
    self_plot_pos <- ggplot(trial_plot_pos, aes(x = value, y = reorder(id, -rank), col = lfcs, group = lfcs))+
      geom_point(size=1.75)+
      scale_color_manual(values = c("coral2", "lightseagreen"))+
      xlim(floor(min(trial_plot_pos$value)), max(trial_plot_pos$value))+
      theme(text = element_text(size=12), legend.position = "none", panel.grid.major.y = element_line(colour="black"), 
            panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
            axis.line.x = element_blank(), axis.line.y = element_blank(), axis.text.x= element_text(),
            axis.ticks.x = element_line(), axis.ticks.y = element_blank(), axis.title.y = element_blank())+
      xlab("LFC")+
      geom_vline(xintercept=0, linetype="dashed", color = "red")
    suppressMessages(ggsave(path = output.directory, filename = paste0("self_enrichment_pos.", plot.format), plot = self_plot_pos, device = plot.format))
    print("Self enrichment plot (positive selection) saved in output directory.")
    
    ## Control LFC plot for negative heatmap
    simple_top_neg <- data.frame(id = top_mg_neg$id, rank = 1:length(top_mg_neg$id))
    control_merge_neg <- merge(x=simple_top_neg, y=sub_control_mg , by = "id")
    control_merge_neg <- merge(x=control_merge_neg, y=all_LFC , by = "id")
    control_merge_neg <- control_merge_neg[order(control_merge_neg$rank),]
    trial_plot_neg <- control_merge_neg %>% pivot_longer(cols = c(LFC, expLFCMeans), names_to = "lfcs")  #pivot plot
    self_plot_neg <- ggplot(trial_plot_neg, aes(x = value, y = reorder(id, -rank), col = lfcs, group = lfcs))+
      geom_point(size=1.75)+
      scale_color_manual(values = c("coral2", "lightseagreen"))+
      xlim(floor(min(trial_plot_neg$value)), max(trial_plot_neg$value))+
      theme(text = element_text(size=12), legend.position = "none", panel.grid.major.y = element_line(colour="black"), 
            panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
            axis.line.x = element_blank(), axis.line.y = element_blank(), axis.text.x= element_text(),
            axis.ticks.x = element_line(), axis.ticks.y = element_blank(), axis.title.y = element_blank())+
      xlab("LFC")+
      geom_vline(xintercept=0, linetype="dashed", color = "red")
    suppressMessages(ggsave(path = output.directory, filename = paste0("self_enrichment_neg.", plot.format), plot = self_plot_neg, device = plot.format))
    print("Self enrichment plot (negative selection) saved in output directory.")
  }
  
  ## REACTOME PATHWAY ANALYSIS
  suppressPackageStartupMessages(library(org.Hs.eg.db)) # Library is loaded here to avoid overlapping function errors
  
  perc_num_pos <- round(nrow(sorted_whole_mg_pos)*0.01)  # Obtain the top 1 %.
  pre_go_pos <- head(sorted_whole_mg_pos, perc_num_pos)
  go_genes_pos <- pre_go_pos$id
  suppressMessages(go_pos <- select(org.Hs.eg.db,
              keys = go_genes_pos,
              columns=c("ENTREZID", "SYMBOL"),
              keytype="SYMBOL"))
  go_pos <- na.omit(go_pos)
  pathways_pos <- enrichPathway(as.data.frame(go_pos)$ENTREZID)
  pathways_pos@result$Description <- gsub("Homo sapiens\r: ", "", as.character(pathways_pos@result$Description))
  pathways_pos@result <- pathways_pos@result[order(-pathways_pos@result$Count),]
  write.csv(pathways_pos@result, paste0(output.directory, "/ReactomePA_pos.csv"), row.names = FALSE)
  pathways_pos_plot <- ggplot(head(pathways_pos@result,10), aes(x = Count, y = reorder(Description, Count)))+
    geom_point(aes(size=Count, colour=p.adjust))+
    scale_size(range = c(4,12))+
    theme(panel.background = element_blank(), axis.line.y = element_line(color="black"),
          axis.line.x = element_line(color="black"))+
    scale_color_gradient(low = "springgreen4", high = "chocolate1")+
    xlim(min(head(pathways_pos@result$Count, 10)), max(head(pathways_pos@result$Count, 10)))
  suppressMessages(ggsave(width = 10, path = output.directory, filename = paste0("ReactomePA_pos.", plot.format), plot = pathways_pos_plot, device = plot.format))
  
  perc_num_neg <- round(nrow(sorted_whole_mg_neg)*0.01)  # Obtain the top 1 %.
  pre_go_neg <- head(sorted_whole_mg_neg, perc_num_neg)
  go_genes_neg <- pre_go_neg$id
  suppressMessages(go_neg <- select(org.Hs.eg.db,
                                    keys = go_genes_neg,
                                    columns=c("ENTREZID", "SYMBOL"),
                                    keytype="SYMBOL"))
  go_neg <- na.omit(go_neg)
  pathways_neg <- enrichPathway(as.data.frame(go_neg)$ENTREZID)
  pathways_neg@result$Description <- gsub("Homo sapiens\r: ", "", as.character(pathways_neg@result$Description))
  pathways_neg@result <- pathways_neg@result[order(-pathways_neg@result$Count),]
  write.csv(pathways_neg@result, paste0(output.directory, "/ReactomePA_neg.csv"), row.names = FALSE)
  pathways_neg_plot <- ggplot(head(pathways_neg@result,10), aes(x = Count, y = reorder(Description, Count)))+
    geom_point(aes(size=Count, colour=p.adjust))+
    scale_size(range = c(4,12))+
    theme(panel.background = element_blank(), axis.line.y = element_line(color="black"),
          axis.line.x = element_line(color="black"))+
    scale_color_gradient(low = "springgreen4", high = "chocolate1")+
    xlim(min(head(pathways_neg@result$Count, 10)), max(head(pathways_neg@result$Count, 10)))
  suppressMessages(ggsave(width = 10, path = output.directory, filename = paste0("ReactomePA_neg.", plot.format), plot = pathways_neg_plot, device = plot.format))
  
  #suppressMessages(ggsave(path = output.directory, filename = paste0("go_pos.", plot.format), plot = dotplot(pathways_pos), device = plot.format))
  detach("package:org.Hs.eg.db", unload=TRUE)
  print("Reactome Pathway Analysis completed")
  
  # GO ENRICHMENT ANALYSIS 
  suppressMessages(library(clusterProfiler))
  go_bp_pos <- enrichGO(as.data.frame(go_pos)$ENTREZID, 'org.Hs.eg.db', ont="BP", pvalueCutoff=0.01)
  go_bp_pos_plot <- ggplot(head(go_bp_pos@result,10), aes(x = Count, y = reorder(Description, Count)))+
    geom_point(aes(size=Count, colour=p.adjust))+
    scale_size(range = c(4,12)) +
    theme(panel.background = element_blank(), axis.line.y = element_line(color="black"),
          axis.line.x = element_line(color="black"))+
    scale_color_gradient(low = "springgreen4", high = "chocolate1")+
    xlim(min(go_bp_pos@result$Count), max(go_bp_pos@result$Count))
  suppressMessages(ggsave(width = 10, path = output.directory, filename = paste0("GO_BP_pos.", plot.format), plot = go_bp_pos_plot, device = plot.format))
}


# sgrna analysis
sgrna_analysis <- function(sginput_files_txt){
  
  # Sub df with columns sgRNA, Gene and LFC.
  sgrna_shortener <- function(sginput_files){
    num <- 1
    sgsub_input_files <- c()
    for (i in sginput_files){
      sgsub_input_files[[num]] <- data.frame(i$sgrna,i$Gene, i$LFC)
      colnames(sgsub_input_files[[num]]) <- c("sgRNA", "Gene", "LFC")
      num <- num + 1
    }
    return(sgsub_input_files)
  }
  
  
  ## Appends a column with the experiment number.
  col_appender <- function(sgsub_input_files){
    num <- 1
    for (i in sgsub_input_files){
      nrows <- nrow(i)
      sgsub_input_files[[num]]$exp_number <- rep(c(num), times=nrows)
      num <- num + 1   
    }
    return(sgsub_input_files)
  }
  
  
  # Execute functions
  sgsub_input_files <- sgrna_shortener(sginput_files_txt)
  sgsub_input_files_coln <- col_appender(sgsub_input_files)
  sgbound <- bind_rows(as.vector(sgsub_input_files_coln)) # Binds all sgRNA sub dfs
  
  # Create boxplot for sgRNAs' LFCs
  sgboxplot <- ggplot(sgbound, aes(x = factor(exp_number), y = LFC))+
    geom_boxplot(aes(colour=factor(exp_number)), outlier.shape = NA)+
    geom_jitter(position = position_jitter(width = 0.25) , size=0.005, alpha = 0.05, aes(colour=factor(exp_number)))+
    xlab("Test experiments")+
    ylab("sgRNA Log2 Fold Change (LFC)")+
    theme(legend.position = "None", axis.text.x = element_text(size = 10))
  suppressMessages(ggsave(path = output.directory, filename = paste0("sgrnas_boxplot.", plot.format), plot = sgboxplot, device = plot.format))
  print("sgRNAs boxplot saved in output directory.")
}


gene_analysis(input_files_txt, control_file)
if(!is.null(opt$sgrna.input.directory)){sgrna_analysis(sginput_files_txt)}
