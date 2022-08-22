## Welcome to MaGplotR
  
options(warn=-1)

## Load libraries
suppressPackageStartupMessages(library(stringr))
print(str_glue(""))
print(str_glue("               ***Welcome to MaGplotR***"))
print(str_glue(""))
print(str_glue("                        v"))
print(str_glue("O       o O       o O   |   o O       o O       o O       o"))
print(str_glue("| O   o | | O   o | | O | o | | O   o | | O   o | | O   o |"))
print(str_glue("| | O | | | | O | | | | | | | | | O | | | | O | | | | O | |"))
print(str_glue("| o   O | | o   O | | o | O | | o   O | | o   O | | o   O |"))
print(str_glue("o       O o       O o   |   O o       O o       O o       O"))
print(str_glue("                        ^"))
print(str_glue(""))
print(str_glue("- Loading libraries..."))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ReactomePA))
suppressMessages(library(clusterProfiler))
suppressMessages(library(matrixStats))
#suppressPackageStartupMessages(library(dplyr))
#suppressPackageStartupMessages(library(tidyr))


#OPTPARSE
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
                     opt_str = c("-s", "--selection"), 
                     type = "character",
                     default = 'pos',
                     dest = 'selection',
                     help="A string: pos for positive selection, or neg for negative selection."
)


parser <- add_option(parser, 
                     opt_str = c("-r", "--sgrna-inputdirectory"), 
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


parser <- add_option(parser, 
                     opt_str = c("-g", "--ontology"), 
                     type = "character",
                     dest = 'gene.ontology',
                     help="A string: BP, MF, CC."
)


parser <- add_option(parser, 
                     opt_str = c("-b", "--blind"), 
                     type = "character",
                     dest = 'colour.blind',
                     help="A string: 'y' or 'n'. 'n' as default."
)


opt = parse_args(parser)


first_dir <- getwd()

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


# Select pos or neg selection
selection = opt$selection


# Select Gene Ontology option
if(!is.null(opt$gene.ontology)){gene_ontology <- opt$gene.ontology}


if(!is.null(opt$colour.blind)){
  col_blind <- opt$colour.blind
} else {
  col_blind <- "n"
}


# Read MaGeCK files
setwd(input_dir_path)  # Set wd to read gene summary files.
MaGeCK_files <- list.files(pattern="*.txt")
input_files_txt = lapply(MaGeCK_files, read.delim)
print(str_glue("- MaGeCK gene summary data detected."))
setwd(first_dir)
if(is.na(opt$control.file)){ 
  control_file <- NULL
  print(str_glue("- No control file detected."))
} else {
  control_file <- read.delim(file.path(control_file_path))
  print(str_glue("- Control file detected."))
}




# Read sgRNA summary files.
if(!is.null(opt$sgrna.input.directory)){
  print(str_glue("- MaGeCK sgRNA data detected."))
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
  print(str_glue("- Analyzing MaGeCK gene summary data..."))
  sub_input_files <- genes_shortener(input_files_txt)
  melted_sub_input_files <- melter(sub_input_files)
  bound <- bind_rows(as.vector(melted_sub_input_files)) # Binds all id_lfc dfs
  
  ## 1. Boxplot with control if supplied
  if (is.null(control_file)){
    boxplot <- ggplot(bound, aes(variable, y = LFC))+
      geom_boxplot(aes(colour=variable),
                   outlier.shape = 1, outlier.size = .75, outlier.stroke = 0.75)+
      geom_jitter(position = position_jitter(width = 0.25) , size = 0.005, alpha = 0.05, aes(colour=variable))+
      xlab("Test experiments")+
      ylab("Genes Log2 Fold Change (LFC)")+
      theme(legend.position = "None", axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
            panel.background = element_blank(), axis.line = element_line(),
            panel.grid.minor.y = element_line(colour = "grey"))
  } else {
    sub_control_mg <- data.frame(id = control_file$id, LFC = control_file$pos.lfc)
    sub_control_mg2 <- data.frame(Control = control_file$id, LFC = control_file$pos.lfc)
    melted_control_mg <- reshape2::melt(sub_control_mg2, id.vars = 2)
    all_data_boxplot_mg <- bind_rows(bound, melted_control_mg)
    boxplot <- ggplot(all_data_boxplot_mg, aes(variable, y= LFC))+
      geom_boxplot(aes(colour=variable),
                   outlier.shape = 1, outlier.size = .75, outlier.stroke = 0.75)+
      geom_jitter(position = position_jitter(width = 0.25) , size = 0.005, alpha = 0.05, aes(colour=variable))+
      xlab("Test experiments")+
      ylab("Genes Log2 Fold Change (LFC)")+
      theme(legend.position = "None", axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
            panel.background = element_blank(), axis.line = element_line(),
            panel.grid.minor.y = element_line(colour = "grey"))
  }
  suppressMessages(ggsave(path = output.directory, filename = paste0("genes_boxplot.", plot.format),
                          plot = boxplot, device = plot.format))
  print(str_glue("- Genes boxplot saved in output directory."))
  
  ## 2. Heatmaps + record csv files (pos and neg).
  ## These 2 lines already existed in the code for the control df purpose.
  sub_mg_LFC_files <- id_LFC_maker(input_files_txt)
  merged_mg_LFC <- sub_mg_LFC_files %>% reduce(inner_join, by = "id")  # Merge all dfs by id
  #######################################################
  if (selection == "pos"){
    ## Preparation for pos
    sub_mg_rank_files_pos <- id_rank_maker_pos(input_files_txt)
    merged_mg_pos <- sub_mg_rank_files_pos %>% reduce(inner_join, by = "id")  # Merge all dfs by id
    ## Add the mean of all exp ranks to each gene as a new column.
    rank_data_mg_pos <- select(merged_mg_pos, !matches("id"))
    merged_mg_pos_x <- merged_mg_pos  # Replicate df to add more cols
    merged_mg_pos_x$RankMeans <- rowMeans(rank_data_mg_pos)  # Add RankMeans
    merged_mg_pos_x$RankSD <- rowSds(as.matrix(rank_data_mg_pos))  # Add Variance of the rank scores
    sorted_whole_mg_pos <- merged_mg_pos_x[order(merged_mg_pos_x$RankMeans),]  # Reorder df by RankMeans
    top_mg_pos <- head(sorted_whole_mg_pos, as.numeric(top_cutoff))  # Head top 25
    ## Preparation for LFC heatmap:
    ordered_vector_pos <- top_mg_pos$id  # Vector with the ordered list of genes (by RankMeans)
    ordered_LFC_heatmap_df_pos <- merged_mg_LFC[match(ordered_vector_pos, merged_mg_LFC$id),]  # Arranges LFCs df by vector
    ordered_LFC_heatmap_df_pos$RankMeans <- top_mg_pos$RankMeans  # Adds RankMeans col
    melted_LFC_pos <- melt(ordered_LFC_heatmap_df_pos,id.vars = c("id", "RankMeans"))  # Preparation for heatmap
    melted_mg_pos <- melt(top_mg_pos,id.vars = c("id", "RankMeans", "RankSD"))  # Preparation for heatmap
    ## Set limits to LFCs. All above 10 and below -5 is considered as max and min respectively.
    melted_LFC_pos$value[melted_LFC_pos$value > 10] <- 10  # All LFCs above 10 are now 10.
    melted_LFC_pos$value[melted_LFC_pos$value < -5] <- -5  # All LFCs below -5 are now -5.
    
    if (col_blind == "y"){
      heatmap_mg_pos <- ggplot(melted_LFC_pos, aes(x=variable, y=reorder(id, -RankMeans), fill=value))+
        geom_tile(colour="black", size=.2)+
        ggtitle("Gene LFC")+
        geom_text(size=1.5, aes(label = round(melted_mg_pos$value, 1)))+
        theme(panel.background = element_blank(),
              legend.position = "left",
              axis.text.x = element_text(angle = 45, hjust = 1),
              axis.title.x = element_blank(),
              axis.title.y = element_blank())+
        scale_fill_distiller(palette = "RdYlBu", direction = 1, name="LFC",
                             limits = c(-5, 10), labels=c("< -5", "0", "5", "> 10"),
                             values = c(0, 0.25, 0.4,1))
    } else {
      heatmap_mg_pos <- ggplot(melted_LFC_pos, aes(x=variable, y=reorder(id, -RankMeans), fill=value))+
        geom_tile(colour="black", size=.2)+
        ggtitle("Gene LFC")+
        geom_text(size=1.5, aes(label = round(melted_mg_pos$value, 1)))+
        theme(panel.background = element_blank(),
              legend.position = "left",
              axis.text.x = element_text(angle = 45, hjust = 1),
              axis.title.x = element_blank(),
              axis.title.y = element_blank())+
        scale_fill_distiller(palette = "RdYlGn", direction = 1, name="LFC",
                             limits = c(-5, 10), labels=c("< -5", "0", "5", "> 10"),
                             values = c(0, 0.25, 0.45,1))
    }
    #scale_fill_distiller(palette = "YlGn", direction = 1, name="LFC",
    #                    limits = c(-5, 10), labels=c("< -5", "0", "5", "> 10"))
    #scale_fill_distiller(palette = 4, direction = 1, name="LFC",
    #                     limits = c(min(melted_LFC$value), max(melted_LFC$value)))
    #scale_fill_distiller(palette = 5, limits = c(min(melted_LFC$value)-2, max(melted_LFC$value)+2), direction = 1)
  } else {
    ## Preparation for neg
    sub_mg_rank_files_neg <- id_rank_maker_neg(input_files_txt)
    merged_mg_neg <- sub_mg_rank_files_neg %>% reduce(inner_join, by = "id")  # Merge all dfs by id
    ## Add the mean of all exp ranks to each gene as a new column.
    rank_data_mg_neg <- select(merged_mg_neg, !matches("id"))
    merged_mg_neg_x <- merged_mg_neg  # Replicate df to add more cols
    merged_mg_neg_x$RankMeans <- rowMeans(rank_data_mg_neg)  # Add RankMeans
    merged_mg_neg_x$RankSD <- rowSds(as.matrix(rank_data_mg_neg))  # Add Variance of the rank scores
    sorted_whole_mg_neg <- merged_mg_neg_x[order(merged_mg_neg_x$RankMeans),]  # Reorder df by RankMeans
    top_mg_neg <- head(sorted_whole_mg_neg, as.numeric(top_cutoff))  # Head top 25
    ## Preparation for LFC heatmap:
    ordered_vector_neg <- top_mg_neg$id  # Vector with the ordered list of genes (by RankMeans)
    ordered_LFC_heatmap_df_neg <- merged_mg_LFC[match(ordered_vector_neg, merged_mg_LFC$id),]  # Arranges LFCs df by vector
    ordered_LFC_heatmap_df_neg$RankMeans <- top_mg_neg$RankMeans  # Adds RankMeans col
    melted_LFC_neg <- melt(ordered_LFC_heatmap_df_neg,id.vars = c("id", "RankMeans"))  # Preparation for heatmap
    melted_mg_neg <- melt(top_mg_neg,id.vars = c("id", "RankMeans", "RankSD"))  # Preparation for heatmap
    ## Set limits to LFCs. All above 10 and below -5 is considered as max and min respectively.
    melted_LFC_neg$value[melted_LFC_neg$value > 10] <- 10  # All LFCs above 10 are now 10.
    melted_LFC_neg$value[melted_LFC_neg$value < -5] <- -5  # All LFCs below -5 are now -5.
    if (col_blind == "y"){
      heatmap_mg_neg <- ggplot(melted_LFC_neg, aes(x=variable, y=reorder(id, -RankMeans), fill=value))+
        geom_tile(colour="black", size=.2)+
        ggtitle("Gene LFC")+
        geom_text(size=1.5, aes(label = round(melted_mg_neg$value, 1)))+
        theme(panel.background = element_blank(),
              legend.position = "left",
              axis.text.x = element_text(angle = 45, hjust = 1),
              axis.title.x = element_blank(),
              axis.title.y = element_blank())+
        scale_fill_distiller(palette = "RdYlBu", direction = 1, name="LFC",
                             limits = c(-5, 10), labels=c("< -5", "0", "5", "> 10"),
                             values = c(0, 0.25, 0.4,1))
    } else {
      heatmap_mg_neg <- ggplot(melted_LFC_neg, aes(x=variable, y=reorder(id, -RankMeans), fill=value))+
        geom_tile(colour="black", size=.2)+
        ggtitle("Gene LFC")+
        geom_text(size=1.5, aes(label = round(melted_mg_neg$value, 1)))+
        theme(panel.background = element_blank(),
              legend.position = "left",
              axis.text.x = element_text(angle = 45, hjust = 1),
              axis.title.x = element_blank(),
              axis.title.y = element_blank())+
        scale_fill_distiller(palette = "RdYlGn", direction = 1, name="LFC",
                             limits = c(-5, 10), labels=c("< -5", "0", "5", "> 10"),
                             values = c(0, 0.25, 0.45,1))
    }
  }
  
  # Save pos or neg
  if (selection == "pos"){
    write.csv(sorted_whole_mg_pos, paste0(output.directory, "/MaGeCK_ranked_genes_pos.csv"), row.names = FALSE)
    print(str_glue("- Positive selection ranked genes (.csv file) saved in output directory."))
    suppressMessages(ggsave(path = output.directory, filename = paste0("heatmap_mageck_pos.", plot.format),
                            plot = heatmap_mg_pos, device = plot.format))
    print(str_glue("- Heatmap (positive selection) saved in output directory."))
  } else {
    write.csv(sorted_whole_mg_neg, paste0(output.directory, "/MaGeCK_ranked_genes_neg.csv"), row.names = FALSE)
    print(str_glue("- Negative selection ranked genes (.csv file) saved in output directory."))
    suppressMessages(ggsave(path = output.directory, filename = paste0("heatmap_mageck_neg.", plot.format),
                            plot = heatmap_mg_neg, device = plot.format))
    print(str_glue("- Heatmap (negative selection) saved in output directory."))
  }
  
  ## Control plots for heatmaps if control file is given
  if(is.null(control_file)){
    control_file <- NULL
  } else {
    ## Add the mean of all exp LFCs to each gene as a new column.
    LFC_data_mg <- select(merged_mg_LFC, !matches("id"))
    merged_mg_LFC$LFCMeans <- rowMeans(LFC_data_mg)
    all_LFC <- data.frame(id = merged_mg_LFC$id, expLFCMeans = merged_mg_LFC$LFCMeans)
    
    if (selection == "pos"){
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
    } else {
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
    }
    
    # Save pos or neg
    if (selection == "pos"){
      suppressMessages(ggsave(path = output.directory, filename = paste0("self_enrichment_pos.", plot.format),
                              plot = self_plot_pos, device = plot.format))
      print(str_glue("- Self enrichment plot (positive selection) saved in output directory."))
    } else {
      suppressMessages(ggsave(path = output.directory, filename = paste0("self_enrichment_neg.", plot.format),
                              plot = self_plot_neg, device = plot.format))
      print(str_glue("- Self enrichment plot (negative selection) saved in output directory."))
    }
  }
  
  ## Gene Expression
  ## Tries to find expression file in current dir, home and Downloads. If found, proceed.
  if (file.exists("./MaGplotR/proteinatlas_gene_ex.tsv") == TRUE){
    ex_table <- read.delim("./MaGplotR/proteinatlas_gene_ex.tsv", check.names = F) ##### Give path to file
  } else if (file.exists("~/MaGplotR/proteinatlas_gene_ex.tsv") == TRUE) {
    ex_table <- read.delim("~/MaGplotR/proteinatlas_gene_ex.tsv", check.names = F)
  } else if (file.exists("~/Downloads/MaGplotR/proteinatlas_gene_ex.tsv") == TRUE) {
    ex_table <- read.delim("~/Downloads/MaGplotR/proteinatlas_gene_ex.tsv", check.names = F)
  } else {
    print(str_glue("- Gene expression file not found"))
  }
  ## If file is found, variable is created and then creates the plot.
  if (exists("ex_table")){
    if (selection == "pos") {
      ordered_ex_table <- na.omit(ex_table[match(ordered_vector_pos, ex_table$Gene),])
    } else {
      ordered_ex_table <- na.omit(ex_table[match(ordered_vector_neg, ex_table$Gene),])
    }
    short_ex_table <- ordered_ex_table[,c("Gene", "A-431", "A549", "Daudi", "HAP1", "HEK 293", "HeLa", "JURKAT", "K-562")]
    short_ex_table2 <- melt(short_ex_table,id.vars = "Gene")
    ordered_vector_ex <- short_ex_table$Gene
    
    plot_ex <- ggplot(short_ex_table2, aes(x=variable, y=factor(Gene, level = rev(ordered_vector_ex)),
                                           size=value, colour = value > 0))+
      geom_point(alpha=0.6, stroke = 1)+
      scale_size_continuous("Expression (nTPM)", limits = c(0, max(short_ex_table2$value)), 
                            range = c(1.5,5))+
      scale_colour_discrete("Expression > 0")+
      theme(panel.background = element_blank())+
      xlab("Cell line")+
      ylab("Gene")
    
    suppressMessages(ggsave(width = 10, path = output.directory, filename = paste0("Expression_plot.", plot.format),
                            plot = plot_ex, device = plot.format))
    print(str_glue("- Gene expression plot saved in output directory."))
  }
  
  
  ## REACTOME PATHWAY ANALYSIS
  suppressPackageStartupMessages(library(org.Hs.eg.db))  # Library is loaded here to avoid overlapping function errors
  
  if (selection == "pos"){
    # Reactome pos
    perc_num_pos <- round(nrow(sorted_whole_mg_pos)*0.01)  # Obtain the top 1 %.
    pre_go_pos <- head(sorted_whole_mg_pos, perc_num_pos)
    go_genes_pos <- pre_go_pos$id
    suppressMessages(go_pos <- select(org.Hs.eg.db,
                                      keys = go_genes_pos,
                                      columns=c("ENTREZID", "SYMBOL"),
                                      keytype="SYMBOL"))
    go_pos <- na.omit(go_pos)
    pathways_pos <- enrichPathway(as.data.frame(go_pos)$ENTREZID, organism = "human", readable = TRUE)
    pathways_pos@result$Description <- gsub("Homo sapiens\r: ", "", as.character(pathways_pos@result$Description))
    pathways_pos@result <- pathways_pos@result[order(-pathways_pos@result$Count),]
    pathways_pos_plot <- ggplot(head(pathways_pos@result,10), aes(x = 0, y = reorder(Description, Count)))+
      geom_point(stat = "identity", aes(size=Count, colour=p.adjust))+
      scale_size(range = c(4,13))+
      theme(panel.background = element_blank(), axis.line.y = element_line(color="black"),
            axis.line.x = element_line(color="black"), axis.text.x = element_blank(),
            axis.ticks.x = element_blank())+
      ylab("Pathway")+
      xlab("Count")+
      scale_color_distiller(palette = "RdPu", values = c(0, 0.25, 1))
    #OrRd #RdPu #YlGn
  } else {
    # Reactome neg
    perc_num_neg <- round(nrow(sorted_whole_mg_neg)*0.01)  # Obtain the top 1 %.
    pre_go_neg <- head(sorted_whole_mg_neg, perc_num_neg)
    go_genes_neg <- pre_go_neg$id
    suppressMessages(go_neg <- select(org.Hs.eg.db,
                                      keys = go_genes_neg,
                                      columns=c("ENTREZID", "SYMBOL"),
                                      keytype="SYMBOL"))
    go_neg <- na.omit(go_neg)
    pathways_neg <- enrichPathway(as.data.frame(go_neg)$ENTREZID, organism = "human", readable = TRUE)
    pathways_neg@result$Description <- gsub("Homo sapiens\r: ", "", as.character(pathways_neg@result$Description))
    pathways_neg@result <- pathways_neg@result[order(-pathways_neg@result$Count),]
    pathways_neg_plot <- ggplot(head(pathways_neg@result,10), aes(x = 0, y = reorder(Description, Count)))+
      geom_point(stat = "identity", aes(size=Count, colour=p.adjust))+
      scale_size(range = c(4,13))+
      theme(panel.background = element_blank(), axis.line.y = element_line(color="black"),
            axis.line.x = element_line(color="black"), axis.text.x = element_blank(),
            axis.ticks.x = element_blank())+
      ylab("Pathway")+
      xlab("Count")+
      scale_color_distiller(palette = "RdPu", values = c(0, 0.25, 1))
  }
  
  # Save pos or neg
  if(selection == "pos"){
    suppressMessages(ggsave(width = 10, path = output.directory, filename = paste0("ReactomePA_pos.", plot.format),
                            plot = pathways_pos_plot, device = plot.format))
    write.csv(pathways_pos@result, paste0(output.directory, "/ReactomePA_pos.csv"), row.names = FALSE)
    # Print most enriched pathway and genes.
    ids_vector <- str_split(pathways_pos@result[pathways_pos@result[1,1], "geneID"], "/")[[1]]
    names_str <- paste(ids_vector, collapse = ", ")
    print(str_glue(
      "- Most enriched pathway (pos) is: '", pathways_pos@result[pathways_pos@result[1,1], "Description"], "' with genes: ", names_str, ".")
    )
  } else {
    suppressMessages(ggsave(width = 10, path = output.directory, filename = paste0("ReactomePA_neg.", plot.format),
                            plot = pathways_neg_plot, device = plot.format))
    write.csv(pathways_neg@result, paste0(output.directory, "/ReactomePA_neg.csv"), row.names = FALSE)
    # Print most enriched pathway and genes.
    ids_vector <- str_split(pathways_neg@result[pathways_neg@result[1,1], "geneID"], "/")[[1]]
    names_str <- paste(ids_vector, collapse = ", ")
    print(str_glue(
      "- Most enriched pathway (neg) is: '", pathways_neg@result[pathways_neg@result[1,1], "Description"], "' with genes: ", names_str, ".")
    )
  }
  print(str_glue("- Reactome Pathway Analysis completed."))

  ## Clustering if > 2 exps.
  if (length(MaGeCK_files) > 2) {
    #### CLUSTERING by clusterProfiler. Top 1 %
    # Pos or neg
    if (selection == "pos"){
      num <- 1
      for (i in input_files_txt){
        input_files_txt[[num]] <- head(input_files_txt[[num]], perc_num_pos)
        num <- num + 1
      }
      # Select just gene ids
      num <- 1
      for (i in input_files_txt){
        input_files_txt[[num]] <- input_files_txt[[num]]$id
        num <- num + 1
      }
      # Annotation
      num <- 1
      sym_ids_vector <- c()
      library(org.Hs.eg.db)
      print(str_glue("- Clustering in progress..."))
      for (i in input_files_txt){
        suppressMessages(sym_ids_vector[[num]] <- select(org.Hs.eg.db,
                                                         keys = input_files_txt[[num]],
                                                         columns=c("ENTREZID", "SYMBOL"),
                                                         keytype="SYMBOL"))
        num <- num + 1
      }
      # Select ENTREZID
      num <- 1
      for (i in sym_ids_vector){
        sym_ids_vector[[num]] <- sym_ids_vector[[num]]$ENTREZID
        num <- num + 1
      }
      # Rename vector indexes
      num <- 1
      for (i in sym_ids_vector){
        names(sym_ids_vector)[[num]] <- num
        num <- num + 1
      }
      # Create the object
      cluster_pos <- compareCluster(geneClusters = sym_ids_vector, fun = "enrichPathway")
      if (!is.null(cluster_pos)){
        cluster_pos@compareClusterResult$Description <- gsub("Homo sapiens\r: ", "", as.character(cluster_pos@compareClusterResult$Description))
        cluster_pos@keytype <- "enrichPathway"
        cluster_pos@readable <- FALSE
        cluster_pos_readable <- setReadable(cluster_pos, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
        suppressMessages(ggsave(height=32, width=32, units=c("cm"), path = output.directory, filename = "cluster.pdf",
                                plot = suppressMessages(cnetplot(cluster_pos_readable)), device = "pdf"))
        print(str_glue("- Clustering completed."))
      } else {
        print(str_glue("- No clusters found."))
      }
      detach("package:org.Hs.eg.db", unload=TRUE)
    } else {
      # neg
      num <- 1
      for (i in input_files_txt){
        input_files_txt[[num]] <- input_files_txt[[num]][order(input_files_txt[[num]]$neg.rank),]
        input_files_txt[[num]] <- head(input_files_txt[[num]], perc_num_neg)
        num <- num + 1
      }
      # Select just gene ids
      num <- 1
      for (i in input_files_txt){
        input_files_txt[[num]] <- input_files_txt[[num]]$id
        num <- num + 1
      }
      # Annotation
      num <- 1
      sym_ids_vector <- c()
      library(org.Hs.eg.db)
      print(str_glue("- Clustering in progress..."))
      for (i in input_files_txt){
        suppressMessages(sym_ids_vector[[num]] <- select(org.Hs.eg.db,
                                                         keys = input_files_txt[[num]],
                                                         columns=c("ENTREZID", "SYMBOL"),
                                                         keytype="SYMBOL"))
        num <- num + 1
      }
      # Select ENTREZID
      num <- 1
      for (i in sym_ids_vector){
        sym_ids_vector[[num]] <- sym_ids_vector[[num]]$ENTREZID
        num <- num + 1
      }
      # Rename vector indexes
      num <- 1
      for (i in sym_ids_vector){
        names(sym_ids_vector)[[num]] <- num
        num <- num + 1
      }
      # Create the object
      cluster_neg <- compareCluster(geneClusters = sym_ids_vector, fun = "enrichPathway")
      if (!is.null(cluster_neg)){
        cluster_neg@compareClusterResult$Description <- gsub("Homo sapiens\r: ", "", as.character(cluster_neg@compareClusterResult$Description))
        cluster_neg@keytype <- "enrichPathway"
        cluster_neg@readable <- FALSE
        cluster_neg_readable <- setReadable(cluster_neg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
        suppressMessages(ggsave(height=32, width=32, units=c("cm"), path = output.directory,
                                filename = "cluster.pdf", plot = suppressMessages(cnetplot(cluster_neg_readable)), device = "pdf"))
        print(str_glue("- Clustering completed."))
      } else {
        print(str_glue("- No clusters found."))
      }
      detach("package:org.Hs.eg.db", unload=TRUE)
    }
    #### End of clustering
  }
  
  
  # GO ENRICHMENT ANALYSIS - only if user chooses
  if(!is.null(opt$gene.ontology)){
    if(selection == "pos"){
      suppressMessages(go_bp_pos <- enrichGO(as.data.frame(go_pos)$ENTREZID, 'org.Hs.eg.db', ont = gene_ontology, pvalueCutoff=0.01))
      go_bp_pos_plot <- ggplot(head(go_bp_pos@result,10), aes(x = Count, y = reorder(Description, Count)))+
        geom_point(aes(size=Count, colour=p.adjust))+
        scale_size(range = c(4,12)) +
        theme(panel.background = element_blank(), axis.line.y = element_line(color="black"),
              axis.line.x = element_line(color="black"))+
        scale_color_gradient(low = "springgreen4", high = "chocolate1")+
        xlim(min(go_bp_pos@result$Count), max(go_bp_pos@result$Count))
      suppressMessages(ggsave(width = 10, path = output.directory, filename = paste0("GO_", gene_ontology, "_pos.", plot.format), plot = go_bp_pos_plot, device = plot.format))
      print(str_glue("- GO Analysis completed."))
    } else {
      suppressMessages(go_bp_neg <- enrichGO(as.data.frame(go_neg)$ENTREZID, 'org.Hs.eg.db', ont = gene_ontology, pvalueCutoff=0.01))
      go_bp_neg_plot <- ggplot(head(go_bp_neg@result,10), aes(x = Count, y = reorder(Description, Count)))+
        geom_point(aes(size=Count, colour=p.adjust))+
        scale_size(range = c(4,12)) +
        theme(panel.background = element_blank(), axis.line.y = element_line(color="black"),
              axis.line.x = element_line(color="black"))+
        scale_color_gradient(low = "springgreen4", high = "chocolate1")+
        xlim(min(go_bp_neg@result$Count), max(go_bp_neg@result$Count))
      suppressMessages(ggsave(width = 10, path = output.directory, filename = paste0("GO_", gene_ontology, "_neg.", plot.format), plot = go_bp_neg_plot, device = plot.format))
      print(str_glue("- GO Analysis completed."))
    }
  }
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
  sgbound <- bind_rows(as.vector(sgsub_input_files_coln))  # Binds all sgRNA sub dfs
  
  # Create boxplot for sgRNAs' LFCs
  sgboxplot <- ggplot(sgbound, aes(x = factor(exp_number), y = LFC))+
    geom_boxplot(aes(colour=factor(exp_number)), outlier.shape = NA)+
    geom_jitter(position = position_jitter(width = 0.25) , size=0.005, alpha = 0.05, aes(colour=factor(exp_number)))+
    xlab("Test experiments")+
    ylab("sgRNA Log2 Fold Change (LFC)")+
    theme(legend.position = "None", axis.text.x = element_text(size = 10))
  suppressMessages(ggsave(path = output.directory, filename = paste0("sgrnas_boxplot.", plot.format), plot = sgboxplot, device = plot.format))
  print(str_glue("- sgRNAs boxplot saved in output directory."))
}


gene_analysis(input_files_txt, control_file)
if(!is.null(opt$sgrna.input.directory)){sgrna_analysis(sginput_files_txt)}
