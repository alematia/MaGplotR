setwd("/home/pox/ScrispR/output_dir")
dir.create("volcano_plots")
dir.create("control_exp_plots")


library(ggplot2)

example_file <- read.delim(file=file.choose())
example_file2 <- read.delim(file=file.choose())
View(example_file)
View(example_file2)

#volcanoplot
volcanoplot <- ggplot(example_file, aes(x = pos.lfc , y = -log2(pos.p.value)))+
  geom_point(size=0.5)

volcanoplot

#control_exp_plot
merged_df <- merge(example_file, example_file2, by = "id")
View(merged_df)
control_exp_plot <- ggplot(merged_df, aes(x = pos.lfc.x, y = pos.lfc.y))+
  geom_point(size=.75, color = "grey20", alpha = 0.8)+
  geom_vline(xintercept=0, linetype="dashed", color = "red")+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  geom_rug(colour="indianred2", alpha=.2, size=0.4)

control_exp_plot
