# MaGplotR
#### A software for Multiple Genetic Screens data visualization

*MaGplotR* produces visualization plots for MaGeCK screen data (RRA analysis). <br/>
This software is designed to run from a linux command line, and it is written in R.<br/>
*MaGplotR* directly uses `gene_summary.txt` files generated by the **test** function from MaGeCK software. Output plots give information of quality control of the screen data and highligths the best hits of the multiple screen experiments.<br/>
Additionally, a control experiment `gene_summary.txt` file can be used to compare the different experiments to a control self-enrichment experiment (or any condition used as control). Also, `sgrna_summary.txt` files can be analyzed for QC as optional arguments.

[Citation](#citation)<br/>
[Motivation](#motivation)<br/>
[Dependencies](#dependencies)<br/>
[Usage](#usage)<br/>



## Citation
If you use this software please cite:<br/>
[PREPRINT] / [PUBLICATION]<br/>


## Motivation
This software was developed to simplify the analysis of MaGeCK screen data from multiple experiments through elegant visualization. Sometimes, several screen experiments are performed testing multiple conditions or with a number of replicates, and it is desired to perform comparisons between these experiments. We put special focus introducing the possibility to add a control experiment to observe the behaviour of the hits in this experiment in a visually quantitative way that may help discarding false positives (i.e. self enrichment of potential hits in control experiment). Two main values are extracted from MaGeCK test summary files; LFC and rank, presenting them in a way that simplifies the whole set of experiments analysis.
Heatmap representation gives a rapid view of how top hits behave in the set of experiements, making very easy to detect if a particular hit had a poorer score in a experiment. The number of hits shown in heatmap plots can be adjusted when running the software from the command line (see options below). Also, the intermediate files to create the heatmap plots are generated and saved in the output directory so the user can check the ranking score of every gene.


## Dependencies
You can either install these packages from R terminal or by running the installation script from Linux terminal.<br/>

- Install from R terminal:<br/>
```r
install.packages("optparse")
install.packages("ggplot2")
install.packages("tidyr")
install.packages("reshape2")
install.packages("dplyr")
install.packages("stringr")
install.packages("tidyverse")
```
- Install using installation script:<br/>
Download `installation_magplotr.R` and run the following command in Linux terminal:
```bash
Rscript installation_magplotr.R
```



## Usage
#### Recommendations:
Create a directory and copy there all the gene summary files you want to analyze. This will be the input directory (-i). These files names usually end with ".gene_summary.txt" as MaGeCK output. MaGplotR uses the experiment name from the filename (i.e. for "exp7.gene_summary.txt" file, "exp7" will be displayed in the plots).<br/>
Optionally, do the same thing with sgRNA summary files (-g). If you have a control gene summary file to analyze (-c), we recommend to put it in the parent directory for simplicity.<br/>
Next, open a terminal in linux and change directory where the MaGplotR.R is located. Then type the commands as in the examples. Make sure you give the complete path to the file/folder as argument (you can easily do this copying the file/folder and clicking paste in the terminal). You can also provide an output folder, where files will be saved, otherwise, input directory is default. You can also choose different output formats and select top cutoff.<br/>

Examples:
```bash
Rscript MaGplotR.R -i path_to_results_directory
```
```bash
Rscript MaGplotR.R -i path_to_results_directory -c path_to_control_file -t 50 -p png -o path_to_output_directory -g path_to_sgRNA_input_directory
```
### Options:
Mandatory arguments:<br/>
-i (input directory): **path** to the folder where gene summary files (test files) are located. All files in this folder will be taken as input.<br/>

Optional arguments:<br/>
-c: (control file): **path** to the control file (no control as default). <br/>
-o: (output directory): **path** to an existing folder where output files will be saved (input directory default).<br/>
-g: (sgRNA input directory): **path** to an existing folder where sgRNA summary files are saved.<br/>
-t: (top cutoff): number of hits to be shown in heatmaps. 25 is default.<br/>
-p: (plot format): just write one among these (png is default): png , pdf, ps, jpeg, tiff, bmp.<br/>


## Output plots and files
Boxplot<br/>
![Screenshot from 2022-04-26 11-57-20](https://user-images.githubusercontent.com/95416488/165274882-804becbc-c209-4eef-ae05-66c687627aaa.png)

Heatmap with control<br/>
Heatmap represents gene ranks in each experiment. Control plot shows the LFC of control (cyan) and the mean LFC of all experiments (red) for each gene.<br/>
![heatmap_and_control](https://user-images.githubusercontent.com/95416488/165736171-025fa334-7f20-4196-ae27-9ba433f86435.jpg)


