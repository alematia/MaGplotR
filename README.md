# MaGplotR
![LICENSE](https://img.shields.io/badge/license-MIT-green)
![R-VERSION](https://img.shields.io/badge/R-%204.2.0-blue)
![MaGplotR-VERSION](https://img.shields.io/badge/release-v0.3.0-orange)
#### A software for Multiple Genetic Screens data visualization

![Screenshot from 2022-05-13 13-02-56](https://user-images.githubusercontent.com/95416488/168270389-73f1f6e9-dee3-468a-ae2c-611c599a8aa9.png)


>*MaGplotR* produces visualization plots for MaGeCK screen data (RRA analysis). <br/>
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
Heatmap representation gives a rapid view of how top hits behave in the set of experiements, making very easy to detect if a particular hit had a poorer score in a experiment. The number of hits shown in heatmap plots can be adjusted when running the software from the command line (see options below). Also, the intermediate files to create the heatmap plots are generated and saved in the output directory so the user can check the ranking score of every gene.<br/>
As extra features, the top 1 % of positive and negative hits are analyzed with [Reactome Pathway Analysis](https://bioconductor.org/packages/release/bioc/html/ReactomePA.html) and [clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html) generating plots for the most enriched pathways and clusters among the data sets. Optionally, GO analysis is available (see options below). These analysis are run using [Genome wide annotation for Human](https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html).


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
install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("ReactomePA")
BiocManager::install("clusterProfiler")
```
- Install using installation script:<br/>
Download `installation_magplotr.R` and run the following command in Linux terminal:
```bash
Rscript installation_magplotr.R
```



## Usage
#### Recommendations:
Create a directory and copy there all the gene summary files you want to analyze. This will be the input directory (-i). These files names usually end with ".gene_summary.txt" as MaGeCK output. MaGplotR uses the experiment name from the filename (i.e. for "exp7.gene_summary.txt" file, "exp7" will be displayed in the plots).<br/>
Optionally, do the same thing with sgRNA summary files (-r). If you have a control gene summary file to analyze (-c), we recommend to put it in the parent directory for simplicity.<br/>
Next, open a terminal in linux and change directory where the MaGplotR.R is located. Then type the commands as in the examples. Make sure you give the  path to the file/folder as argument, either the complete path, or from current directory: ```./input_directory/``` . You can also provide an output folder, where files will be saved, otherwise, input directory is default. You can also choose whether to analyze by positive or negative selction, run downstream GO enrichment analysis, select different output formats and top cutoff.<br/>

Examples:
```bash
Rscript MaGplotR.R -i path_to_results_directory/
```
```bash
Rscript MaGplotR.R -i path_to_results_directory/ -c path_to_control_file -s neg -t 50 -p png -o path_to_output_directory/ -r path_to_sgRNA_input_directory -g MF
```
### Options:
Mandatory arguments:<br/>
-i (input directory): **path** to the folder where gene summary files (test files) are located. All files in this folder will be taken as input.<br/>

Optional arguments:<br/>
-c: (control file): **path** to the control file (no control as default). <br/>
-o: (output directory): **path** to an existing folder where output files will be saved (input directory default).<br/>
-s: (selection): write pos or neg, for positive or negative selection (pos is default).<br/>
-r: (sgRNA input directory): **path** to an existing folder where sgRNA summary files are saved.<br/>
-t: (top cutoff): number of hits to be shown in heatmaps. 25 is default.<br/>
-p: (plot format): just write one among these (png is default): png , pdf, ps, jpeg, tiff, bmp.<br/>
-g: (GO type of terms): write BP, MF or CC (no GO analysis as default).<br/>


## Output plots and files
#### Boxplot<br/>
Representation of all gene LFCs in each experiment (and control if supplied). Gives a quick view of selection / scattering for every experiment.
![genes_boxplot](https://user-images.githubusercontent.com/95416488/172137508-eac5ca52-4b8b-4362-a672-c6c89ad4f35c.png)


#### Heatmap with control<br/>
Heatmap represents gene ranks in each experiment. Control plot shows the LFC of control (cyan) and the mean LFC of all experiments (red) for each gene.<br/>
![heatmap_and_control](https://user-images.githubusercontent.com/95416488/165736171-025fa334-7f20-4196-ae27-9ba433f86435.jpg)

#### Reactome Pathway Analysis<br/>
Reactome PA of top 1 % hits<br/>
![Screenshot from 2022-05-04 15-23-43](https://user-images.githubusercontent.com/95416488/166690308-e08bb1cd-734a-43e0-9278-154074f44b03.png)

#### Cluster plot<br/>
![Screenshot from 2022-05-12 14-46-15](https://user-images.githubusercontent.com/95416488/168256030-fb922cac-a18c-40df-8be2-5eda7cffe121.png)


#### Example of terminal display:<br/>
![Screenshot from 2022-06-06 11-46-42](https://user-images.githubusercontent.com/95416488/172137707-6d300a97-5fb1-4e31-a77b-a8341b8790d4.png)

