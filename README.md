# MaGplotR
![LICENSE](https://img.shields.io/badge/license-MIT-green)
![R-VERSION](https://img.shields.io/badge/R-%204.2.0-blue)
![MaGplotR-VERSION](https://img.shields.io/badge/release-v1.0.1-orange)
#### A software for Multiple Genetic Screens data visualization

![Screenshot from 2022-05-13 13-02-56](https://user-images.githubusercontent.com/95416488/168270389-73f1f6e9-dee3-468a-ae2c-611c599a8aa9.png)


>*MaGplotR* produces visualization plots for MaGeCK screen data (RRA analysis). <br/>
This software is designed to run from a Linux, MacOS or Windows command line, and it is written in R.<br/>
*MaGplotR* directly uses `gene_summary.txt` files generated by the **test** function from MaGeCK software. Output plots provide information on quality control of the screen data and highligths the best hits of the multiple screen experiments.<br/>
Additionally, a control experiment `gene_summary.txt` file can be used to compare the different experiments to a control self-enrichment experiment (or any condition used as control). Also, `sgrna_summary.txt` files can be analyzed for QC as optional arguments.

[Citation](#citation)<br/>
[Motivation](#motivation)<br/>
[Dependencies](#dependencies)<br/>
[Usage](#usage)<br/>



## Citation
If you use this software please cite:<br/>
<br/>
**MaGplotR: a software for the analysis and visualization of multiple MaGeCK screen datasets through aggregation**<br/>
Alejandro Matia, Maria M. Lorenzo, Duo Peng<br/>
bioRxiv 2023.01.12.523725; doi: https://doi.org/10.1101/2023.01.12.523725<br/>


## Motivation
This software was developed to simplify the analysis of MaGeCK screen data outputs from multiple experiments through elegant visualization. Sometimes, several screen experiments are performed testing multiple conditions or with a number of replicates, and it is desired to perform comparisons between these experiments. We put special focus introducing the possibility of adding a control experiment to observe the behaviour of the hits in this experiment in a visually quantitative way that may help discarding false positives (i.e. self enrichment of potential hits in control experiment). Two main values are extracted from MaGeCK test summary files; LFC and rank, presenting them in a way that simplifies the whole set of experiments analysis.
Heatmap representation gives a rapid view of how top hits behave in the set of experiements, making very easy to detect if a particular hit had a poorer score in a experiment. The number of hits shown in heatmap plots can be adjusted when running the software from the command line (see options below). Also, the intermediate files to create the heatmap plots are generated and saved in the output directory so the user can check the ranking score of every gene. To assay the feasibility of hits, a plot representing the expression of the hit genes in the most relevant cell lines is also generated with data from [The Human Protein Atlas](https://www.proteinatlas.org/about/download).<br/>

As extra features, the top 1 % of positive and negative hits are analyzed with [Reactome Pathway Analysis](https://bioconductor.org/packages/release/bioc/html/ReactomePA.html) and [clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html) generating plots for the most enriched pathways and clusters among the data sets. Optionally, GO analysis is available (see options below). These analysis are run using [Genome wide annotation for Human](https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html).


## Getting started
```bash
git clone https://github.com/alematia/MaGplotR.git
```

## Dependencies
You can either install these packages from R or by running the installation script from Linux/MacOS terminal.<br/>

- Install from R. Execute the following commands in R:<br/>
```r
install.packages("tidyverse")
install.packages("optparse")
install.packages("reshape2")
install.packages("matrixStats")
install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("ReactomePA")
BiocManager::install("clusterProfiler")
```
- Install using installation script:<br/>
Download `installation_magplotr.R` and run the following command in Linux or MacOS terminal:
```bash
Rscript installation_magplotr.R
```
*Note: some of these packages like tidyverse may have additional dependencies: install libcurl, openssl, libxml-2.0, libfontconfig1-dev and libfreetype6-dev libraries. Example for debian based systems:

```bash
sudo apt install libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev libfreetype6-dev libtiff5-dev libharfbuzz-dev libfribidi-dev
```


## Usage
#### Recommendations:
1) Install R and its dependencies as shown above.<br/>
2) Clone or download the repository and put the MaGplotR folder either in your Downloads or your Home directory.<br/>
3) Put all the gene summary files you want to analyze in a designated folder. Do not put your control file in this folder. This will be the input directory (-i). These files names usually end with ".gene_summary.txt" as MaGeCK output. MaGplotR uses the experiment name from the filename (i.e. for "exp7.gene_summary.txt" file, "exp7" will be displayed in the plots).<br/>
4) For simplicity, we recommend you to put the input folder, the control file, and the output folder (optionally create it first) in the MaGplotR folder where the program files are stored.<br/>
5) Next, open a terminal in Linux, MacOS or Windows and type the commands as in the examples. Make sure you give the path to the file/folder as argument, either the complete path, or from current directory: ```./input_directory/``` . You can also provide an output folder, where files will be saved, otherwise, input directory is default.
6) WINDOWS: if you are using the tool from a Windows terminal, make sure to assign R to the PATH first like this: ```$env:Path += ';C:\Program Files\R\R-4.2.1\bin\x64\' ```
<br/>

Examples:
```bash
Rscript MaGplotR.R -i example_test_files/
```
```bash
Rscript MaGplotR.R -i path_to_results_directory/ -o path_to_output_directory/ -c path_to_control_file -s neg -t 50 -p png -r path_to_sgRNA_input_directory -g BP,MF -b y
```
### Options:
Mandatory arguments:<br/>
-i (input directory): **path** to the folder where gene summary files (test files) are located. All files in this folder will be taken as input.<br/>

Optional arguments:<br/>
-c: (control file): **path** to the control file (no control as default). <br/>
-o: (output directory): **path** to an existing folder where output files will be saved (input directory default).<br/>
-s: (selection): write pos or neg, for positive or negative selection (pos is default), i.e.: `-s neg`.<br/>
-r: (sgRNA input directory): **path** to an existing folder where sgRNA summary files are saved.<br/>
-t: (top cutoff): number of hits to be shown in heatmaps (25 is default), i.e.: `-t 50`.<br/>
-x: (threshold): top % of hits to be used as for Pathway and Gene Ontology analysis. 1 % is default. i.e.: `-x 1.5`.<br/>
-p: (plot format): just write one among these (png is default): png , pdf, ps, jpeg, tiff, bmp, i.e.: `-p pdf`.<br/>
-g: (GO categories): write BP, MF or CC (no GO analysis as default) i.e.: `-g BP`. Also write several parameters at once: i.e.: `-g BP,MF,CC`.<br/>
-b: (colour blind): write y or n (no is default), i.e.: `-b y`.<br/>


## Output plots and files
#### Boxplot<br/>
Representation of all gene LFCs in each experiment (and control if supplied). Gives a quick view of selection / scattering for every experiment.
<p align="center">
<img src="https://user-images.githubusercontent.com/95416488/208246015-28e49484-bed7-404f-8f47-246fbe34d701.png" width="650" height="650" align="center">
</p>

#### PCA plot<br/>
PCA dimensionality reduction of the rank scores of every experiment. Groups experiments by similarity.
<p align="center">
<img src="https://github.com/alematia/MaGplotR/assets/95416488/b0e0803d-092b-47a4-8e3d-7c773de1be31.png" width="650" height="650" align="center">
</p>

#### Heatmap with control<br/>
Heatmap cells are filled with each gene's LFC. Numbers inside the cells are gene ranks in each experiment. Control plot shows the LFC of control (cyan) and the mean LFC of all experiments (red) for each gene.<br/>

<p float="left">
  <img src="https://user-images.githubusercontent.com/95416488/210894756-0dd52b21-4b78-4300-895e-65293ee4fd3a.png" width="700" height="420" align="center">
  <img src="https://user-images.githubusercontent.com/95416488/236638451-19dfab28-738b-4789-8471-f9e1633aa420.png" width="100" height="70" align="center">
 </p>

#### Colorblind heatmap<br/>
When using the option ```-b y```

<p float="left">
  <img src="https://user-images.githubusercontent.com/95416488/210895190-5a1208d0-c7f3-48a5-b454-d7ee1e9d7acc.png" width="700" height="420" align="center">
  <img src="https://user-images.githubusercontent.com/95416488/236638451-19dfab28-738b-4789-8471-f9e1633aa420.png" width="100" height="70" align="center">
 </p>


#### Expression plot<br/>
Is only generated if MaGplotR folder is located in Home or Downloads directory, using the expression file that is inside of it.
<p align="center">
<img src="https://user-images.githubusercontent.com/95416488/210894845-83a55a75-cc65-45b0-9bb7-72c3d197527e.png" width="464" height="650" align="center">
</p>

#### Reactome Pathway Analysis<br/>
By using the top 1 % of hits, a Reactome Pathway Analysis is performed. The plot represents the number of genes (Count) that belong to a pathway.<br/>
![ReactomePA_pos](https://user-images.githubusercontent.com/95416488/210894920-c16d4531-98c6-4356-aae5-7ca203932773.png)


#### Cluster plot<br/>
Is only generated when the number of experiments (input) is > 2. Genes used for clustering are the top % genes chosen by user.
![Screenshot from 2022-12-19 13-49-51](https://user-images.githubusercontent.com/95416488/208430073-2fa6c062-4d20-43d1-860a-3ddad6139383.png)


#### Example of terminal display:<br/>
![Screenshot from 2022-12-19 13-46-19](https://user-images.githubusercontent.com/95416488/208429514-11102732-432f-425e-af1d-450d348037ec.png)

