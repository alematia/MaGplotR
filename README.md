# ScrispR
A software for Multiple Genetic Screens data visualization

## Dependencies
You can either install these packages from R terminal or by running the installation script from linux terminal.<br/>

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
Download `installation_scrispr.R` and run the following command in Linux terminal:
```bash
Rscript installation_scrispr.R
```



### USAGE

Example:
```bash
Rscript ScripsR.R -i path_to_results_directory -f mg -c path_to_control_file -p png -o path_to_output_directory
```

Mandatory arguments:<br/>
-i (input directory): **path** to the folder where all the test files are located. <br/>
-f (format): specifies the input format. Just write mg for MaGeCK gene_summary files or sb for ScreenBEAM files.<br/>

Optional arguments:<br/>
-c: (control file): **path** to the control file (no control as default). <br/>
-p: (plot format): just write one among these (png is default): png , pdf, ps, jpeg, tiff, bmp.<br/>
-o: (output directory): **path** to an existing folder where output files will be saved (input directory default). <br/>

