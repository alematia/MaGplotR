# ScrispR | MaG-scRipt
A software for Multiple Genetic Screens data visualization

[Dependencies](#dependencies)<br/>
[Usage](#usage)<br/>

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
Rscript ScripsR.R -i path_to_results_directory -c path_to_control_file -p png -o path_to_output_directory -s path_to_sgRNA_input_directory
```

Mandatory arguments:<br/>
-i (input directory): **path** to the folder where gene summary files (test files) are located. <br/>

Optional arguments:<br/>
-c: (control file): **path** to the control file (no control as default). <br/>
-o: (output directory): **path** to an existing folder where output files will be saved (input directory default).<br/>
-g: (sgRNA input directory): **path** to an existing folder where sgRNA summary files are saved.<br/>
-p: (plot format): just write one among these (png is default): png , pdf, ps, jpeg, tiff, bmp.<br/>

