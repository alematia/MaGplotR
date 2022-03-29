# ScrispR
A software for Genetic Screens data visualization

# USAGE

Example:
Rscript ScripsR.R -i results_directory/ -f mg -c path_to_control_file -p png -o output_directory


Mandatory arguments:<br/>
-i (input): must be a path to the folder where all the test files are located. <br/>
-f (format): specifies the input format. Just write mg for MaGeCK gene_summary files or sb for ScreenBEAM files.<br/>

Optional arguments:<br/>
-c: <br/>
-p: (plot format): Just write one among these (png is default): png , pdf, ps, jpeg, tiff, bmp.<br/>
-o: <br/>

