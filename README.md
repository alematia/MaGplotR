# ScrispR
A software for Genetic Screens data visualization

# USAGE

Example:<br/>
Rscript ScripsR.R -i results_directory/ -f mg -c path_to_control_file -p png -o output_directory/


Mandatory arguments:<br/>
-i (input directory): **path** to the folder where all the test files are located. <br/>
-f (format): specifies the input format. Just write mg for MaGeCK gene_summary files or sb for ScreenBEAM files.<br/>

Optional arguments:<br/>
-c: (control file): **path** to the control file (no control as default). <br/>
-p: (plot format): just write one among these (png is default): png , pdf, ps, jpeg, tiff, bmp.<br/>
-o: (output directory): **path** to an existing folder where output files will be saved (input directory default). <br/>

