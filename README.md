# *** CONCATABOMINATION PIPELINE ***

1) TO INSTALL:

You need a unix-based system (Mac OS X or Linux). You will also need to have perl installed.
Type in "perl -v" in the command line (or in terminal window in Mac) to find out the version you have. If you don't have perl you can get it from:
http://www.perl.org/get.html

PIPELINE INSTALLATION:
In a terminal window (or command line in Linux) type:
tar -xvzf concatabomination-pipeline-v4.1-parallel.tgz
cd concatabomination-pipeline-v4.1-parallel
sh make.sh


For visualisation of instability of taxa, you will need Cytoscape. 
Download and install Cytoscape in the following link:
http://www.cytoscape.org/download.html

Download version 2.8.3, the latest of the 2.x Series, for your corresponding Operating System. 
(this version is requested because the tutorial is written to use this version of Cytoscape). If you are familiarised with Cytoscape, then give it a go with a newer version.

Go to step XX to follow tutorial on how to open and visualize the result file in Cytoscape 2.8.3.


2) TESTING
Once finished with the installation, test it with one of the attached test-datasets (which should be in the same folder as the "concatabomination_pipeline_v4.1_parallel"):
type:
concatabomination_pipeline_v4.1_parallel -i Gauthier_1986.coding.nex

The screen will print the following:

Running concatabomination pipeline, starting at step 1 and ending at step 12
starting Step 3
starting Step 4
starting Step 5
starting Step 6 (needed for steps 8 onwards)
starting Step 7 (needed for steps 8 onwards)
starting Steps 8 & 9
starting Step 10
starting Step 11
starting Step 12

And produce several output files. 

If you had no error messages, now you are ready to go!!
IMPORTANT: Don't forget that you need your input files to be in simple non-interleaved Nexus format, like the example files.


3) QUICK TUTORIAL

To start the pipeline and see the different options, type:
concatabomination_pipeline_v4.1_parallel -h

This will print the following to screen:

OPTIONS:
                
        -i      <input file>    
                simplified nexus format with names of taxa in quotes
                
        -n      <number>                
                Which PerlEQ output to filter (table including non-informative sites = 2; or excluding them = 1) Default = 2

        -s      <number>                
                optional set to the First step of the pipline to run, default = 1
                
        -e      <number>                
                optional set to the Last step of the pipline to run, default = 12
                
        -f      <filtered|non-filtered> 
                "filtered" deletes potential taxonomic equivalents prior to concatabomination (Default). 
                "non-filtered" includes all taxa in concatabomination.
                
        -k      <TRUE|FALSE>                    
                keeptmp. Specifies whether to keep temporary files. By default set to FALSE.
                
        -b      <TRUE|FALSE>                    
                binary coding matrix. If set to "TRUE", assumes coding matrices only include "1", "0", "-" or "?". 
                Uses fast method to calculate compatibility. Default = "FALSE";
                
        -p      <number>                
                the number of processors to use at the compatibility step. Default = 1;
                Maximum processors detected on this system = $MAXPROC
        
        -h
                Displays this message


The usage is:
concatabomination_pipeline_v4_parallel -i inputfile.nex options


