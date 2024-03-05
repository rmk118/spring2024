----------------------- DOWNLOAD LINKS ------------------------
//Raw Reads Links//

Raw reads for this study are available on the European Nucleotide Database (ENA) under project PRJEB60884 (https://www.ebi.ac.uk/ena/browser/view/PRJEB60884 -> link at time of uploading the README file).

The reads from King were downloaded from ENA (PRJEB50679) and metadata was downloaded form https://doi.org/10.6084/m9.fgshare.19453889.v1. 
King, N.G., Moore, P.J., Thorpe, J.M. et al. Consistency and Variation in the Kelp Microbiota: Patterns of Bacterial Community Structure Across Spatial Scales. Microb Ecol 85, 1265–1275 (2023). https://doi.org/10.1007/s00248-022-02038-0


// Tide height observations //

Downloaded from the Governement of Canada Website (Station 07735) https://www.tides.gc.ca/en/stations/7735 




----------------------- ABOUT THE FILES ------------------------
// CODE - all .R files //
Scripts should be run in order from 1-4 to reproduce the full analysis pipeline. If you only want to run a script to make a figure, find that script (4-script_name.R) and use the phyloseq objects (.RDS files) uploaded for your convinience. 

Script 1 is the dada2 pipeline from "Callahan, B., McMurdie, P., Rosen, M. et al. DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods 13, 581–583 (2016). https://doi.org/10.1038/nmeth.3869" . See dada2 tutorial here https://benjjneb.github.io/dada2/tutorial.html 

Script 2 is assigning the taxonomy 

Script 3 if filtering and rarefying the data

Scripts 4 are the scripts used foe data analysis and making the figures for the manuscript. The title of the script matches they type of analysis.


// .CSV FILES //
These files include the unrarefied taxonomy and otu tables along with the study metadata (bacterial, kelp condition, logger data). 


// .RDS FILES //
These are the files used for bacterial data analysis (scripts 4-script_name.RDS)