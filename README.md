Overview:

This repo contains the datasets and code to reproduce a preliminary figure in “Identification of Aberrant B Cell Activity in Ulcerative Colitis Patients” by Jeff Hsiao and Constantine Tzouanas for 20.440 at MIT. In order to study the contributions of B cells in ulcerative colitis (UC) development, we will be analyzing the transcriptional programs of these cells in health and disease based on published datasets of single-cell RNA sequencing (Smillie et al., 2019). B cells are of particular interest due to their association with multiple autoimmune and chronic inflammatory diseases. The code provided is implemented in R and incorporates Seurat v3 (Stuart et al., 2019). Currently, it performs dataset quality control, dimensional reduction for PCA, k-nearest neighbors clustering, dimensional reduction for UMPA visualization, and health vs. UC dataset differentiation. Future work involves identifying the B cell populations via marker genes to distinguish the cell types in each cluster. We will then compare the transcriptional programs of B cells between healthy and UC patients.

Citations:
1.	Smillie, C.S., et al. "Intra- and Inter-cellular Rewiring of the Human Colon during Ulcerative Colitis," Cell, 178(3), 714-730 (2019). DOI: 10.1016/j.cell.2019.06.029.
2.	Stuart, T., et al. "Comprehensive Integration of Single-Cell Data," Cell, 177(7), 1888-1902 (2019). DOI: 10.1016/j.cell.2019.05.031.

Data:

Smillie et al. obtained 68 biopsies from colonoscopy examinations of 12 healthy individuals and 18 UC patients. They then performed single-cell RNA sequencing using the 10x platform. These samples led to them generating 366,650 high-quality single-cell transcriptomes. The data containing immune cells from their work can be found at the Broad Institute Single-Cell Portal (https://singlecell.broadinstitute.org/single_cell/study/SCP259/intra-and-inter-cellular-rewiring-of-the-human-colon-during-ulcerative-colitis#study-summary). The file entitled “gene_sorted_Imm.matrix.mtx” should be downloaded into the “RawData_Folder” folder in the repo prior to running the code due to file size limitations on GitHub. To process the data, we perform the functionalities of the code detailed in the “Overview” section.


Folder Structure:

This repo contains four files (two hidden) and two folders. The .gitignore file is a hidden file for git to ignore file types specified within when changes are pushed to or pulled from the repo. The .gitignore file here is based on https://github.com/github/gitignore for R. The second hidden file is .gitattributes, specifying the file types that require git-lfs; herein, it’s specified to do so for *.mtx files and *.tsv files. The other two files with this repo are this README (summarizing the project and contents of this repo) and the code (implemented in R for data analysis). The “RawData_Folder” folder contains raw data downloaded as specified in the “Data” section of this README. The “ProcessedData_Folder” contains output data from the .R code. The “Plot_Folder” is the output location of the figures from the .R code.

Installation:

To reproduce the figure in “Plot_Folder”, please first download the .R code as well as all the folders contained within this repo. Please note that the code and folders should all be in the same overall folder. The file entitled “gene_sorted_Imm.matrix.mtx” (https://singlecell.broadinstitute.org/single_cell/study/SCP259/intra-and-inter-cellular-rewiring-of-the-human-colon-during-ulcerative-colitis#study-summary) should then be downloaded into “RawData_Folder.” In order to run the code, Seurat, data.table, ggplot2, gridExtra, gtable, plotly, scales, spam, dplyr, and patchwork should be installed. Finally, run the .R code by a method of your choice. If any library listed is not already installed, this can be done by "install.packages("<the package's name>")" in an R session command line.
