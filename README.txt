A shiny app software (created by R) to compare statistical analysis of gene expression values of normal vs abnormal sex-chromosomed sample types. 
The main task of this project is to develop a tool application using R language programming which
enables the user to visualize and thoroughly investigate the different analyses results all together
at the same time, save or compare between them and finely customize display settings in order to
serve his purpose and help him to extract relevant biological outcomes. This software application
is valid for any type of data set by slight adjustment of titles according to the data set. This could
help not just the researchers but also the normal users for an easier result survey.
The software designed performs all the statistical analysis tests on the provided data set simultaneously
and display the produced results in the display area in the form of tabs (a tab for each). These
tests (discussed below) involves the comparison of a pair of sample types at a time and includes
t-SNE analysis, PCA analysis and linear fitting with p-value calculation which is compared to log
fold change to determine the differentially expressed genes. The results are displayed in the form
of plots, heatmaps and tables. The differentially expressed genes determined from the previous
analyses then undergoes a gene ontology analysis to determine the most significant GO terms from
which biological findings are obtained. Then from the four different sample types the effect on gene
expression of the two abnormal samples is compared with the two controls.

The files attached are:
1- App.R : the front-end application software
2- Appcode.R : the back-end application operating code (statistical analyses and plotting functions)
3- Preprocessing.R : the preprocessing code for the raw illumina data (using lumi package)
4- Manhattan and GO annotation.R : code for adding some gene annotations for manhattan plot and gene ontology
