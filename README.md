# Glanders Paper Scripts

The scripts in the repository are for the paper titled "Integrative genomic and transcriptomic analysis provides a more rational approach for improvement of glanders serodiagnosis".

All scripts are provided as is, for visibility purposes only. No attempt has been made to optimise the code.

Here is a short description of the contents of this repository:
- OperonAnalysis.py is used to combine the results from OperonSEQer with the annotation (GFF) files in order to obtain more information about the predicted operons. Once complete, operons were screened for presence of disruptions by structural variants. Afterwards, promoter regions were predicted for these operons. SNPs and INDELs found in these regions were identified to create a list of potentially disrupted promoters. 
- The VariantCircosPlots folder contains bash scripts and configuration files to run circos. These scripts and configurations files were used to generate the circos heatmap figures showing the distribution of different variant types across the genomes. These files were adapted from https://github.com/mldmort/Mortazavi2021_B6B10. In order to generate the figures, run `bash 1-cmd_prepare_data.sh` first to generate the required files. Then run `bash 2-cmd_run_circos.sh` to create the figures.
