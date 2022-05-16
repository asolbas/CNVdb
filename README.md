# CNVdb
Statistical analysis, annotation and functional enrichment analsysis of the IIS-FJD CNV database of allele frequencies. 

### Usage
These scripts were developed to analyse the CNV database of the Instituto de Investigación Sanitaria Fundación Jiménez Díaz (IIS-FJD). They are specific for this database and therefore they will not work for CNV databases that do not share the same structure. 

### Roadmap
The files contained in this repository are: 
- __StatisticalAnalysis.R__: statistical and exploratory data analysis of the database. 
- __DensityPlot.R__: tools for the data visualization of the distribution of the genomic regions contained in the database. 
- __split_matrix.py__: divides of the data base of CNVs in two -one for duplications and another for deletions-. 
- __run_AnnotSV.sh__: bash script to annotate the regions of the CNV data base. 
- __genes.R__: statistical analysis of the genes annotated int the database. 
- __GO_v2.R, KEGG_v2.R, InterPro_v2.R and HPO_v2__: functional annotation of the genes mapped to the regions of the CNV database and enrichment analysis. 

### Contact

Ana Solbas Casajús - a.solbas@alumnos.upm.es

Project Link: https://github.com/asolbas/CNVdb
