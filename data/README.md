# Population genetics and invasion history of the European Starling across Aotearoa New Zealand

[https://doi.org/10.5061/dryad.6djh9w1bd](https://doi.org/10.5061/dryad.6djh9w1bd)

## Description of the data and file structure

**Data Files:**

starling_noduprel_qual_miss_filt.recode.vcf - final filtered genetic data file.
Metadata_NZ_AU_UK_BE_ReplicatesSibRemoved2.csv - metadata file containing coordinate and location information for all samples in the final data file.

**Scripts:**
01_sequencecleaning_code.pdf: Trimming data from all the the sequence data batches. 
02_SNPcalling_code.pdf: Calling SNPs using BCFtools. 
03_SNPfiltering_code.pdf: Code for the SNP filtering methods. 
04_BioinformaticSexing_code.pdf: Profiling the sex of individuals using sex chromosome heterozygosity levels. 
05_Popgen_Selection_code.pdf: Baypass analysis and plots for assessing signals of selection within invasive lineages. 
06_SFS_code.7z: Code for the folded site frequency spectrum fitering and plots. 
07_jaccard_diversity_code.ipynb: Genetic pairwise dissimilarity using Jaccard's.
08_MMRR_code.R: Isolation by distance and environment analysis, performed using MMRR r package.
09_PSMC_code.ipynb: demographic analysis using PSMC on whole genome resequencing data.
10_Stairwayplot2_code.pdf: demographic analysis using StairwayPlot2 on the site frequency spectrum data taken from the reduced representation sequencing data.

### Files and variables

#### File: Metadata_NZ_AU_UK_BE_ReplicatesSibRemoved2.csv

**Description:** 

##### Variables

* id:Individuals ID matching VCF file ID
* pop: Sampling locations for individuals.
* pop2: Population location groupings for individuals used in the manuscript.
* Con: Continent code for the population.
* lat: Latitude of capture site.
* lon: Latitude of capture site.
* loc: More detailed discription of capture/colleciton/source location.

#### File: starling_noduprel_qual_miss_filt.recode.vcf

**Description:** 
Final genetic variant file after filtering.

**Access information:**

The raw DArT-seq data have been deposited under BioProject accession no. PRJNA1164936 (MRL samples) and PRJNA1168657 (AUK, WHA, WEL, and CAN samples) in the NCBI BioProject database ([https://www.ncbi.nlm.nih.gov/bioproject/](https://www.ncbi.nlm.nih.gov/bioproject/)). The whole genome resequencing data of the three MRL individuals used in the psmc analysis, along with an additional three individuals from MRL individuals that were not presented in any analysis in this manuscript, have been deposited under BioProject accession no. PRJNA1165315.