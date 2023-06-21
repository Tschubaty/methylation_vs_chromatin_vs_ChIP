#################################################################
##  CpG Islands according to HG38
##  
##  input: cpgIslandExt.hg38.bed
##  wget -qO- http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cpgIslandExt.txt.gz
##
##  output:      
##  v01 - 04.03.2021
##  Author: Daniel Batyrev 777634015
#################################################################

#Clear R working environment 
#rm(list=ls())

this.dir <- dirname(rstudioapi::getSourceEditorContext()$path)
#this.dir <- "/ems/elsc-labs/meshorer-e/daniel.batyrev/Encode3/"

setwd(this.dir)

# field	example	SQL type	info	description
# bin	585	smallint(6)	range	Indexing field to speed chromosome range queries.
# chrom	chr1	varchar(255)	values	Reference sequence chromosome or scaffold
# chromStart	28735	int(10) unsigned	range	Start position in chromosome
# chromEnd	29737	int(10) unsigned	range	End position in chromosome
# name	CpG: 111	varchar(255)	values	CpG Island
# length	1002	int(10) unsigned	range	Island Length
# cpgNum	111	int(10) unsigned	range	Number of CpGs in island
# gcNum	731	int(10) unsigned	range	Number of C and G in island
# perCpg	22.2	float	range	Percentage of island that is CpG
# perGc	73	float	range	Percentage of island that is C or G
# obsExp	0.85	float	range	Ratio of observed(cpgNum) to expected(numC*numG/length) CpG in island
CpG_bed_colum_names <- c("chrom",	"chromStart",	"chromEnd",	"name",	"length",	"cpgNum",	"gcNum",	"perCpg",	"perGc",	"obsExp")
CpG_file_name <- "cpgIslandExt.hg38.bed"
island_bed <- read.csv(file = CpG_file_name,
               header = FALSE,
               sep = "\t",
               col.names = CpG_bed_colum_names)
