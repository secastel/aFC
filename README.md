
# aFC
**a**llelic **F**old **C**hange

Calculates allelic fold change (aFC) using standard input files for [fastQTL](http://fastqtl.sourceforge.net).

Method developed by [Pejman Mohmammadi](mailto:pmohammadi@nygenome.org), software by [Stephane E. Castel](mailto:scastel@nygenome.org) both in the [Lappalainen Lab](http://tllab.org) at the New York Genome Center and Columbia University Department of Systems Biology.

Runs on Python 2.7.x and has the following dependencies: [pandas](http://pandas.pydata.org), [statsmodels](http://statsmodels.sourceforge.net), [scikits-bootstrap](https://github.com/cgevans/scikits-bootstrap), [NumPy](http://www.numpy.org).

#Usage
Requires tabix indexed gzip compressed VCF file containing genotypes and BED file containing phenotypes, identical to the inputs of [fastQTL](http://fastqtl.sourceforge.net), and a list of QTL to calculate aFC for. If provided, covariates will be regressed out of the phenotype values. Outputs the aFC and corresponding 95% confidence interval for each input QTL.

#Arguments
##Required
* **--vcf** - Tabix indexed and gzipped VCF file containing sample genotypes. See [fastQTL](http://fastqtl.sourceforge.net/pages/format_vcf.html) for format details.
* **--pheno** - Tabix indexed and gzipped BED file containing sample phenotypes. See [fastQTL](http://fastqtl.sourceforge.net/pages/format_bed.html) for format details.
* **--qtl** - File containing QTL to calculate allelic fold change for. Should contain tab separated columns 'pid' with phenotype (gene) IDs and 'sid' with SNP IDs. Optionally can include the columns 'sid_chr' and 'sid_pos', which will facilitate tabix retrieval of genotypes, greatly reducing runtime.
* **--geno** - Which field in VCF to use as the genotype. By default 'GT' = genotype. Setting to 'DS' will use dosage rounded to the nearest integer (IE 1.75 = 2 = 1|1).
* **--chr** - Limit to a specific chromosome.
* **--log_xform** - The data has been log transformed (1/0). If so, please set --log_base.
* **--o** - Output file.

#Optional
* **--cov** _()_ - Covariates file. See [fastQTL](http://fastqtl.sourceforge.net/pages/format_cov.html) for format details.
* **--matrix_o** _()_ - Output the raw data matrix used to calculate aFC for each QTL into the specific folder.
* **--boot** _(100)_ - Number of bootstraps to perform for effect size confidence interval. Can be set to 0 to skip confidence interval calculation, which will greatly reduce runtimes.
* **--ecap** _(log2(100))_ - Absolute aFC cap in log2.
* **--log_base** _(2)_ - Base of log applied to data. If other than 2, data will be converted to log2.

#Output File
sid	pid	sid_chr	sid_pos	log2_esize	log2_esize_lower	log2_esize_upper

* 1 - **sid** - Variant ID.
* 2 - **pid** - Phenotype (gene) ID.
* 3 - **log2_aFC** - allelic Fold Change in log2.
* 4 - **log2_aFC_lower** - Lower estimate of 95% confidenace interval of log2(aFC).
* 5 - **log2_aFC_upper** - Upper estimate of 95% confidenace interval of log2(aFC).
