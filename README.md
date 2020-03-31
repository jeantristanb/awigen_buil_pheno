# build_genotype_awigen
build and transform data for AWIGEN project

## what need :
R : library : kableExtra knitr scales
## how to used 
Rscript ./launchAnalyse_all.r [param]
## parameters
###parametre :
 * `pheno_file` : file phenotyep
  * `sep` : separator : , ; tab space
  * `val_sex` : valeur for females in sex
  * `head_id` : head of id must be 2 separate by comma [FID,IID]
  * `head_sex` : head of sex, if not, no analyse by sexa
  * `head_site` : head of sex, all individuals are assigned to all
 * `params` : file of param, see section file params
 * `tr_var` : file contains transformation 
 * `fam_file`: fam file to compare individuals common individual in pheno and sex
 * `dir_out` : dir out ['default ./']
 * `output` : dir out ['default output']

###file params : `params`

### file transformation : `tr_var`
## old version 
  * contains old version used previously in awigen project
