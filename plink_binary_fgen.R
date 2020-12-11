library(data.table)
library(magrittr)
library(dplyr)
install.packages("tidyverse")
library(tidyverse)


#STEP 1: Create family id for each individual 
f_name <- "new_six_fam_eff"
snps_n <- 100

require(data.table)
require(magrittr)
require(tidyverse)
print(paste0("File header is: ", f_name))
print(paste0("Number of snps to spike effect in: ", snps_n))
# Is it always 100 SNP effects? If so, why?


#STEP 1.1: Make sure base plink files have bed, bim, and fam name before continuing
if(is.na(f_name)) { 
  stop("Error: Need a base plink file (bed, bim, bam, fam)")
}
#STEP 1.2: Make sure we have at least 1 snp for phenotype creation
if(snps_n < 1 | is.na(snps_n)) { 
  stop("Error: Need at least one snp to pull from plink to make pheno")
}

# Creates effect distribution for number of snps and writes it out
set.seed(414)
#eff <- runif(snps_n, min = -.2, max = .2)
## Null first
eff <- c(rep(log(1.5), 60), rep(log(.5), 40))

print("Effects created")
print(eff)
write.table(eff, file = "output/eff.txt", sep = "\t", col.names = F, 
            row.names = F,
            quote = F)

# Call to take trans file and make into plink- make sure string call is correct
bed_call <- paste0("plink -tfile output/trans --make-bed --out ", "output/",f_name)
print(bed_call)
#TRY to call printed file 
try(system(bed_call, intern = F))

