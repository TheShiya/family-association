######## PARAMETERS #############
# Odds ration/ effect size : 1.5 for 60 sleected SNPs
# base prevelance for random trait : 0.1
# family effect: 
# individual effect:
##################################
library(data.table)
library(magrittr)
library(dplyr)
library(tidyverse)
#Creates family id



f_name <- "new_six_fam_eff"
snps_n <- 100

require(data.table)
require(magrittr)
require(tidyverse)

print(paste0("File header is: ", f_name))
print(paste0("Number of snps to spike effect in: ", snps_n))


if(is.na(f_name)) { 
  stop("Error: Need a base plink file (bed, bim, bam, fam)")
}

if(snps_n < 1 | is.na(snps_n)) { 
  stop("Error: Need at least one snp to pull from plink to make pheno")
}

# Creates effect distribution for number of snps and writes it out
set.seed(414)
#eff <- runif(snps_n, min = -.2, max = .2)
## Null first
# Pick 100 snps from genome and test the effect of these snps in GWAS, despite noise of family pods of relnatedness of error, see effect that is picked up. 
# there values of 100 are the causal effect that are going to be given to snps selected at random. We want to see consisnent positive an negative effect on trait
#  
eff <- c(rep(log(1.5), 60), rep(log(.5), 40))

# log(1.5), log(0.5)

print("Effects created")
print(eff)
write.table(eff, file = "output/effect.txt", sep = "\t", col.names = F, 
            row.names = F,
            quote = F)

# Call to take trans file and make into plink
bed_call <- paste0("plink -tfile output/trans --make-bed --out ", "output/",f_name)
print(bed_call)
try(system(bed_call, intern = F))

# Call to get mafs for each allele
# we may not want to sample just ant random snp, as some snps may just occur in one person and not affet the rest of the chosen population. THerefore we assess
# the minor allele frequencey which provides the second second most comon allele in the given population- we will use this to smple snps
maf_call <- paste0("plink -bfile ", "output/", f_name," --freq --out ", "output/",f_name)
print(maf_call)
try(system(maf_call, intern = F))

# reads maf, selects 100 snps above .4 frequency randomly and writes it out- arange by allele frequency 
maf <- fread(paste0("output/",f_name, ".frq")) %>% arrange(desc(MAF)) %>% as_tibble()
set.seed(414)
###### This will need be MAF .05
# These are selected snps, with MAFs >= .09 & MAFs <= .102 
snps_af <- maf %>% filter(MAF >= .090 & MAF <= .102) %>% sample_n(snps_n) %>% 
  select(SNP, MAF) %>% arrange(SNP)
snps <- snps_af$SNP
print("Snps chosen from plink file")
print(snps)
# write out snp files
write.table(snps, file = "output/selected_snps.txt", sep = "\t", col.names = F, 
            row.names = F,
            quote = F)


# LD scores - find null that are not in LD with effect snps lower than .1
# creation of ld file which will have linkage disequilibruim score which is how much is each snp is correlected with each oher snp
# trailing snpa in manhattan plot may be in ld with effect snps. Filter these snps at one point. This is just testing for snp correlation.
ld_call<-paste0("plink --bfile output/",f_name,
                " --ld-snp-list output/my_snps.txt --out output/", 
                f_name, " --r2 --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0")
print(ld_call)

ld<-fread(paste0("output/", f_name,".ld")) %>% as_tibble()
ld
# chromosone snp is on, base pair for snp, and at the end is correlation btw snp
# correlation coefficent of snp
# if in perfect correlation, the r^2 would have a value of 1
# if not as frequent, then now r^2 would be low- we dont end up using much of this but just for record keeping it is there!
snps_no_ld<-ld %>% filter(SNP_A %in% snps) %>% arrange(SNP_A) %>% filter(R2 <= .1)

write.table(snps_no_ld, file = "output/snps_not_in_ld_with_selected.txt", 
            sep = "\t", col.names = T, 
            row.names = F,
            quote = F)

ld_check<-fread("output/snps_not_in_ld_with_selected.txt")

#Writes out a merged snp and effect size file
snp_eff_write <- cbind(snps, eff) %>% as_tibble() %>%  mutate(eff = as.numeric(eff)) 
write.table(snp_eff_write, file = "output/snp_eff.txt", sep = "\t", col.names = F, 
            row.names = F,
            quote = F)

print("snps and associated effects")
print(snp_eff_write)

# Subsets based on the snps that were wrote out earlier in plink
# taking plink files created earlier and extracting selected snps from file - essentially just those snps by the ~16k individual
sub_call <- paste0("plink -bfile ", "output/", f_name, 
                   " --extract output/my_snps.txt --recodeAD --out ", "output/",f_name)
sub_call
try(system(sub_call, intern = F))
# Reads the plink file and arranges by FID for the genetypes
# saving the extracted snps into .raw file 
rs <- fread(paste0("output/", f_name, ".raw")) %>% arrange(FID) %>% as_tibble()
rs
# cormatting: IID for individual code, PAT, MAT, SEX, addetive genetic effect (_C), and dominamt genetic effect
# addetive - homozygous recessive - trait does not show up - 0 for less dominant allele
# heterozygous - one allele that is dominant one that is recessive, then you will have a 1
# homozygous dominant - 2

# dominant model- _HET- snps with dominant allele 
#formatting corretly 
pc<-c("P1", "P2", "C1", "C2", "C3", "C4", "C5","C6")
pc_order<-rep(pc, length(unique(rs$FID)))
pc_correct <- c("P1", "C1", "C2", "C3", "C4", "C5", "C6", "P2")

# Correct order for families
rs_ordered<-rs %>% mutate(pc_rel = pc_order) %>% 
  group_by(FID) %>% slice(match(pc_correct, pc_rel)) %>% select(-pc_rel)


# Take additive model takes specifically the additive subsetted genotypes and 
# writes them out

# Selecing only addative model snps 
rs1 <- rs_ordered %>% select(names(rs)[!str_detect(names(rs), "HET")]) %>% 
  ungroup(FID) %>% 
  select(-(FID:PHENOTYPE))

# selecting only dominant model- anything with a 1 or a 2 is a 1
rs2<-sapply(rs1, function(x) ifelse(x==2,1,x)) %>% as_tibble()

# writing out different genetic effects
write.table(rs1, file = "output/genos_add.txt", sep = "\t", col.names = T, 
            row.names = F,
            quote = F)

write.table(rs2, file = "output/genos_dom.txt", sep = "\t", col.names = T, 
            row.names = F,
            quote = F)

# Reads in the previously written out files
rs1 <- fread("output/genos_add.txt")
snp_eff <- fread("output/snp_eff.txt")
# Applies effects as multiple of allele for each individual (to be added into pheno)

#######################################################################
# This fixes previous bug where snp effects were out of order when applied to
# genos
# get rid of _C to properly align with SNP names compared to SNP effects
names(rs1) <- substr(names(rs1),1,nchar(names(rs1))-2)
# set by snp_order
setcolorder(rs1, as.character(snp_eff$V1))
##################################################
# Double check the names are okay
stopifnot(snp_eff$V1 == names(rs1))

# snp effect matrix creation- mutiplying the effect the effects to entirety of genotype matrix
# for snps that have causal effect
snps_eff_mat <- as.matrix(rs1) %*% diag(snp_eff$V2)

# Keep zeroes 0-->0, 1-->1 2-->1
#snp_eff_mat_dom <-#

# sumation of effects
# effects are summed row wise, causual phenotype built out of genotype- addetive model
r1_snp_eff <- rs %>% select(IID) %>% bind_cols(rowSums(snps_eff_mat) %>% 
                                                 as_tibble())
print("These are the snp add effects per individual")
print(r1_snp_eff)

# write out
write.table(r1_snp_eff , file = "output/iid_snp_eff.txt", sep = "\t", col.names = T, 
            row.names = F,
            quote = F)





# This creates correct family ids for everyone
fid<-NULL
q<-1
fid_all<-NULL
fam_seq<-seq(from=1, to=2000, by=10)


# Added a second chunk for the second generation of children
# the tfam of course needs to be modified
for (q in 1:200){
  
  start_fid <- fam_seq[q]
  end_fid <- fam_seq[q]+9
  
  fid<-c(rep(start_fid:end_fid, each=2), rep(start_fid:end_fid, times=2), 
         rep(start_fid:end_fid, times=2), rep(start_fid:end_fid, times=2))
  
  fid_all<-c(fid_all, fid)
}


#fid_all[1:40]
#phenotye

library(mvtnorm)

# base matrix - has to be mirrored from start to end
# starting with parent and followed by 6 children and ending with spouse
# coeficcent of relatedness with yourself is 1/2 = 0.5
                               # your children is 0.5/2 = 0.25
                               # your spouce is 0/2 = 0.
# matrix of relatedness- multivarient normal distibution bc we use coevairent matrix to influence in multiple random normals for individuals
# coevareience matrix- when adding noise in family effect, we create a shared enviornmental impact between family members- genes + enviornment = phenotype

eight_fam_mat <- matrix(
  c(
    0.5, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0,
    0.25, 0.5, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25,
    0.25, 0.25, 0.5, 0.25, 0.25, 0.25, 0.25, 0.25,
    0.25, 0.25, 0.25, 0.5, 0.25, 0.25, 0.25, 0.25,
    0.25, 0.25, 0.25, 0.25, 0.5, 0.25, 0.25, 0.25,
    0.25, 0.25, 0.25, 0.25, 0.25, 0.5, 0.25, 0.25,
    0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.5, 0.25,
    0, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.5
  ),
  byrow = TRUE,
  ncol = 8
)
beta_fam <- NULL # family effect 
beta_eff <- NULL # casual effect

i <- 0
j <- 0

# This is going to be changed by family size
# Select the appropriate matrix with the family

trait_l <- list()
trait_l_bin <- list()
#Bump up to 50 and test
#Create 

# This creates one binary pheno and one continuous for each family
for (j in 1:1){
  for (i in 1:2000) {
    beta <-
      rmvnorm(n = 1,
              mean = c(0, 0, 0, 0, 0, 0, 0, 0),
              sigma = eight_fam_mat)[1, ] #family effect
    
    beta_fam <- c(beta_fam, beta)
    
  }
  
  print(length(beta_fam))
  
  # normal error 
  epsilon <- rep(rnorm(n = 8000)) #individual error
  
  # base prevelance for random trait 
  prev <- .1 # Prevalence
  
  # function of prevelence used to create binary trait 
  alpha <- log(prev/(1-prev))
  
  # total effect for family and casual effect
  beta_tot <- (beta_fam + rowSums(snps_eff_mat)) # Snps plus family effect
  
  # propability function and converting to exponential to later sample from it 
  trait_l_bin[[j]] <- exp(alpha + beta_tot) / (1 + exp(alpha + beta_tot))
  
  # these should not be data.frames make sure they are vectors
  # continuous trait in list format that is as a function of the beta family, total effect of rows and casual effect, and random noise
  trait_l[[j]] <- beta_tot + epsilon 
  
  beta_fam <- NULL
}

# worth doing another random seed for the purpose of random sampling
seed(917)
samp_func <- function(X) (sample(0:1, size=1, prob=c(1-X,X)))

#normally we would have some random noise with this but for the binary one, because we have sampling prodedure based on probability, we do not need to 
# manually account for random noise (epsilon). 
# take binary probabiliy and samples 
traits_bin<-apply(trait_l_bin %>% bind_cols(), MARGIN=c(1,2), samp_func) %>% 
  as_tibble()
# binary phenotype
names(traits_bin) <- "V11"

# continuous
traits_cont<-bind_cols(trait_l)
names(traits_cont) <- "V1"


read.table(paste0("output/",f_name,".fam")) %>% rename(
  "id" = "V2",
  "disease" = "V1",
  "trait" = "V3",
  "age" = "V4",
  "sex" = "V5",
  "num" = "V6"
) %>% mutate(num=c(1:16000)) %>% 
  mutate(
    family_id = fid_all,
    pop_id = rep(1:200, each = 80)
  ) %>%
  arrange(family_id) %>% mutate(pc_rel = pc_order) %>% 
  group_by(family_id) %>% slice(match(pc_correct, pc_rel)) %>% 
  select(-pc_rel) %>%
  # This step is important because the order of the traits assumes this family 
  # order
  bind_cols(traits_cont, traits_bin) %>%
  arrange(num) %>% 
  as_tibble() %>% 
  select(id, family_id, pop_id, V1, V11) %>% 
  rename("FID"="family_id", "IID" = "id") -> pheno2


# This can be used in plink
write.table(pheno2, paste0("output/",f_name,"_phenotypes_structured.txt"), sep="\t", quote=F, 
            row.names = F, col.names = T)
# This is clean for gemma
write.table(pheno2 %>% select(V1, V11), paste0("output/",f_name,"_phenotypes_structured_one.txt"), sep="\t",
            col.names = F, row.names = F, quote = F)

#### Check with linear regression in R and V1 being pheno with X as the genotype
#### which is the plink file with the 0,1,2


############### New, this subsets families from 8 to 4-8 and then subsets those
############### to families with at least two affectd individuals


# This takes the most recent fam file
fam_file <- fread(paste0("output/",f_name,".fam")) %>% as_tibble()

fam_file %>% group_by(V1) %>% arrange((V1))

fam_temp <- 0
fam_siz <- 0
sub_set_fams <- list()

# Selects children based on US household proportions subsets eventually
# via plink

##########################################################
# Sample how to get a fixed seed to repeat random samples
# Change seed value to see it work
# set.seed(415)
# vec3 <- c()
# for(u in 1:2000){
# vec3[u] <- sample(2:6, size=1, prob=c(0.58,	0.26,	0.10,	0.06, 0.06))
# }
##########################################################

  
for(i in 1:2000){
  fam_temp <- fam_file[fam_file$V1==i,]
 
  fam_size <- sample(2:6, size=1, prob=c(0.58,	0.26,	0.10,	0.06, 0.06))
  
  fam_temp <- fam_temp[1:(2+fam_size),] # ordered 
  
  sub_set_fams[[i]] <- cbind(fam_temp$V1, fam_temp$V2) 
}

# Binds into dataframe so families are subsetted
sub_set_list<-do.call(rbind,sub_set_fams)

write.table(sub_set_list, paste0("output/",f_name,"_sample_fams.txt"), sep="\t",
            col.names = F, row.names = F, quote = F)

sub_set_list <- fread(paste0("output/",f_name,"_sample_fams.txt"))

print(sub_set_list)

print(paste0("Number of individuals after sampling families: ", length(sub_set_list$V2)))

sub_fams <- paste0("plink -bfile ", "output/", f_name, 
                   " --keep output/", f_name,"_sample_fams.txt --out ", "output/",f_name, "_subbed",
                   " --make-bed")
print(sub_fams)


# This takes those subsetted families and then subsets based on two or more
# cases for my binary phenotype

affected_fids <- pheno2 %>% filter(IID %in% sub_set_list$V2) %>% group_by(FID) %>% 
  select(FID, V11) %>% summarise(affected = sum(V11)) %>% filter(affected >= 2)

selected_fams<-pheno2 %>% filter(FID %in% affected_fids$FID)

# This is for plink the subsetted pheno file
write.table(selected_fams, paste0("output/",f_name,
            "_phenotypes_structured_selected_fams.txt"), sep="\t", quote=F, 
            row.names = F, col.names = T)
# This is clean for gemma the subsetted pheno file of fams with >= 2 affected
write.table(selected_fams %>% select(V1, V11), paste0("output/",f_name,
            "_phenotypes_structured_one_selected_fams.txt"), sep="\t",
            col.names = F, row.names = F, quote = F)

# This subsets the relevant families to pull genotype data
write.table(selected_fams %>% select(FID, IID), 
            paste0("output/",f_name,"_sample_affected_fams.txt"), sep="\t",
            col.names = F, row.names = F, quote = F)

sub_set_affected <- fread(paste0("output/",f_name,
                                 "_phenotypes_structured_selected_fams.txt"))

print(paste0("Number of families after sampling >=2 cases per fam: ", 
             length(unique(selected_fams$FID))))

print(paste0("Number of individuals after sampling >=2 cases per fam: ", 
             length(selected_fams$FID)))

#selected_fams

sub_affected <- paste0("plink -bfile ", "output/", f_name, 
                   " --keep output/", f_name,"_sample_affected_fams.txt --out ", "output/",f_name, "_affected",
                   " --make-bed")
print(sub_affected)

# s_test <- c()
# 
# for(i in 1:1000){
# s_test[i]<-sample(2:6, size=1)
# }

# Bed file, fam file, gemma file, pheno (null)

# Create grm and run gwas

gemma_grm <- paste0("gemma -bfile output/",f_name,"_affected -gk 2 -o ",f_name,"_affected_grm")

print(gemma_grm)

gemma_assoc <- paste0("gemma -bfile output/",f_name,
                      "_affected -p output/",f_name,
       "_phenotypes_structured_one_selected_fams.txt -k output/",f_name,
       "_affected_grm.sXX.txt -lmm 4 -n 1 -o ", f_name)

print(gemma_assoc)

library(qqman)
asso_test<-paste0("output/",f_name,".assoc.txt")
gm <- fread(asso_test) %>% as_tibble()
my_snps <- fread("output/snp_eff.txt")

qq(gm$p_wald)
names(gm)
manhattan(gm, chr = "chr", bp = "ps", p = "p_wald", snp = "rs", highlight = my_snps$V1)


gm %>% filter(rs %in% my_snps$V1) %>% select(rs, af, beta, se) %>% 
  arrange(rs) %>% left_join(my_snps, by=c("rs"="V1")) %>% 
  rename("GivenEff"="V2") %>% arrange(beta) %>% print(n=Inf)


pheno <- fread(paste0("output/",f_name,"_phenotypes_structured_selected_fams.txt"))

length(unique(pheno$FID))

fams<-fread(paste0("output/",f_name,"_sample_affected_fams.txt"))

unique(fams$V1) %>% length()






