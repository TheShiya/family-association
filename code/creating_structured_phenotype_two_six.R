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
eff <- c(rep(log(1.5), 60), rep(log(.5), 40))

# log(1.5), log(0.5)



print("Effects created")
print(eff)
write.table(eff, file = "output/eff.txt", sep = "\t", col.names = F, 
            row.names = F,
            quote = F)

# Call to take trans file and make into plink
bed_call <- paste0("plink -tfile output/trans --make-bed --out ", "output/",f_name)
print(bed_call)
try(system(bed_call, intern = F))

# Call to get mafs for each allele
maf_call <- paste0("plink -bfile ", "output/", f_name," --freq --out ", "output/",f_name)
print(maf_call)
try(system(maf_call, intern = F))

# reads maf, selects 100 snps above .4 frequency randomly and writes it out
maf <- fread(paste0("output/",f_name, ".frq")) %>% arrange(desc(MAF)) %>% as_tibble()
set.seed(414)
###### This will need be MAF .05
snps_af <- maf %>% filter(MAF >= .090 & MAF <= .102) %>% sample_n(snps_n) %>% 
  select(SNP, MAF) %>% arrange(SNP)
snps <- snps_af$SNP
print("Snps chosen from plink file")
print(snps)
write.table(snps, file = "output/my_snps.txt", sep = "\t", col.names = F, 
            row.names = F,
            quote = F)


# LD scores - find null that are not in LD with effect snps lower than .1
ld_call<-paste0("plink --bfile output/",f_name,
                " --ld-snp-list output/my_snps.txt --out output/", 
                f_name, " --r2 --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0")
print(ld_call)

ld<-fread(paste0("output/", f_name,".ld")) %>% as_tibble()

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
sub_call <- paste0("plink -bfile ", "output/", f_name, 
                   " --extract output/my_snps.txt --recodeAD --out ", "output/",f_name)
sub_call
try(system(sub_call, intern = F))
# Reads the plink file and arranges by FID for the genetypes

rs <- fread(paste0("output/", f_name, ".raw")) %>% arrange(FID) %>% as_tibble()

pc<-c("P1", "P2", "C1", "C2", "C3", "C4", "C5","C6")
pc_order<-rep(pc, length(unique(rs$FID)))
pc_correct <- c("P1", "C1", "C2", "C3", "C4", "C5", "C6", "P2")

# Correct order for families
rs_ordered<-rs %>% mutate(pc_rel = pc_order) %>% 
  group_by(FID) %>% slice(match(pc_correct, pc_rel)) %>% select(-pc_rel)


# Take additive model takes specifically the additive subsetted genotypes and 
# writes them out

# This is for additive model
rs1 <- rs_ordered %>% select(names(rs)[!str_detect(names(rs), "HET")]) %>% 
  ungroup(FID) %>% 
  select(-(FID:PHENOTYPE))

# This creates dominant model
rs2<-sapply(rs1, function(x) ifelse(x==2,1,x)) %>% as_tibble()

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
# This fixes previous bug where snp effects were out pf order when applied to
# genos
names(rs1) <- substr(names(rs1),1,nchar(names(rs1))-2)
setcolorder(rs1, as.character(snp_eff$V1))
##################################################
# Double check the names are okay
stopifnot(snp_eff$V1 == names(rs1))

snps_eff_mat <- as.matrix(rs1) %*% diag(snp_eff$V2)

# Keep zeroes 0-->0, 1-->1 2-->1
#snp_eff_mat_dom <-#

r1_snp_eff <- rs %>% select(IID) %>% bind_cols(rowSums(snps_eff_mat) %>% 
                                                 as_tibble())
print("These are the snp add effects per individual")
print(r1_snp_eff)

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


# four_fam_mat <- matrix(
#   c(
#     0.5, 0, 0.25, 0.25,
#     0, 0.5, 0.25, 0.25,
#     0.25, 0.25, 0.5, 0.25,
#     0.25, 0.25, 0.25, 0.5
#   ),
#   byrow = TRUE,
#   ncol = 4
# )
# 
# five_fam_mat <- matrix(
#   c(
#     .5, 0.25, 0.25, 0.25, 0,
#     0.25, 0.5, 0.25, 0.25, 0.25,
#     0.25, 0.25, 0.5, 0.25, 0.25,
#     0.25, 0.25, 0.25, 0.5, 0.25,
#     0, 0.25, 0.25, 0.25, 0.25, 0.5
#   ),
#   byrow = TRUE,
#   ncol = 5
# )
# 
# six_fam_mat <- matrix(
#   c(
#     0.5, 0.25, 0.25, 0.25, 0.25, 0,
#     0.25, 0.5, 0.25, 0.25, 0.25, 0.25,
#     0.25, 0.25, 0.5, 0.25, 0.25, 0.25,
#     0.25, 0.25, 0.25, 0.5, 0.25, 0.25,
#     0.25, 0.25, 0.25, 0.25, 0.5, 0.25,
#     0, 0.25, 0.25, 0.25, 0.25, 0.5,
#   ),
#   byrow = TRUE,
#   ncol = 6
# )
# 
# seven_fam_mat <- matrix(
#   c(
#     0.5, 0.25, 0.25, 0.25, 0.25, 0.25, 0,
#     0.25, 0.5, 0.25, 0.25, 0.25, 0.25, 0.25,
#     0.25, 0.25, 0.5, 0.25, 0.25, 0.25, 0.25,
#     0.25, 0.25, 0.25, 0.5, 0.25, 0.25, 0.25,
#     0.25, 0.25, 0.25, 0.25, 0.5, 0.25, 0.25,
#     0.25, 0.25, 0.25, 0.25, 0.25, 0.5, 0.25,
#     0, 0.25, 0.25, 0.25, 0.25, .25, 0.5
#   ),
#   byrow = TRUE,
#   ncol = 7
# )

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
beta_fam <- NULL
beta_eff <- NULL

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
  
  epsilon <- rep(rnorm(n = 8000)) #individual error
  
  zeta <- rep(rnorm(200, sd=.25), each = 40) #sub population effect
  
  prev <- .1 # Prevalence
  
  alpha <- log(prev/(1-prev))
  
  
  beta_tot <- (beta_fam + rowSums(snps_eff_mat)) # Snps plus family effect
  
  
  trait_l_bin[[j]] <- exp(alpha + beta_tot) / (1 + exp(alpha + beta_tot))
  
  # these should not be data.frames make sure they are vectors
  trait_l[[j]] <- beta_tot + epsilon + zeta*0
  
  beta_fam <- NULL
}



samp_func <- function(X) (sample(0:1, size=1, prob=c(1-X,X)))


traits_bin<-apply(trait_l_bin %>% bind_cols(), MARGIN=c(1,2), samp_func) %>% 
  as_tibble()
names(traits_bin) <- "V11"

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
#rename("IID"="id","FID"="family_id", "PID"="pop_id") -> pheno

# mutate(ind_error=epsilon) %>%
# arrange(family_id) %>% mutate(fam_effect=fam_eff_tot) %>%
# arrange(pop_id) %>% mutate(pop_effect=rep(Z, each=40)) %>%
# mutate(trait=(pop_effect+fam_effect+ind_error))


# Reread in genotype data and look at IDs



# lm_runs <- list()
# 
# for (i in names(rs1)){
#   lm_runs[i]<-lm(pheno2$V1 ~ rs1[[i]])$coef[[2]]
# }



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






