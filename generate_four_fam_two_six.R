#Adds needed packages

library(dplyr)
library(magrittr)
library(stringr)
library(data.table)

#Cretaing segmented children

#Options function allows user to set and examine global options which affects the way R computes and displays results 
#Scripen is a penalty to be applied when deciding to print numeric values in fixed or exponential notation
#Positive values bias towards fixed and negative towards scientific notation: fixed notation will be preferred unless it is more than ‘scipen’ digits wider.
options(scipen=999)

# strrep function: Repeat the character strings in a character vector a given number of times
seq3 <- strrep("40 ", 200)
#print(paste0(seq3))
#trim whitespace from the end of the string
seq4 <- str_trim(seq3)
#print(paste0(seq4))

returns an absolute date-time value
start <- Sys.time()

#Gerenates a static row id seed to start at 414 each time
set.seed(414)
row_id = read.table("new_id_row.txt")$V1
# [1/3] Note: "new_id_row.txt" is essentially Joe's original .vec file but now is inverted, meaning the pairs are still randomly generated 
# [2/3] but now the pairs are flipped. In theory, this should not affect results since child creation pairs are still random and remain static. 
# [3/3] The script I used to do this can be found here 

# [1/2] Before proceeding, you will need to download U-Chicago's ms program for enerating samples underneutral models. Documentation and download can 
# [2/2] be found here: https://uchicago.app.box.com/s/l3e5uf13tikfjm7e1il1eujitlsjdx13
path_ms = "/Users/samanthafigeuredo/Downloads/self-replication/"
paste0("gcc -o ", path_ms,"ms ",path_ms,"ms.c ", path_ms,"streec.c ",
       path_ms,"rand1t.c -lm")


paste0("gcc -o ", path_ms,"ms ",path_ms,"ms.c ", path_ms,"streec.c ",
       path_ms,"rand1t.c -lm")


#Function from within software for msdir
source(paste0(path_ms,"readms.output.R"))
#
# Expect command line args at the end. 
args = commandArgs(trailingOnly = TRUE)
# Extract and cast as numeric from character


bash_batch_ind <- 1


batch_begin <- seq(from = 0, to = 9875, by = 125)

start_sim <- batch_begin[bash_batch_ind] + 1

end_sim <- batch_begin[bash_batch_ind] + 125

print(paste0("This is batch ", bash_batch_ind))



k<-0


seeds <- list()

start_time <- Sys.time()

for (k in start_sim:end_sim) {
  #Randomness within this function some seed
  row_id_gc1<-unlist(c(replicate(20,as.list(sample(1:0, 2)))))
  row_id_gc2<-unlist(c(replicate(20,as.list(sample(1:0, 2)))))
  
  #Added so siblings are not related
  
  #######################
  haps1 <- sample(0:1, 1)
  haps2 <- sample(0:1, 1)
  haps3 <- sample(0:1, 1)
  haps4 <- sample(0:1, 1)
  haps5 <- sample(0:1, 1)
  haps6 <- sample(0:1, 1)
  ########################

  try(system(paste0("gcc -o ", path_ms,"ms ",path_ms,"ms.c ", path_ms,"streec.c ",
                  path_ms,"rand1t.c -lm"), intern=FALSE, wait=TRUE))
  
  temp_out<-try(system(paste0(path_ms,"ms 8000 1 -t 5.0 -r 0.1 10000 -I 200 ", 
		 seq4," 4000 -s 100"), intern=TRUE))

	seeds[[k+1]] <- temp_out[2]


  #Reads and interprets generate output
  msout_d <- read.ms.output(temp_out)
  
  
  #Convert to table and filter based upon row_id
  pop <- msout_d$gametes[[1]] %>% as.data.frame() %>% slice(1:8000) %>% 
    as_tibble()
  
  pop$subpop = rep(1:200, times = 1, each = 40)
    
  pop$id = rep(1:4000, times = 1, each = 2)
    
  pop$subpop_id = rep(1:20, times = 200, each = 2)
    
  pop$row_id = row_id
    
  pop <- pop %>% select(id, subpop, subpop_id, row_id, 1:100)
  
  popl <- split.data.frame(pop, pop$subpop)



  ##############################################################################
  
  sub1 <- 0
  sub2 <- 0
  sub3 <- 0
  sub4 <- 0
  sub5 <- 0
  sub6 <- 0
  i <- 0
  # See if removing this works
  #set.seed(414)
  #Maybe add this outside of the loop
 
 
  # Fix row ids for sibs and grandkids
  
  # Vary kids by some number set 
  # sample some number 2-6 and assign that to founders
  # Associate an FID with the children
  
  for (i in 1:length(popl)) {
  
    first_gen <- NULL
    second_gen <- NULL
    
    sub1 <- popl[[i]] %>% filter(row_id == haps1)
    
    sub1$subpop_id <- rep(21:30, times = 1, each = 2)
    
    sub2 <- popl[[i]] %>% filter(row_id == haps2)
    
    sub2$subpop_id <- rep(31:40, times = 1, each = 2)
    
    sub3 <- popl[[i]] %>% filter(row_id == haps3)
    
    sub3$subpop_id <- rep(41:50, times = 1, each = 2)
    
    sub4 <- popl[[i]] %>% filter(row_id == haps4)
    
    sub4$subpop_id <- rep(51:60, times = 1, each = 2)
    
    sub5 <- popl[[i]] %>% filter(row_id == haps5)
    
    sub5$subpop_id <- rep(61:70, times = 1, each = 2)
    
    sub6 <- popl[[i]] %>% filter(row_id == haps6)
    
    sub6$subpop_id <- rep(71:80, times = 1, each = 2)
    

    first_gen <- bind_rows(sub1, sub2, sub3, sub4, sub5, sub6)
    
    # sub3 <- first_gen %>% filter(subpop_id >= 21) %>% mutate(row_id = row_id_gc1) %>% 
    #   filter(row_id==haps3)
    # 
    # sub3$subpop_id <- rep(41:50, times = 1, each = 2)
    # 
    # sub4 <- first_gen %>% filter(subpop_id >= 21) %>% mutate(row_id = row_id_gc2) %>% 
    #   filter(row_id==haps4)
    # 
    # sub4$subpop_id <- rep(51:60, times = 1, each = 2)
    # 
    # second_gen <- bind_rows(sub3, sub4)
    
    popl[[i]] <- bind_rows(list(popl[[i]], first_gen))
    
  }
  
  
  
 
  
  
  
  
  
  
  
  ##############################################################################
  #Need to create columns and row names and make sure it works
  
 
  geno_l <- list()
  t1 <- 0
  j <- 0
  
  for (j in 1:length(popl)) {
    t1 <- popl[[j]] %>% select(subpop, subpop_id, V1:V100) %>%
      
      group_by(subpop, subpop_id) %>%
      
      summarise_all(list(sum)) %>%
      
      ungroup() %>%
      
      arrange((subpop_id))
    
    geno_l[[j]] <- t1
  }
  
  ##############################################################################
  
  genotype_pop <- bind_rows(geno_l)
  
  genotype_pop$id = paste0("id_",genotype_pop$subpop_id, "_", genotype_pop$subpop) 

  genotype_pop <- genotype_pop %>% select(id, V1:V100)

	options(scipen=999)
  
  geno.colnames <- paste0("rs", as.character((k - 1) * 100 + 1:100))
    
  df<-genotype_pop %>% as.data.frame(stringsAsFactors=F) %>% t() %>% as_tibble()
  
  names(df) <- df[1, ]
  
  df <- df[-1, ]
  
  df$rsid <- geno.colnames
  
  #These are the actual positions
  df1 <- df %>% mutate(r1=rep(1,100), r3=rep(0,100), num=(k - 1) * 100 + 1:100) %>% 
    select(r1,rsid,r3,num,id_1_1:id_80_200)
  
  df1[5:length(df1)]<-sapply(df1[5:length(df1)], function(x) gsub("0", "A A", x))
  df1[5:length(df1)]<-sapply(df1[5:length(df1)], function(x) gsub("1", "A C", x))
  df1[5:length(df1)]<-sapply(df1[5:length(df1)], function(x) gsub("2", "C C", x))
  
  fwrite(df1, file = paste0("output/", "out_", k, ".txt"), sep="\t", 
         quote=FALSE, col.names = FALSE)
}


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


trans.tfam <- genotype_pop %>% mutate(family_id = fid_all, individual_id = id, 
                                      paternal_id = rep(0, length(df1)-4), 
                                      maternal_id=rep(0,length(df1)-4), 
                                      sex=rep(1, length(df1)-4), 
                                      phenotype=rep(1,length(df1)-4)) %>%
              select(family_id, individual_id, paternal_id, maternal_id, sex, 
                     phenotype)

fwrite(trans.tfam, file="output/trans.tfam", sep="\t", quote=FALSE, 
       col.names = FALSE)



seeds<-unlist(seeds) %>% as.data.frame()

fwrite(seeds, file=paste0("output/",bash_batch_ind,"_seeds.txt"), 
       sep="\t")



end_time = Sys.time()


print(end_time-start_time)
