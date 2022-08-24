# USE: Run bootstrap analysis HCM/DCM paper style
# SCOPE: HCM/DCM 
# When: 2/16  
# Who: McNally Lab (LLP)
# NOTE: designed to function on the Lab's server (check library path and  so on)
#       Look for "AD HOC" that's where there's where a somewhat arbitrary choice is made
#        please don't assume it's the only place.
# Usually the opposite order of paths is better, but in this case not because
# the wrong version of dependencies would be otherwise loaded. 
# It seems like a bug of this machine installation
.libPaths( c( "/data2/R/x86_64-pc-linux-gnu-library/4.1",.libPaths()) )


library(tidyverse)
library(magrittr) # to make available sometools that tidyverse masks
library(tools) # mostly md5
library(readxl)

# Let's check env --commented is where we were when this was run/written
getRversion()
# [1] ‘4.1.1’
getwd()
# [1] "/data2/Bootstrap"
# read the preprocessed data -------------------------------------------------
var <- read_rds(paste0('Output/','var.rds'))
ancestry_data <- read_rds(paste0('Output/','ancestry_data.rds'))
# Note that the md5sum and time stamps of the files are in "0_process_raw_input.R"
# Run some checks
var %>% str()
nrow(var)
#[1] 3074320
#[1] 4314824
var %>% distinct(Chr)
# "chr1" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19"  "chr2" "chr20" "chr21" "chr22"  "chr3"  "chr4"  "chr5"  "chr6"  "chr7"  "chr8"  "chr9"  "chrX" 
var %>% summarize(across(c(Id,Chr,phen,Gene,AF,ENCODE_BL,RF), ~n_distinct(.x)))
# Id Chr phen  Gene AF ENCODE_BL RF
# 1 264  23    8 20372  2         1  1
# List of IDs, in this case all of them
ID <- var %>% distinct(Id) %>% .$Id 
length(ID)
# [1] 264
var <- droplevels(var)

# Set seed for repeatability.
set.seed(3)

#Test function to Create a random sample of Ids with replacement should be about .632 and I am not going to say why
sample(ID,size = length(ID), replace=T) %>% n_distinct() / length(ID) 

# plot the frequencies vs AF to see if anything looks weird.
var %>% ggplot()+geom_density(aes(x=gnomAD_AF_G, col=as.factor(AF)))


# Check MTCH2 --------------------------------------------------------------------
MTCH2_length <- 262*3
tmp <- var %>% filter(Gene=='MTCH2') %>% left_join(ancestry_data, by='Id') 
nrow(tmp)
# [1] 133
# [1] 543
tmp <- droplevels(tmp)
tmp %>% 
   group_by(Loc,AF, ANCESTRY,gnomAD_AF_G, gnomAD_AF_AFR_G, gnomAD_AF_NFE_G, gnomAD_AF_EAS_G) %>% 
  summarize(n_var = n()) %>% pivot_wider(names_from = ANCESTRY, values_from = n_var, values_fill =0)

# All the ones who have a triple of variants, zero are HOM
tmp %>% filter(between(Loc, 47660294,47660301)) %>% select(Id, Loc) %>% group_by(Id) %>% summarise(n= n()) %>% select(n) %>% table
# 0 
tmp %>% filter(between(Loc, 47660294,47660301)) %>% filter(AF != .5) %>% nrow
#0
# Plot density for MTCH2 and the rest
coeff <- 2.5
ggplot() + 
  geom_density(data = var, aes(x=gnomAD_AF_G)) + 
  #geom_density(data = var %>% filter(Gene == 'MTCH2'), aes(x=gnomAD_AF_G), adjust = 1 )+
  geom_histogram(data = var %>% filter(Gene == 'MTCH2'), aes(x=gnomAD_AF_G, y = ..count../sum(..count..)*coeff), fill = 'mediumpurple1', alpha = .4 )+
  xlab('gnomAD allele frequency')+
  scale_y_continuous(
    # Features of the first axis
    name = "density",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*coeff,name="Rescaled Counts")
  )+
  theme_bw()

# Create allele counts and number of chromosome variables
var <- var %>% mutate(AC=2*AF,NC=2)

# Add Phenotype information so we have a single and consistent df
var <- var %>% left_join(ancestry_data %>% transmute(Id=IID,gnomad_ancestry,Sex=GENETIC.SEX),by='Id')
#Correction for actual ploidy and weird calls on chrX for males AD HOC fix====================
var[var$Sex=='m' & var$Chr=='chrX', ]$AC <- 1
var[var$Sex=='m' & var$Chr=='chrX', ]$NC <- 1

create_phen <- function(phen_in){
  # function used to transform available phenotype into a phenotype useful for the analysis
  phen_out <- case_when(
    # In this case we are just returning one value because we are not comparing anything
    phen_in == 'HCM' ~ 'HCM',
    phen_in == 'DCM' ~ 'DCM',
    T ~ 'CM'
  )
  phen_out
}


# Check for Genes "duplication" and make AD HOC fix===========================
var %>% distinct(Gene,Chr) %>% group_by(Gene) %>% filter(n() > 1)
# Gene  Chr  
# 1 LSP1  chr11
# 2 LSP1  chr13
# 3 CKS1B chr1 
# 4 CKS1B chr5 
#var[var$Gene == 'LSP1',]$Chr <- 'chr11'
#var[var$Gene == 'CKS1B',]$Chr <- 'chr1'

# Create dataframes for bootstrap------------------------------------------------------------------
# FILTERING first frequency and phenotype
var.bs <- var %>% 
  filter( 
    between(gnomAD_AF_AFR_G, left = .25, right = .5) | 
    between(gnomAD_AF_NFE_G, left = .25, right = .5) | 
    between(gnomAD_AF_EAS_G, left = .25, right = .5) 
    ) %>% # Select frequency cutoff that matches any of the ancestries
  filter(phen %in% c('HCM','DCM')) %>% 
  transmute(Id,type.of.CM=create_phen(phen),gnomad_ancestry,Sex,Gene,Chr,Loc,Ref,Alt,AC,NC, gnomAD_AF_AFR_G,gnomAD_AF_NFE_G, gnomAD_AF_EAS_G)
# Check if there are any double entries, like when variants are reported both as high and moderate by snpEFF
nrow(var.bs)
#[1] 1207673 # .25-.5 only HCM and DCM
# OLD to remove 
# [1] 4314824 # all vars
# [1] 1391748 # .1-.25 modifier region worked
# [1] 1293278 # .1-.25 only HCM and DCM
# [1] 1892385 # .25-.5 only HCM and DCM
#[1] 1864629 
# some sanity checks
# Check how many near fixated in at least one pop we are left with
var.bs %>% filter(near(gnomAD_AF_NFE_G,1,.1) | near(gnomAD_AF_AFR_G,1,.1) | near(gnomAD_AF_EAS_G,1,.1) ) %>% nrow()
# [1] 6116 .1 .24
# [1] [1] 47042 .25-.5 when keeping only HCM/DCM 

# Check how many rare in at least one pop we are left with
var.bs %>% filter(near(gnomAD_AF_NFE_G,0,.01) | near(gnomAD_AF_AFR_G,0,.01) | near(gnomAD_AF_EAS_G,0,.01) ) %>% nrow()
#[1] 34026 .25-.5 only HCM & DCM

# [1] 163429
# [1] 153466 only HCM DCM

var.bs %>% group_by(Id,Chr,Loc,Ref,Alt) %>% filter(n() > 1) %>% nrow()
# [1] 0
# Plot density of variants by variant frequency by ancestry used
var.bs %>% ggplot() +
  geom_density(aes(gnomAD_AF_AFR_G), col='orange') +
  annotate(geom='text',x=.1, y=.3, col='orange', label='AFR')+
  geom_density(aes(gnomAD_AF_NFE_G),col='blue') + 
  annotate(geom='text',x=.3, y=.3, col='blue', label='NFE')+
  geom_density(aes(gnomAD_AF_EAS_G),col='red')+
  annotate(geom='text',x=.5, y=.3, col='red', label='EAS')+
  theme_bw()

ancestry_data.bs <- left_join(var.bs %>% group_by(Id) %>% summarize(n_var = n()) %>% ungroup(), ancestry_data %>% select(-IID), by=c('Id'))
# Plot the variant counts to find a 'reasonable NFE set, 
vline <- 5500
ancestry_data.bs %>% ggplot() + geom_histogram(aes(x=n_var, col = ANCESTRY), binwidth =100) + geom_vline(xintercept = vline, col='green')
# FILTER:: SubSample selection by ancestry =====================================
ancestry_data.bs <- ancestry_data.bs %>% 
#  filter( ANCESTRY == 'EUR' & n_var <= vline) # Keep only EUR removing the high var count EUR
#  filter( ANCESTRY == 'AFR') # Keep only EUR removing the high var count EUR
   filter(!(ANCESTRY == 'EUR' & n_var > vline)) %>%  # Remove the high variant count EUR
   filter(!(ANCESTRY == 'ASN')) %>%  # Remove Asians because we have too few
   filter(!(ANCESTRY == 'HIS')) # Remove HIS for whom we have no good reference AFs
  
nrow(ancestry_data.bs)
# [1] 142 keep only EUR, 30 for AFR only
# [1] D 172 C -ASN and with only HCM and DCM
# 181 C B - ASN
# 186 B remove HIS and high variant EUR
# 264 A All option A
# Plot nr of variant with genomic group coloring after sample filtering
ancestry_data.bs %>%  ggplot() + geom_histogram(aes(x=n_var, col = ANCESTRY), binwidth =50) #+ geom_vline(xintercept = vline, col='green')
# Look for outliers
ancestry_data.bs %>% filter(! near(n_var, mean(n_var), 5*10^3))
# FILTER Remove outliers ==============================================
ancestry_data.bs <- ancestry_data.bs %>% filter(near(n_var, mean(n_var), 5*10^3))
# Phenotype filtering 
var.bs <- var.bs %>% semi_join(ancestry_data.bs, by = 'Id')
nrow(ancestry_data.bs)
#[1] 262 A
#[1] 186 B
#[1] 172 D no EAS, only HCM + DCM
#[1] 142 E only EUR with less than vline vars, F 30 for AfR

# Plot nr variants by gnomAD probability by variant (so each variant appears only once)
# Normalization isn't correct like this and the plots aren't that comparable as a consequence.
var.bs %>% distinct(Chr, Loc,Ref,Alt,gnomAD_AF_AFR_G,gnomAD_AF_EAS_G,gnomAD_AF_NFE_G) %>% 
  mutate(gnomAD_AF_AFR_G = pmax(gnomAD_AF_AFR_G, 10^-4)) %>%
  mutate(gnomAD_AF_EAS_G = pmax(gnomAD_AF_EAS_G, 10^-4)) %>%
  mutate(gnomAD_AF_NFE_G = pmax(gnomAD_AF_NFE_G, 10^-4)) %>%
  ggplot() +
  geom_density(aes(gnomAD_AF_AFR_G), col='orange') +
  annotate(geom='text',x=.1, y=.3, col='orange', label='AFR')+
  geom_density(aes(gnomAD_AF_NFE_G),col='blue') + 
  annotate(geom='text',x=.3, y=.3, col='blue', label='NFE')+
  geom_density(aes(gnomAD_AF_EAS_G),col='red')+
  annotate(geom='text',x=.5, y=.3, col='red', label='EAS')+
  scale_y_continuous(trans='log')+
#  scale_x_continuous(trans='log')+
  theme_bw()

var.bs <- var.bs %>% semi_join(ancestry_data.bs, by = 'Id')
nrow(var.bs)
#[1] 811053 .25-.5 D
# Old
#[1] 881270  D  .1-.25 only EUR and AFR and only HCM and DCM
#[1] 2866848 C
#[1] 2943579 B  952189 for .1-.25. modifiers
#[1] 4278193 A
#[1] 666606 E
#[1] 144447 F 

#Utility dataframes
ID <- var.bs %>% distinct(Id) %>% .$Id 
length(ID)
# [1] 30  F
# [1] 142 E
# [1] 172 D
# [1] 262 A
# [1] 186 B
# [1] 181 C 
var.bs <- droplevels(var.bs)

chrs <- var.bs %>% distinct(Chr) %>% .$Chr
length(chrs)
# [1] 23
pheno.4 <- var.bs %>% distinct(Id, Sex, gnomad_ancestry, type.of.CM) %>% droplevels #simplifies phenotype 
dim(pheno.4)
#[1] 30  # F
#[1] 142 # E
#[1] 172   4 # D
#[1] 181   4 # C
#[1] 186   4 # B
# [1] 262   4 # A
geneByChr <- var.bs %>% distinct(Gene,Chr) #List of genes with the chromosome they are on
dim(geneByChr)
#[1] 5867    2
#[1] 5873    2 E
# [1] 9313    2 D [1] 6011    2 .25-.5
# [1] 22057     2 A
# [1] 21381     2 B 9404    2 modifier
# [1] 21244     2
# Check if there's any weird "duplication"
geneByChr[geneByChr$Gene %in% geneByChr[duplicated(geneByChr$Gene),]$Gene,]


computeEffectiveAF <- function(B.db,ID.B){
  #Function to compute the statistic for a dataframe of variants. it assumes that the bootstrap sample is "collapsed" as 
  # subjects by weight, where the weight is the nr of time a subject appears in the sample with replacement.
  tmp <- B.db %>% mutate(wAN=AC*w,wNC=NC*w) #Multiply the nr of alleles and chromosomes by the weight in the B sample
  # sum the weighted allele numbers by variant race and condition
  #NOTE: Assumes variants appear only once per subject, which is clearly true unless, for example
  #      there are SNPEFF High and Moderate variants that are the same variant and should have been filtered out
  tmp <- tmp %>% group_by(Chr,Loc,Ref,Alt, gnomAD_AF_AFR_G,gnomAD_AF_NFE_G, gnomAD_AF_EAS_G,Gene,type.of.CM,gnomad_ancestry) %>%
    summarise(tvAN=sum(wAN)) %>% ungroup 
  
  
  #Outer product of *all* chr in full db times all ids (with weight) in subsample
  idByChr <- merge(ID.B  %>% unique, data.frame(Chr=chrs))
  #      (to avoid missing chr because there are no variants for a specific subsample)
  #      If a chr has no variants in the original set, it would never be observed anyway
  #      because it can't appear in a bootstrap sample   
  #Contingency table of weighted nr of chromosomes in type.of.CM outer Race (wNC), 
  #considering male or female for chrX
  count.B <- idByChr %>% 
    group_by(type.of.CM,Sex,Chr,gnomad_ancestry) %>% 
    summarise(nSubj=sum(w)) %>% spread(Sex,nSubj,fill=0) %>% #Nr of males and females counting repeat samples
    mutate(wNC=f*2 +m*ifelse(Chr=='chrX',1,2)) %>% select(-m,-f) %>% ungroup # if chrX, count males as 1*weight
  #Merge it to the variant table so each variant for each HCM by race group has the total number of chromosomes or possible alleles
  tmp <- merge(count.B ,tmp)
  
  # Marginalize the count.B table, to obtain nr of chromosomes by condition, irrespective of ancestry 
  count.B <- count.B %>% group_by(type.of.CM,Chr) %>% summarise(totNC=1.0*sum(wNC)) %>% select(type.of.CM,Chr,totNC) %>% ungroup
  count.B <- merge(count.B,geneByChr,by='Chr') %>% select(-Chr) #nr of possible alleles by disease and gene
  
  rm(idByChr) #Clean up some memory in case we are tight
  
  #Produce the excess variants for this allele, given gnomAD, race and nr of subjects. Just sum the allele frequencies found
  ##NOTE: alternative approach that account for weighted averaging across races. The alleles are summer across variants and race as alleles 
  ##Compute difference in Allele counts between cohort and gnomAD, not differences in frequencies
  tmp <- tmp %>% mutate( excAN = 1.0*tvAN - wNC*
                           case_when(
                             gnomad_ancestry=='NFE' ~ gnomAD_AF_NFE_G,
                             gnomad_ancestry=='AFR' ~ gnomAD_AF_AFR_G,
                             gnomad_ancestry=='EAS' ~ gnomAD_AF_EAS_G
                           )
                  ) 
  ##Sum all the excess allele counts across race and variant alike
  tmp <- tmp %>% group_by(Gene,type.of.CM) %>% summarise( CexcAN=sum(excAN) ) %>% ungroup
  
  ## Divide by the total nr of chromosomes  
  tmp <- merge( tmp, count.B) 
  tmp <- tmp %>% mutate(EAF = CexcAN/totNC) %>% select(Gene,type.of.CM,EAF)
  return(tmp)
  
}

#Bootstrap loop
nBoot <- 5000 
#Create matrix with all combinations of genes and type.of.CM possible from the data
geneAF <- merge(geneByChr %>% select(Gene) %>% unique, pheno.4 %>% select(type.of.CM) %>% unique ) %>% arrange(Gene)
dim(geneAF)
# [1] 6011    2 D .25-.5 # 12022     2 for HCM vs DCM
# 22057     2     A
# [1] 21381     2 B 9404    2
#[1] 21244     2  C

gc()

#Bootstrap loop
options(dplyr.summarise.inform = FALSE) # So the screen isn't packed with garbage
for (i in seq(nBoot)){
  #Create weights for the various IDs depending upon resampling frequency
  #list of resampled IDs with how many times they were resampled
  ID.B <- data.frame(Id=sample(ID,size = length(ID), replace=T)) %>% group_by(Id) %>% summarize(w=n()) #
  #Add Race & Gender 
  ID.B <- ID.B %>% left_join(pheno.4 %>% select(Id,type.of.CM,gnomad_ancestry,Sex),by='Id')
  #Filter for Ids in bootstrap sample
  var.B <- var.bs %>% filter(Id %in% ID.B$Id) 
  #Add bootstrap weight to the genomic data 
  var.B <- var.B %>% left_join(ID.B %>% select(Id,w), by='Id') 
  B.res <- computeEffectiveAF(var.B,ID.B) # Compute statistic
  B.res %>% select(Gene,type.of.CM,EAF)
  geneAF <- left_join(geneAF,B.res,by = c("Gene", "type.of.CM"),suffix=c("",toString(i)))
  if ( i %% 100 == 0 ) {
    print(paste("Done with Sample",i))
    gc() # remove changes of crashing
  }
}

geneAF.CM <- geneAF
geneAF.HCM_DCM <- geneAF
# Save bootstrap raw data
saveRDS(geneAF.CM,'geneAF.CM.rds')
md5sum('geneAF.CM.rds')
# geneAF.CM.rds 
# "f93535d88f038218187e9dd4e0335db0" 
saveRDS(geneAF.HCM_DCM,'geneAF.HCM_DCM.rds')
md5sum('geneAF.HCM_DCM.rds')
# geneAF.HCM_DCM.rds 
# "3d31dee77c2374a6fe732b8f4d16d3de" 
saveRDS(geneAF,'geneAF.HCM_DCM.EUR.rds')
md5sum('geneAF.HCM_DCM.EUR.rds')
# geneAF.HCM_DCM.EUR.rds 
# "8cd7fd5ea2fc7255d019859ef52d24a0" 
saveRDS(geneAF,'geneAF.HCM_DCM.AFR.rds')
md5sum('geneAF.HCM_DCM.AFR.rds')
# geneAF.HCM_DCM.AFR.rds 
# "1a22511e38ac9c6c3b6d9ed4665d2b82" 

# tail size for the empirical distribution (p-value/2 divided by correction factor)
FDR    <- .01 # .20 is for 20%
ntests <- var.bs %>% select(Gene) %>% unique %>% nrow #nr of observed genes: we are not testing genes for which we have no data
# most stringent value is 1, one can go up until the ith CI contains 0. If the first n pass at the cut with i=1, then the next set
# can be set a i=n+1
i <- 1
tail <- (i*FDR/ntests)/2  # Maximal correction for FDR is i=1

# Single phenotype test
#Convert dataframe into a matrix of numbers
geneAF[,1:2] -> Cnames
tmp1 <- as.matrix(geneAF[,3:ncol(geneAF)])
n <- 1
tests <- data.frame(gene=character(),LB=numeric(),UB=numeric())
for (i in 1:(nrow(geneAF))){
  gene<-Cnames[i,1]
  CM1 <-Cnames[i,2]
  x1 <-tmp1[i,]
  x1[is.na(x1)] <- 0
  ICM1<-quantile(x1,c(c(tail,1-tail)))
  if ( gene == 'MTCH2') { print(paste('gene=',gene, 'LB=',ICM1[[1]],'UB=',ICM1[[2]]))}
  if ( ICM1[[1]] > 0 | ICM1[[2]] < 0){
   tests <- rbind( tests,
                   data.frame(gene=gene, LB=ICM1[[1]],UB=ICM1[[2]]))
    n <- n+1
  }
  
}
# Option D "gene= MTCH2 LB= 0.083027702647949 UB= 0.536675636552803" out of 184 genes




# For comparing groups  ==========================================================
FDR=.01 # .20 is for 20%
ntests <- 3* (var.bs %>% select(Gene) %>% unique %>% nrow) # nr of observed genes: running 3 tests for each gene
# most stringent value is 1, one can go up until the ith CI contains 0. If the first n pass at the cut with i=1, then the next set
# can be set a i=n+1
i <- 23
tail <- (i*FDR/ntests)/2  # Maximal correction for FDR is i=1

geneAF[,1:2] -> Cnames
tmp1 <- as.matrix(geneAF[,3:ncol(geneAF)])
n <- 1
#tests <- data.frame(gene=character(),LB=numeric(),UB=numeric())
for (i in 1:(nrow(geneAF)/2)){
  gene1<-Cnames[i*2-1,1]
  CM1 <-Cnames[i*2-1,2]
  gene2<-Cnames[i*2,1]
  CM2 <-Cnames[i*2,2]
  if(gene1 != gene2){ print(paste("ERROR:",gene2,gene1))}
  x1 <-tmp1[i*2-1,]
  x1[is.na(x1)] <- 0
  x2 <-tmp1[i*2,]
  x2[is.na(x2)] <- 0
  x12 <- x2-x1
  ICM1<-quantile(x1,c(c(tail,1-tail)))
  ICM2<-quantile(x2,c(c(tail,1-tail)))
  DIFF<-quantile(x1-x2,c(c(tail,1-tail)))
  if ( gene1 == 'MTCH2') {
    cat(paste(
      gene2," DELTA: [",sprintf("%6.4f",DIFF[[1]]),"-",sprintf("%6.4f",DIFF[[2]]),"],",
      CM2,": [",sprintf("%6.4f",ICM2[[1]]),"-",sprintf("%6.4f",ICM2[[2]]),"],",
      CM1,": [",sprintf("%6.4f",ICM1[[1]]),"-",sprintf("%6.4f",ICM1[[2]]),"]\n",
      sep=""))
  }
  if (DIFF[[1]] > 0 | DIFF[[2]] < 0 |  ( ICM1[[1]] > 0  & ICM2[[1]] > 0) | (ICM1[[2]] < 0 & ICM2[[2]] < 0)){
    #if (ICM2[[1]] > 0 | ICM2[[2]] < 0){
    cat(paste(
      gene2," DELTA: [",sprintf("%6.4f",DIFF[[1]]),"-",sprintf("%6.4f",DIFF[[2]]),"],",
      CM2,": [",sprintf("%6.4f",ICM2[[1]]),"-",sprintf("%6.4f",ICM2[[2]]),"],",
      CM1,": [",sprintf("%6.4f",ICM1[[1]]),"-",sprintf("%6.4f",ICM1[[2]]),"]\n",
      sep=""))
    n <- n + 1
  }
  
}



