# USE: Prepare the data to Run bootstrap analysis HCM/DCM paper style
# SCOPE: HCM/DCM 
# When: 2/09  
# Who: McNally Lab (LLP)
# NOTE: designed to function on the Lab's server (check library path and  so on)
#.      script filters out variables that failed gnomAD in some way, see below
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

#This is the list of variants that failed one or more filters from gnomAD
# Produce by SDK using a python script
file <- "/data/Genomics_tools/Databases/gnomAD/download/vcf/genomes/all_non_PASS_variants_wSegdupLCRDecoy.txt"
# there should be about 24 mil variants, produced by awk selecting for the filter field being RF
#file <- "/data/Genomics_tools/Databases/gnomAD/download/vcf/genomes/RF_variants.txt"
# "b2c77de7f2761fa469a9620ce30fe4c7"
# [1] "2022-02-07 21:10:56 CST"
md5sum(file)
#/data/Genomics_tools/Databases/gnomAD/download/vcf/genomes/all_non_PASS_variants_wSegdupLCRDecoy.txt 
#"cb90c338174f275bd1fd8cb774f29b4a"
file.info(file)$mtime
#[1] "2022-03-01 21:46:07 CST"
fail_pass_list_gnomAD <- read.csv(file,header = F,sep = '\t')
names(fail_pass_list_gnomAD) <- c('Chr', 'Loc', 'RS', 'Ref', 'Alt','QUAL','RF','segdup_lcr' )
fail_pass_list_gnomAD <- fail_pass_list_gnomAD %>% mutate(Chr = paste0('chr',Chr)) 
dim(fail_pass_list_gnomAD)
# [1] 58621369        8
# [1] 24025998        7
fail_pass_list_gnomAD %>% group_by(RF) %>% summarize(n())
# A tibble: 8 × 2
# RF                        `n()`
# 1 AC0                     1935314
# 2 AC0;InbreedingCoeff        1904
# 3 AC0;InbreedingCoeff;RF    10595
# 4 AC0;RF                  5695542
# 5 InbreedingCoeff          116851
# 6 InbreedingCoeff;RF       167427
# 7 PASS                   26667738
# 8 RF                     24025998
fail_pass_list_gnomAD %>% filter(RF == 'PASS') %>% group_by(segdup_lcr) %>% summarize(n())
# All variants have problems and should be removed from the variant list

# Read the variant_list -----------------------------------------------------------
#This is the dataframe made from VariantAnalysisBySubject.pl. 
#The only modification I made to it is that it didn't filter by any gene list, only removed subjects
file <- "/data2/reannotate_Cardiomyopathy_Genomes/Cohort/SS600_Genomes/all_genes.remove_excluded.dataframe"
md5sum(file)
#/data2/reannotate_Cardiomyopathy_Genomes/Cohort/SS600_Genomes/all_genes.remove_excluded.dataframe 
# "1beabfd3ec2e38a058ca5bdee29c9ae0" 
file.info(file)$mtime
# [1] "2022-02-15 10:48:15 CST"
tmp <- read.csv(file,header = T,sep = '\t')
dim(tmp)
# [1] 4723184      23
# Find variants that are not flagged for RF=====================================
var <- tmp %>% anti_join(fail_pass_list_gnomAD, by=c('Chr','Loc','Ref', 'Alt')) 
var$RF = F
dim(var)
# [1] 3212310      24
# [1] 4483417      24
# Find the flagged variants
tmp.1 <- tmp %>%  semi_join(fail_pass_list_gnomAD, by=c('Chr','Loc','Ref', 'Alt')) 
dim(tmp.1)
# [1] 1510874      23
# [1] 239767     23
# Check that the numbers add
nrow(tmp.1) + nrow(var) == nrow(tmp)
# T
tmp.1$RF = T
gc()

var <- rbind(var,tmp.1)
var <- var %>% arrange(Id,Chr,Loc,Ref,Alt)
tmp <- tmp %>% arrange(Id,Chr,Loc,Ref,Alt)

nrow(var) == nrow(tmp)
#T
all.equal(var %>% select(-RF), tmp)
# T 
rm(tmp, tmp.1,fail_pass_list_gnomAD)
# Clean up some memory
gc() 

# Eliminate duplicated variants that arrived as both H and M effects ======================
# Keep 'H'
var %>% group_by(Id, Chr,Loc,Ref,Alt) %>% filter(n() > 1) %>% nrow()
# [1] 34606
tmp <- var %>% group_by(Id, Chr,Loc,Ref,Alt) %>% filter(n() > 1 & Eff == 'M') 
nrow(tmp)
#[1] 17303
var <- var %>% anti_join(tmp)
nrow(var)
#[1] 4705881
rm(tmp)
gc()

# turns out that MTCH2 lost on average more variants that the rest
tbl <- var %>% transmute(Gene == 'MTCH2',RF) %>% table() 
tbl %>% summary()
tbl %>% prop.table(1)
# RF
# Gene == "MTCH2"     FALSE      TRUE
# FALSE 0.6797129 0.3202871
# TRUE  0.1619976 0.8380024
# nr of gene
var%>% select(Id,Gene) %>%  summarise(across(.cols=everything(),~n_distinct(.))) # Id  Gene
# 1 264 22663

# Read the demographics data -----------------------------------------------------------
#Here is the HCM/DCM demographics table (this does not include plink generated ancestry)
file <- '/data2/reannotate_Cardiomyopathy_Genomes/Cohort/HCM_DCM_demographics_table.csv'
md5sum(file)
#/data2/reannotate_Cardiomyopathy_Genomes/Cohort/HCM_DCM_demographics_table.csv 
#"50779be70975af90b453d0bd015313e8" 
file.info(file)$mtime
#[1] "2020-05-15 18:18:48 CDT"
subj_data <- read.csv(file,header = T,sep = ',')
subj_data  %>%  summarise(across(.cols=everything(), ~n_distinct(.x))) 
# ID patID Race.Genetic Sex HOM.HET.score
# 1 317   317           16   2           317
# Fix ID, first checking for consistency
subj_data %>% transmute(pre = stringr::str_sub(ID, start=1, end=5)) %>% distinct()
# pre
# 1 SS600
subj_data %>% mutate(Id = stringr::str_sub(ID, start=6)) %>% distinct(Id) %>%nrow()
# 317 
subj_data <-subj_data %>% mutate(Id = stringr::str_sub(ID, start=6)) 
# Find subjects in variant file but not in subject files
missing_demo <- unique(var$Id)[!unique(var$Id) %in% stringr::str_sub(subj_data$ID,start = 6)]
length(missing_demo)
# 71
# Cluster Race data (seems to be more self-assigned race than ancestry)
subj_data <- subj_data %>% mutate(
  Race_or_Genetic = case_when(
    Race.Genetic %in% c('Caucasian','CAUCASIAN','Caucasian (French descent)','EA') ~ 'EA',
    Race.Genetic %in% c('African American','AA')                                   ~ 'AA',
    Race.Genetic %in% c('Hispanic','HISPANIC','Caucasian/Hispanic',
                        'Hispanic(more than one race)','HA')                       ~ 'HA',
    T                                                                              ~ 'OA')
)
    
#Ancestry clustering from running plink on all of our genomes (Including NUgene and SDY).
#This doesn't have every person because we've gotten new genomes since running this, 
#and we did not include genomes run on a machine other than xten.
file <- '/data2/ancestry_run/HCM_DCM_full/HCM_DCM_race_class_coding_1_2_wSex.csv'
md5sum(file)
#/data2/ancestry_run/HCM_DCM_full/HCM_DCM_race_class_coding_1_2_wSex.csv 
#"ccb7b6f1d16c703151c82f76073ed307" 
file.info(file)$mtime
# [1] "2022-02-15 17:04:59 CST"
ancestry_data <- read.csv(file,header = T,sep = ',',colClasses = c('factor','numeric','numeric','integer','factor','numeric','factor'))
ancestry_data <- ancestry_data %>% mutate(Id = IID) 
ancestry_data %>% nrow()
# 411
# We have everyone
unique(var$Id)[!unique(var$Id) %in% ancestry_data$Id]
#integer(0)
ancestry_data %>% summarise(across(.cols=everything(),~n_distinct(.x)))
#IID  C1  C2 CLUST ANCESTRY HOM.HET.SCORE GENETIC.SEX  Id
#1 411 410 410     4        4           411           2 411
ancestry_data %>% summarise(across(.cols=everything(),~sum(is.na(.x))))
# IID C1 C2 CLUST ANCESTRY HOM.HET.SCORE GENETIC.SEX Id
# 1   0  0  0     0        0             0           0  0
tmp <- ancestry_data %>% inner_join(subj_data, by = 'Id') 
# Check if the sex variable seems to work, particularly focusing on the values around 1. 
# We're good
tmp %>% ggplot() + geom_histogram(aes(x=HOM.HET.score, col=Sex), binwidth = .5)
tmp %>% filter(near(HOM.HET.score, 1, .25))  %>% select(HOM.HET.score,GENETIC.SEX,Sex, ID)
# Compare the ancestry from one file with the other, of looks fine since one doesn't have ASN and it's cratered
tmp %>% select(Race_or_Genetic, ANCESTRY) %>% table()
# ANCESTRY
# Race_or_Genetic AFR ASN EUR HIS
# AA  47   0   0   3
# EA   0   1 163   5
# HA   0   8   2  20
# OA  13   4  41  10
# Visually compare the two sets of calls, they are reasonably concordant, yellow dots are subjects with uncertain sex estimation
tmp %>% ggplot() + geom_point(aes(x=C1,y=C2,col= ANCESTRY), size = 4)+
  geom_point(aes(x=C1,y=C2,shape= Race_or_Genetic), size = 2) +
  geom_point(data = tmp %>% filter(near(HOM.HET.SCORE, 1, .1)), aes(x=C1,y=C2), size = 1, col='yellow')
rm(tmp)
# We have AFR, NFE and EAS in gnomAD, so we collapse the ancestral data in 3 categories
ancestry_data <- ancestry_data %>% mutate(gnomad_ancestry = case_when(
      ANCESTRY %in% c('EUR','HIS') ~ 'NFE',
      ANCESTRY %in% c('AFR') ~ 'AFR',
      ANCESTRY %in% c('ASN') ~ 'EAS',
      T ~ 'very wrong'
        )
)
ancestry_data %>% group_by(gnomad_ancestry) %>% summarize(n())
# gnomad_ancestry `n()`
# 1 AFR                64
# 2 EAS                15
# 3 NFE               332
tbl <- ancestry_data %>% select(gnomad_ancestry, GENETIC.SEX) %>% table(., useNA = 'ifany')
summary(tbl)
# Number of cases in table: 411 
# Number of factors: 2 
# Test for independence of all factors:
#   Chisq = 0.9297, df = 2, p-value = 0.6282
rm(tbl)
# Check if there are other uncertaint sex scores, same as above, we are good
ancestry_data %>% filter(near(HOM.HET.SCORE, 1, .1))  
# Split out failed RF, process the variants =====================================
nrow(var)
# [1] 4705881
var_rf <- var %>% filter(RF)
nrow(var_rf)
# [1] 1507658
var <- var %>% filter(!RF)
nrow(var)
# [1] 3198223

# Check missing data to filter out garbage ---------------------------------------------
var %>% summarise(across(.cols=everything(),~sum(is.na(.x))))
# Id Eff phen Gene Chr Loc Ref Alt PP2 AF GERP GTEx  CLNSIG CLNSIGCONF DCM_FREQ HCM_FREQ NUGENOMES_FREQ gnomAD_AF_G gnomAD_AF_AFR_G gnomAD_AF_NFE_G
# 1  0   0    0    0   0   0   0   0   0  0    0    0 2702120    3188828        0        0              0       36083           36083           36083
# gnomAD_AF_EAS_G gnomAD_COV_G ENCODE_BL RF
# 1           36083         9904   3196133  0
# Check the variants that have <NA> for depth in gnomAD and less than 8, eliminate brutally ====================
var_no_depth <- var %>% filter(is.na(gnomAD_COV_G) | gnomAD_COV_G < 8) 
var <- var %>% filter( !(is.na(gnomAD_COV_G) | gnomAD_COV_G < 8) )
nrow(var)
#[1] 3161964

# No missing data in gnomAD
var %>% summarise(across(.cols=everything(),~sum(is.na(.x))))
# Id Eff phen Gene Chr Loc Ref Alt PP2 AF GERP GTEx  CLNSIG CLNSIGCONF DCM_FREQ HCM_FREQ NUGENOMES_FREQ gnomAD_AF_G gnomAD_AF_AFR_G gnomAD_AF_NFE_G
# 1  0   0    0    0   0   0   0   0   0  0    0    0 2666021    3152569        0        0              0           0               0               0
# gnomAD_AF_EAS_G gnomAD_COV_G ENCODE_BL RF
# 1               0            0   3160927  0

# Identify compound heterozygous and fix ============================================
tmp <- var %>% filter(if_any(starts_with('gnomAD') | c('Alt','AF'), ~ stringr::str_detect(.x, pattern = ","))) 
nrow(tmp)
# [1] 13587
# Select the complement of var relative to tmp (the non comp het)
var <- var %>% anti_join(tmp)
# Split variants and add them back to var
tmp_start <- tmp %>% mutate(across(starts_with('gnomAD') | c(Alt,AF), ~stringr::str_extract(.x,pattern='^[^,]*'))) 
tmp_end <- tmp %>% mutate(across(starts_with('gnomAD') | c(Alt,AF), ~stringr::str_extract(.x,pattern='[^,]*$'))) 
var <- rbind(var, tmp_start, tmp_end)
rm(tmp, tmp_start, tmp_end)

var <- var %>% mutate(across(.cols=starts_with('gnomAD') | c(AF), ~as.numeric(.x))) 
dim(var)
# [1] 3175551      24
# No NA left, we are good to go.
# var %>% summarise(across(.cols=everything(),~sum(is.na(.x))))
# Id Eff phen Gene Chr Loc Ref Alt PP2 AF GERP GTEx  CLNSIG CLNSIGCONF DCM_FREQ HCM_FREQ NUGENOMES_FREQ gnomAD_AF_G gnomAD_AF_AFR_G gnomAD_AF_NFE_G
# 1  0   0    0    0   0   0   0   0   0  0    0    0 2679608    3166156        0        0              0           0               0               0
# gnomAD_AF_EAS_G gnomAD_COV_G ENCODE_BL RF
# 1               0            0   3174488  0
# Encode blacklist =====================================
var %>% summarize(across(c(Gene, Eff,ENCODE_BL,RF), ~n_distinct(.x)))
#Gene Eff ENCODE_BL RF
#1 20385   2         5  1
var %>% select(ENCODE_BL) %>% table(.,useNA = 'ifany')
#(GAATG)n     centromeric_repeat Low_mappability_island       Satellite_repeat                   <NA> 
#  674                     97                    257                     35                3174488 
# Remove all non NA and remove the variable
var <- var %>% filter(is.na(ENCODE_BL))
var <- var %>% mutate(Id = as.factor(Id))
nrow(var)
#[1] 3174488

# Run some stats over the vars ====================================
var %>% summarize(across(c(Gene, Id), ~n_distinct(.x)))
var %>% group_by(Id) %>% summarize(n=n()) %>% summarize(min=min(n), median=median(n),mean=mean(n),max=max(n))
#min median   mean   max
#8115 11292. 12025. 21274
# plot nr of variants by race
var %>% group_by(Id) %>% summarize(n_var=n()) %>% left_join(ancestry_data, by=c('Id' = 'IID')) %>% 
  ggplot()+geom_density(aes(x=n_var, col=ANCESTRY))
# Just plot the SNPs
var %>% filter(str_length(Ref) == 1 & str_length(Alt) == 1) %>% group_by(Id) %>% summarize(n_var=n()) %>% left_join(ancestry_data, by=c('Id' = 'IID')) %>% 
  ggplot()+geom_density(aes(x=n_var, col=ANCESTRY))
# Check substitutions
var %>% filter(str_length(Ref) == 1 & str_length(Alt) == 1) %>% select(Ref,Alt) %>% table()
# Alt
# Ref      A      C      G      T
# A      0 110405 473225  69903
# C 117071      0 174065 503853
# G 546109 177726      0 112725
# T  69365 480800 108163      0
# Count indels and SNPs
var %>% transmute(SNP = str_length(Ref) == 1 & str_length(Alt) == 1) %>% table() 
# FALSE    TRUE 
# 231078 2943410 
# Check mismatched frequencies
var  %>% filter(near(gnomAD_AF_G, 0, .0001)) %>% transmute(NU_FREQ = as.numeric(stringr::str_extract(NUGENOMES_FREQ,pattern='^[^\\/]*'))) %>% ggplot()  + geom_density(aes(NU_FREQ), adjust=.05, fill='purple') +
   geom_vline(xintercept = 50)
#
var %>% filter( as.numeric(stringr::str_extract(NUGENOMES_FREQ,pattern='^[^\\/]*')) < 5 ) %>% ggplot() + geom_density(aes(gnomAD_AF_G), fill='purple') 
nrow(var)
# [1] 3174488
# Eliminate mismatched variants, very likely sequencing artifacts
var <- var %>% filter( 
    !(
      as.numeric(stringr::str_extract(NUGENOMES_FREQ,pattern='^[^\\/]*')) > 50 &
#      as.numeric(stringr::str_extract(DCM_FREQ,pattern='^[^\\/]*')) > 30 &
#        as.numeric(stringr::str_extract(HCM_FREQ,pattern='^[^\\/]*')) > 30 &
       near(gnomAD_AF_G, 0, .0001)
    )
  ) 
nrow(var)
#[1] 3074320
var %>% transmute(SNP = str_length(Ref) == 1 & str_length(Alt) == 1) %>% table() 
# .
# FALSE    TRUE 
# 206983 2867337 

# Plot PCs of ancestry with self reported race when available and unusual HOM/HET ratios for the subjects that we have in
# the variant file
tmp <- ancestry_data %>% filter( IID %in% (var %>% distinct(Id))$Id) %>% left_join(subj_data %>% transmute(IID = Id, Race_or_Genetic), by='IID')
ggplot() + geom_point(data = tmp, aes(x=C1,y=C2,col = ANCESTRY), size = 4) +
  geom_point(data=subset(tmp, !is.na(Race_or_Genetic)), aes(x=C1,y=C2,shape= Race_or_Genetic), size = 2) +
  geom_point(data = subset(tmp, near(HOM.HET.SCORE, 1, .1)), aes(x=C1,y=C2), size = 1, col='yellow')


# Write out the processed files
dir.create('Output')
file <- paste0('Output/','var.rds')
write_rds(var,file)
md5sum(file)
# Output/var.rds 
# "996db523f51c6db5fff2f02038d06454" 
# Output/var.rds 
# ""75ab3d7ee6bd774dfb6b302b1a55dc7c" 
file.info(file)$mtime
# [1] "2022-03-02 13:02:25 CST"
# "2022-02-18 13:15:25 CST"
file <- paste0('Output/','ancestry_data.rds')
write_rds(ancestry_data,file)
md5sum(file)
# Output/ancestry_data.rds 
# "3357362b61c59656b4d23b3834304d33" 
file.info(file)$mtime
#[1] "2022-03-02 13:03:08 CST"

# we don't locally store subj_data as it seems useless


