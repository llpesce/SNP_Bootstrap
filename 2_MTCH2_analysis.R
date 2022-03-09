# USE: Complete analysis of Julie's data on MTCH2
# SCOPE: HCM/DCM + gene + variant 
# When: 2/09  
# Who: McNally Lab (LLP & SDK)
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

# var is the set of all the variants that were left after clean up and filtering for artifacts
# var.bs is the dataframe used in the bootstrap analysis in file "1_run_bootstrap.R" so should be .25-.1
dim(var)
# [1] 3074320      28
var %>% summarise(across(everything(),~n_distinct(.)))
# Id Eff phen  Gene Chr    Loc  Ref  Alt  PP2 AF GERP  GTEx CLNSIG CLNSIGCONF DCM_FREQ HCM_FREQ NUGENOMES_FREQ gnomAD_AF_G gnomAD_AF_AFR_G gnomAD_AF_NFE_G
# 1 264   2    8 20372  23 157112 2643 7553 1003  2 4272 20015     37        423      171      138            901       60557           39009           39621
# gnomAD_AF_EAS_G gnomAD_COV_G ENCODE_BL RF AC NC gnomad_ancestry Sex
# 1           12191        15198         1  1  2  2               3   2

dim(var.bs)
# [1] 811053     14
var.bs %>% summarise(across(everything(),~n_distinct(.)))
# Id type.of.CM gnomad_ancestry Sex Gene Chr  Loc Ref Alt AC NC gnomAD_AF_AFR_G gnomAD_AF_NFE_G gnomAD_AF_EAS_G
# 1 172          2               2   2 6011  23 9697 266 212  2  2            9445            9594            6492

# There is only one variant in the data selected by the bootstrap
var.bs %>% filter(Gene == 'MTCH2') %>% distinct(Chr, Loc, Ref, Alt)
# Chr      Loc Ref Alt
# 1 chr11 47640429   G   C

# AF comparison with NU GENE, hard coded numbers from Sam
var.bs %>% filter(Loc == 47640429 & Chr == 'chr11') %>% summarize(sum(AC), sum(NC))
# sum(AC) sum(NC)
# 1      77     138
x <- c(77, 441)
n <- c(344, 1544)
prop.test(x, n, p = NULL, alternative = "two.sided",
          correct = TRUE)
# 
# 2-sample test for equality of proportions with continuity correction
# 
# data:  x out of n
# X-squared = 5.0882, df = 1, p-value = 0.02409
# alternative hypothesis: two.sided
# 95 percent confidence interval:
#   -0.1130366 -0.0105325
# sample estimates:
#   prop 1    prop 2 
# 0.2238372 0.2856218 

# Finding the IDs by genotype 
tmp <- left_join(var.bs %>% distinct(Id) , 
                 var.bs %>% filter(Loc == 47640429 & Chr == 'chr11')  %>% 
                   transmute(Id, G = if_else(AC == '2', 'HOM', 'HET')) 
)  %>% 
  mutate(G = replace_na(G,replace = 'REF'))
tmp %>% select(G) %>% table
# HET HOM REF 
# 61   8 103 
tmp <- tmp %>% left_join(subj_data %>% select(Id,patID ), by= 'Id') 
write_csv(tmp, 'subjects_MTCH2_var.csv')
md5sum('subjects_MTCH2_var.csv')
# subjects_MTCH2_var.csv 
# "07925ef11a0d9d0bab5e96138c4e5304" 
subj_MTCH2 <- tmp
rm(tmp)

tmp <- readxl::read_excel(path = '/data2/Cohort_Patient_Data_MASTERUPDATE090418.xlsx')
dim(tmp)

tmp <- inner_join(subj_MTCH2, tmp, by=c('Id' = 'ID'))

# Check what variables we have =====
tmp %>% summarize(across(everything(), function(x) {sum(is.na(x)) }))
names(tmp) <- make.names(names(tmp))
tmp$age.1st.present <- as.numeric(tmp$age.1st.present)
# Check if we can predict standard cardiac vas --------------------
for (var in c('age.1st.present','IVSd_BSA','EF', 'LVIDd_BSA', 'LVPWd_BSA')){
  print(paste('FITTING',var))
  form <- paste(var,"~ type.of.CM + G1")
  print('Additive Model')
  tmp %>% mutate(G1 = case_when(G == 'REF' ~ 0, G == 'HET' ~ 1, T ~ 2 ))  %$%
    lm(eval(parse(text=form))) %>% summary %>% print()
  print('Dominant Model')
  tmp %>% mutate(G1 = case_when(G == 'REF' ~ "REF",  T ~ 'VAR' )) %$% 
    lm(eval(parse(text=form))) %>% summary %>% print()
  print('')
  print('')
}

tmp.1 <- tmp %>% select(Id,G, IVSd_BSA,EF, LVIDd_BSA,LVPWd_BSA) %>% filter(across(everything(), ~!is.na(.))) 

pc <- prcomp(x=as.matrix(tmp.1[,3:6]), scale=T, center = T)

lm(pc$x[,'PC1'] ~ (tmp.1 %>% mutate(G1 = case_when(G == 'REF' ~ 0, G == 'HET' ~ 1, T ~ 2 )))$G1) %>% summary()
lm(pc$x[,'PC2'] ~ (tmp.1 %>% mutate(G1 = case_when(G == 'REF' ~ 0, G == 'HET' ~ 1, T ~ 2 )))$G1) %>% summary()
lm(pc$x[,'PC3'] ~ (tmp.1 %>% mutate(G1 = case_when(G == 'REF' ~ 0, G == 'HET' ~ 1, T ~ 2 )))$G1) %>% summary()
lm(pc$x[,'PC4'] ~ (tmp.1 %>% mutate(G1 = case_when(G == 'REF' ~ 0, G == 'HET' ~ 1, T ~ 2 )))$G1) %>% summary()


pc$x %>% as.data.frame() %>% cbind(.,tmp.1) %>% ggplot() + geom_point(aes(PC1,PC2, col = G))
pc$x %>% as.data.frame() %>% cbind(.,tmp.1) %>% ggplot() + geom_point(aes(PC2,PC3, col = G))
