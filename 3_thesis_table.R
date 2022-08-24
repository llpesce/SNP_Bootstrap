# when August 24th 2022
# Who  McNally lab (LLP)
# What making a table1 for the 172 subjects in the bootstrap set in Julie's thesis
# Note, it uses the boostrap data from produced by the processing input and running the boostrap files

library('tableone')

# The boostrap data has ancestry meant to be matched with the gnomAD categories, we need the original ones
# that are in the ancestry dataframe
tbl <- var.bs %>% distinct(Id,type.of.CM,Sex,gnomad_ancestry) %>% 
  inner_join(ancestry_data.bs, by = 'Id') %>% droplevels()
nrow(tbl)
# [1] 172

#table 1
listVars <- c("ANCESTRY", "type.of.CM", "Sex")
catVars <- listVars
exact <- catVars
table1 <- CreateTableOne(vars = listVars, data = tbl,factorVars = listVars)# , catVars)
print(table1, quote = TRUE, noSpaces = TRUE,showAllLevels = TRUE)
# ""
# ""                "level" "Overall"   
# "n"              ""      "172"       
# "ANCESTRY (%)"   "AFR"   "30 (17.4)" 
# ""               "EUR"   "142 (82.6)"
# "type.of.CM (%)" "DCM"   "105 (61.0)"
# ""               "HCM"   "67 (39.0)" 
# "Sex (%)"        "f"     "77 (44.8)" 
# ""               "m"     "95 (55.2)" 
