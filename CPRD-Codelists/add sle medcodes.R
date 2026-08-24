## Thijs' code for searching through medcode lists for a specific condition

library(aurum)
library(tidyverse)


## load all CPRD Aurum medcodes
medcodelist <- read.delim("C:/Users/tj358/OneDrive - University of Exeter/CPRD/Aurum codelists/CPRDAurumMedical.txt", sep = "\t", quote = "", header=T)

medcodelist_sle <- medcodelist %>% filter(grepl("systemic lupus|disseminated lupus|lupus erythematosus|lupus nephritis|cerebral lupus", Term, ignore.case = T))

medcodelist_sle <- medcodelist_sle %>% filter(!grepl("drug|activity|history|cutaneous|discoid|cell|limited|test|nail|melanosis|wart|nodul|tumid|factor|chillblain|chilblain|migrans|mucous|oral|lichen|rash|derm|cheilitis|calcinosis", Term, ignore.case = T))

#sense-check all lupus entries we wouldn't be capturing:
check <- medcodelist %>% filter(grepl("lupus", Term, ignore.case = T)) %>% filter(!Term %in% medcodelist_sle$Term)

write.table(medcodelist_sle, file = "C:/Users/tj358/OneDrive - University of Exeter/CPRD/Aurum codelists/medcodes/exeter_medcodelist_sle.tsv", sep= "\t", quote=FALSE, row.names=FALSE) 
