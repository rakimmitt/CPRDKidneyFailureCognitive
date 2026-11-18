library(aurum)
library(tidyverse)
library(dplyr)

## load CPRD medcodes
 medcodelist <- readr::read_tsv("C:\\Users\\rk535\\OneDrive\\1 - PhD\\Data Science\\CPRD\\202406_EMISMedicalDictionary.txt", 
                          quote = "", 
                          col_names=T,
                          col_types = readr::cols(
                          MedCodeId = readr::col_character(),
                          SnomedCTDescriptionId = readr::col_character(),
                          SnomedCTConceptId = readr::col_character()
                                )
                          )

#include anything matching these strings:
medcodelist_highldl <- medcodelist %>% filter(grepl("hypercholesterolaemia|hyperlipidaemia|hypercholesterolemia|hyperlipidemia|high cholesterol|cholesterol high|cholesterol very high|high LDL|lipids abnormal|lipids high|dyslipid|low-density lipoproteinaemia|low density lipoproteinaemia|low-density lipoproteinemia|low density lipoproteinemia|familial xanthomatosis|xanthomatosis, familial|hyperbetalipoproteinaemia|lipidaemia|hyperlipoproteinaemia", Term, ignore.case = T))

#exclude entries if matching any of these strings:
medcodelist_highldl <- medcodelist_highldl %>% filter(!grepl("high density lipoprotein below|hypolipidaemia|retinalis|provision|uncertain significance|not detected|CGI|genotyping|drugs|criteria|lung function|behaviour|assessment|test|screen|electrop|gallbladder|FH|family history", Term, ignore.case = T))

View(medcodelist_highldl)

readr::write_tsv(medcodelist_high_LDL, "C:\\Users\\rk535\\OneDrive\\1 - PhD\\Data Science\\CPRD\\Github clone\\CPRDKidneyFailureCognitive\\CPRD-Codelists\\Medcodes\\exeter_medcodelist_highldl.txt")