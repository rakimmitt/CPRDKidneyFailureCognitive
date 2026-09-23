library(aurum)
library(tidyverse)
library(dplyr)

## load CPRD medcodes
medcodelist <- read.delim("C:\\Users\\rk535\\OneDrive\\1 - PhD\\Data Science\\CPRD\\202406_EMISMedicalDictionary.txt", 
                          sep = "\t", 
                          quote = "", 
                          header=T,  
                          col_types = readr::cols(
                          MedCodeId = readr::col_character(),
                          SnomedCTDescriptionId = readr::col_character(),
                          SnomedCTConceptId = readr::col_character()
                                )
                          )

#include anything matching these strings:
medcodelist_dementia <- medcodelist %>% filter(grepl("dement|alzheimer|pick'?s? disease|binswang", Term, ignore.case = T))

#exclude entries if matching any of these strings:
medcodelist_dementia <- medcodelist_dementia %>% filter(!grepl("menthol|debridement|evidement|assessment|review|monitor|screen|protocol|therapy|quality|indicator|tool|manag|awareness|support|risk|staff group|nursing|reference|detection|remission|measurement|provision|information|education|prescribing|understanding|planning|pillar|PDS|question|link worker|registered|ten item|dementia test|leaflet|referral|refer to|advisor|signpost|society|plan|wellbeing|score|carer|healthcare|family history|FH:|declined|rating|scale|relative|occupant|motorcycle|pedestrian|cyclist|driver|collision|clothing|arizona|pickwick|pickle|nail|finger|heart|dental|vein|ice-pick|ice pick|picking|pickpocket|pickerel|pickup|pick-up|pick up|picks up|picker|pick's cell|picket|middle ear pick|picky|infant|fiber|nodule|adenoma|burkholder|between meals|ralston|butterfly|falling|suspected|antidementia|pseudodementia|not superimposed|nieman|glial|neurofibrillary|pick cell|glial cell", Term, ignore.case = T))

#substance related things I think we should remove - but can review and reconsider
medcodelist_dementia <- medcodelist_dementia %>% filter(!grepl("Mental & behav dis|Mental and behav dis|hallucinogen|cannabi|opioids|hypntcs|cocaine|stimulant|inhalant|psych subs|sedative|hypnotic", Term, ignore.case = T))

##  --> this creates a large list with a lot of extra entries explicitly describing dementia which are not in the alldementia.txt codelist



## comparing lists:

exeter_medcodelist_alldementia <- read.delim("C:/Users/rk535/OneDrive/1 - PhD/Data Science/CPRD/Github clone/CPRDKidneyFailureCognitive/CPRD-Codelists/Medcodes/exeter_medcodelist_alldementia.txt")

# Retain one row per medcodeid in each codelist
tj_codes <- medcodelist_dementia %>%
  mutate(MedCodeId = as.character(MedCodeId)) %>%
  distinct(MedCodeId, .keep_all = TRUE)

exeter_codes <- exeter_medcodelist_alldementia %>%
  mutate(MedCodeId = as.character(MedCodeId)) %>%
  distinct(MedCodeId, .keep_all = TRUE)

# Codes present in both codelists
codes_in_both <- tj_codes %>%
  semi_join(exeter_codes, by = "MedCodeId")


# Codes in medcodelist_dementia but not in the Exeter codelist
codes_only_in_tj <- tj_codes %>%
  anti_join(exeter_codes, by = "MedCodeId")


# Optional: codes in the Exeter codelist but not medcodelist_dementia
codes_only_in_original <- exeter_codes %>%
  anti_join(tj_codes, by = "MedCodeId")

# Summary

overlap_summary <- tibble(
  comparison = c(
    "Unique codes in medcodelist_dementia",
    "Unique codes in exeter_medcodelist_alldementia",
    "Codes present in both",
    "Codes only in medcodelist_dementia",
    "Codes only in exeter_medcodelist_alldementia"
  ),
  number_of_codes = c(
    nrow(tj_codes),
    nrow(exeter_codes),
    nrow(codes_in_both),
    nrow(codes_only_in_tj),
    nrow(codes_only_in_original)
  )
)

overlap_summary

View(codes_only_in_tj)
write.table(codes_only_in_tj, file = "C:\\Users\\rk535\\Downloads\\codes_only_in_tj", sep= "\t", quote=FALSE, row.names=FALSE) 

## Previous code


## search for relevant terms a|b|c (a OR b OR c)
medcodelist_renalaccessinfection <- medcodelist %>% filter(grepl("terms|here", Term, ignore.case = T))

## search for relevant terms using AND logic between lines (and OR logic within lines)
medcodelist_renalaccessinfection <- medcodelist %>%
  filter(
    grepl("infect|sepsis", Term, ignore.case = TRUE) &
    grepl("dialysis|venous|cannula|catheter|central", Term, ignore.case = TRUE)
  )

## inspect medcodelist in observation order

medcodelist_renalaccessinfection %>%
  arrange(Observations) %>%
  View()

## Filter out irrelevant terms

medcodelist_renalaccessinfection <- medcodelist_renalaccessinfection %>% filter(!grepl("urin|nerv|bladd|vir", Term, ignore.case = T))

#sense-check all renal access infection entries we wouldn't be capturing:
check <- medcodelist %>% filter(grepl("access", Term, ignore.case = T)) %>% filter(!Term %in% medcodelist_renalaccessinfection$Term)

write.table(medcodelist_renalaccessinfection, file = "C:\\Users\\Robert\\OneDrive\\1 - PhD\\Data Science\\CPRD\\Codelists\\medcodes\\exeter_medcodelist_renalaccessinfection.txt", sep= "\t", quote=FALSE, row.names=FALSE) 

# delirium

medcodelist <- medcodelist %>%
  mutate(across(
    c(MedCodeId, SnomedCTDescriptionId, SnomedCTConceptId),
    as.character
  ))

medcodelist_delirium <- medcodelist %>% filter(grepl("deliri|confus", Term, ignore.case = T))#
medcodelist_delirium <- medcodelist_delirium %>% filter(!grepl("tremens|dis due|disease due|dis multi|kors|wernick|sexual|identity|withdrawal|withdr|wthdr|arousal|epilep|seizure|ictal|migraine|cannabi|thc|pcp|CRB|CAM|assessment|monitoring|CURB|ACVPU|urea|respiratory rate|concussion|trauma|left-right|dissociative|tetrameres|cottus|lactobacillus|peziza|weissella|taenia|rutland|coryne|risk|level|care|management|education|emergence|exhaustion|O/E - confused|O/E - mentally confused|mentally confused|stimulant induced|stimulant-induced|delirium caused by|delirium due|drug induced delirium|drug-induced delirium|opioid delirium|cocaine induced|cocaine-induced|intoxication|hallucinogen induced|amphetamine|phencycl|chronic confusion|reactive|Confused speech content|not confused|Onset of confusion|Confused speech|synthetic cath|Psychoactive substance-induced|cocaine delirium|inhalant delirium|inhalant induced|opioid induced|hallucinogen delirium|Orientation confused|psychogenic confusion|senile confusion|intermittent confusion", Term, ignore.case=T))
medcodelist_delirium <- medcodelist_delirium %>% filter(!grepl("^Confusion$|^Confused$|^Confusional state$|\\[D\\]Confusion", Term, ignore.case=F))

View(medcodelist_delirium)

