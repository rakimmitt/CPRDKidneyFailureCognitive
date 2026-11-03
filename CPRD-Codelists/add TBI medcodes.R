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
medcodelist_tbi <- medcodelist %>% filter(grepl("traumatic brain|concussion|intracranial injury|brain contusion|contusion of brain|traumatic subdural|traumatic subarachnoid|traumatic extradural", Term, ignore.case = T))

#exclude entries if matching any of these strings:
medcodelist_tbi <- medcodelist_tbi %>% filter(!grepl("labyrinthine|tooth|questionnaire|nontraumatic|spinal|non-traumatic|software|without mention intracranial|no intracranial|without mention of intracranial|without intracranial", Term, ignore.case = T))

View(medcodelist_tbi)

save(medcodelist_tbi, file = paste0("C:\\Users\\rk535\\OneDrive\\1 - PhD\\Data Science\\CPRD\\Github clone\\CPRDKidneyFailureCognitive\\CPRD-Codelists\\Medcodes\\exeter_medcodelist_tbi.txt"))
