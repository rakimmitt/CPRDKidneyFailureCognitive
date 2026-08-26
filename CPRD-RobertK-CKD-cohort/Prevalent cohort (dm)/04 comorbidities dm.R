############################################################################################

#Setup
library(tidyverse)
library(aurum)
library(EHRBiomarkr)
rm(list=ls())

cprd = CPRDData$new(cprdEnv = "diabetes-jun2024",cprdConf = "C:\\Users\\rk535\\OneDrive\\1 - PhD\\Data Science\\CPRD\\.aurum.yaml")
codesets = cprd$codesets()
codes_2024 = codesets$getAllCodeSetVersion(v = "01/06/2024")

# Load Robert's additional codelists

# Use a new version name whenever the contents of your codelists change
custom_version <- "rk_2026-08-25"

# Replace this with the exact path copied from File Explorer
codelist_root <- paste0(
  "C:/Users/rk535/OneDrive/1 - PhD/Data Science/CPRD/Github clone/CPRDKidneyFailureCognitive/CPRD-Codelists"
)

custom_codelist_directories <- c(
  file.path(codelist_root, "Medcodes"),
  file.path(codelist_root, "ICD10"),
  file.path(codelist_root, "OPCS4")
)

missing_directories <- custom_codelist_directories[
  !dir.exists(custom_codelist_directories)
]

if (length(missing_directories) > 0) {
  stop(
    "The following codelist directories were not found: ",
    paste(missing_directories, collapse = ", ")
  )
}

# Reads compatible .txt files and loads them into the analysis database
codesets$loadAll(
  paths = custom_codelist_directories,
  version = custom_version
)

# Obtain database-backed versions of the newly loaded codelists
custom_codes <- codesets$getAllCodeSetVersion(
  version = custom_version
)

# Add new codelists to codes_2024.
# If a name already exists, combine the standard and custom lists.
for (code_name in names(custom_codes)) {

  if (code_name %in% names(codes_2024)) {

    codes_2024[[code_name]] <- dplyr::union(
      codes_2024[[code_name]],
      custom_codes[[code_name]]
    )

  } else {

    codes_2024[[code_name]] <- custom_codes[[code_name]]

  }
}

analysis_prefix = "ckd"

####################################

comorbids <- c("acutepancreatitis",
               "af",
               "angina",
               "anxiety_disorders",
               "asthma",
               "bph",
               "bronchiectasis",
               "chronicpancreatitis",
               "ckd5_code",
               "cld",
               "copd",
               "cysticfibrosis",
               "dementia",
               "diabeticnephropathy",
               "dka",
               "falls",
               "fh_premature_cvd", #family history of premature CVD - for QRISK2
               "frailty_simple",
               "haem_cancer",
               "heartfailure",
               "hosp_cause_majoramputation",
               "hosp_cause_minoramputation",
               "hypertension",
               "ihd", #ischaemic heart disease
               "incident_mi",
               "incident_stroke",
               "lowerlimbfracture",
               "micturition_control",
               "myocardialinfarction",
               #"neuropathy",
               "osteoporosis",
               "otherneuroconditions",
               "pad", #peripheral arterial disease
               "pulmonaryfibrosis",
               "pulmonaryhypertension",
               "severe_retinopathy",
               "non_severe_retinopathy",
               "revasc", #revascularisation procedure
               "rheumatoidarthritis",
               "solid_cancer",
               "solidorgantransplant",
               "stroke",
               "tia",  #transient ischaemic attack
               "ukpds_photocoagulation",
               "unstableangina",
               "urinary_frequency",
               #"vitreoushemorrhage",
               "volume_depletion",
               "genital_infection",
               "genital_infection_nonspec",
               "aav",
               "adpkd",
               "alport",
               "antigbm",
               "fabry",
               "fsgs",
               "gn_nos",
               "haemorrhagicstroke",
               "igan",
               "incident_stroke",
               "incident_mi",
               "ischaemicstroke",
               "kfdeath",
               "mcd",
               "membranous",
               "mpgn",
               "other_pkd",
               "sle",
               "GIinfection",
               "LRTI",
               "URTI",
               "UTI",
               "acutecholecystitis",
               "acutesinusitis",
               "boneinfection",
               "candidiasis",
               "cellulitis",
               "covid",
               "endocarditis",
               "eyeinfection",
               "infectiveotitisexterna",
               "influenza",
               "jointinfection",
               "meningitis",
               "otherfungalinfection",
               "otherskininfection",
               "pneumonia",
               "sepsis",
               "surgicalsiteinfection",
               "tuberculosis"
)

custom_comorbids <- c(
               "all dementia",
               "alzheimers",
               "ckd5_noKRT",
               "ckd5",
               "delirium",
               "haemodialysis",
               "mci",
               "peritoneal_dialysis",
               "renalaccessinfection",
               "transplant",
               "vascular_dementia",
               "uti",
               "skininfection",
               "respiratorytractinfection"
)

comorbids <- unique(c(comorbids, custom_comorbids))

############################################################################################

analysis = cprd$analysis("all_patid")


for (i in comorbids) {
  
  if (length(codes_2024[[i]]) > 0) {
    print(paste("making", i, "medcode table"))
    
    raw_tablename <- paste0("raw_", i, "_medcodes")
    
    data <- cprd$tables$observation %>%
  inner_join(codes_2024[[i]], by = "medcodeid") %>%
  analysis$cached(raw_tablename, indexes=c("patid", "obsdate"))
    
    assign(raw_tablename, data)
    
  }
  
  if (length(codes_2024[[paste0("icd10_", i)]]) > 0 & i!="hypertension") {
    print(paste("making", i, "ICD10 code table"))
    
    raw_tablename <- paste0("raw_", i, "_icd10")
    
    data <- cprd$tables$hesDiagnosisEpi %>%
  inner_join(
    codes_2024[[paste0("icd10_", i)]],
    sql_on = "LHS.ICD LIKE CONCAT(icd10,'%')" ) %>%
      analysis$cached(raw_tablename, indexes=c("patid", "epistart"))
    
    assign(raw_tablename, data)
  }
  
  if (length(codes_2024[[paste0("opcs4_", i)]]) > 0) {
    print(paste("making", i, "OPCS4 code table"))
    
    raw_tablename <- paste0("raw_", i, "_opcs4")
    
    data <- cprd$tables$hesProceduresEpi %>%
      inner_join(
        codes_2024[[paste0("opcs4_", i)]],
        by = c("OPCS" = "opcs4")
      ) %>%
      analysis$cached(raw_tablename, indexes=c("patid", "evdate"))
    
    assign(raw_tablename, data)
    
  }
  
}


## Add to beginning of list so don't have to remake interim tables when add new comorbidity to end of above list
# Make new primary cause hospitalisation for heart failure, incident MI, and incident stroke comorbidities

raw_primary_hhf_icd10 <- raw_heartfailure_icd10 %>%
  filter(d_order==1) %>%
  analysis$cached("raw_primary_hhf_icd10", indexes=c("patid", "epistart"))

raw_primary_incident_mi_icd10 <- raw_incident_mi_icd10 %>%
  filter(d_order==1) %>%
  analysis$cached("raw_primary_incident_mi_icd10", indexes=c("patid", "epistart"))

raw_primary_incident_stroke_icd10 <- raw_incident_stroke_icd10 %>%
  filter(d_order==1) %>%
  analysis$cached("raw_primary_incident_stroke_icd10", indexes=c("patid", "epistart"))


## Add to beginning of list so don't have to remake interim tables when add new comorbidity to end of above list
comorbids <- c("primary_hhf", "primary_incident_mi", "primary_incident_stroke", comorbids)


# Separate frailty by severity into three different categories
## Add to beginning of list so don't have to remake tables when add new comorbidity to end of above list
raw_frailty_mild_medcodes <- raw_frailty_simple_medcodes %>% filter(frailty_simple_cat=="Mild")
raw_frailty_moderate_medcodes <- raw_frailty_simple_medcodes %>% filter(frailty_simple_cat=="Moderate")
raw_frailty_severe_medcodes <- raw_frailty_simple_medcodes %>% filter(frailty_simple_cat=="Severe")
comorbids <- setdiff(comorbids, "frailty_simple")
comorbids <- c("frailty_mild", "frailty_moderate", "frailty_severe", comorbids)


############################################################################################

## Get index date

analysis = cprd$analysis(analysis_prefix)

# get dates at 6 month intervals
dates <- seq(from = as.Date("2019-03-01"),
             to   = as.Date("2024-03-01"),
             by   = "6 months")

date_strings <- format(dates, "%Y-%m-%d")


for (d in date_strings) {
  
  index_date <- as.Date(d)
  print(d)
  
  ## Clean comorbidity data and combine with index date
  
  for (i in comorbids) {
    
    print(paste("merging index date ", d, " with", i, "code occurrences"))
    
    index_date_merge_tablename <- paste0(d, "_full_", i, "_merge")
    
    medcode_tablename <- paste0("raw_", i, "_medcodes")
    icd10_tablename <- paste0("raw_", i, "_icd10")
    opcs4_tablename <- paste0("raw_", i, "_opcs4")
    
    
    if (exists(medcode_tablename)) {
      
      medcodes <- get(medcode_tablename) %>%
        select(patid, date=obsdate, code=medcodeid) %>%
        mutate(source="gp")
      
    }
    
    if (exists(icd10_tablename)) {
      
      icd10_codes <- get(icd10_tablename) %>%
        select(patid, date=epistart, code=ICD) %>%
        mutate(source="hes")
      
    }
    
    if (exists(opcs4_tablename)) {
      
      opcs4_codes <- get(opcs4_tablename) %>%
        select(patid, date=evdate, code=OPCS) %>%
        mutate(source="hes")
      
    }
    
    
    if (exists("medcodes")) {
      
      all_codes <- medcodes
      rm(medcodes)
      
      if(exists("icd10_codes")) {
        all_codes <- all_codes %>%
          union_all(icd10_codes)
        rm(icd10_codes)
      }
      
      if(exists("opcs4_codes")) {
        all_codes <- all_codes %>%
          union_all(opcs4_codes)
        rm(opcs4_codes)
      }
    }
    
    else if (exists("icd10_codes")) {
      
      all_codes <- icd10_codes
      rm(icd10_codes)
      
      if(exists("opcs4_codes")) {
        all_codes <- all_codes %>%
          union_all(opcs4_codes)
        rm(opcs4_codes)
      }
    }
    
    else if(exists("opcs4_codes")) {
      
      all_codes <- opcs4_codes
      rm(opcs4_codes)
    }
    
    all_codes_clean <- all_codes %>%
      inner_join(cprd$tables$validDateLookup, by="patid") %>%
      filter(date>=min_dob & ((source=="gp" & date<=gp_end_date) | (source=="hes" & date<=as.Date("2023-03-31")))) %>%
      select(patid, date, source, code)
    
    rm(all_codes)
    
    data <- all_codes_clean %>%
      mutate(datediff=datediff(date, index_date)) %>%
      analysis$cached(index_date_merge_tablename, index="patid")
    
    rm(all_codes_clean)
    
    assign(index_date_merge_tablename, data)
    
    rm(data)
    
  }
  
  
  ############################################################################################
  
  # Find earliest pre-index date, latest pre-index date and first post-index date dates
  
  comorbidities <- cprd$tables$patient %>%
    select(patid)
  
  for (i in comorbids) {
    
    print(paste("working out pre- and post-index date code occurrences for ", i, " at ", d))
    
    index_date_merge_tablename <- paste0(d, "_full_", i, "_merge")
    interim_comorbidity_table <- paste0(d, "_comorbidities_im_", i)
    pre_index_date_earliest_date_variable <- paste0("pre_index_date_earliest_", i)
    pre_index_date_latest_date_variable <- paste0("pre_index_date_latest_", i)
    pre_index_date_variable <- paste0("pre_index_date_", i)
    post_index_date_date_variable <- paste0("post_index_date_first_", i)
    
    pre_index_date <- get(index_date_merge_tablename) %>%
      filter(date<=index_date) %>%
      group_by(patid) %>%
      summarise({{pre_index_date_earliest_date_variable}}:=min(date, na.rm=TRUE),
                {{pre_index_date_latest_date_variable}}:=max(date, na.rm=TRUE)) %>%
      ungroup()
    
    post_index_date <- get(index_date_merge_tablename) %>%
      filter(date>index_date) %>%
      group_by(patid) %>%
      summarise({{post_index_date_date_variable}}:=min(date, na.rm=TRUE)) %>%
      ungroup()
    
    pre_index_date_earliest_date_variable <- as.symbol(pre_index_date_earliest_date_variable)
    
    comorbidities <- comorbidities %>%
      left_join(pre_index_date, by="patid") %>%
      mutate({{pre_index_date_variable}} := !is.na(!!pre_index_date_earliest_date_variable))
      left_join(post_index_date, by="patid") %>%
      analysis$cached(interim_comorbidity_table, unique_indexes="patid")
  }
  
  comorbidities <- comorbidities %>% analysis$cached(paste0(d, "_comorbidities"), unique_indexes="patid")
  
}

