############################################################################################

#Setup
library(tidyverse)
library(aurum)
library(EHRBiomarkr)
rm(list=ls())

cprd = CPRDData$new(cprdEnv = "diabetes-jun2024",cprdConf = "C:\\Users\\rk535\\OneDrive\\1 - PhD\\Data Science\\CPRD\\.aurum.yaml")
codesets = cprd$codesets()
codes_2024 = codesets$getAllCodeSetVersion(v = "01/06/2024")

# Read dementia-medication product-code lists locally

dementia_meds_directory <- paste0(
  "C:/Users/rk535/OneDrive/1 - PhD/Data Science/CPRD/Github clone/CPRDKidneyFailureCognitive/CPRD-Codelists/Prodcodes/Dementia medications"
)

if (!dir.exists(dementia_meds_directory)) {
  stop(
    "The dementia-medication codelist directory was not found: ",
    dementia_meds_directory
  )
}

dementia_medication_files <- list.files(
  dementia_meds_directory,
  pattern = "\\.txt$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(dementia_medication_files) == 0) {
  stop(
    "No .txt codelist files were found in: ",
    dementia_meds_directory
  )
}

read_dementia_medication_file <- function(file) {

  codelist <- readr::read_tsv(
    file,
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE,
    progress = FALSE
  ) %>%
    rename_with(stringr::str_to_lower)

  if (!"prodcodeid" %in% names(codelist)) {
    stop(
      "The following file does not contain a prodcodeid column: ",
      file
    )
  }

  codelist %>%
    transmute(
      prodcodeid = stringr::str_remove_all(
        prodcodeid,
        "[^0-9]"
      )
    ) %>%
    filter(
      !is.na(prodcodeid),
      prodcodeid != ""
    ) %>%
    mutate(
      prodcodeid = bit64::as.integer64(prodcodeid)
    )
}

dementia_medication_codes <- purrr::map_dfr(
  dementia_medication_files,
  read_dementia_medication_file
) %>%
  distinct(prodcodeid)

if (nrow(dementia_medication_codes) == 0) {
  stop("The dementia-medication codelist contains no valid product codes")
}

print(
  paste(
    nrow(dementia_medication_codes),
    "unique dementia-medication product codes loaded"
  )
)

analysis_prefix = "ckd"

############################################################################################


# Define medications

meds <- c("ace_inhibitors",
          "arb",
          "statins",
          "finerenone",
          "beta_blockers",
          "ca_channel_blockers",
          "thiazide_diuretics",
          "loop_diuretics",
          "ksparing_diuretics",
          "definite_genital_infection_meds",
          "topical_candidal_meds",
          "dementia_medications")


############################################################################################

# Pull out clean script instances

analysis <- cprd$analysis("all_patid")

for (i in meds) {

  print(i)

  raw_tablename <- paste0("raw_", i, "_prodcodes")

  if (i == "dementia_medications") {

    # This codelist is stored locally, so copy=TRUE is required
    code_list <- dementia_medication_codes
    copy_code_list <- TRUE

  } else {

    # All other lists are database-backed standard codelists
    code_list <- codes_2024[[i]]
    copy_code_list <- FALSE
  }

  if (is.null(code_list)) {
    stop("No product-code codelist was found for: ", i)
  }

  data <- cprd$tables$drugIssue %>%
    inner_join(
      code_list,
      by = "prodcodeid",
      copy = copy_code_list
    ) %>%
    select(
      patid,
      date = issuedate
    ) %>%
    analysis$cached(
      raw_tablename,
      indexes = c("patid", "date")
    )

  assign(raw_tablename, data)
}

clean_oha_prodcodes <- clean_oha_prodcodes %>% analysis$cached("clean_oha_prodcodes")

drugclasses <- c("Acarbose", "DPP4", "GIPGLP1", "Glinide", "GLP1", "MFN", "SGLT2", "SU", "TZD")

for (dc in drugclasses) {
  
  prodcode_list_name = paste0("clean_", dc, "_prodcodes")
  
  prodcode_list <- clean_oha_prodcodes %>%
    filter(drug_class_1 == dc | drug_class_2 == dc) %>%
    select(patid, date)
  
  assign(prodcode_list_name, prodcode_list)
  
  
}

raw_dementia_medications_prodcodes %>%
  summarise(
    prescription_issues = n(),
    patients = n_distinct(patid)
  ) %>%
  collect()

clean_insulin_prodcodes <- clean_insulin_prodcodes %>% analysis$cached("clean_insulin_prodcodes")

clean_insulin_prodcodes <- clean_insulin_prodcodes %>%
  select(patid, date) %>%
  union_all(clean_oha_prodcodes %>%
              filter(drug_class_1=="INS" | drug_class_2=="INS") %>%
              select(patid, date))


meds <- c(meds, drugclasses, "insulin")

############################################################################################

analysis = cprd$analysis(analysis_prefix)


# get dates at 6 month intervals
dates <- seq(from = as.Date("2019-03-01"),
             to   = as.Date("2024-03-01"),
             by   = "6 months")

date_strings <- format(dates, "%Y-%m-%d")


for (d in date_strings) {
  print(d)
  index_date <- as.Date(d)
  
  
  for (i in meds) {
    
    print(i)
    
    if (i %in% drugclasses | i=="insulin") {
      
      clean_tablename <- paste0("clean_", i, "_prodcodes")
      index_date_merge_tablename <- paste0(d, "_full_", i, "_merge")
      
      data <- get(clean_tablename) %>%
        mutate(datediff=datediff(date, index_date)) %>%
        analysis$cached(index_date_merge_tablename, indexes="patid")
      
      assign(index_date_merge_tablename, data)
      
      rm(data)
      
    } else {
    
    raw_tablename <- paste0("raw_", i, "_prodcodes")
    index_date_merge_tablename <- paste0(d, "_full_", i, "_merge")
    
    data <- get(raw_tablename) %>%
      
      inner_join(cprd$tables$validDateLookup, by="patid") %>%
      # filter(date>=min_dob & date<=gp_ons_end_date) %>%
      filter(date>=min_dob & date<=gp_end_date) %>%
      select(patid, date) %>%
      
      mutate(datediff=datediff(date, index_date)) %>%
      
      analysis$cached(index_date_merge_tablename, indexes="patid")
    
    assign(index_date_merge_tablename, data)
    
    rm(data)
    }
  }
  
  
  ############################################################################################
  
  # Find earliest pre-index date, latest pre-index date and first post-index date dates
  
  
  medications <- cprd$tables$patient %>%
    select(patid)
  
  
  for (i in meds) {
    
    print(paste("working out pre- and post- index date code occurrences for", i, " at ", d))
    
    index_date_merge_tablename <- paste0(d, "_full_", i, "_merge")
    interim_medications_table <- paste0(d, "_meds_im_", i)
    pre_index_date_earliest_date_variable <- paste0("pre_index_date_earliest_", i, "")
    pre_index_date_latest_date_variable <- paste0("pre_index_date_latest_", i, "")
    post_index_date_date_variable <- paste0("post_index_date_first_", i, "")
    
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
    
    medications <- medications %>%
      left_join(pre_index_date, by="patid") %>%
      left_join(post_index_date, by="patid") %>%
      analysis$cached(interim_medications_table, unique_indexes="patid")
    
  }
  
  
  # Cache final version
  
  medications <- medications %>% analysis$cached(paste0(d, "_medications"), unique_indexes="patid")
  
}
