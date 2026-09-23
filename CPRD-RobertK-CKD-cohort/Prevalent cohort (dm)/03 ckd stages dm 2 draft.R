############################################################################################

# Setup
library(tidyverse)
library(aurum)
library(EHRBiomarkr)
rm(list=ls())

cprd = CPRDData$new(cprdEnv = "diabetes-jun2024",cprdConf = "C:\\Users\\rk535\\OneDrive\\1 - PhD\\Data Science\\CPRD\\.aurum.yaml")

codesets = cprd$codesets()
codes = codesets$getAllCodeSetVersion(v = "01/06/2024")

analysis_prefix <- "rk_ckd"
source_analysis <- cprd$analysis("all_patid")

clean_egfr_medcodes <- source_analysis$cached(
  name = "clean_egfr_medcodes"
)
analysis <- cprd$analysis(analysis_prefix)

clean_egfr_medcodes %>% count()

############################################################################################

# Convert eGFR to CKD stage

ckd_stages_from_all_egfr <- clean_egfr_medcodes %>%
  rename(egfr = testvalue) %>%
  mutate(ckd_stage=ifelse(egfr<15, "stage_5",
                          ifelse(egfr<30, "stage_4",
                                 ifelse(egfr<45, "stage_3b",
                                        ifelse(egfr<60, "stage_3a",
                                               ifelse(egfr<90, "stage_2",
                                                      ifelse(egfr>=90, "stage_1", NA)))))))

################################################################################################################################

# Only keep CKD stages if >1 consecutive test with the same stage, and if time between earliest and latest consecutive test with same stage are >=90 days apart

## For each patient:
### A) Define period from current test until next test as having the ckd_stage of current test
### B) Join together consecutive periods with the same ckd_stage
### C) If period contains >1 test, and there is >=90 days between the first and last test in the period, it is 'confirmed'

### A) Define period from current test until next test as having the ckd_stage of current test

#### Add in row labelling within each patient's values + max number of rows for each patient

ckd_stages_from_algorithm <- ckd_stages_from_all_egfr %>%
  group_by(patid) %>%
  dbplyr::window_order(date) %>%
  mutate(patid_row_id=row_number()) %>%
  mutate(patid_total_rows=max(patid_row_id, na.rm=TRUE)) %>%
  ungroup()

#### For rows where there is a next test, use this as end date; for last row, use start date as end date

ckd_stages_from_algorithm <- ckd_stages_from_algorithm %>%
  mutate(next_row=patid_row_id+1) %>%
  left_join(ckd_stages_from_algorithm, by=c("patid","next_row"="patid_row_id")) %>%
  mutate(ckd_start=date.x,
         ckd_end=if_else(is.na(date.y),date.x,date.y),
         ckd_stage=ckd_stage.x,
         egfr=egfr.x) %>%
  select(patid, patid_row_id, ckd_stage, ckd_start, ckd_end, egfr)

### B) Join together consecutive periods with the same ckd_stage

ckd_stages_from_algorithm <- ckd_stages_from_algorithm %>%
  group_by(patid, ckd_stage) %>%
  dbplyr::window_order(patid, ckd_stage, patid_row_id) %>%
  mutate(lead_var=lead(ckd_start),
         cummax_var=cummax(ckd_end)) %>%
  mutate(compare=cumsum(lead_var>cummax_var)) %>%
  mutate(indx=ifelse(row_number()==1, 0L, lag(compare))) %>%
  ungroup() %>%
  group_by(patid, ckd_stage ,indx) %>%
  summarise(first_test_date=min(ckd_start,na.rm=TRUE),
            last_test_date=max(ckd_start,na.rm=TRUE),
            maximum_date=max(ckd_end,na.rm=TRUE),
            test_count=max(patid_row_id, na.rm=TRUE)-min(patid_row_id, na.rm=TRUE)+1) %>%
  ungroup() %>%
  analysis$cached("ckd_stages_from_algorithm_interim_1",indexes=c("patid", "ckd_stage", "test_count", "first_test_date", "last_test_date"))

ckd_stages_from_algorithm %>% count()
#40,051,589

ckd_stages_from_algorithm %>% summarise(total=sum(test_count, na.rm=TRUE))
#total number of tests: 112,440,442 as above

### C) Remove periods with 1 reading, or with multiple readings but <90 days between first and last test, and cache

ckd_stages_from_algorithm <- ckd_stages_from_algorithm %>%
  filter(test_count>1 & datediff(last_test_date, first_test_date)>=90) %>%
  analysis$cached("ckd_stages_from_algorithm_interim_2",indexes=c("patid","ckd_stage","first_test_date"))

ckd_stages_from_algorithm %>% count()
#17,886,174

################################################################################################################################

# Combine with CKD5 medcodes/ICD10/OPCS4 codes

codelist_root <- paste0(
  "C:/Users/rk535/OneDrive/1 - PhD/Data Science/CPRD/",
  "Github clone/CPRDKidneyFailureCognitive/CPRD-Codelists"
)

read_ckd5_list <- function(folder, wanted_name, code_column) {
  directory <- file.path(codelist_root, folder)
  if (!dir.exists(directory)) stop("Codelist directory missing: ", directory)

  files <- list.files(directory, pattern = "\\.txt$", ignore.case = TRUE,
                      recursive = TRUE, full.names = TRUE)
  names_in_files <- files %>%
    basename() %>%
    tools::file_path_sans_ext() %>%
    str_to_lower() %>%
    str_remove("^exeter_medcodelist_") %>%
    str_remove("^exeter_")
  selected <- files[names_in_files == wanted_name]

  if (length(selected) == 0L) return(NULL)
  if (length(selected) > 1L) stop("Multiple codelists found for: ", wanted_name)

  codes <- readr::read_tsv(
    selected, col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE, progress = FALSE
  ) %>%
    rename_with(str_to_lower) %>%
    select(all_of(code_column)) %>%
    filter(!is.na(.data[[code_column]]), .data[[code_column]] != "") %>%
    distinct()

  if (nrow(codes) == 0L) stop("Empty codelist: ", selected)
  codes
}

ckd5_medcodes <- read_ckd5_list("Medcodes", "ckd5", "medcodeid")
ckd5_icd10 <- read_ckd5_list("ICD10", "icd10_ckd5", "icd10")
ckd5_opcs4 <- read_ckd5_list("OPCS4", "opcs4_ckd5", "opcs4")

if (is.null(ckd5_medcodes) && is.null(ckd5_icd10) && is.null(ckd5_opcs4)) {
  stop("No custom ckd5 codelist was found in any coding system.")
}

## Clean, find earliest date per person, and re-cache

# Find patient records containing the CKD5 codes.
# Only include coding systems for which a codelist was found.

ckd5_records <- list()

# Medcodes
if (!is.null(ckd5_medcodes)) {

  ckd5_records[["gp"]] <- cprd$tables$observation %>%
    inner_join(
      ckd5_medcodes,
      by = "medcodeid",
      copy = TRUE
    ) %>%
    select(patid, date = obsdate) %>%
    mutate(source = "gp")
}

# ICD10
if (!is.null(ckd5_icd10)) {

  ckd5_records[["icd10"]] <- cprd$tables$hesDiagnosisEpi %>%
    inner_join(
      ckd5_icd10,
      sql_on = "LHS.ICD LIKE CONCAT(RHS.icd10, '%')",
      copy = TRUE
    ) %>%
    select(patid, date = epistart) %>%
    mutate(source = "hes")
}

# OPCS4
if (!is.null(ckd5_opcs4)) {

  ckd5_records[["opcs4"]] <- cprd$tables$hesProceduresEpi %>%
    inner_join(
      ckd5_opcs4,
      by = c("OPCS" = "opcs4"),
      copy = TRUE
    ) %>%
    select(patid, date = evdate) %>%
    mutate(source = "hes")
}

# Combine the available patient-record tables
all_ckd5_records <- purrr::reduce(
  ckd5_records,
  dplyr::union_all
)

# Apply your existing date rules and find the earliest CKD5 record
earliest_clean_ckd5 <- all_ckd5_records %>%
  inner_join(
    cprd$tables$validDateLookup,
    by = "patid"
  ) %>%
  filter(
    date >= min_dob &
      (
        (source == "gp" & date <= gp_end_date) |
        (source == "hes" &
           (is.na(gp_end_date) | date <= gp_end_date))
      )
  ) %>%
  group_by(patid) %>%
  summarise(
    first_test_date = min(date, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  analysis$cached(
    name = "earliest_clean_ckd5",
    indexes = c("patid", "first_test_date")
  )

## Combine CKD5 and other codes

ckd_stages_from_algorithm <- ckd_stages_from_algorithm %>%
  select(patid, ckd_stage, first_test_date) %>%
  union_all(earliest_clean_ckd5 %>% mutate(ckd_stage="stage_5")) %>%
  analysis$cached("ckd_stages_from_algorithm_interim_3",indexes=c("patid","ckd_stage","first_test_date"))

ckd_stages_from_algorithm %>% count()        
#12,130,677

################################################################################################################################

# Define date of onset for each stage

## For each person, define date of onset of each stage (earliest incident) - assume no returning to less severe stages

ckd_stages_from_algorithm <- ckd_stages_from_algorithm %>%
  group_by(patid, ckd_stage) %>%
  summarise(ckd_stage_start=min(first_test_date, na.rm=TRUE)) %>% 
  ungroup()

## Remove where start date of less severe stage is later than start date of more severe stage
### Reshape wide first

ckd_stages_from_algorithm <- ckd_stages_from_algorithm %>%
  pivot_wider(id_cols=patid,
              names_from=ckd_stage,
              values_from=ckd_stage_start) %>%
  mutate(stage_1=ifelse(!is.na(stage_1) & !is.na(stage_2) & stage_1>stage_2, NA, stage_1),
         stage_1=ifelse(!is.na(stage_1) & !is.na(stage_3a) & stage_1>stage_3a, NA, stage_1),
         stage_1=ifelse(!is.na(stage_1) & !is.na(stage_3b) & stage_1>stage_3b, NA, stage_1),
         stage_1=ifelse(!is.na(stage_1) & !is.na(stage_4) & stage_1>stage_4, NA, stage_1),
         stage_1=ifelse(!is.na(stage_1) & !is.na(stage_5) & stage_1>stage_5, NA, stage_1),
         stage_2=ifelse(!is.na(stage_2) & !is.na(stage_3a) & stage_2>stage_3a, NA, stage_2),
         stage_2=ifelse(!is.na(stage_2) & !is.na(stage_3b) & stage_2>stage_3b, NA, stage_2),
         stage_2=ifelse(!is.na(stage_2) & !is.na(stage_4) & stage_2>stage_4, NA, stage_2),
         stage_2=ifelse(!is.na(stage_2) & !is.na(stage_5) & stage_2>stage_5, NA, stage_2),
         stage_3a=ifelse(!is.na(stage_3a) & !is.na(stage_3b) & stage_3a>stage_3b, NA, stage_3a),
         stage_3a=ifelse(!is.na(stage_3a) & !is.na(stage_4) & stage_3a>stage_4, NA, stage_3a),
         stage_3a=ifelse(!is.na(stage_3a) & !is.na(stage_5) & stage_3a>stage_5, NA, stage_3a),
         stage_3b=ifelse(!is.na(stage_3b) & !is.na(stage_4) & stage_3b>stage_4, NA, stage_3b),
         stage_3b=ifelse(!is.na(stage_3b) & !is.na(stage_5) & stage_3b>stage_5, NA, stage_3b),
         stage_4=ifelse(!is.na(stage_4) & !is.na(stage_5) & stage_4>stage_5, NA, stage_4)) %>%
  analysis$cached("ckd_stages_from_algorithm_interim_4", unique_indexes="patid")

ckd_stages_from_algorithm %>% count()        
#8,466,065

######################################################################################

# Save the completed staging table

ckd_stages <- ckd_stages_from_algorithm %>%
  analysis$cached(
    name = "stages",
    unique_indexes = "patid"
  )

# Inspect the result

ckd_stages %>% count()
ckd_stages %>% head()

