#Setup
library(tidyverse)
library(aurum)
library(EHRBiomarkr)
library(tidyverse)
library(MatchIt)
rm(list=ls())

cprd = CPRDData$new(cprdEnv = "diabetes-jun2024",cprdConf = "C:\\Users\\rk535\\OneDrive\\1 - PhD\\Data Science\\CPRD\\.aurum.yaml")
codesets = cprd$codesets()
codes = codesets$getAllCodeSetVersion(v = "01/06/2024")

analysis = cprd$analysis("rk_ckd")

# Settings
max_controls <- 4L
max_age_gap_years <- 50
set.seed(123)

# 1. Load the required columns; both groups must have HES linkage
columns <- c(
  "patid", "pracid", "dob", "gender", "regstartdate",
  "gp_end_date", "hes_end_date", "with_hes"
)

cases <- advanced_ckd_cohort %>%
  filter(with_hes == 1) %>%
  select(all_of(c(columns, "index_date"))) %>%
  collect()

controls <- non_ckd_cohort %>%
  filter(with_hes == 1) %>%
  select(all_of(columns)) %>%
  collect()

# Preserve patient IDs as text and ensure dates have Date class
prepare <- function(x) {
  x %>%
    mutate(
      across(c(patid, pracid, gender), as.character),
      across(
        any_of(c("dob", "regstartdate", "gp_end_date",
                 "hes_end_date", "index_date")),
        as.Date
      )
    ) %>%
    arrange(patid)
}

cases <- prepare(cases)
controls <- prepare(controls)

# Each patient must appear once, and cases cannot also be controls
stopifnot(
  !anyNA(cases$patid), !anyNA(controls$patid),
  !anyDuplicated(cases$patid), !anyDuplicated(controls$patid),
  !any(cases$patid %in% controls$patid)
)

# Controls with missing matching/eligibility information cannot be matched
controls <- controls %>%
  drop_na(patid, pracid, dob, gender,
          regstartdate, gp_end_date, hes_end_date)

empty_pairs <- tibble(
  matched_case_patid = character(),
  patid = character()
)

# 2. Match within one practice
match_practice <- function(ca, co) {

  # Cases with missing information remain in the final case table,
  # but cannot receive controls.
  ca <- ca %>% drop_na(pracid, dob, gender, index_date)

  if (nrow(ca) == 0L || is.null(co) || nrow(co) == 0L) {
    return(empty_pairs)
  }

  # Rows = cases; columns = controls.
  # All dates are inclusive. Age gap uses 365.25 days per year.
  eligible <-
    outer(as.numeric(ca$index_date),
          as.numeric(co$regstartdate), `>=`) &
    outer(as.numeric(ca$index_date),
          as.numeric(co$gp_end_date), `<=`) &
    outer(as.numeric(ca$index_date),
          as.numeric(co$hes_end_date), `<=`) &
    abs(outer(as.numeric(ca$dob), as.numeric(co$dob), `-`)) <=
      max_age_gap_years * 365.25

  # Remove patients with no possible pairing from this matching call only
  keep_cases <- rowSums(eligible) > 0L
  ca <- ca[keep_cases, , drop = FALSE]
  eligible <- eligible[keep_cases, , drop = FALSE]

  if (nrow(ca) == 0L) return(empty_pairs)

  keep_controls <- colSums(eligible) > 0L
  co <- co[keep_controls, , drop = FALSE]
  eligible <- eligible[, keep_controls, drop = FALSE]

  dat <- bind_rows(
    mutate(ca, is_case = 1L),
    mutate(co, is_case = 0L)
  ) %>%
    mutate(
      dob_days = as.numeric(dob),
      gender = factor(gender)
    ) %>%
    as.data.frame()

  rownames(dat) <- dat$patid

  # Omit variables that are constant within this practice
  variables <- c("dob_days", "gender")
  variables <- variables[
    vapply(dat[variables], function(x) n_distinct(x) > 1L, logical(1))
  ]

  formula <- if (length(variables)) {
    reformulate(variables, response = "is_case")
  } else {
    is_case ~ 1
  }

  # A single possible pair needs no distance ranking
  distances <- if (nrow(ca) == 1L && nrow(co) == 1L) {
    matrix(0, 1L, 1L)
  } else if (length(variables)) {
    MatchIt::mahalanobis_dist(formula, data = dat)
  } else {
    matrix(0, nrow(ca), nrow(co))
  }

  dimnames(distances) <- list(ca$patid, co$patid)
  stopifnot(all(is.finite(distances)))

  # Prohibit pairs that fail the date or age-gap requirements
  distances[!eligible] <- Inf

  fit <- matchit(
    formula,
    data = dat,
    method = "nearest",
    distance = distances,
    ratio = min(max_controls, nrow(co)),
    replace = FALSE,
    m.order = "random"
  )

  # Extract actual matches; NA entries are unfilled control slots
  mm <- fit$match.matrix

  tibble(
    matched_case_patid = rep(rownames(mm), each = ncol(mm)),
    patid = as.vector(t(mm))
  ) %>%
    filter(!is.na(patid))
}

# 3. Run separately within each practice: this guarantees exact pracid
case_groups <- split(cases, cases$pracid)
control_groups <- split(controls, controls$pracid)

pairs <- lapply(names(case_groups), function(practice) {
  match_practice(case_groups[[practice]], control_groups[[practice]])
}) %>%
  bind_rows(empty_pairs)

# 4. Create the matched-control table with the requested case linkage
matched_controls <- pairs %>%
  left_join(controls, by = "patid") %>%
  left_join(
    cases %>%
      transmute(
        matched_case_patid = patid,
        matched_case_index_date = index_date,
        index_date = index_date
      ),
    by = "matched_case_patid"
  )

# Retain ALL HES-linked cases, including cases with zero controls
matched_cases <- cases %>%
  left_join(
    pairs %>% count(matched_case_patid, name = "n_controls"),
    by = c("patid" = "matched_case_patid")
  ) %>%
  mutate(n_controls = coalesce(n_controls, 0L))

# Optional combined table: case rows carry their own ID and index date
matched_cohort <- bind_rows(
  matched_cases %>%
    mutate(
      is_case = 1L,
      matched_case_patid = patid,
      matched_case_index_date = index_date
    ),
  matched_controls %>% mutate(is_case = 0L)
)

analysis = cprd$analysis("rk_ckd")

matched_cohort %>%
  select(patid, is_case, matched_case_patid, matched_case_index_date,
         index_date, dob, gender, pracid, regstartdate, gp_end_date,
         hes_end_date) %>%
  analysis$cached("matched_cohort", unique_indexes="patid",
                  indexes=c("is_case", "matched_case_patid", "index_date"))

# 5. Check control reuse and display the number of controls per case
stopifnot(
  !anyDuplicated(matched_controls$patid),
  all(matched_cases$n_controls <= max_controls),
  nrow(matched_cases) == nrow(cases)
)

matched_cases %>%
  count(n_controls, name = "number_of_cases") %>%
  print(n = Inf)

