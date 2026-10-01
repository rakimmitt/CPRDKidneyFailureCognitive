#Setup
library(tidyverse)
library(aurum)
library(EHRBiomarkr)
library(tidyverse)
library(MatchIt)
rm(list=ls())

cprd = CPRDData$new(cprdEnv = "nondiabetes-jun2024",cprdConf = "C:\\Users\\rk535\\OneDrive\\1 - PhD\\Data Science\\CPRD\\.aurum.yaml")
codesets = cprd$codesets()
codes = codesets$getAllCodeSetVersion(v = "01/06/2024")

analysis = cprd$analysis("rk_ckd")

advanced_ckd_cohort <- advanced_ckd_cohort %>%
  analysis$cached("advanced_ckd_cohort")
  advanced_ckd_cohort %>% count()


  non_ckd_cohort <- non_ckd_cohort %>%
  analysis$cached("non_ckd_cohort")
  non_ckd_cohort %>% count()

# Settings
max_controls <- 4L
max_age_gap_years <- 5
set.seed(123)

# 1. Load the required columns; both groups must have HES linkage
columns <- c(
  "patid", "pracid", "dob", "gender", "regstartdate",
  "gp_end_date", "hes_end_date", "with_hes", "ethnicity_5cat", "imd_decile"
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

    # Include missing ethnicity as an explicit category.
    ethnicity_5cat = coalesce(
      as.character(ethnicity_5cat), "Missing"
    ),

    # Preserve the original IMD, including NA.
    # na_if() also handles "Missing" if previously assigned.
    imd_decile = as.numeric(
    na_if(as.character(imd_decile), "Missing")
   ),

    # Additional variables used only for matching.
    imd_missing = as.integer(is.na(imd_decile)),
    imd_for_matching = coalesce(imd_decile, 5.5),

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
      gender = factor(gender),
      ethnicity_5cat = factor(ethnicity_5cat)
    ) %>%
    as.data.frame()

  rownames(dat) <- dat$patid

  # Omit variables that are constant within this practice
  variables <- c("dob_days", "gender", "ethnicity_5cat", "imd_for_matching", "imd_missing")
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
         index_date, dob, gender, ethnicity_5cat, imd_decile, pracid, regstartdate, gp_end_date,
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





########

########

# Checking quality of matches

# Combine matched cases and controls
# NB because of the variable ratio matching, there will be imbalances in overall age and gender distributions
# Pair-wise matching should be tighter, and we can summarise weighted matching characteristics

comparison <- bind_rows(
  matched_cases %>%
    filter(n_controls > 0) %>%
    transmute(
      group = "Cases",
      dob, gender,
      index_date = index_date
    ),

  matched_controls %>%
    transmute(
      group = "Controls",
      dob, gender,
      index_date = matched_case_index_date
    )
) %>%
  mutate(
    age_at_index = as.numeric(index_date - dob) / 365.25
  )

# Summarise age at index
comparison %>%
  group_by(group) %>%
  summarise(
    n = n(),
    missing_age = sum(is.na(age_at_index)),
    mean_age = mean(age_at_index, na.rm = TRUE),
    sd_age = sd(age_at_index, na.rm = TRUE),
    median_age = median(age_at_index, na.rm = TRUE),
    youngest = min(age_at_index, na.rm = TRUE),
    oldest = max(age_at_index, na.rm = TRUE)
  ) %>%
  print(width = Inf)

# Gender counts and percentages, including missing values
comparison %>%
  count(group, gender) %>%
  group_by(group) %>%
  mutate(percent = round(100 * n / sum(n), 1)) %>%
  ungroup() %>%
  print(n = Inf)

# Pair-wise comparison

pair_comparison <- matched_controls %>%
  select(patid, matched_case_patid, dob, gender) %>%
  left_join(
    matched_cases %>%
      select(
        matched_case_patid = patid,
        case_dob = dob,
        case_gender = gender
      ),
    by = "matched_case_patid"
  ) %>%
  mutate(
    age_gap_years = abs(as.numeric(dob - case_dob)) / 365.25,
    same_gender = gender == case_gender
  )

pair_comparison %>%
  summarise(
    n_pairs = n(),
    mean_age_gap = mean(age_gap_years, na.rm = TRUE),
    median_age_gap = median(age_gap_years, na.rm = TRUE),
    p95_age_gap = quantile(age_gap_years, 0.95, na.rm = TRUE),
    largest_age_gap = max(age_gap_years, na.rm = TRUE),
    percent_same_gender = 100 * mean(same_gender, na.rm = TRUE)
  ) %>%
  print(width = Inf)

# Weighted comparison matching

comparison_weighted <- bind_rows(
  matched_cases %>%
    filter(n_controls > 0) %>%
    transmute(
      group = "Cases",
      age_at_index = as.numeric(
        as.Date(index_date) - as.Date(dob)
      ) / 365.25,
      weight = 1
    ),

  matched_controls %>%
    left_join(
      matched_cases %>%
        select(matched_case_patid = patid, n_controls),
      by = "matched_case_patid"
    ) %>%
    transmute(
      group = "Controls",
      age_at_index = as.numeric(
        as.Date(matched_case_index_date) - as.Date(dob)
      ) / 365.25,
      weight = 1 / n_controls
    )
)

comparison_weighted %>%
  group_by(group) %>%
  summarise(
    n = n(),
    mean_age = weighted.mean(age_at_index, weight, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print(width = Inf)

## Weighted SMD for age

balance_data <- comparison_weighted %>%
  mutate(is_case = as.integer(group == "Cases"))

age_balance <- cobalt::bal.tab(
  is_case ~ age_at_index,
  data = balance_data,
  weights = balance_data$weight,
  method = "weighting",
  estimand = "ATT",
  s.d.denom = "treated",   # Cases are the "treated" group here
  continuous = "std",
  un = TRUE,
  thresholds = c(m = 0.1)
)

print(age_balance)

# NB diff. adj is the weighted SMD, which is the relevant metric for this variable. The unweighted SMD is also reported for reference.




  # Weight gender matching

  gender_comparison <- bind_rows(
  matched_cases %>%
    filter(n_controls > 0) %>%
    transmute(
      group = "Cases",
      gender = as.character(gender),
      weight = 1
    ),

  matched_controls %>%
    left_join(
      matched_cases %>%
        select(matched_case_patid = patid, n_controls),
      by = "matched_case_patid"
    ) %>%
    transmute(
      group = "Controls",
      gender = as.character(gender),
      weight = 1 / n_controls
    )
)

# Weighted percentages within each group.

# Missing gender, if present, is displayed as its own category.
gender_comparison %>%
  group_by(group, gender) %>%
  summarise(
    n = n(),
    weighted_n = sum(weight),
    .groups = "drop"
  ) %>%
  group_by(group) %>%
  mutate(weighted_percent = 100 * weighted_n / sum(weighted_n)) %>%
  ungroup() %>%
  print(n = Inf, width = Inf)


# Weighted gender SMD

unique(gender_comparison$gender)

gender_category <- "1"

gender_balance <- gender_comparison %>%
  filter(!is.na(gender)) %>%
  group_by(group) %>%
  summarise(
    proportion = weighted.mean(gender == gender_category, weight),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = group, values_from = proportion) %>%
  transmute(
    case_percent = 100 * Cases,
    control_percent = 100 * Controls,
    difference_percentage_points = 100 * (Cases - Controls),
    weighted_smd = if_else(
      Cases > 0 & Cases < 1,
      (Cases - Controls) / sqrt(Cases * (1 - Cases)),
      NA_real_
    ),
    absolute_smd = abs(weighted_smd)
  )

print(gender_balance, width = Inf)