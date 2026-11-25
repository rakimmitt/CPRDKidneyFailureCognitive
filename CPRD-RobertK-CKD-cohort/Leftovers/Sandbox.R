
diabetes_2024 %>%
    filter(is_case == 1) %>%
    summarise(follow_up_years = mean(follow_up_years, na.rm = TRUE))

    diabetes_2024 %>%
    summarise(follow_up_years = mean(follow_up_years, na.rm = TRUE))



    names(diabetes_2024)[
  grepl("smok|hypert|htn", names(diabetes_2024), ignore.case = TRUE)
]

diabetes_2024 %>% count(gender)
diabetes_2024 %>% count(ethnicity_5cat)

###

column_name <- "preegfr""

diabetes_2024 %>%
filter(is_case == 0) %>%
  summarise(
    n_observed = sum(!is.na(.data[[column_name]])),
    n_missing = sum(is.na(.data[[column_name]])),
    mean = mean(.data[[column_name]], na.rm = TRUE),
    sd = sd(.data[[column_name]], na.rm = TRUE),
    median = median(.data[[column_name]], na.rm = TRUE),
    iqr = IQR(.data[[column_name]], na.rm = TRUE),
    max = max(.data[[column_name]], na.rm = TRUE),
    min = min(.data[[column_name]], na.rm = TRUE)
  )

###

diabetes_summary %>%
    filter(is_case == 0) %>%
    summarise(
        n_nockd = sum(is.na(index_date_ckd_stage)),
        n_ckd1 = sum(index_date_ckd_stage == 1),
        n_ckd2 = sum(index_date_ckd_stage == 2),
        n_ckd3 = sum(index_date_ckd_stage == 3),
        n_ckd4 = sum(index_date_ckd_stage == 4),
        n_ckd5 = sum(index_date_ckd_stage == 5)
    )

diabetes_summary %>%
  count(is_case, index_date_ckd_stage, name = "n") %>%
  group_by(is_case) %>%
  mutate(
    percentage = round(100 * n / sum(n), 1)
  ) %>%
  ungroup() %>%
  arrange(is_case, index_date_ckd_stage) %>%
  print(n = Inf)


diabetes_summary <- diabetes_summary %>%
  mutate(
    time_to_dementia = as.numeric(
      difftime(
        post_index_date_first_alldementia,
        index_date,
        units = "days"
      ) / 365.25
    )
  )


   diabetes_2024 %>%
    filter(is_case == 0) %>%
    summarise(dm_dur_all = mean(dm_dur_all, na.rm = TRUE))


dementia_prediction = glm(post_index_alldementia ~ index_date_age + is_case + follow_up_years + dm_dur_all, diabetes_summary, family='binomial')
summary(dementia_prediction)
sumdementia_prediction = summary(dementia_prediction)


#### Look at competing risk

competing_data <- diabetes_summary %>%
  mutate(
    across(
      c(index_date, earliest_all_dementia, death_date,
        gp_end_date, hes_end_date),
      as.Date
    ),

    # Require both follow-up end dates to be known.
    censor_date = pmin(gp_end_date, hes_end_date),

    # Treat these predictors as categorical.
    is_case = factor(is_case, levels = c(0, 1)),
    gender = factor(gender)
  ) %>%
  filter(
    !is.na(index_date),
    !is.na(censor_date),

    # Exclude prevalent dementia, including diagnoses on index date.
    is.na(earliest_all_dementia) |
      earliest_all_dementia > index_date,

    # Exclude deaths on/before index date.
    is.na(death_date) | death_date > index_date,

    censor_date > index_date
  ) %>%
  mutate(
    # censor_date is non-missing, so every row has an end date.
    end_date = pmin(
      earliest_all_dementia,
      death_date,
      censor_date,
      na.rm = TRUE
    ),

    event_status = case_when(
      !is.na(earliest_all_dementia) &
        earliest_all_dementia == end_date ~ 1L,

      !is.na(death_date) &
        death_date == end_date ~ 2L,

      TRUE ~ 0L
    ),

    follow_up_years = as.numeric(end_date - index_date) / 365.25
  )

# Explicit complete-case analysis for the four predictors.
model_data <- competing_data %>%
  filter(
    if_all(
      all_of(c("index_date_age", "is_case", "dm_dur_all", "gender")),
      ~ !is.na(.x)
    )
  )

# Check exclusions.
tibble(
  original_records = nrow(diabetes_summary),
  eligible_follow_up_records = nrow(competing_data),
  complete_case_records = nrow(model_data)
)

# Check outcomes within each group.
model_data %>%
  count(is_case, event_status, name = "n")

  dementia_cox <- coxph(
  Surv(follow_up_years, event_status == 1) ~
    index_date_age + is_case + dm_dur_all + gender,
  data = model_data,
  cluster = patid,
  robust = TRUE,
  ties = "efron",
  x = TRUE
)

summary(dementia_cox)

cox_summary <- summary(dementia_cox)

cox_results <- tibble(
  predictor = rownames(cox_summary$coefficients),
  hazard_ratio = cox_summary$conf.int[, "exp(coef)"],
  lower_95_CI = cox_summary$conf.int[, "lower .95"],
  upper_95_CI = cox_summary$conf.int[, "upper .95"],
  p_value = cox_summary$coefficients[, "Pr(>|z|)"]
)

cox_results %>%
  print(width = Inf)

ph_check <- cox.zph(dementia_cox)
print(ph_check)
plot(ph_check)


load("C:\\Users\\rk535\\OneDrive\\1 - PhD\\Data Science\\CPRD\\Github clone\\CPRDKidneyFailureCognitive\\CPRD-Codelists\\Medcodes\\exeter_medcodelist_tbi.txt")

View(medcodelist_tbi)

readr::write_tsv(medcodelist_tbi, "C:\\Users\\rk535\\OneDrive\\1 - PhD\\Data Science\\CPRD\\Github clone\\CPRDKidneyFailureCognitive\\CPRD-Codelists\\Medcodes\\exeter_medcodelist_tbi.txt")





library(tidyverse)

ggplot(mtcars, aes(x = wt, y = mpg, colour = factor(cyl))) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~cyl, scales = "free_x") +
  labs(title = "Scatter plot of MPG vs Weight",
       x = "Weight (1000 lbs)",
       y = "Miles per Gallon") +
  theme_minimal()

glimpse(penguins)

ggplot(penguins,
       aes(x = flipper_len, y = body_mass)) +
  geom_point(aes(colour = species)) +
  geom_smooth(method = "lm") +
  facet_wrap(~island, scales = "fixed")

ggplot(penguins,
       aes(x = body_mass, bins = 10)) +
       geom_histogram()

ggplot(penguins,
        aes(x = species, y = body_mass)) +  
        geom_boxplot()  

ggplot(penguins, 
       aes(x = fct_infreq(species))) +
       geom_bar()


penguins %>%
  filter(!is.na(sex)) %>%
  mutate(sex = fct_recode(sex, "Male" = "male", "Female" = "female")) %>%
ggplot() +
aes(x = species, fill = sex) +
       geom_bar(position = "dodge") +
       labs(title = "Count of penguins by species and sex",
            x = "Species",
            y = "Count",
            fill = "Sex") +
            theme_minimal()

penguins |>
     ggplot() +
     aes(x = flipper_len, y = body_mass, colour = species) +
     geom_point() +
     geom_smooth(method = "lm")


# [55] "qrisk2_smoking_cat"                               
# [56] "qrisk2_smoking_cat_uncoded"                       
# [57] "regstartdate"                                     
# [58] "smoking_cat" 

Sys.setenv(QUARTO_PATH = paths[file.exists(paths)][1])
quarto::quarto_version()
quarto::quarto_render("baseline_table.qmd")

dementia_concordance <- with(
  combined_2024,
  table(
    `Post-index dementia code` = factor(
      !is.na(post_index_date_first_alldementia),
      levels = c(FALSE, TRUE),
      labels = c("Absent", "Present")
    ),
    `Post-index dementia medication` = factor(
      !is.na(post_index_date_first_dementia_medications),
      levels = c(FALSE, TRUE),
      labels = c("Absent", "Present")
    )
  )
)

dementia_concordance