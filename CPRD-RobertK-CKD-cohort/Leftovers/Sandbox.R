
diabetes_2024 %>%
    filter(is_case == 0) %>%
    summarise(follow_up_years)