
# perform logistic regression ---------------------------------------------
library(tidyverse)

df <- read_csv("/data/foci/01_foci_quant.csv") %>%
  filter(time == 20) %>%
  mutate(dose = dose %>% str_remove_all("dose_") %>% as.numeric())

# subset only 0 and 10
doselist_20_subset <- df %>% filter(dose == 10 | dose == 0)

doselist_20_subset <- doselist_20_subset %>% mutate(state = case_when(dose == 0 ~ "inactivated",
                                                                      dose == 10 ~ "activated")) %>% mutate(state = as.factor(state))

doselist_20_subset$state <- doselist_20_subset$state %>% relevel("inactivated")

# perform logistic regression ---------------------------------------------

logistic <- glm(state ~ foci, data = doselist_20_subset, family = "binomial")

# predict for other dose points
doselist_20$probs <- predict(logistic,
        newdata = doselist_20,
        type = "response")

doselist_20 <- doselist_20 %>% mutate(state = case_when(probs >= 0.5 ~ "activated",
                                         probs < 0.5 ~ "inactivated"))

# get percentage of activated cells
results_df <-doselist_20 %>%
  group_by(dose) %>%
  select(state) %>%
  table() %>%
  as.data.frame.matrix() %>%
  mutate(percentage_activated = activated/(activated+inactivated)*100) %>% mutate(dose = doselist_20 %>% pull(dose) %>% unique())



