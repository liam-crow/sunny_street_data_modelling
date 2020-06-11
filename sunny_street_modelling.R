
#### Initial set up ####
library(readr)
library(snakecase)
library(dplyr)
library(ggplot2)
library(ggthemes)

theme_set(theme_minimal())

#### Load in data ####

# import excel file
patients <- read_csv("Sunny Street Patients.csv", guess_max = 4000)
# clean up names, lowercase and underscores
names(patients) <- to_snake_case(names(patients))

# unwanted columns, and filter out bad data
patients_select <- patients %>% 
    select(-contains("_active_ind"), -contains("_date")) %>% 
    filter(age != 999) #removes 3% of the data points

# clear up NA's
patients_select[is.na(patients_select)] <- 0

# remove columns that only have 1 unique value
check_unique <- function(data) {length(unique(data)) != 1}
patients_unique <- Filter(check_unique, patients_select)

dim(patients)
dim(patients_unique)
#removed 54 rows and 368 columns


ggplot(patients_unique, aes(x = age, after_stat(count), fill = alcohol_status_name)) + 
    geom_density(alpha = 0.5) + 
    ggtitle("Density of Ages by Drinking Status") +
    xlab("Age") + ylab("") +
    labs(fill = "Alcohol Status")

patients_unique %>% 
    group_by(alcohol_status_name) %>% count()

ggplot(patients_unique, aes(x = age, after_stat(count), fill = gender_name)) + 
    geom_density(alpha = 0.5) + 
    # geom_rug() +
    ggtitle("Density of Ages by Gender") +
    xlab("Age") + ylab("") +
    labs(fill = "Gender")

patients_unique %>% 
    group_by(gender_name) %>% count()

ggplot(patients_unique, aes(x = age, after_stat(count), fill = as.character(diag_chronic_disease_ind))) + 
    geom_density(alpha = 0.5, position = 'fill') + 
    # geom_rug() +
    ggtitle("Proportion of Chronic Disease Diag by Age") +
    xlab("Age") + ylab("") +
    labs(fill = "Chronic Illness")

#### select the most common diagnoses 
min_ratio <- function(data) {
    table_data <- table(data)
    min(table_data)/sum(table_data)*100
}
column_data <- sapply(patients_unique, min_ratio)

most_common_diag <- tibble::enframe(column_data) %>% 
    arrange(-value) %>% 
    # mutate(row_no = row_number())
    filter(
        grepl('diag', name)
    ) %>% 
    filter(
        row_number() %in% c(1:10)
    ) %>% pull(name)


patients_unique %>% 
    group_by(indigenous_status_name) %>% 
    count()

# combine Indigenous populations into the same group,
patients_model <- patients_unique %>% 
    select(
        quarter_name, patient_id, patient_link_id, alcohol_status_name, alcohol_consumption_frequency_name,
        alcohol_consumption_amount, allergy_count, age, gender_name, indigenous_status_name, sexuality_name,
        pensioner_ind, most_common_diag
    ) %>% 
    mutate(
        pensioner_ind = as.character(pensioner_ind),
        indigenous_status_name = gsub("Torres Strait Islander but not Aboriginal origin", "Indigenous", indigenous_status_name),
        indigenous_status_name = gsub("Aboriginal but not Torres Strait Islander origin", "Indigenous", indigenous_status_name),
        indigenous_status_name = gsub("Both Aboriginal and Torres Strait Islander origin", "Indigenous", indigenous_status_name)
    )

library(brms)

sapply(patients_model %>% select(most_common_diag), min_ratio)

formula_text <- " ~ age + gender_name + indigenous_status_name + quarter_name + alcohol_consumption_frequency_name + 
        + allergy_count + pensioner_ind"

# "diag_chronic_disease_ind" 
# "diag_mental_health_ind"   
# "diag_depression_ind"     
# "diag_anxiety_ind"         
# "diag_drug_abuse_ind"      
# "diag_respiratory_ind"    
# "diag_asthma_ind"          
# "diag_cardiovascular_ind"  
# "diag_musculoskeletal_ind"
# "diag_schizophrenia_ind"  

#### diag_chronic_disease_ind ####

null_f <- bf(paste0('diag_chronic_disease_ind', formula_text))

null_model <- brm(
    formula = null_f,
    data    = patients_model,
    family  = "bernoulli", #logistic
    iter    = 1000, #mcmc, how many times do we run it
    warmup  = 300,
    cores   = 4, #how many cores to run on?
    chains  = 4, #how many chains do we want to run? 
                 #(normally make it equal to cores, running parallel)
    control = list(max_treedepth = 10, adapt_delta = 0.8) #increase to improve acc of model
                                                          #sig increase in time taken
)

library(rstan)
library(bayesplot)

# show histograms of the posterior distributions
mcmc_plot(null_model, type = "hist")

# plot some diagnostics of the sampler
mcmc_plot(null_model, type = "neff")
mcmc_plot(null_model, type = "rhat")

# plot some diagnostics specific to the NUTS sampler
mcmc_plot(null_model, type = "nuts_acceptance")
mcmc_plot(null_model, type = "nuts_divergence")

plot(null_model)

posterior <- extract(null_model, inc_warmup = TRUE, permuted = FALSE)
ppc_dens_overlay(null_model)

color_scheme_set("mix-blue-pink")
p <- mcmc_trace(posterior,  pars = c("mu", "tau"), n_warmup = 300,
                facet_args = list(nrow = 2, labeller = label_parsed))

mcmc_plot(null_model)
plot(conditional_effects(null_model)) 
plot(conditional_effects(null_model, effects = "age:gender_name"))
plot(conditional_effects(null_model, effects = "age:indigenous_status_name"))
plot(conditional_effects(null_model, effects = "gender_name:indigenous_status_name"))

#### diag_mental_health_ind ####

null_f_2 <- bf(paste0('diag_mental_health_ind', formula_text))

null_model_2 <- brm(
    formula = null_f_2,
    data    = patients_model,
    family  = "bernoulli",
    iter    = 1000,
    cores   = 4,
    chains  = 4,
    control = list(max_treedepth = 10, adapt_delta = 0.8)
)

mcmc_plot(null_model_2)
plot(conditional_effects(null_model_2)) 
plot(conditional_effects(null_model_2, effects = "age:gender_name"))
plot(conditional_effects(null_model_2, effects = "age:indigenous_status_name"))
plot(conditional_effects(null_model_2, effects = "gender_name:indigenous_status_name"))

#### diag_mental_health_ind ####

null_f_3 <- bf(paste0('diag_depression_ind', formula_text))

null_model_3 <- brm(
    formula = null_f_3,
    data    = patients_model,
    family  = "bernoulli",
    iter    = 1000,
    cores   = 4,
    chains  = 4,
    control = list(max_treedepth = 10, adapt_delta = 0.8)
)

mcmc_plot(null_model_3)
plot(conditional_effects(null_model_3)) 
plot(conditional_effects(null_model_3, effects = "age:gender_name"))
plot(conditional_effects(null_model_3, effects = "age:indigenous_status_name"))
plot(conditional_effects(null_model_3, effects = "gender_name:indigenous_status_name"))

#### diag_drug_abuse_ind ####

null_f_4 <- bf(paste0('diag_drug_abuse_ind', formula_text))

null_model_4 <- brm(
    formula = null_f_4,
    data    = patients_model,
    family  = "bernoulli",
    iter    = 1500,
    cores   = 4,
    chains  = 4,
    control = list(max_treedepth = 10, adapt_delta = 0.9)
)

mcmc_plot(null_model_4)
mcmc_trace(null_model_4)
plot(conditional_effects(null_model_4)) 
plot(conditional_effects(null_model_4, effects = "age:gender_name"))
plot(conditional_effects(null_model_4, effects = "age:indigenous_status_name"))
plot(conditional_effects(null_model_4, effects = "gender_name:indigenous_status_name"))


