# code and analysis for Dale Fisher textbook chapter on field epidemiology
# cholera scenario



# Load packages -----------------------------------------------------------
pacman::p_load(
     sitrep,
     rio,
     here,
     janitor,
     lubridate,
     randomNames,
     flextable,
     tidyverse
)


# Generate cholera data ---------------------------------------------------
chol_raw <- sitrep::gen_data("cholera")


# Prepare cholera data for textbook use -----------------------------------
chol <- chol_raw %>% 
     # clean names
     clean_names() %>% 
     
     # rename manual
     rename(
          case_id = case_number,
          date_admit = date_of_consultation_admission,
          date_onset = date_of_onset,
          town = patient_origin,
     ) %>% 
     
     # add patient names
     mutate(name = randomNames::randomNames(n = nrow(.),
                                            gender = sex)) %>% 
     # add age groups
     mutate(age_grp = age_categories(age_years,
                                       breakers = c(0,5,10,15,20,30,40,50,60))) %>% 
     
     # recode unknown gender to NA
     mutate(sex = na_if(sex, "U")) %>% 
     
     # recode towns to showcase "waves" of epidemic by town
     mutate(town = case_when(
          # first cases reported in C
          date_onset < ymd("20180120") ~ "C",
          
          # then cases reported from other jurisdictions
          date_onset < ymd("20180215") ~ sample(c("A", "B", "C"), nrow(.), replace = TRUE, prob = c(.1, .2, .7)),
          
          # area C declining
          date_onset < ymd("20180315") ~ sample(c("A", "B", "C"), nrow(.), replace = TRUE, prob = c(.3, .3, .4)),
          
          # area B most
          date_onset < ymd("20180415") ~ sample(c("A", "B", "C"), nrow(.), replace = TRUE, prob = c(.3, .6, .1)),

          # area A increasing
          date_onset < ymd("20180515") ~ sample(c("A", "B", "C"), nrow(.), replace = TRUE, prob = c(.6, .3, .14)),
          
          TRUE ~ "A")) %>% 
     
     # re-code case_id
     mutate(case_id = )
          
     # recode pregnant away from MSF acronyms      
     mutate(
          pregnant = na_if(pregnant, "NA"),
          pregnant = recode(pregnant,
                 "Y" = "Yes",
                 "N" = "No",
                 "W" = "No")) %>% 
     
     # assign new symptoms/risk factor columns based on probabilities
     mutate(
          # Most vaccinated are towns A and B. Vax data has more NA than other columns
          prior_vax = case_when(
               town == "A" ~ sample(c("Yes", "No", NA), nrow(.), replace = TRUE, prob = c(.4, .4, .2)),
               town == "B" ~ sample(c("Yes", "No", NA), nrow(.), replace = TRUE, prob = c(.7, .2, .1)),
               town == "C" ~ sample(c("Yes", "No", NA), nrow(.), replace = TRUE, prob = c(.2, .7, .1)),
          ),
          # worst delays to care and dehydration seen from town C
          dehy_at_admit = case_when(
               town == "A" ~ sample(c("Severe", "Moderate", NA), nrow(.), replace = TRUE, prob = c(.4, .4, .05)),
               town == "B" ~ sample(c("Severe", "Moderate", NA), nrow(.), replace = TRUE, prob = c(.3, .6, .05)),
               town == "C" ~ sample(c("Severe", "Moderate", NA), nrow(.), replace = TRUE, prob = c(.8, .1, .05)),
          ),
          # Very few people have bloody stools generally
          stool_bloody = case_when(
               town == "A" ~ sample(c("Yes", "No", NA), nrow(.), replace = TRUE, prob = c(.1, .8, .05)),
               town == "B" ~ sample(c("Yes", "No", NA), nrow(.), replace = TRUE, prob = c(.2, .7, .05)),
               town == "C" ~ sample(c("Yes", "No", NA), nrow(.), replace = TRUE, prob = c(.1, .8, .05)),
          )
     ) %>% 
     
     # re-calculate risk of death, based on risk and protective factors
     mutate(
          # risk factors
          death_chance = 0,
          death_chance = ifelse(age_grp %in% c("0-4", "5-10", "60+"), death_chance + 0.15, death_chance), # vulnerable ages
          death_chance = ifelse(dehy_at_admit %in% c("Severe"), death_chance + 0.40, death_chance), # dehydration
          death_chance = ifelse(dehy_at_admit %in% c("Moderate"), death_chance + 0.20, death_chance), # dehydration
          death_chance = ifelse(!is.na(date_onset) & date_onset < ymd("20180301"), death_chance + 0.40, death_chance),  # first wave
          #death_chance = ifelse(!is.na(loose_stool) & loose_stool == "Yes", death_chance + 0.20, death_chance), # diarrhea
          death_chance = ifelse(!is.na(town) & town == "C", death_chance + 0.20, death_chance), # remote, poor infrastructure, epicentre
          death_chance = ifelse(!is.na(pregnant) & pregnant == "Yes", death_chance + 0.15, death_chance), # complications
          death_chance = ifelse(!is.na(oedema) & oedema == "Yes", death_chance + 0.20, death_chance),

          # protective factors
          death_chance = ifelse(!is.na(prior_vax) & prior_vax == "Yes", death_chance - 0.80, death_chance),    # prior vaccination
          death_chance = ifelse(!is.na(stool_bloody) & stool_bloody == "Yes", death_chance - 0.20, death_chance), # likely not cholera
          death_chance = ifelse(!is.na(date_onset) & date_onset > ymd("20180415"), death_chance - 0.40, death_chance),  # ending protection
          
          ) %>%
     
     # normalise death chance to between 1 and 0
     mutate(death_norm = (death_chance-min(death_chance))/(max(death_chance)-min(death_chance))) %>% 
     
     # create binary death column
     mutate(outcome = ifelse(death_norm > 0.7, "Death", "Alive")) %>% 
     
     # re-arrange columns
     select(case_id,
            #name,
            town, sex, age_years, age_grp, date_onset,
            #date_admit,  # dates of admit appear to be before dates of onset?
            prior_vax, dehy_at_admit, oedema, pregnant, stool_bloody,
            outcome)
     



# Outputs -----------------------------------------------------------------

# head - first 5 cases
chol %>% 
     arrange(date_onset) %>% 
     head() %>% 
     qflextable() %>% 
     fontsize(size = 8, part = "all") %>%
     line_spacing(space = 0.5) %>% 
     theme_box() %>% 
     set_caption(caption = "caption")

# Epicurve to show waves of epidemic
ggplot(data = chol, mapping = aes(x = date_onset, fill = town))+
     geom_histogram()

# Epicurve showing outcomes improving as response infrastructure improves
ggplot(data = chol, mapping = aes(x = date_onset, fill = outcome))+
     geom_histogram()

# show CFR declining over time
chol %>% 
     group_by(week = floor_date(date_onset, "week")) %>% 
     summarise(CFR = sum(outcome == "Death", na.rm=T) / n()) %>% 
     ggplot(mapping = aes(x = week, y = CFR))+
     geom_line(size = 2)+
          scale_y_continuous(labels = scales::percent)+
          scale_x_date(labels = scales::label_date_short())+
          theme_minimal(base_size = 14)

