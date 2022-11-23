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
     tsibble,
     tidyverse
)

pacman::p_load(rio,          # File import
               here,         # File locator
               tidyverse,    # data management + ggplot2 graphics
               tsibble,      # handle time series datasets
               slider,       # for calculating moving averages
               imputeTS,     # for filling in missing values
               feasts,       # for time series decomposition and autocorrelation
               forecast,     # fit sin and cosin terms to data (note: must load after feasts)
               trending,     # fit and assess models 
               tmaptools,    # for getting geocoordinates (lon/lat) based on place names
               ecmwfr,       # for interacting with copernicus sateliate CDS API
               stars,        # for reading in .nc (climate data) files
               units,        # for defining units of measurement (climate data)
               yardstick,    # for looking at model accuracy
               surveillance  # for aberration detection
)


# Generate cholera data ---------------------------------------------------
chol_raw <- sitrep::gen_data("cholera")



# Prepare cholera data for textbook use -----------------------------------

# Generate fake ages in realistic distribution
dat <- tribble(
     ~min, ~max, ~prop,
     0, 5, 0.2,
     5, 10, 0.2,
     10, 20, 0.2,
     20, 30, 0.2,
     30, 40, 0.1,
     40, 50, 0.1,
     50, 75, 0.05)
rows <- sample(nrow(dat), 300, replace=TRUE, prob=dat$prop)
ages <- round(dat$min[rows] + runif(300) * (dat$max[rows] - dat$min[rows]))
hist(ages)

# clean dataset
chol <- chol_raw %>% 
     # clean names
     clean_names() %>% 
     
     # rename manual
     rename(
          case_id = case_number,
          date_admit = date_of_consultation_admission,
          date_onset = date_of_onset,
          town = patient_origin,
          preg = pregnant,
          rdt = cholera_rdt_result
     ) %>% 
     
     # add patient names
     mutate(name = randomNames::randomNames(n = nrow(.),
                                            gender = sex)) %>% 
     
     # adjust ages
     mutate(age_years = ages) %>% 
     
     # add age groups
     mutate(age_grp = age_categories(age_years,
                                       breakers = c(0,5,10,15,20,30,40,50,60))) %>% 
     
     # recode unknown gender to NA
     mutate(sex = na_if(sex, "U")) %>% 
     
     # recode rdt result
     mutate(rdt = case_when(
          rdt == "ND" ~ NA_character_,
          rdt == "N"  ~ "Neg",
          rdt == "P2" ~ "Neg",
          rdt == "I"  ~ "Neg",
          TRUE        ~ "Pos"
          
     )) %>% 
     
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
     #mutate(case_id = )
          
     # recode pregnant away from MSF acronyms      
     mutate(
          preg = na_if(preg, "NA"),
          preg = recode(preg,
                 "Y" = "Yes",
                 "N" = "No",
                 "W" = "No")) %>% 
     
     # recode sex away from MSF acronyms      
     mutate(sex = as.character(sex), 
            sex = na_if(sex, " U"),
            sex = na_if(sex, "NA")) %>% 
     
     # assign new symptoms/risk factor columns based on probabilities
     mutate(
          # Most vaccinated are towns A and B. Vax data has more NA than other columns
          prior_vax = case_when(
               town == "A" ~ sample(c("Yes", "No", NA), nrow(.), replace = TRUE, prob = c(.4, .4, .2)),
               town == "B" ~ sample(c("Yes", "No", NA), nrow(.), replace = TRUE, prob = c(.7, .2, .1)),
               town == "C" ~ sample(c("Yes", "No", NA), nrow(.), replace = TRUE, prob = c(.2, .7, .1)),
          ),
          # worst delays to care and dehydration seen from town C
          dehy = case_when(
               town == "A" ~ sample(c("Severe", "Moderate", NA), nrow(.), replace = TRUE, prob = c(.4, .4, .05)),
               town == "B" ~ sample(c("Severe", "Moderate", NA), nrow(.), replace = TRUE, prob = c(.3, .6, .05)),
               town == "C" ~ sample(c("Severe", "Moderate", NA), nrow(.), replace = TRUE, prob = c(.8, .1, .05)),
          ),
          # Very few people have bloody stools generally
          blood_stool = case_when(
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
          death_chance = ifelse(dehy %in% c("Severe"), death_chance + 0.40, death_chance), # dehydration
          death_chance = ifelse(dehy %in% c("Moderate"), death_chance + 0.20, death_chance), # dehydration
          death_chance = ifelse(!is.na(date_onset) & date_onset < ymd("20180301"), death_chance + 0.40, death_chance),  # first wave
          #death_chance = ifelse(!is.na(loose_stool) & loose_stool == "Yes", death_chance + 0.20, death_chance), # diarrhea
          death_chance = ifelse(!is.na(town) & town == "C", death_chance + 0.20, death_chance), # remote, poor infrastructure, epicentre
          death_chance = ifelse(!is.na(preg) & preg == "Yes", death_chance + 0.15, death_chance), # complications
          death_chance = ifelse(!is.na(oedema) & oedema == "Yes", death_chance + 0.20, death_chance),

          # protective factors
          death_chance = ifelse(!is.na(prior_vax) & prior_vax == "Yes", death_chance - 0.80, death_chance),    # prior vaccination
          death_chance = ifelse(!is.na(blood_stool) & blood_stool == "Yes", death_chance - 0.20, death_chance), # likely not cholera
          death_chance = ifelse(!is.na(date_onset) & date_onset > ymd("20180415"), death_chance - 0.40, death_chance),  # ending protection
          
          ) %>%
     
     # normalise death chance to between 1 and 0
     mutate(death_norm = (death_chance-min(death_chance))/(max(death_chance)-min(death_chance))) %>% 
     
     # create binary death column
     mutate(outcome = ifelse(death_norm > 0.7, "Death", "Alive")) %>% 
     
     # re-arrange columns
     select(case_id,
            name,
            town, sex, age_years, 
            age_grp,  # for use later
            date_onset,
            date_admit,  # dates of admit appear to be before dates of onset?
            prior_vax, dehy, oedema, preg, blood_stool,
            rdt,
            outcome)
     


# Generate historical data ------------------------------------------------

counts <- data.frame(
     date = seq.Date(ymd("20130101"), ymd("20180401"), "day"),
     cases = sample(0:10, 1917, replace = TRUE)) %>% 
     mutate(month = month(date)) %>% 
     mutate(cases = ifelse(month %in% c(1,2,3), cases * 1.5, cases)) %>% 
     mutate(cases = ifelse(date > ymd("20180101"), cases * 1.5, cases)) %>% 
     group_by(epiweek = floor_date(date, "week")) %>% 
     summarise(cases = sum(cases,na.rm=T))


ggplot()+
     geom_line(data = counts, mapping = aes(x = epiweek, y = cases))+
     labs(
          title = "Weekly reports of diarrhoea, 2013-current",
          x = "Week",
          y = "Weekly cases",
          caption = "")
# 
# # time series
# ## define start date (when observations began)
# start_date <- min(counts$epiweek)
# 
# ## define a cut-off week (end of baseline, start of prediction period)
# cut_off <- yearweek("2018-01-01")
# 
# ## define the last date interested in (i.e. end of prediction)
# end_date <- yearweek("2018-04-30")
# 
# ## find how many weeks in period (year) of interest
# num_weeks <- as.numeric(end_date - cut_off)
# 
# 
# ## add in missing weeks till end of year 
# counts <- counts %>%
#      ## group by the region
#      group_by_key() %>%
#      ## for each group add rows from the highest epiweek to the end of year
#      group_modify(~add_row(.,
#                            epiweek = seq(max(.$epiweek) + 1, 
#                                          end_date,
#                                          by = 1)))
# 
# counts <- counts %>% 
#      mutate(
#           ## combine fourier terms for weeks prior to  and after 2010 cut-off date
#           ## (nb. 2011 fourier terms are predicted)
#           fourier = rbind(
#                ## get fourier terms for previous years
#                fourier(
#                     ## only keep the rows before 2011
#                     filter(counts, 
#                            epiweek <= cut_off), 
#                     ## include one set of sin cos terms 
#                     K = 1
#                ), 
#                ## predict the fourier terms for 2011 (using baseline data)
#                fourier(
#                     ## only keep the rows before 2011
#                     filter(counts, 
#                            epiweek <= cut_off),
#                     ## include one set of sin cos terms 
#                     K = 1, 
#                     ## predict 52 weeks ahead
#                     h = num_weeks
#                )
#           )
#      )








# Outputs -----------------------------------------------------------------

# head - first 5 cases
chol %>% 
     select(case_id,
            #name,
            town, sex, age_years, 
            #age_grp,  # for use later
            date_onset,
            #date_admit,  # dates of admit appear to be before dates of onset?
            prior_vax, dehy,
            #oedema,
            preg, blood_stool,
            rdt,
            outcome) %>% 
     arrange(date_onset) %>% 
     head() %>% 
     qflextable() %>% 
     fontsize(size = 8, part = "all") %>%
     line_spacing(space = 1) %>% 
     theme_box() %>% 
     padding(padding = 1, part = "all") %>% 
     fit_to_width(max_width = 6)
     #autofit()

# Demographic pyramid
age_pyramid(chol,
            age_group = age_grp,
            split_by = sex,
            proportional = TRUE,
            show_midpoint = FALSE)+
     labs(title = "Age and sex of cholera cases",
          subtitle = str_glue("Total of {nrow(chol)} cases as of {format(max(chol$date_admit, na.rm=T), '%d %b, %Y')}"),
          fill = "Sex",
          y = "Proportion of total",
          x = "Age group")+
     theme_minimal(base_size = 8)+
     theme(
          legend.position = "top"
     )

ggsave("demographic_pyramid.png", width = 3, height = 3)




# Epicurve to show waves of epidemic
epiweeks <- seq.Date(floor_date(min(chol$date_onset, na.rm=T), "week"),
                     ceiling_date(max(chol$date_onset, na.rm=T), "week"),
                     "week")

ggplot(data = chol, mapping = aes(x = date_onset, fill = town))+
     geom_histogram(breaks = epiweeks, closed = "left",
                    color = "black")+
     scale_fill_brewer(type = "qual")+
     scale_x_date(
          date_breaks = "2 weeks",
          labels = scales::label_date_short())+
     scale_y_continuous(
          breaks = seq(0, 40, 5)
     )+
     facet_wrap(~town, ncol = 1)+
     theme_minimal()+
     labs(
          title = "Epidemic curve by town",
          subtitle = "300 cases as of 30 April, 2018",
          x = "Week of symptom onset",
          y = "Weekly cases",
          fill = "Town"
     )+
     theme(legend.position = "none")

ggsave("epicurve by town.png", width = 3, height = 6)
     

# Epicurve showing outcomes improving as response infrastructure improves
ggplot(data = chol, mapping = aes(x = date_onset, fill = outcome))+
     geom_histogram()+
     labs()+
     facet_wrap(~town)




# show CFR declining over time
chol %>% 
     group_by(week = floor_date(date_onset, "week")) %>% 
     summarise(CFR = sum(outcome == "Death", na.rm=T) / n()) %>% 
     ggplot(mapping = aes(x = week, y = CFR))+
     geom_line(size = 2)+
          scale_y_continuous(labels = scales::percent)+
          scale_x_date(labels = scales::label_date_short())+
          theme_minimal(base_size = 14)

