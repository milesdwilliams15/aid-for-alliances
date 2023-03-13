## Main Analysis Script ##
## Started: 7/14/2022   ##
## Updated: 2/21/2023   ##

# Note: After being rejected from ISQ I'm doing a few things:
# 1. dropping the stuff about aid channels
# 2. trying out a generalized diff-in-diff design
# 3. I'll also check out some interactions with donor gdp


# setup -------------------------------------------------------------------
rm(list = ls())
library(tidyverse) # because tidyverse is life
library(here)      # makes getting files easier
library(estimatr)  # for linear models with robust inference
library(coolorrr)
set_theme()
set_palette()

# some helpers ------------------------------------------------------------

de_mean  <- function(x) x - mean(x, na.rm=T)
re_scale <- function(x) (x - mean(x, na.rm=T)) / sd(x, na.rm=T)


# data --------------------------------------------------------------------

dt <- read.csv(
  here("01_data/cleaned_data.csv")
)
dt <- dt %>%
  group_by(
    donor, year
  ) %>%
  mutate(
    wgi_stand = re_scale(wgi_rank)
  ) %>%
  ungroup %>%
  mutate(
    trt_oda = ifelse(trt_oda <= 0, 0, trt_oda),
    p_bypass = (total_oda - gov_oda) / total_oda,
    p_public = 1 - p_bypass
  )
dt <- dt %>% # lag alliance variables and covariates
  group_by(dyad) %>%
  mutate(
    across(
      c(total_oda, wgi_stand, defense, nonagg, income, pop, fh_total, trade, fdi, usmil),
      ~ lag(.x, order_by = year),
      .names = "{.col}_lag"
    )
  ) %>% drop_na(defense) %>%
  mutate(
    dyad = as.factor(dyad)
  )

# estimate models ---------------------------------------------------------

base_form <- asinh(total_oda) ~ 
  # controls
  asinh(income_lag) + asinh(pop_lag) + asinh(disaster) + 
  civilwar+ fh_total_lag + asinh(dist) + asinh(trade_lag) +
  asinh(fdi_lag) + wgi_stand_lag + asinh(usmil) 

forms <- list(
  update(base_form, . ~ - . + nonagg),
  update(base_form, . ~ - . + defense),
  update(base_form, . ~ - . + nonagg + defense),
  update(base_form, . ~ - . + nonagg + defense + nonagg:defense),
  update(base_form, . ~ - . + asinh(total_oda_lag) + nonagg),
  update(base_form, . ~ - . + asinh(total_oda_lag) + defense),
  update(base_form, . ~ - . + asinh(total_oda_lag) + nonagg + defense),
  update(base_form, . ~ - . + asinh(total_oda_lag) + nonagg + defense + nonagg:defense),
  update(base_form, . ~ . + nonagg),
  update(base_form, . ~ . + defense),
  update(base_form, . ~ . + nonagg + defense),
  update(base_form, . ~ . + nonagg + defense + nonagg:defense),
  update(base_form, . ~ . + asinh(total_oda_lag) + nonagg),
  update(base_form, . ~ . + asinh(total_oda_lag) + defense),
  update(base_form, . ~ . + asinh(total_oda_lag) + nonagg + defense),
  update(base_form, . ~ . + asinh(total_oda_lag) + nonagg + defense + nonagg:defense)
)

## Via OLS ----

main_ols_fits <- forms %>%
  map(
    ~ lm_robust(
      .x,
      data = dt,
      fixed_effects = ~ dyad + year,
      clusters = dyad,
      se_type = 'stata'
    )
  )
save(
  main_ols_fits,
  file = here::here("04_regression_output",
                    "main_ols_fits.R")
)
# ## Via PPML ----
# 
# ppml <- function(..., clusters = NULL) {
#   cat("\nWorking...")
#   fit <- glm(..., family = quasipoisson)
#   if(is.null(clusters)) {
#     vcv <- sandwich::vcovHC(fit, type = 'HC1')
#   } else {
#     vcv <- sandwich::vcovCL(fit, clusters = clusters)
#   }
#   return(
#     list(
#       fit = fit,
#       vcv = vcv
#     )
#   )
#   cat("\nDone!\n")
# }
# main_ppml_fits <- forms %>%
#   map(
#     ~ ppml(
#       update(.x, ~ . + dyad + year),
#       data = dt,
#       clusters = dt$dyad
#     )
#   )


# regression tables -------------------------------------------------------

library(texreg)
cmap <- list(
  'nonagg' = 'Nonaggression',
  'defense' = 'Defense',
  'nonagg:defense' = 'Interaction'
)
screenreg( # only FEs
  main_ols_fits[1:4], include.ci=F,
  custom.coef.map = cmap,
  stars = c(0.1,0.05,0.01,0.001)
)
screenreg( # FEs + DV lag 
  main_ols_fits[5:8], include.ci=F,
  custom.coef.map = cmap,
  stars = c(0.1,0.05,0.01,0.001)
)
screenreg( # FEs + controls
  main_ols_fits[9:12], include.ci=F,
  custom.coef.map = cmap,
  stars = c(0.1,0.05,0.01,0.001)
)
screenreg( # FEs + lag + controls
  main_ols_fits[13:16], include.ci=F,
  custom.coef.map = cmap,
  stars = c(0.1,0.05,0.01,0.001)
)

# stylized facts ----------------------------------------------------------
dt %>%
  mutate(
    alliance = (nonagg | defense)+0
  ) %>%
  group_by(
    dyad
  ) %>%
  mutate(
    ally_switch = c(NA, diff(alliance, na.rm=T)),
    aid_change = c(NA, diff(total_oda))
  ) %>%
  ungroup %>%
  filter(ally_switch!=0) %>%
  group_by(ally_switch) %>%
  summarize(
    increase = #sum(aid_change > 0),
      median(aid_change[aid_change>0]),
    decrease = #sum(aid_change < 0)
      median(aid_change[aid_change<0])
  )

dt %>%
  group_by(year) %>%
  summarize(
    allies = sum(nonagg | defense) ,
    defense = sum(defense) ,
    nonagg = sum(nonagg) 
  )

dt %>%
  group_by(donor) %>%
  mutate(
    have_def = any(defense)
  ) %>%
  ungroup %>%
  filter(have_def) %>%
  count(donor)


# jack-knife --------------------------------------------------------------

# I need to make a helper that residualizes the data by dyad
# and year so I can do the leave-one-donor out validation
demean <- function(form, data) {
  mm <- model.frame(update(form, ~ . + dyad + year + donor), data)
  mm %>%
    select(-c(dyad, year, donor)) %>%
    mutate(
      across(
        everything(),
        ~ resid(lm(.x ~ mm$dyad + mm$year))
      ),
      donor = mm$donor
    )
}
# I only want to run this once bc it takes forever!
if(!exists(here::here("01_data", "jackdata.csv"))) {
  jackforms <- list(
    base_form,
    update(base_form, ~ . + nonagg),
    update(base_form, ~ . + defense),
    update(base_form, ~ . + nonagg + defense),
    update(base_form, ~ . + total_oda_lag),
    update(base_form, ~ . + total_oda_lag + 
             nonagg),
    update(base_form, ~ . + total_oda_lag +
             defense),
    update(base_form, ~ . + total_oda_lag +
             nonagg + defense)
  )
  jackdata <- demean(jackforms[[8]], dt)
  # whew! that took a while...
  names(jackdata) <- names(jackdata) %>%
    str_remove("asinh") %>%
    gsub("[()]", "", .)
  write_csv(jackdata, here::here("01_data", "jackdata.csv"))
} else {
  jackdata <- read_csv(here::here("01_data", "jackdata.csv"))
}

# make formulas
base_form <- total_oda ~ 
  # controls
  income_lag + pop_lag + disaster + 
  civilwar+ fh_total_lag + dist + trade_lag +
  fdi_lag + wgi_stand_lag + usmil 

jackforms <- list(
  base_form,
  update(base_form, ~ . + nonagg),
  update(base_form, ~ . + defense),
  update(base_form, ~ . + nonagg + defense),
  update(base_form, ~ . + total_oda_lag),
  update(base_form, ~ . + total_oda_lag + 
           nonagg),
  update(base_form, ~ . + total_oda_lag +
           defense),
  update(base_form, ~ . + total_oda_lag +
           nonagg + defense)
)

donors <- dt$donor %>%
  unique()

donors %>%
  map_dfr(
    ~ {
      drop_donor <- .x
      jackforms %>%
        map_dfr(
          ~ {
            lm(
              .x,
              data = filter(jackdata, donor != drop_donor)
            ) -> fit
            pred_data <- filter(jackdata, donor == drop_donor)
            predict(
              fit, newdata = pred_data
            ) -> pred
            r <- sqrt(mean((pred_data$total_oda - pred)^2))
            tibble(
              R2 = r
            ) %>%
              mutate(
                donor = drop_donor,
                fit = list(fit)
              )
    })
  }
) -> jack_val

jack_val %>%
  mutate(
    model = rep(
      c("+ 0",
        "+ Nonagg.",
        "+ Defense",
        "+ Nonagg. + Defense"),
      len = n()
    ),
    design = rep(
      c("D-in-D", "lagged-DV"),
      each = 4
    ) %>% rep(len = n())
  ) %>%
  group_by(donor, design) %>%
  mutate(
    R2change = R2 - R2[model == "+ 0"]
  ) %>%
  ungroup -> jack_val

jack_val %>%
  filter(R2change != 0) %>%
  group_by(
    model, design
  ) %>%
  summarize(
    mean = mean(R2change),
    median = median(R2change),
    min = min(R2change),
    max = max(R2change),
    improved = paste0(round(100 * mean(R2change < 0), 2), "%")
  ) %>%
  ggplot() +
  aes(y = model) +
  facet_wrap(~ design) +
  geom_errorbarh(
    aes(xmin = min,
        xmax = max),
    height = 0.25
  ) +
  geom_point(
    aes(x = mean),
    color = "darkred"
  ) +
  geom_point(
    aes(x = median),
    color = "darkblue",
    shape = 12
  ) +
  geom_text(
    aes(x = (mean + max)/2,
        label = improved),
    vjust = -1
  ) +
  labs(
    x = "Difference in RMSE",
    y = NULL
  )

set_palette(
  binary = c("darkblue",
             "lightblue"),
  from_coolors = F
)
jack_val %>%
  arrange(R2change) %>%
  filter(model =="+ Nonagg. + Defense") %>%
  select(donor, R2change, design) %>%
  ggplot() +
  aes(x = R2change,
      y = reorder(donor, -R2change),
      shape = design,
      color = design) +
  geom_point() +
  geom_vline(
    xintercept = 0
  ) +
  ggpal(
    type = "binary"
  ) +
  facet_wrap(~ ifelse(
    R2change < 0, "Predictions Improved",
    "Predictions Worsened"
  ),
  scales = "free_y") +
  labs(
    x = "Change in RMSE",
    y = NULL,
    color = NULL,
    shape = NULL
  ) +
  theme(
    legend.position = "right",
    legend.direction = "vertical"
  )
ggsave(
  here::here("03_figures",
             "crosval.png"),
  dpi = 500,
  height = 6,
  width = 6
)

dt %>%
  group_by(donor) %>%
  summarize(
    defense = sum(defense==1)+0,
    nonagg = sum(nonagg==1)+0
  ) %>%
  right_join(
    jack_val, by = "donor"
  ) %>%
  group_by(defense, nonagg) %>%
  filter(model == "+ Nonagg. + Defense") %>%
  summarize(
    improved = 100 * mean(R2change < 0),
    donors = paste0(unique(donor), collapse = ", ")
  ) %>%
  ungroup %>%
  data.frame %>%
  stargazer::stargazer(
    summary = F
  )

# sensitivity analysis with sensemakr -------------------------------------

library(sensemakr)
dd_def_sen <- sensemakr(
  jackforms[[4]], data = jackdata,
  treatment = "defense",
  benchmark_covariates = c("income_lag",
                           "pop_lag"),
  kd = 1:3
)
dd_non_sen <- sensemakr(
  jackforms[[4]], data = jackdata,
  treatment = "nonagg",
  benchmark_covariates = c("income_lag",
                           "pop_lag"),
  kd = 1:3
)
ld_def_sen <- sensemakr(
  jackforms[[8]], data = jackdata,
  treatment = "defense",
  benchmark_covariates = c("income_lag",
                           "pop_lag"),
  kd = 1:3
)
ld_non_sen <- sensemakr(
  jackforms[[8]], data = jackdata,
  treatment = "nonagg",
  benchmark_covariates = c("income_lag",
                           "pop_lag"),
  kd = 1:3
)

summary(dd_def_sen)
summary(dd_non_sen)
summary(ld_def_sen)
summary(ld_non_sen)

ovb_minimal_reporting(dd_def_sen)
ovb_minimal_reporting(dd_non_sen)
ovb_minimal_reporting(ld_def_sen)
ovb_minimal_reporting(ld_non_sen)
