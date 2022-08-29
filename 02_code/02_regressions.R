## Main Analysis Script ##
## Started: 7/14/2022   ##
## Updated: 7/14/2022   ##


# setup -------------------------------------------------------------------
rm(list = ls())
library(tidyverse) # because tidyvers is life
library(here)      # makes getting files easier
library(estimatr)  # for linear models with robust inference
library(AER)       # for tobit
library(censReg)   # for tobit with random intercepts
library(mgcv)      # for logit with random intercepts


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
      c(defense, nonagg, income, pop, fh_total, trade, fdi, usmil),
      ~ lag(.x, order_by = year)
    )
  ) %>% drop_na(defense)

# estimate models ---------------------------------------------------------

base_form <- ~ nonagg + defense + # outcomes of interest
  # controls
  asinh(income) + asinh(pop) + asinh(disaster) + 
  civilwar + fh_total + asinh(dist) + asinh(trade) +
  asinh(fdi) + asinh(usmil) + colony

main_forms <- list(
  update(base_form, asinh(total_oda) ~ .),
  update(base_form, asinh(gov_oda) ~ .),
  update(base_form, asinh(total_oda - gov_oda) ~ .),
  update(base_form, asinh(ngo_oda) ~ .),
  # update(base_form, asinh(ppp_oda) ~ .), # no public private partnerships
  update(base_form, asinh(mlt_oda) ~ .),
  update(base_form, asinh(trt_oda) ~ .),
  update(base_form, asinh(otr_oda) ~ .)
)
bi_forms <- list(
  update(base_form, total_oda > 0 ~ .),
  update(base_form, gov_oda > 0 ~ .),
  update(base_form, (total_oda - gov_oda) > 0 ~ .),
  update(base_form, ngo_oda > 0 ~ .),
  update(base_form, mlt_oda > 0 ~ .),
  update(base_form, trt_oda > 0 ~ .),
  update(base_form, otr_oda > 0 ~ .)
)

# intr_forms <- main_forms %>%
#   map(
#     ~ update(.x, ~ . + wgi_stand:nonagg + wgi_stand:defense)
#   )
  

## Via Truncated OLS ----

main_ols_fits <- main_forms %>%
  map(
    ~ lm_robust(
      .x,
      data = dt,
      fixed_effects = ~ donor + year,
      clusters = dyad,
      se_type = 'stata'
    )
  )
# intr_ols_fits <- intr_forms %>%
#   map(
#     ~ lm_robust(
#       .x,
#       data = dt,
#       fixed_effects = ~ donor + year,
#       clusters = dyad,
#       se_type = 'stata'
#     )
#   )

## Via Tobit ----

main_tfe_fits <- main_forms %>%
  map(
    ~ {
      tob <- tobit(
        update(.x, ~ . + donor + as.factor(year)),
        data = dt 
      )
      tob$var <- vcovCL(tob, cluster = dt$dyad, type = 'HC1')
      tob
    }
  )
# intr_tfe_fits <- intr_forms %>%
#   map(
#     ~ {
#       tob <- tobit(
#         update(.x, ~ . + donor + as.factor(year)),
#         data = dt
#       )
#       tob$var <- vcovCL(tob, cluster = dt$dyad, type = 'HC1')
#       tob
#     }
#   )


## Via Mixed Effects Tobit ----

main_tme_fits <- main_forms %>%
  map(
    ~ censReg(
      update(.x, ~ . + donor + poly(year, 3)),
      data = plm::pdata.frame(dt, index = 'dyad'),
      method = 'BHHH',
      nGHQ = 4
    )
  )


## Via Mixed Effects Logit ----
main_lme_fits <- bi_forms[1:5] %>%
  map(
    ~ gam(
      update(.x, ~ . + donor + poly(year, 3) + s(id, bs = 're')),
      data = dt %>% mutate(id = as.numeric(as.factor(dyad))),
      family = binomial,
      method = 'REML'
    )
  )
# intr_tme_fits <- intr_forms %>%
#   map(
#     ~ censReg(
#       update(.x, ~ . + donor + poly(year, 3)),
#       data = plm::pdata.frame(dt, index = 'dyad'),
#       method = 'BHHH',
#       nGHQ = 4
#     )
#  )

## Via PPML ----

ppml <- function(..., clusters = NULL) {
  fit <- glm(..., family = quasipoisson)
  if(is.null(clusters)) {
    vcv <- sandwich::vcovHC(fit, type = 'HC0')
  } else {
    vcv <- sandwich::vcovCL(fit, clusters = clusters)
  }
  return(
    list(
      fit = fit,
      vcv = vcv
    )
  )
}
main_ppml_fits <- list(
  update(base_form, total_oda ~ .),
  update(base_form, gov_oda ~ .),
  update(base_form, total_oda - gov_oda ~ .),
  update(base_form, ngo_oda ~ .),
  # update(base_form, asinh(ppp_oda) ~ .), # no public private partnerships
  update(base_form, mlt_oda ~ .),
  update(base_form, trt_oda ~ .),
  update(base_form, otr_oda ~ .)
) %>%
  map(
    ~ ppml(
      update(.x, ~ . + donor + poly(year, 3)),
      data = dt,
      clusters = dt$dyad
    )
  )

main_ppml_fits[1:3] %>%
  map_dfr(
    ~ coeftest(
      .x$fit, vcov. = .x$vcv
    ) %>%
      estimatr::tidy() %>%
      filter(term %in% c('nonagg', 'defense'))
  )
# estimate logged-proportional odds model via mixed effects logit

# lme_fit <- gam(
#   update(base_form, p_public ~ . + s(id, bs = 're')),
#   data = dt%>% mutate(id = as.numeric(as.factor(dyad))),
#   family = binomial,
#   method = 'REML'
# )



# regression tables -------------------------------------------------------

library(texreg)
cmap <- list(
  'nonagg' = 'Nonaggression',
  'defense' = 'Defense'
)
model_names <- c(
  'Total',
  'Gov-to-Gov',
  'Bypass',
  'NGOs',
  'MOs',
  'TRT',
  'Other'
)
screenreg( # OLS fits by total aid and public vs. bypass channels
  main_ols_fits[1:3], include.ci=F,
  custom.coef.map = cmap,
  custom.model.names = model_names[1:3],
  stars = c(0.1,0.05,0.01,0.001)
)
screenreg( # OLS fits by bypass channel type
  main_ols_fits[-c(1:3)], include.ci = F,
  custom.coef.map = cmap,
  custom.model.names = model_names[-c(1:3)]
)
# screenreg(
#   intr_ols_fits, include.ci=F,
#   custom.coef.map = cmap
# )
screenreg( # Tobit fits by total aid and public vs. bypass channels
  main_tfe_fits[1:3], include.ci=F,
  custom.coef.map = cmap,
  custom.model.names = model_names[1:3]
)
screenreg( # Tobit fits by bypass channel type
  main_tfe_fits[-c(1:3)], include.ci = F,
  custom.coef.map = cmap,
  custom.model.names = model_names[-c(1:3)]
)
# screenreg(
#   intr_tfe_fits,
#   custom.coef.map = cmap
# )
screenreg( # ME Tobit fits by total aid and public vs. bypass channels
  main_tme_fits[1:3], include.ci=F,
  custom.coef.map = cmap,
  custom.model.names = model_names[1:3]
)
screenreg( # ME Tobit fits by bypass channel type (only first two)
  main_tme_fits[4:7], include.ci = F,
  custom.coef.map = cmap,
  custom.model.names = model_names[4:7]
)
screenreg( # ME logit fits for bypass channels that won't converge with ME Tobit
  main_lme_fits[4:7], include.ci = F,
  custom.coef.map = cmap,
  custom.model.names = model_names[4:7]
)
# screenreg(
#   intr_tme_fits,
#   custom.coef.map = cmap
# )
# screenreg(
#   lme_fit,
#   custom.coef.map = cmap
# )

# coefplots ---------------------------------------------------------------

tidy_fits1 <- main_tfe_fits %>%
  map_dfr(
    ~ coeftest(.x) %>%
      estimatr::tidy()
  ) %>%
  filter(
    term %in% c('nonagg', 'defense')
  ) %>%
  mutate(
    ODA = rep(
      model_names,
      each = 2
    ),
    method = 'Tobit (Rob. SEs)'
  )
tidy_fits2 <- main_tme_fits %>%
  map_dfr(
    ~ coeftest(.x) %>%
      estimatr::tidy()
  ) %>%
  filter(
    term %in% c('nonagg', 'defense')
  ) %>%
  mutate(
    ODA = rep(
      model_names,
      each = 2
    ),
    method = 'Tobit (Dyadic MEs)'
  )
tidy_fits3 <- main_lme_fits %>%
  map_dfr(
    ~ coeftest(.x) %>%
      estimatr::tidy()
  ) %>%
  filter(
    term %in% c('nonagg', 'defense')
  ) 
tidy_fits4 <- main_ppml_fits %>%
  map_dfr(
    ~ coeftest(.x$fit, vcov. = .x$vcv) %>%
      estimatr::tidy()
  ) %>%
  filter(
    term %in% c('nonagg', 'defense')
  )
tidy_fits3$ODA <- rep(model_names[1:5], each = 2)
tidy_fits3$method <- "Logit (Dyadic MEs)"
tidy_fits4$ODA <- rep(model_names, each = 2)
tidy_fits4$method <- "PPML (Robust SEs)"
  # mutate(
  #   ODA = rep(
  #     model_names,
  #     each = 2
  #   ),
  #   method = 'Logit (Dyadic MEs)'
  # )
tidy_fits <- bind_rows(tidy_fits1, tidy_fits2)


library(patchwork)
p1 <- tidy_fits2 %>%
  filter(ODA=='Total') %>%
  ggplot() +
  aes(
    x = estimate,
    xmin = estimate - 1.96 * std.error,
    xmax = estimate + 1.96 * std.error,
    y = term,
    label = round(estimate, 2)
  ) +
  geom_point() +
  geom_errorbarh(
    height = 0
  ) +
  geom_text(
    vjust = -1
  ) +
  geom_vline(
    xintercept = 0,
    lty = 2
  ) +
  scale_y_discrete(
    labels = c(
      'Defense',
      'Nonaggression'
    )
  ) +
  labs(
    x = 'Estimate with 95% CI',
    y = NULL,
    title = 'Total ODA'
  ) +
  theme_light()

p2 <- ggplot(tidy_fits2 %>%
         filter(ODA %in% c('Gov-to-Gov', 'Bypass'))) +
  aes(
    x = estimate,
    xmin = estimate - 1.96 * std.error,
    xmax = estimate + 1.96 * std.error,
    y = term,
    color = ODA,
    label = round(estimate, 2)
  ) +
  geom_point(
    position = ggstance::position_dodgev(-.5)
  ) +
  geom_errorbarh(
    position = ggstance::position_dodgev(-.5),
    height = 0
  ) +
  geom_text(
    show.legend = F,
    vjust = -1,
    position = ggstance::position_dodgev(-.5)
  ) +
  geom_vline(
    xintercept = 0,
    lty = 2
  ) +
  labs(
    x = 'Estimate with 95% CIs',
    y = NULL,
    color = 'Channel',
    title = 'ODA by delivery channel'
  ) +
  scale_y_discrete(
    labels = c('Defense', 'Nonaggression')
  ) +
  theme_light() +
  theme(
    legend.position = c(0.2,0.7),
    legend.background = element_rect(
      color = 'black'
    )
  )
p1 + p2
ggsave(
    here("03_figures/main_effects.png"),
    dpi = 1200,
    height = 3,
    width = 10
  )

p1 <- tidy_fits4 %>%
  filter(ODA=='Total') %>%
  ggplot() +
  aes(
    x = estimate,
    xmin = estimate - 1.96 * std.error,
    xmax = estimate + 1.96 * std.error,
    y = term,
    label = round(estimate, 2)
  ) +
  geom_point() +
  geom_errorbarh(
    height = 0
  ) +
  geom_text(
    vjust = -1
  ) +
  geom_vline(
    xintercept = 0,
    lty = 2
  ) +
  scale_y_discrete(
    labels = c(
      'Defense',
      'Nonaggression'
    )
  ) +
  labs(
    x = 'Estimate with 95% CI',
    y = NULL,
    title = 'Total ODA'
  ) +
  theme_light()

p2 <- ggplot(tidy_fits4 %>%
               filter(ODA %in% c('Gov-to-Gov', 'Bypass'))) +
  aes(
    x = estimate,
    xmin = estimate - 1.96 * std.error,
    xmax = estimate + 1.96 * std.error,
    y = term,
    color = ODA,
    label = round(estimate, 2)
  ) +
  geom_point(
    position = ggstance::position_dodgev(-.5)
  ) +
  geom_errorbarh(
    position = ggstance::position_dodgev(-.5),
    height = 0
  ) +
  geom_text(
    show.legend = F,
    vjust = -1,
    position = ggstance::position_dodgev(-.5)
  ) +
  geom_vline(
    xintercept = 0,
    lty = 2
  ) +
  labs(
    x = 'Estimate with 95% CIs',
    y = NULL,
    color = 'Channel',
    title = 'ODA by delivery channel'
  ) +
  scale_y_discrete(
    labels = c('Defense', 'Nonaggression')
  ) +
  theme_light() +
  theme(
    legend.position = c(0.2,0.7),
    legend.background = element_rect(
      color = 'black'
    )
  )
p1 + p2
ggsave(
  here("03_figures/main_ppml_effects.png"),
  dpi = 1200,
  height = 3,
  width = 10
)

p1 <- tidy_fits3 %>%
  filter(ODA=='Total') %>%
  ggplot() +
  aes(
    x = estimate,
    xmin = estimate - 1.96 * std.error,
    xmax = estimate + 1.96 * std.error,
    y = term,
    label = round(estimate, 2)
  ) +
  geom_point() +
  geom_errorbarh(
    height = 0
  ) +
  geom_text(
    vjust = -1
  ) +
  geom_vline(
    xintercept = 0,
    lty = 2
  ) +
  scale_y_discrete(
    labels = c(
      'Defense',
      'Nonaggression'
    )
  ) +
  labs(
    x = 'Estimate with 95% CI',
    y = NULL,
    title = 'Total ODA'
  ) +
  theme_light()

p2 <- ggplot(tidy_fits3 %>%
               filter(ODA %in% c('Gov-to-Gov', 'Bypass'))) +
  aes(
    x = estimate,
    xmin = estimate - 1.96 * std.error,
    xmax = estimate + 1.96 * std.error,
    y = term,
    color = ODA,
    label = round(estimate, 2)
  ) +
  geom_point(
    position = ggstance::position_dodgev(-.5)
  ) +
  geom_errorbarh(
    position = ggstance::position_dodgev(-.5),
    height = 0
  ) +
  geom_text(
    show.legend = F,
    vjust = -1,
    position = ggstance::position_dodgev(-.5)
  ) +
  geom_vline(
    xintercept = 0,
    lty = 2
  ) +
  labs(
    x = 'Estimate with 95% CIs',
    y = NULL,
    color = 'Channel',
    title = 'ODA by delivery channel'
  ) +
  scale_y_discrete(
    labels = c('Defense', 'Nonaggression')
  ) +
  theme_light() +
  theme(
    legend.position = c(0.2,0.7),
    legend.background = element_rect(
      color = 'black'
    )
  )
p1 + p2
ggsave(
  here("03_figures/main_logit_effects.png"),
  dpi = 1200,
  height = 3,
  width = 10
)

`%nin%` <- Negate(`%in%`)
p1 <- ggplot(tidy_fits1 %>%
         filter(ODA %nin% c('Total','Gov-to-Gov', 'Bypass'))) +
  aes(
    x = estimate,
    xmin = estimate - 1.96 * std.error,
    xmax = estimate + 1.96 * std.error,
    y = term,
    color = ODA
  ) +
  geom_point(
    position = ggstance::position_dodgev(-.5)
  ) +
  geom_errorbarh(
    position = ggstance::position_dodgev(-.5),
    height = 0
  ) +
  geom_vline(
    xintercept = 0,
    lty = 2
  ) +
  labs(
    x = 'Estimate with 95% CIs',
    y = NULL,
    color = 'Bypass Channel',
    title = 'Classic Tobit with cluster-robust SEs'
  ) +
  scale_y_discrete(
    labels = c('Defense', 'Nonaggression')
  ) +
  theme_light() +
  theme(
    legend.position = c(0.2,0.65),
    legend.background = element_rect(
      color = 'black'
    )
  )
p2 <- ggplot(tidy_fits3 %>%
               filter(ODA %nin% c('Total','Gov-to-Gov', 'Bypass'))) +
  aes(
    x = estimate,
    xmin = estimate - 1.96 * std.error,
    xmax = estimate + 1.96 * std.error,
    y = term,
    color = ODA
  ) +
  geom_point(
    position = ggstance::position_dodgev(-.5)
  ) +
  geom_errorbarh(
    position = ggstance::position_dodgev(-.5),
    height = 0
  ) +
  geom_vline(
    xintercept = 0,
    lty = 2
  ) +
  labs(
    x = 'Estimate with 95% CIs',
    y = NULL,
    color = 'Bypass Channel',
    title = 'Mixed effects logistic regression'
  ) +
  scale_y_discrete(
    labels = c('Defense', 'Nonaggression')
  ) +
  theme_light() +
  theme(
    legend.position = 'none'
  )

p1 + p2
ggsave(
  here("03_figures/bypass_effects.png"),
  dpi = 1200,
  height = 3,
  width = 10
)

# save the main regression results
save(
  main_ols_fits,
  file = here("04_regression_output/main_ols_fits.R")
)
save(
  main_tfe_fits,
  file = here("04_regression_output/main_tfe_fits.R")
)
save(
  main_tme_fits,
  file = here("04_regression_output/main_tme_fits.R")
)
save(
  main_lme_fits,
  file = here("04_regression_output/main_lme_fits.R")
)
save(
  main_ppml_fits,
  file = here("04_regression_output/main_ppml_fits.R")
)



