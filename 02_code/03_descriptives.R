## Descriptive Analysis Script ##


# setup -------------------------------------------------------------------
rm(list = ls())
library(tidyverse) # because tidyverse is life
library(here)      # makes navigating files easier

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

# average aid by allies and non allies ------------------------------------

theme_set(theme_light())
dt %>%
  filter(year!=2005) %>%
  mutate(
    ally = ifelse(
      nonagg | defense, "Alliance", "No Alliance"
    )
  ) %>%
  group_by(
    year, donor
  ) %>%
  mutate(
    aid_share = total_oda / sum(total_oda)
  ) %>%
  group_by(
    year, ally
  ) %>%
  summarize(
    aid = median(aid_share)
  ) %>%
  ggplot() +
  aes(
    x = year,
    y = aid,
    color = ally
  ) +
  geom_line(
    size = .75
  ) +
  scale_color_manual(
    values = c(
      "Alliance" = "indianred3",
      "No Alliance" = "royalblue"
    )
  ) + 
  scale_y_continuous(
    labels = scales::percent
  ) +
  labs(
    x = NULL,
    y = "Median Bilateral ODA Share",
    color = NULL
  ) +
  theme(
    legend.position = 'top'
  )
ggsave(
  here("03_figures/median_aid_by_alliance.png"),
  height = 3,
  width = 6
)

# aid by alliance for the US ----------------------------------------------

dt %>%
  filter(
    donor %in% c("USA", "JPN"),
    year != 2005
  ) %>%
  mutate(
    ally = ifelse(defense | nonagg, "Alliance", "No Alliance")
  ) %>%
  group_by(
    donor, year
  ) %>%
  group_by(
    donor, year, ally
  ) %>%
  summarize(
    aid = median(total_oda)
  ) %>%
  ggplot() +
  aes(
    x = year,
    y = aid,
    color = ally
  ) +
  facet_wrap(
    ~ donor,
    scales = "free_y"
  ) +
  geom_line(
    size = 0.75
  ) +
  scale_x_continuous(
    breaks = c(2006, 2014)
  ) +
  scale_y_continuous(
    labels = scales::comma
  ) +
  scale_color_manual(
    values = c(
      "Alliance" = "indianred3",
      "No Alliance" = "royalblue"
    )
  ) +
  labs(
    x = NULL,
    y = "Median ODA\n(commitments in mil. $)",
    color = NULL
  ) +
  theme(
    legend.position = 'top',
    axis.text = element_text(
      hjust = c(0, 1)
    )
  )
ggsave(
  here("03_figures/us_jp_aid.png"),
  height = 3,
  width = 6
)

# no Israel ---------------------------------------------------------------

dt %>%
  filter(
    donor %in% c("USA", "JPN"),
    recipient != "ISR",
    year != 2005
  ) %>%
  mutate(
    ally = ifelse(defense | nonagg, "Alliance", "No Alliance")
  ) %>%
  group_by(
    donor, year
  ) %>%
  group_by(
    donor, year, ally
  ) %>%
  summarize(
    aid = median(total_oda)
  ) %>%
  ggplot() +
  aes(
    x = year,
    y = aid,
    color = ally
  ) +
  facet_wrap(
    ~ donor,
    scales = "free_y"
  ) +
  geom_line(
    size = 0.75
  ) +
  scale_x_continuous(
    breaks = c(2006, 2014)
  ) +
  scale_y_continuous(
    labels = scales::comma
  ) +
  scale_color_manual(
    values = c(
      "Alliance" = "indianred3",
      "No Alliance" = "royalblue"
    )
  ) +
  labs(
    x = NULL,
    y = "Median ODA\n(commitments in mil. $)",
    color = NULL
  ) +
  theme(
    legend.position = 'top',
    axis.text = element_text(
      hjust = c(0, 1)
    )
  )
ggsave(
  here("03_figures/us_jp_aid_no_isr.png"),
  height = 3,
  width = 6
)

# alliances types by donor ------------------------------------------------

dt %>%
  filter(
    year > 2005,
    donor %in% c("USA", "JPN")
  ) %>%
  group_by(
    donor, year
  ) %>%
  summarise(
    Defense = mean(defense),
    Nonaggression = mean(nonagg)
  ) %>%
  pivot_longer(
    Defense:Nonaggression
  ) %>%
  ggplot() +
  aes(
    x = year,
    y = value,
    color = name
  ) +
  facet_wrap(
    ~ donor
  ) +
  geom_line(
    size = 0.75
  ) +
  scale_color_manual(
    values = c(
      "Defense" = "indianred3",
      "Nonaggression" = "royalblue"
    )
  ) +
  scale_x_continuous(
    breaks = c(2006, 2014)
  ) +
  scale_y_continuous(
    labels = scales::percent
  ) +
  labs(
    x = NULL,
    y = "Recipients with Alliance Type",
    color = NULL
  ) +
  theme(
    legend.position = 'top',
    axis.text = element_text(
      hjust = c(0, 1)
    )
  )
ggsave(
  here("03_figures/us_jp_alliances.png"),
  height = 3,
  width = 6
)

# look at Pakistan --------------------------------------------------------

dt %>%
  filter(
    recipient == "PAK",
    donor == "USA"
  ) %>%
  ggplot() +
  aes(year, total_oda) +
  geom_line(
    size = .75
  ) +
  geom_vline(
    aes(xintercept = min(year[nonagg==1]))
  ) +
  scale_x_continuous(
    breaks = seq(2005, 2014, by = 2)
  ) +
  scale_y_continuous(
    labels = scales::comma
  ) +
  labs(
    x = NULL,
    y = "Bilateral ODA from US\n(commitments in mil. $)"
  )
ggsave(
  here("03_figures/us_pak_aid.png"),
  height = 3,
  width = 6
)

# aid to albania ----------------------------------------------------------

dt %>%
  filter(
    recipient == "ALB"
  ) %>%
  group_by(
    donor
  ) %>%
  mutate(
    nato = max(defense)
  ) %>%
  ungroup %>%
  ggplot() +
  aes(
    x = year, 
    y = total_oda,
    color = ifelse(nato==1, 'NATO', 'non-NATO')
  ) +
  stat_summary(
    fun = ~ sinh(mean(asinh(.x))), geom = "line",
    size = 0.75
  ) +
  geom_vline(
    aes(xintercept = min(year[defense==1]))
  ) +
  scale_x_continuous(
    breaks = seq(2005, 2014, by = 2)
  ) +
  labs(
    x = NULL,
    y = "Mean ODA\n(commitments in mil. $)",
    color = NULL
  ) +
  scale_color_manual(
    values = c(
      "NATO" = "royalblue",
      "non-NATO" = "indianred3"
    )
  )
ggsave(
  here("03_figures/alb_nato_aid.png"),
  height = 3,
  width = 6
)
