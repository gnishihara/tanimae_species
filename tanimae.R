# NMDS analysis of seaweed taxa in
# Arikawa Bay, Nagasaki
# Data: Shinichiro Tanimae
# 2021 / 06 / 27

library(vegan)
library(tidyverse)
library(lubridate)
library(gnnlab)
library(ggrepel)
library(showtext)
library(magick)

# REFERENCE
# https://jkzorz.github.io/2020/04/04/NMDS-extras.html

font_add_google("Noto Sans", "notosans")
showtext_auto()

# Load dataset
fname = "~/Lab_Data/tanimaes/share_files/seaweed_sp_data_210421-210623.csv"
tanimae = read_csv(fname) |>   mutate(station = factor(station))
tanimae = tanimae |>   rename(g = class)

# Some species have not been identified to the genus level
t1 = tanimae |>
  group_by(date, station, g) |> 
  summarise(n = length(g))

ggplot(t1) + 
  geom_col(aes(y = g, x = n, fill = station),
           position = position_dodge2()) +
  facet_grid(rows = vars(date))

# Shannon Diversity of the taxa
t2 = tanimae |> 
  group_by(date, station, g) |> 
  summarise(n = length(g)) |> 
  ungroup()|> pivot_wider(names_from = g,
                  values_from= n,
                  values_fill = 0) |> 
  print()


# NMDS 
t3 = t2 |> group_nest(date) |> 
  mutate(shannon = map(data, function(x) {
    z = x |> pull(station)
    s = x |> dplyr::select(-station) |> diversity()
    tibble(station = z, shannon = s)
  })) |> 
  mutate(bray = map(data, function(x) {
    # Pair wise dissimilarity index (Bray-Curtis) 
    # Range is from 0 (similar) to 1 (dissimilar)
    z = x |> pull(station)
    s = x |> dplyr::select(-station) |> vegdist("bray")
    tibble(station = z, bray = list(s))
  })) 

M = t3 |> select(date, data) |> slice(1) |> unnest(data) |> select(-date)
# First try, k = 2 but the stress was > 0.05.
# Since stress was > 0.05, set k = 3
nmdsout = M |> select(-station) |> as.matrix() |> metaMDS(k = 2)

nmds_scores = scores(nmdsout) |> as_tibble() |> 
  mutate(station = 1:7) |> 
  mutate(location = ifelse(station %in% c(3,5,6, 7), 
                           "Outside", "Inside")) |> 
  mutate(station = factor(station, 
                          levels = 1:7, 
                          labels = str_glue("St. {1:7}"))) 

# Plot of NMDS1 and NMDS2 for all times.
ggplot(nmds_scores) + 
  geom_point(aes(x = NMDS1, y = NMDS2, color = location)) +
  geom_text_repel(aes(x = NMDS1, y = NMDS2, 
                      color = location, label = station)) +
  scale_color_manual(values = viridis::viridis(3))  +
  scale_x_continuous(parse(text = "NMDS[1]")) +
  scale_y_continuous(parse(text = "NMDS[2]")) +
  guides(color = guide_legend(override.aes = list(label = ""))) +
  ggtitle("Family level NMDS plot") +
  theme_gray(base_family = "notosans") +
  theme(legend.position = c(0,0),
        legend.justification = c(0,0),
        legend.background = element_blank(),
        legend.title = element_blank())

## Save figure as pdf and png.
wh = aseries(7)
fname = "tanimae.pdf"
ggsave(fname, width = wh[2], height = wh[1], units = "mm")

image_read_pdf(fname) |> 
  image_resize("600x") |> 
  image_write(str_replace(fname, "pdf", "png"))

#################################################################################
## Try looking at monthly data.

t5 = t3 |> 
  mutate(nmds = map(data, function(x) {
    x |> select(-station) |> as.matrix() |> metaMDS(k = 2)
  })) |> 
  mutate(nmds = map(nmds, function(x) {
    s = scores(x) |> as_tibble()
    s |> mutate(station = 1:7) |> 
      mutate(station = factor(station))
  })) |> 
  unnest(nmds) 

t5 |> 
  mutate(month = month(date)) |> 
  mutate(month = month.abb[month]) |> 
  mutate(location = ifelse(station %in% c(3,5,6, 7), "Outside", "Inside")) |> 
  ggplot() + 
  geom_point(aes(x = NMDS1, y = NMDS2, color = location),size = 2) +
  geom_text_repel(aes(x = NMDS1, y = NMDS2, label = station, color = location), size = 10) +
  scale_color_manual(values = viridis::viridis(3)) +
  guides(color = guide_legend(override.aes = list(label = ""))) +
  ggtitle("Family level NMDS plot") +
  theme_gray(base_family = "notosans") +
  theme(legend.position = c(0,0),
        legend.justification = c(0,0),
        legend.background = element_blank(),
        legend.title = element_blank()) +
  facet_wrap("date")


