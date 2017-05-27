required <- c("nandb", "filesstrings", "stringr", "tidyverse",
              "magrittr", "EBImage", "rcompanion")
for (r in required) if (!require(r, character.only = TRUE)) install.packages(r)
invisible(sapply(required, library, character.only = TRUE))

SimulateBleaching <- function(arr3d, bleaching.trace) {
  d <- dim(arr3d)
  n.frames <- d[3]
  stopifnot(bleaching.trace[1] == 1, n.frames == length(bleaching.trace),
            all(bleaching.trace > 0))
  frame.totals <- apply(arr3d, 3, sum)
  bleached.frame.totals <- round(frame.totals * bleaching.trace)
  frames.to.add <- bleached.frame.totals - frame.totals
  pixels.per.frame <- d[1] * d[2]
  bleached <- arr3d
  for (i in seq_along(frame.totals)) {
    frame.i.indices <- seq((i - 1) * pixels.per.frame + 1, i * pixels.per.frame)
    if (frames.to.add[i] > 0) {
      samp <- sample(frame.i.indices, frames.to.add[i], replace = TRUE)
      for (j in samp) {
        bleached[j] <- bleached[j] + 1
      }
    } else if (frames.to.add[i] < 0) {
      sample.x <- rep(frame.i.indices, bleached[frame.i.indices])
      samp <- sample(sample.x, -frames.to.add[i])
      for (j in samp) {
        bleached[j] <- bleached[j] - 1
      }
    }
  }
  bleached
}

FiveTenMers <- function(sim, brightness) {
  brightness[brightness < 1] <- NA
  ones <- brightness[NthNumber(sim, 1) == 1]
  fives <- brightness[NthNumber(sim, 1) == 5]
  tens <- brightness[NthNumber(sim, 1) == 10]
  fivemers <- expand.grid(ones, fives) %>% {(.[[2]] - 1) / (.[[1]] - 1)}
  tenmers <- expand.grid(ones, tens) %>% {(.[[2]] - 1) / (.[[1]] - 1)}
  tibble(fivemers, tenmers)
}

decays <- tibble(
  frame = seq_len(500),
  linear0.5 = seq(1, 0.5, length.out = 500),
  linear0.75 = seq(1, 0.75, length.out = 500),
  linear0.9 = seq(1, 0.9, length.out = 500),
  linear0.95 = seq(1, 0.95, length.out = 500),
  linear0.99 = seq(1, 0.99, length.out = 500),
  exponential0.5 = exp(log(0.5) * seq(0, 1, length.out = 500)),
  exponential0.75 = exp(log(0.75) * seq(0, 1, length.out = 500)),
  exponential0.9 = exp(log(0.9) * seq(0, 1, length.out = 500)),
  exponential0.95 = exp(log(0.95) * seq(0, 1, length.out = 500)),
  exponential0.99 = exp(log(0.99) * seq(0, 1, length.out = 500)),
  power0.5 = seq_len(500) ^ log(0.5, 500),
  power0.75 = seq_len(500) ^ log(0.75, 500),
  power0.9 = seq_len(500) ^ log(0.9, 500),
  power0.95 = seq_len(500) ^ log(0.95, 500),
  power0.99 = seq_len(500) ^ log(0.99, 500)
)

sim.tifs <- setdiff(list.files(pattern = "^Sim.*tif$"),
                    list.files(pattern = "bleach.*tif$"))
sim.imgs <- map(sim.tifs, ReadImageData)
for (decay in names(decays)[-1]) {
  print(decay)
  bleached.imgs <- map(sim.imgs, ~ SimulateBleaching(., decays[[decay]]))
  map2(bleached.imgs, paste0(BeforeLastDot(sim.tifs),
                             "_bleach=", decay, ".tif"),
       ~ WriteIntImage(.x, .y))
}

decays %>% gather(line, relative_intensity, -1) %>%
  ggplot + aes(frame, relative_intensity, colour = line) + geom_line() +
  theme(legend.title = element_blank())
ggsave("decays.pdf", width = 7, height = 4)

BrightnessTxtFolder(tau = NA)
BrightnessTxtFolder(tau = 10)
BrightnessTxtFolder(tau = "auto", seed = 1)

csvs <- list.files(pattern = "Sim.*csv$")
master <- tibble(
  sim = csvs %>% StrAfterNth("Simulation", 1) %>%
    NthNumber(1, TRUE) %>% paste0("M"),
  replicate = csvs %>% StrAfterNth("rep", 1) %>% NthNumber(1),
  bleaching = map_chr(csvs, function(string) {
    if (!str_detect(string, "bleach")) return ("aaa")
    StrAfterNth(string, "bleach=", 1) %>%
      BeforeLastDot() %>%
      str_split("_") %>%
      unlist %>%
      nth(1)
  }),
  tau = csvs %>% StrAfterNth("tau", 1) %>% NthNumber(1),
  auto_tau = str_detect(csvs, "tau=auto"),
  brightness = csvs %>% lapply(ReadImageTxt) %>% map_dbl(mean, na.rm = TRUE)
) %>% arrange(sim, replicate, bleaching) %>%
  mutate(bleaching = if_else(bleaching == "aaa", "none", bleaching))
write_csv(master, "master.csv")

filter(master, auto_tau, bleaching != "none") %>%
  mutate(bleaching = paste0(str_sub(bleaching, 1, 3),
                            NthNumber(bleaching, 1, T, T, T, T))) %>%
  ggplot + aes(bleaching, tau) +
  geom_jitter(height = 0, width = 0.2) + theme(legend.position = "none")
ggsave("clust.pdf", width = 5, height = 3)

olig <- master %>%
  mutate(tau = ifelse(auto_tau, "auto", tau), auto_tau = NULL,
         bleaching = ifelse(bleaching == "none", "aaa", bleaching)) %>%
  group_by(bleaching, tau) %>%
  summarize(oligomerization = list(FiveTenMers(sim, brightness))) %>%
  ungroup %>%
  mutate(mean_5mer = map_dbl(oligomerization, ~ median(.[[1]], na.rm = TRUE)),
         pm95_5mer = map_dbl(oligomerization, ~ 1.96 * sd(.[[1]], na.rm = TRUE)),
         mean_10mer = map_dbl(oligomerization, ~ median(.[[2]], na.rm = TRUE)),
         pm95_10mer = map_dbl(oligomerization, ~ 1.96 * sd(.[[2]], na.rm = TRUE)),
         oligomerization = NULL,
         bleaching = ifelse(bleaching == "aaa", "none", bleaching))
write_csv(olig, "olig.csv")

set.seed(3)
fold.change.mix <- master %>%
  mutate(tau = ifelse(auto_tau, "auto", tau),
         auto_tau = NULL) %>%
  mutate(tau = ifelse(is.na(tau), "NA", tau)) %>%
  filter(map_lgl(NthNumber(bleaching, decimals = TRUE) %in% c(0.75, 0.95),
                 isTRUE),
         !map_lgl(bleaching, ~ isTRUE(str_detect(., "power")))) %>%
  group_by(tau) %>%
  summarize(fivetenmers = list(FiveTenMers(sim, brightness))) %>%
  mutate(fivemers = map(fivetenmers, 1), tenmers = map(fivetenmers, 2),
         fivetenmers = NULL) %>%
  unnest %>%
  gather(oligomer, fold_change, fivemers:tenmers) %>%
  na.omit %>%
  groupwiseMedian(data = ., var = "fold_change", group = c("oligomer", "tau"),
                  R = 10^5, conf = 0.95)
write_csv(fold.change.mix, "fold_change_mix.csv")

clust <- master %>%
  mutate(tau = ifelse(auto_tau, "auto", tau)) %>%
  group_by(sim, tau, bleaching) %>%
  summarise(mean_brightness = mean(brightness), sd_brightness = sd(brightness)) %>%
  ungroup %>%
  mutate(bleaching_type = NthNonNumeric(bleaching, 1),
         bleaching_extent = as.factor(NthNumber(bleaching, -1, decimals = TRUE)))
write_csv(clust, "clust.csv")

ggplot(clust, aes(bleaching_extent, mean_brightness,
                  colour = tau, shape = bleaching_type)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_point(position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(x = bleaching_extent, y = mean_brightness,
                    ymin = mean_brightness - sd_brightness,
                    ymax = mean_brightness + sd_brightness),
                position = position_dodge(width = 0.8)) +
  scale_y_continuous(breaks = 0:12) +
  guides(col = guide_legend(override.aes = list(shape = 15, size = 10)))
ggsave("clust.pdf", width = 6, height = 4)
