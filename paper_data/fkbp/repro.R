lapply(c("nandb", "filesstrings", "stringr", "tidyverse", "xtable",
         "magrittr", "EBImage", "autothresholdr", "rcompanion", "clipr"),
       library, character.only = TRUE)

FoldChange <- function(drug, brightness) {
  epsilon.before <- brightness[drug == "before"] - 1
  epsilon.after <- brightness[drug != "before"] - 1
  epsilon.after / epsilon.before
}

BrightnessTxtFolder(tau = NA, mst = "MinError")
BrightnessTxtFolder(tau = 10, mst = "MinError")
BrightnessTxtFolder(tau = "auto", mst = "MinError", seed = 0)

csvs <- list.files(pattern = "brightness.*csv$")

master <- tibble(
  well = NthNumber(csvs, 1),
  cell = csvs %>% StrAfterNth("cell") %>% NthNumber,
  drug = ifelse(str_detect(csvs, "before"), "before", "after"),
  concentration = ifelse(str_detect(csvs, "after"),
                         csvs %>%
                           StrBeforeNth("nM") %>%
                           NthNumber(-1, decimals = TRUE),
                         0),
  tau = csvs %>% StrAfterNth("tau", 1) %>% NthNumber(1),
  auto_tau = str_detect(csvs, "tau=auto"),
  brightness = csvs %>% lapply(ReadImageTxt) %>% map_dbl(mean, na.rm = TRUE)
)
write_csv(master, "master.csv")

master %>% mutate(tau = ifelse(auto_tau, "auto", tau), auto_tau = NULL) %>%
  filter(concentration > 0) %>%
  mutate(concentration = as.factor(concentration),
         cell = as.factor(cell)) %>%
  ggplot + aes(tau, brightness, colour = concentration, shape = cell) +
  geom_point(position = position_jitterdodge()) + ylim(0.9, 2.5)+
  guides(col = guide_legend(override.aes = list(shape = 15, size = 5))) +
  scale_y_continuous(breaks = seq(1, 2.1, 0.2))
ggsave("compareconcs.pdf", width = 6, height = 3.5)
fold_changes.sum <- master %>%
  mutate(tau = ifelse(auto_tau, "auto", tau), auto_tau = NULL) %>%
  group_by(tau, well, cell) %>%
  summarise(fold_change = FoldChange(drug, brightness),
            concentration = unique(concentration[concentration != 0])) %>%
  group_by(tau) %>%
  summarise(median_fold_change = median(fold_change),
            mad_fold_change = mad(fold_change))
write_csv(fold_changes.sum, "fold_changes_sum.csv")

fold_changes.cells <- master %>%
  mutate(tau = ifelse(auto_tau, "auto", tau), auto_tau = NULL) %>%
  group_by(well, cell, tau) %>%
  summarise(fold_change = FoldChange(drug, brightness),
            concentration = unique(concentration[concentration != 0])) %>%
  ungroup %>%
  arrange(tau)
write_csv(fold_changes.sum, "fold_changes_cells.csv")
sumstat <- groupwiseMedian(data = filter(fold_changes.cells), var = "fold_change",
                           group = "tau", bca = TRUE, R = 10^5)
fold_changes.cells %>% mutate(concentration = as.factor(concentration),
                              cell = as.factor(cell)) %>%
  ggplot + aes(tau, fold_change) +
  geom_boxplot() + geom_jitter() +
  scale_y_continuous("fold change", breaks = -1:8, limits = c(-1.5, 9))

master %>%
  mutate(tau = ifelse(auto_tau, "auto", tau), auto_tau = NULL) %>%
  ggplot + aes(tau, brightness, colour = before) + geom_jitter()

tifs <- list.files(pattern = "tif$")
print(tifs[i])
ReadImageData(tifs[i]) %>% normalize %>% display(method = "r")
ReadImageData(tifs[i]) %>% mean_stack_thresh("MinE") %>%
  normalize %>% display(method = "r")
i <- i + 1

tifs <- list.files(pattern = "tif$")
ReadImageData(tifs[i]) %>% mean_stack_thresh("Huang") %>% median(na.rm = T)
i <- i + 1