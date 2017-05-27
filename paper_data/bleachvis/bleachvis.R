required <- c("nandb", "gridExtra", "tidyverse", "EBImage")
for (r in required) if (!require(r, character.only = TRUE)) install.packages(r)
invisible(sapply(required, library, character.only = TRUE))

img <- ReadImageData("main_before.tif")
means <- autothresholdr::mean_stack_thresh(img, "MinError") %>%
  apply(3, mean, na.rm = TRUE)
means.tau10 <- autothresholdr::mean_stack_thresh(img, "MinError") %>%
  CorrectForBleaching(10) %>%
  apply(3, mean, na.rm = TRUE)
set.seed(2)
best.tau <- autothresholdr::mean_stack_thresh(img, "MinError") %>% BestTau
means.tauauto <- autothresholdr::mean_stack_thresh(img, "MinError") %>%
  CorrectForBleaching(best.tau) %>%
  apply(3, mean, na.rm = TRUE)
bleach.vis <- tibble(frame = seq_along(means), `no detrend` = means,
                     `10` = means.tau10, auto = means.tauauto) %>%
  gather(tau, mean_intensity, -1)
ggplot(bleach.vis, aes(frame, mean_intensity, colour = tau)) + geom_line()
ggsave(filename = "bleachvis1.pdf", width = 6, height = 4)
bleach.vis %>%
  filter(tau != "no detrend") %>%
  ggplot + aes(frame, mean_intensity, colour = tau) + geom_line() +
  scale_colour_manual(values = scales::hue_pal()(3)[1:2])
ggsave(filename = "bleachvis2.pdf", width = 6, height = 4)
