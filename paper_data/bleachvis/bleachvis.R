img <- ReadImageData("bleachvis.tif")
means <- autothresholdr::mean_stack_thresh(before, "MinError") %>%
  apply(3, mean, na.rm = TRUE)
means.tau10 <- autothresholdr::mean_stack_thresh(before, "MinError") %>%
  CorrectForBleaching(10) %>%
  apply(3, mean, na.rm = TRUE)
means.before.tauauto <- autothresholdr::mean_stack_thresh(before, "MinError") %>%
  CorrectForBleaching(160) %>%
  apply(3, mean, na.rm = TRUE)
bleach.vis <- tibble(frame = seq_along(means.before), `no detrend` = means.before,
                     `10` = means.before.tau10, auto = means.before.tauauto) %>%
  gather(tau, mean_intensity, -1)
ggplot(bleach.vis, aes(frame, mean_intensity, colour = tau)) + geom_line()
ggsave(filename = "bleachvis1.pdf", width = 6, height = 4)
bleach.vis %>%
  filter(tau != "no detrend") %>%
  ggplot + aes(frame, mean_intensity, colour = tau) + geom_line() +
  scale_colour_manual(values = scales::hue_pal()(3)[1:2])
ggsave(filename = "bleachvis2.pdf", width = 6, height = 4)
