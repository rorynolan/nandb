required <- c("nandb", "gridExtra", "tidyverse", "EBImage", "magrittr")
invisible(sapply(required, library, character.only = TRUE))

before <- ReadImageData("cell4_fkbp_before.tif")
after <- ReadImageData("cell4_fkbp_with.tif")
mean.intensity.before <- MeanPillars(before)
mean.intensity.after <- MeanPillars(after)
intensity.before.plot <- MatrixRasterPlot(MedianFilterB(mean.intensity.before),
                        colours = c("black", "white"),
                        na.colour = "red", limits = c(0, 60), clip = TRUE,
                        scale.name = "intensity") +
  theme(legend.position = "none")
intensity.after.plot <- MatrixRasterPlot(MedianFilterB(mean.intensity.after),
                         colours = c("black", "white"),
                         na.colour = "red", limits = c(0, 60), clip = TRUE,
                         scale.name = "intensity")
brightness.before <- ReadImageTxt(
  "cell4_fkbp_before_brightness_frames=500_tau=auto=41_mst=Triangle_filter=NA.csv")
brightness.before.plot <- MatrixRasterPlot(MedianFilterB(brightness.before),
                                        limits = c(0.9, 1.25), clip = T,
                                        scale.name = "brightness",
                                        include.breaks = 1, log.trans = T) +
  theme(legend.position = "none")
brightness.after <- ReadImageTxt(
  "cell4_fkbp_with_brightness_frames=500_tau=auto=125_mst=Triangle_filter=NA.csv")
brightness.after.plot <- MatrixRasterPlot(MedianFilterB(brightness.after),
                                         limits = c(0.9, 1.25), clip = T,
                                         scale.name = "brightness",
                                         include.breaks = 1, log.trans = T)

before.hex <- ArrArrHexPlot(MeanPillars(before),
                         brightness.before, limits = c(1, 500), log.trans = T) +
  xlab("intensity") + ylab("brightness") + ylim(0.75, 1.5) +
  theme(legend.position = "none", axis.title = element_text(size = 16),
        axis.text = element_text(size = 15))
after.hex <- ArrArrHexPlot(MeanPillars(after), brightness.after,
                          limits = c(1, 500), log.trans = T) +
  xlab("intensity") + ylab("brightness") + ylim(0.75, 1.5) +
  theme(legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.position = "right", axis.title = element_text(size = 16),
        axis.text = element_text(size = 15))

density.plot <- tibble(before = as.vector(brightness.before),
                     after = as.vector(brightness.after)) %>%
  gather(line, value, 1:2) %>%
  mutate(line = forcats::fct_relevel(line, c("before", "after"))) %>%
  ggplot + aes(value, colour = line) + geom_density() +
  xlim(0.8, 1.4) + scale_colour_manual(values = c("black", "red")) +
  geom_segment(x = mean(brightness.before, na.rm = TRUE),
               xend = mean(brightness.before, na.rm = TRUE), y = 0, yend = 6,
               colour = "black", linetype = "dashed") +
  geom_segment(x = mean(brightness.after, na.rm = TRUE),
               xend = mean(brightness.after, na.rm = TRUE), y = 0, yend = 6,
               colour = "red", linetype = "dashed") +
  geom_hline(yintercept = 0, colour = "white") + xlab("brightness")

grid.arrange(intensity.before.plot, intensity.after.plot,
             brightness.before.plot, brightness.after.plot,
             before.hex, after.hex,
             density.plot,
             layout_matrix = rbind(c(1, 2),
                                   c(3, 4),
                                   c(5, 6),
                                   c(7, 7)
             )
)
fig <- arrangeGrob(intensity.before.plot, intensity.after.plot,
            brightness.before.plot, brightness.after.plot,
            before.hex, after.hex,
            density.plot,
            layout_matrix = rbind(c(1, 2),
                                  c(3, 4),
                                  c(5, 6),
                                  c(7, 7)
            )
)
ggsave(file = "fig.pdf", fig, width = 2.75)

means.before <- autothresholdr::mean_stack_thresh(before, "tri") %>%
  apply(3, mean, na.rm = TRUE)
means.before.tau10 <- autothresholdr::mean_stack_thresh(before, "tri") %>%
  CorrectForBleaching(10) %>%
  apply(3, mean, na.rm = TRUE)
means.before.tauauto <- autothresholdr::mean_stack_thresh(before, "tri") %>%
  CorrectForBleaching(41) %>%
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
