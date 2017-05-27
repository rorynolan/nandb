required <- c("nandb", "tidyverse")
for (r in required) if (!require(r, character.only = TRUE)) install.packages(r)
invisible(sapply(required, library, character.only = TRUE))

tibble(time = 0:100, original = time, corrected = ExpSmooth(time, 40),
             naive = nandb:::ExpSmoothNaive(time, 40)) %>%
  gather(series, value, -time) %>%
  ggplot + aes(time, value, colour = series) + geom_line()
ggsave("dontbenaive.pdf", height = 3, width = 5)
