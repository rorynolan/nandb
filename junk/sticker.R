pkgs <- c("hexSticker", "tidyverse", "here", "magick")
invisible(lapply(pkgs, library, character.only = TRUE))

sticker(here("junk", "nandb.png"),
        package = "nandb",
        filename = here("junk", "sticker.png"),
        s_width = 0.7,
        s_x = 0.99, s_y = 0.9,
        p_y = 1.63,
        url = "github.com/rorynolan/nandb",
        u_x = 1.13, u_color = "black", u_size = 0.8,
        h_fill = "black", h_color = "grey")
image_read(here("junk", "sticker.png"))
