theme_ptf <- function(base_size = 22,
                      legend.position = "right",
                      panel.grid = TRUE,
                      panel.border = TRUE) {
    thm <- theme_bw(base_size = base_size) +
        theme(
            strip.text = element_text(size = base_size),
            axis.text = element_text(color = "black", size = base_size * 0.8),
            axis.title = element_text(color = "black", size = base_size),
            legend.position = legend.position,
            legend.title = element_text(size = base_size * 0.9),
            legend.text = element_text(size = base_size * 0.8),
            strip.background = element_rect(fill = "grey90", color = NA)
        )

    if (!panel.grid) {
        thm <- thm + theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
        )
    }

    if (!panel.border) {
        thm <- thm + theme(panel.border = element_blank())
    }

    return(thm)
}

text_base_size = 12
