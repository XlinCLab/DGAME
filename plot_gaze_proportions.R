library(dplyr)
library(ggplot2)

plot_gaze_proportions <- function(response_time, median_d_onset, outfile = NULL) {
  plotti1 <- response_time %>%
    filter(!(condition == "no_conflict" & AOI == "aoi_comp")) %>%
    select(Prop, Time, subj, AOI, condition) %>%
    group_by(Time, AOI, subj, condition) %>%
    summarize_all(mean) %>%
    ungroup() %>%
    mutate(condition = case_when(condition == "conflict" ~ 'pair', TRUE ~ 'singleton')) %>%
    ggplot(aes(x = Time, y = Prop, color = AOI, linetype = condition)) +
    geom_smooth(method = 'loess', span = 0.13, lwd = 1, se = TRUE, level = 0.95) +
    theme_light() +
    geom_vline(xintercept = 0, linetype = "dotted") +
    geom_vline(xintercept = 471, linetype = "dotted") +
    geom_vline(xintercept = median_d_onset, linetype = "dotted") +
    #  geom_rect(data = data.frame(condition = "conflict"), aes(xmin = 1100, xmax = 1900, ymin = 0.51, ymax = 0.52), alpha = 1, fill="black", inherit.aes = FALSE,colour="black")+
    theme(
      axis.text = element_text(size = 15),
      axis.title = element_text(size = 15),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.text = element_text(size = 15),
      strip.text.x = element_text(size = 15),
      strip.text.y = element_text(size = 15),
      legend.background = element_rect(fill = alpha("white", 0))
    )

  if (!is.null(outfile)) {
    ggsave(filename = outfile, plot = plotti1, width = 8, height = 6, dpi = 300)
  }

  invisible(plotti1)
}
