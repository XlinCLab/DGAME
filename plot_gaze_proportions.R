library(dplyr)
library(ggplot2)

plot_ti1 <- function(response_time, median_d_onset, outfile = NULL) {
  plotti1 <- response_time %>%
    filter(!(condition == "no_conflict" & AOI == "aoi_comp")) %>%
    select(Prop, Time, subj, AOI, condition) %>%
    group_by(Time, AOI, subj, condition) %>%
    summarize_all(base::mean) %>%
    ungroup() %>%
    mutate(condition = case_when(condition == "conflict" ~ "pair", TRUE ~ "singleton")) %>%
    ggplot(aes(x = Time, y = Prop, color = AOI, linetype = condition)) +
    geom_smooth(method = "loess", span = 0.13, lwd = 1, se = TRUE, level = 0.95) +
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


plot_ti3 <- function(response_time_comp, median_n_offset, outfile = NULL) {
  plotti3 <- response_time_comp %>%
    select(Prop, Time, aoi_fct, subj) %>%
    group_by(Time, aoi_fct, subj) %>%
    summarize_all(base::mean) %>%
    ungroup() %>%
    dplyr::rename(AOI = aoi_fct) %>%
    ggplot(aes(x = Time, y = Prop, color = AOI)) +
    geom_smooth(method = "loess", span = 0.1, lwd = 1, se = TRUE, level = 0.95) +
    theme_light() +
    geom_vline(xintercept = 0, linetype = "dotted") +
    geom_vline(xintercept = median_n_offset, linetype = "dotted") +
    theme(
      axis.text = element_text(size = 15),
      axis.title = element_text(size = 15),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.text = element_text(size = 15),
      strip.text.x = element_text(size = 15),
      strip.text.y = element_text(size = 15),
      legend.background = element_rect(fill = alpha("white", 0))
    ) +
    #annotate("rect",xmin = -2900, xmax = -2200, ymin = 0.049, ymax = 0.05,colour="black",fill = "black")+
    #annotate("rect",xmin = 1400, xmax = 3300, ymin = 0.049, ymax = 0.05,colour="black",fill = "black")+
    scale_color_manual(values=c('red','blue'))

  if (!is.null(outfile)) {
    ggsave(filename = outfile, plot = plotti3, width = 8, height = 6, dpi = 300)
  }

  invisible(plotti3)

}
