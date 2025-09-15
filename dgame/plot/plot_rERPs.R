#!/usr/bin/env Rscript
library(tidyverse)

create_noun_rERP_plots <- function(nouns_data, plot_outdir = NULL) {
    # Load dataframe from file if filepath is passed as input arg
    if (is.character(nouns_data)) {
        nouns_data <- read.csv(nouns_data)
    }
    if (!is.data.frame(nouns_data)) {
        stop("Input must be a data frame or path to a CSV file.")
    }

    # Plot the rERPs time-locked to noun onset by condition and ROI from -200 to 1000ms
    noun_condition_plot <- nouns_data %>%
      drop_na() %>%
      mutate(laterality = fct_relevel(laterality,'left','central','right'),
             saggitality = fct_relevel(saggitality,'prefrontal','frontal','central','posterior','occipital')
      ) %>% 
      filter(time >= -200 & time <= 1000) %>% 
      mutate(time = (time)/1000) %>%
      group_by(saggitality,laterality,time,condition,subject) %>%
      summarize_all(mean) %>% 
      ungroup() %>% 
      select(-subject) %>% 
      group_by(saggitality,laterality,time,condition) %>%
      summarize_all(mean) %>%
      mutate(condition = case_when(condition == "conflict" ~ 'pair',
                                   TRUE ~ 'singleton')) %>% 
      ggplot(aes(x=time,y=data,color =condition,))+
      facet_grid(saggitality~laterality)+
      geom_smooth(method='loess',span = 0.1,lwd = 0.5,level=0.99,se = T)+
      theme_minimal()+
      geom_vline(xintercept = 0.471,linetype="dotted")+
      theme(panel.background = element_rect(fill = 'white', colour = 'grey50'),
            axis.text=element_text(size=16),
            axis.title=element_text(size=16),
            legend.title=element_blank(), 
            legend.text=element_text(size=16),
            strip.text.x = element_text(size = 16),
            strip.text.y = element_text(size = 16),
            legend.position = "bottom",
            legend.direction="horizontal",
            panel.spacing = unit(2,"lines")
      )+
      guides(color = guide_legend(nrow=1))+
      scale_y_reverse()+
      geom_vline(xintercept = 0) +
      geom_hline(yintercept = 0) +
      geom_vline(xintercept = 0,linetype="dotted")+
      scale_color_manual(values=c('red','blue'))


    # Plot nouns by mean fixation time and condition along the saggital axis from -200 to 1000 ms
    noun_fixation_time_plot = nouns_data %>% 
      drop_na() %>% 
      select(data,time,subject,mean_target_fixation,saggitality,condition) %>%
      mutate(
        saggitality = fct_relevel(saggitality,'prefrontal','frontal','central','posterior','occipital')
      ) %>% 
      filter(time >= -200 & time <= 1000) %>% 
      mutate(time = (time)/1000) %>%
      mutate(mean_target_fixation = fct_relevel(mean_target_fixation, "before_noun","during_noun","after_noun")) %>% 
      group_by(time,subject,mean_target_fixation,saggitality,condition) %>%
      summarize_all(mean) %>% 
      ungroup() %>% 
      select(-subject) %>% 
      group_by(time,mean_target_fixation,saggitality,condition) %>%
      summarize_all(mean) %>%
      mutate(condition = case_when(condition == "conflict" ~ 'pair',
                                  TRUE ~ 'singleton')) %>% 
      ggplot(aes(x=time,y=data,color = condition, linetype = mean_target_fixation))+
      facet_grid(saggitality~condition)+
      geom_smooth(method='loess',span = 0.1,lwd = 0.5,level=0.99,se = T)+
      theme_minimal()+
      geom_vline(xintercept = 0.471,linetype="dotted")+
      theme(panel.background = element_rect(fill = 'white', colour = 'grey50'),
            axis.text=element_text(size=16),
            axis.title=element_text(size=16),
            legend.title=element_blank(), 
            legend.text=element_text(size=16),
            strip.text.x = element_text(size = 16),
            strip.text.y = element_text(size = 16),
            legend.position = "bottom",
            legend.direction="horizontal",
            panel.spacing = unit(2,"lines")
      )+
      guides(color = guide_legend(nrow=1))+
      scale_y_reverse()+
      geom_vline(xintercept = 0) +
      geom_hline(yintercept = 0) +
      geom_vline(xintercept = 0,linetype="dotted")+
      scale_color_manual(values=c('red','blue'))

    # Save plots to output directory
    if (!is.null(plot_outdir)) {
        ggsave(
            filename = file.path(plot_outdir, "noun_rERPs_condition.png"),
            plot = noun_condition_plot,
            width = 16,
            height = 12,
            dpi = 300
    )
        ggsave(
            filename = file.path(plot_outdir, "noun_rERPs_fixation-time.png"),
            plot = noun_fixation_time_plot,
            width = 16,
            height = 12,
            dpi = 300
    )
  }

  invisible(noun_condition_plot)
  invisible(noun_fixation_time_plot)
        
}

# Parse command line arguments:
# 1) path to noun rERP data CSV
# 2) path to fixation time rERP data CSV
# 3) path to output directory for plots
args <- commandArgs(trailingOnly = TRUE)
noun_data <- args[1]
fixation_data <- args[2]
plot_outdir <- args[3]
create_noun_rERP_plots(noun_data, plot_outdir)
#create_fixation_rERP_plots(fixation_data, plot_outdir)


# plotti5 =fixation_data %>% 
#   drop_na() %>% 
#   ungroup() %>% 
#   select(data,time,subject,fix_at,condition,channel) %>%
#   filter(channel %in% c("Fpz","Cz","Oz")) %>% 
#   mutate(channel = fct_relevel(channel,'Fpz','Cz','Oz')) %>% 
#    filter(time >= -200 & time <= 1000) %>% 
#   #filter(fix_at == "target" | fix_at == "other") %>% 
#   mutate(time = (time)/1000) %>%
#   group_by(subject,time,fix_at,condition,channel) %>%
#   summarize_all(mean) %>% 
#   ungroup() %>% 
#   select(-subject) %>% 
#   group_by(time,fix_at,condition,channel) %>%
#   summarize_all(mean) %>% 
#   mutate(condition = case_when(condition == "conflict" ~ 'pair',
#                                TRUE ~ 'singleton')) %>% 
#   ggplot(aes(x=time,y=data,color = condition))+
#   facet_grid(channel~fix_at) + 
#   geom_smooth(method='loess',span = 0.15,lwd = 0.5,level=0.99,se = T)+
#   theme_minimal()+
#    theme(panel.background = element_rect(fill = 'white', colour = 'grey50'),
#         axis.text=element_text(size=16),
#         axis.title=element_text(size=16),
#         legend.title=element_blank(), 
#         legend.text=element_text(size=16),
#         strip.text.x = element_text(size = 16),
#         strip.text.y = element_text(size = 16),
#         legend.position = "bottom",
#         legend.direction="horizontal",
#         panel.spacing = unit(2,"lines")
#   )+
#   guides(color = guide_legend(nrow=1))+
#   scale_y_reverse()+
#   geom_vline(xintercept = 0) +
#   geom_hline(yintercept = 0) +
#   geom_vline(xintercept = 0,linetype="dotted")+
#   scale_color_manual(values=c('red','blue'))


# plotti6 =fixation_data %>% 
#   filter(fix_at == "target") %>% 
#   drop_na() %>% 
#   ungroup() %>% 
#   select(data,time,subject,fix_at,condition,laterality,saggitality) %>%
#   mutate(laterality = fct_relevel(laterality,'left','central','right'),
#          saggitality = fct_relevel(saggitality,'prefrontal','frontal','central','posterior','occipital')
#   ) %>% 
#   filter(time >= -200 & time <= 1000) %>% 
#   mutate(time = (time)/1000) %>%
#   group_by(subject,time,fix_at,condition,laterality,saggitality) %>%
#   summarize_all(mean) %>% 
#   ungroup() %>% 
#   select(-subject) %>% 
#   group_by(time,fix_at,condition,laterality,saggitality) %>%
#   summarize_all(mean) %>% 
#   mutate(condition = case_when(condition == "conflict" ~ 'pair',
#                                TRUE ~ 'singleton')) %>% 
#   ggplot(aes(x=time,y=data,color = condition))+
#   facet_grid(saggitality~laterality) + 
#   geom_smooth(method='loess',span = 0.1,lwd = 0.5,level=0.99,se = T)+
#   theme_minimal()+
#   theme(panel.background = element_rect(fill = 'white', colour = 'grey50'),
#         axis.text=element_text(size=16),
#         axis.title=element_text(size=16),
#         legend.title=element_blank(), 
#         legend.text=element_text(size=16),
#         strip.text.x = element_text(size = 16),
#         strip.text.y = element_text(size = 16),
#         legend.position = "bottom",
#         legend.direction="horizontal",
#         panel.spacing = unit(2,"lines")
#   )+
#   guides(color = guide_legend(nrow=1))+
#   scale_y_reverse()+
#   geom_vline(xintercept = 0) +
#   geom_hline(yintercept = 0) +
#   geom_vline(xintercept = 0,linetype="dotted")+
#   scale_color_manual(values=c('red','blue'))



# plotti7 =fixation_data %>% 
#   drop_na() %>% 
#   ungroup() %>% 
#   filter(fix_at == "target") %>%
#    filter(channel %in% c("Fpz","Cz","Oz")) %>% 
#   mutate(channel = fct_relevel(channel,'Fpz','Cz','Oz')) %>% 
#   select(data,time,subject,fix_time,condition,channel) %>%
#   filter(time >= -200 & time <= 1000) %>% 
#   mutate(time = (time)/1000) %>%
#   select(data,time,subject,fix_time,condition,channel) %>%
#   group_by(subject,time,fix_time,condition,channel) %>%
#   summarize_all(mean) %>% 
#   ungroup() %>% 
#   select(-subject) %>% 
#   group_by(time,fix_time,condition,channel) %>%
#   summarize_all(mean) %>% 
#   mutate(fix_time = case_when(fix_time == ">after_noun" ~ "after_noun",
#                               TRUE ~ fix_time)) %>% 
#   mutate(condition = case_when(condition == "conflict" ~ 'pair',
#                                TRUE ~ 'singleton')) %>% 
#   ggplot(aes(x=time,y=data,linetype=fix_time,color = condition))+
#   facet_grid(channel~condition) + 
#   geom_smooth(method='loess',span = 0.1,lwd = 0.5,level=0.99,se = T)+
#   theme_minimal()+
#   theme(
#         axis.text=element_text(size=16),
#         axis.title=element_text(size=16),
#         legend.title=element_blank(), 
#         legend.text=element_text(size=16),
#         strip.text.x = element_text(size = 16),
#         strip.text.y = element_text(size = 16),
#         legend.position = "bottom",
#         legend.direction="horizontal",
#         panel.spacing = unit(2,"lines")
#   )+
#   guides(color = guide_legend(nrow=1))+
#   scale_y_reverse()+
#   geom_vline(xintercept = 0) +
#   geom_hline(yintercept = 0) +
#   geom_vline(xintercept = 0,linetype="dotted")+
#   scale_color_manual(values=c('red','blue'))
