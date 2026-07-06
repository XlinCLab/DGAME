library(data.table)
library(tidyverse)
library(future.apply)
script_path <- rstudioapi::getActiveDocumentContext()$path
dir <- dirname(script_path)
setwd(dir)  #set working directory to script path
# Only set up parallel backend if you actually want to use it
# plan(multisession, workers = parallel::detectCores())
# 1) global options & working directory
options(contrasts = c("contr.treatment", "contr.Poly"), contr.sum.show.levels = TRUE)

# 2) load electrode coordinates once
coords_path = 'r_channel_positions.txt'
coords <- fread(coords_path, col.names = c('channel','lat','sag','z')) %>%
  mutate(channel = as.character(channel))

load_dataset <- function(
    data_dir,
    subject_ids = c('02'),#,'03','05','07','09','10','11','12','14','15','16','17','18','20','22','23','24','25','26','27','28','29','30','31','32','33'),
    mode = c("NOUN","FIX"),
    pattern_template = list(
      NOUN = "%s_.*_unfold_N\\.csv$",
      FIX    = "%s_.*_unfold_FIX\\.csv$"
    ),
    use_parallel = FALSE
) {
  mode <- match.arg(mode)
  
 
  
  # 3) list all CSV files once to avoid repeated directory scans
  all_files <- list.files(pattern = sprintf("_unfold_(FIX|N)\\.csv$"), recursive = TRUE)
  
  results <- vector("list", length(subject_ids))
  names(results) <- subject_ids
  
  # (optional) set up parallel plan if requested
  if (use_parallel) {
    library(future.apply)
    plan(multisession, workers = parallel::detectCores())
  }
  
  # 4) process each subject
  for (s in subject_ids) {
    message("Subject ", s, " (mode=", mode, ")…")
    
    # subset the global file list by subject and mode
    pat <- sprintf(pattern_template[[mode]], s)
    files <- grep(pat, all_files, value = TRUE)
    if (length(files) == 0) {
      warning("No files for subject ", s)
      next
    }
    
    # choose lapply or future_lapply based on use_parallel
    read_fun <- if (use_parallel) future.apply::future_lapply else lapply
    
    # load all CSVs
    tables <- read_fun(files, function(fp) {
      dt <- fread(fp)
      dt[, channel := as.character(channel)]
      dt
    })
    df <- bind_rows(tables)
    rm(tables); gc()
    
    # annotate with coordinates
    df <- df %>% left_join(coords, by = "channel")
    
    # compute laterality & sagittality
    df <- df %>%
      mutate(
        laterality = case_when(
          lat <  0           ~ "left",
          lat >  0           ~ "right",
          TRUE               ~ "central"
        ),
        saggitality = case_when(
          sag >  0 & sag <= 0.0714  ~ "frontal",
          sag >  0.0714             ~ "prefrontal",
          sag <  0 & sag >= -0.0929 ~ "posterior",
          sag < -0.0929             ~ "occipital",
          TRUE                       ~ "central"
        )
      )
    
    # mode-specific processing
    if (mode == "FIX") {
      df <- df %>%
        mutate(
          fix_time = case_when(
            trial_time < 0    ~ "before_noun",
            trial_time >= 0 & trial_time <= .5 ~ "during_noun",
            trial_time >  .471                 ~ "after_noun",
            TRUE                                ~ NA_character_
          )
        ) %>%
        select(-event, -trial_time) %>%
        group_by(time, condition, fix_at, fix_time,
                 saggitality, laterality,
                 channel) %>%
        summarize_all(mean, na.rm = TRUE) %>%
        ungroup()
    } else if(mode == "NOUN") {
      df <- df %>%
        mutate(mean_target_fixation = as.factor(case_when(mean_target_fixation < 0  ~ 'before_noun',
                                                          mean_target_fixation >= 0 & mean_target_fixation <= 0.471 ~ 'during_noun',
                                                          mean_target_fixation > 0.471 ~ 'after_noun'
                                                         )))
    } else {
      df <- df %>% mutate(subject = s)
    } 
    
    results[[s]] <- df
    gc()
  }
  
  # 5) combine & return
  bind_rows(results, .id = "subject")
}


 dat_n   <- load_dataset(
   data_dir    = paste(dir,'/eeg/unfold_out/results/',s,sep=""),
   mode        = "NOUN"
 ) %>% 
   mutate(mean_target_fixation = fct_relevel(mean_target_fixation,"before_noun","during_noun","after_noun" )) 
 dat_fix <- load_dataset(
   data_dir    = paste(dir,'/eeg/unfold_out/results/',s,sep=""),
   mode        = "FIX"
 ) %>% 
   mutate(fix_time = fct_relevel(fix_time,"before_noun","during_noun","after_noun" )) 
#plot the rERPs time-locked to noun onset by condition and ROI from -200 to 1000ms
plot_n_condition =dat_n %>% 
  drop_na() %>% 
  select(data,time,subject,saggitality,laterality,condition) %>%
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
  #rename condition for paper
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

gc()
plot_n_condition
#plot nouns by mean fixation time and condition along the saggital axis from -200 to 1000 ms
plot_n_fix_time =dat_n %>% 
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
  #rename condition for paper
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



gc()
plot_n_fix_time


plotti5 =dat_fix %>% 
  drop_na() %>% 
  ungroup() %>% 
  select(data,time,subject,fix_at,condition,channel) %>%
  filter(channel %in% c("Fpz","Cz","Oz")) %>% 
  mutate(channel = fct_relevel(channel,'Fpz','Cz','Oz')) %>% 
   filter(time >= -200 & time <= 1000) %>% 
  #filter(fix_at == "target" | fix_at == "other") %>% 
  mutate(time = (time)/1000) %>%
  group_by(subject,time,fix_at,condition,channel) %>%
  summarize_all(mean) %>% 
  ungroup() %>% 
  select(-subject) %>% 
  group_by(time,fix_at,condition,channel) %>%
  summarize_all(mean) %>% 
  mutate(condition = case_when(condition == "conflict" ~ 'pair',
                               TRUE ~ 'singleton')) %>% 
  ggplot(aes(x=time,y=data,color = condition))+
  facet_grid(channel~fix_at) + 
  geom_smooth(method='loess',span = 0.15,lwd = 0.5,level=0.99,se = T)+
  theme_minimal()+
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
plotti5
gc()

plotti6 =dat_fix %>% 
  filter(fix_at == "target") %>% 
  drop_na() %>% 
  ungroup() %>% 
  select(data,time,subject,fix_at,condition,laterality,saggitality) %>%
  mutate(laterality = fct_relevel(laterality,'left','central','right'),
         saggitality = fct_relevel(saggitality,'prefrontal','frontal','central','posterior','occipital')
  ) %>% 
  filter(time >= -200 & time <= 1000) %>% 
  mutate(time = (time)/1000) %>%
  group_by(subject,time,fix_at,condition,laterality,saggitality) %>%
  summarize_all(mean) %>% 
  ungroup() %>% 
  select(-subject) %>% 
  group_by(time,fix_at,condition,laterality,saggitality) %>%
  summarize_all(mean) %>% 
  mutate(condition = case_when(condition == "conflict" ~ 'pair',
                               TRUE ~ 'singleton')) %>% 
  ggplot(aes(x=time,y=data,color = condition))+
  facet_grid(saggitality~laterality) + 
  geom_smooth(method='loess',span = 0.1,lwd = 0.5,level=0.99,se = T)+
  theme_minimal()+
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
plotti6
gc()

plotti7 =dat_fix %>% 
  drop_na() %>% 
  ungroup() %>% 
  filter(fix_at == "target") %>%
   filter(channel %in% c("Fpz","Cz","Oz")) %>% 
  mutate(channel = fct_relevel(channel,'Fpz','Cz','Oz')) %>% 
  select(data,time,subject,fix_time,condition,channel) %>%
  filter(time >= -200 & time <= 1000) %>% 
  mutate(time = (time)/1000) %>%
  select(data,time,subject,fix_time,condition,channel) %>%
  group_by(subject,time,fix_time,condition,channel) %>%
  summarize_all(mean) %>% 
  ungroup() %>% 
  select(-subject) %>% 
  group_by(time,fix_time,condition,channel) %>%
  summarize_all(mean) %>% 
  mutate(fix_time = case_when(fix_time == ">after_noun" ~ "after_noun",
                              TRUE ~ fix_time)) %>% 
  mutate(condition = case_when(condition == "conflict" ~ 'pair',
                               TRUE ~ 'singleton')) %>% 
  ggplot(aes(x=time,y=data,linetype=fix_time,color = condition))+
  facet_grid(channel~condition) + 
  geom_smooth(method='loess',span = 0.1,lwd = 0.5,level=0.99,se = T)+
  theme_minimal()+
  theme(
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
plotti7
gc()
