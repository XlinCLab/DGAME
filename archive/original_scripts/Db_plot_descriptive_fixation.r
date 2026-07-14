#get first/last target/comp gaze for trials
library(ggplot2)
library(tidyverse)
library(patchwork)

# FUNCTIONS
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

# median CI

median_cl_boot <- function(x, conf = 0.83) {
  lconf <- (1 - conf)/2
  uconf <- 1 - lconf
  require(boot)
  bmedian <- function(x, ind) median(x[ind])
  bt <- boot(x, bmedian, 1000)
  bb <- boot.ci(bt, type = "perc")
  data.frame(y = median(x), ymin = quantile(bt$t, lconf), ymax = quantile(bt$t, 
                                                                          uconf))
}


script_path <- rstudioapi::getActiveDocumentContext()$path
dir <- dirname(script_path)
setwd(dir)  #set working directory to script path

n_s = '27'
dat = tibble()
for(s in n_s){
  cat(paste(s))
#  setwd(paste('/media/ingmar/ieeg_data/experiments/DGame2/preproc/eyetracking/',s,'/',sep="")) #only for full data
  s_fix = tibble()
  pat = paste('fixations_times_.*_trials.csv')
  
  file_list = list.files(pattern = pat ,recursive=T) 
  for(file in file_list){
    tmp = read.csv(file)
    # Combine the datasets into a single data frame, exception for subject 29
    if(s == '29'){
      if(nrow(tmp > 0)){
        tmp$y_scaled = NA
        tmp = tmp %>% select(-director,-X44)
      }
      s_fix = rbind(s_fix,tmp)
    }
    gc()
  }
  dat = rbind(dat,s_fix)
}

dat = dat %>% 
  mutate(fix_at = as.factor(case_when(
    aoi_target == TRUE ~ "target",
    aoi_otherTarget == TRUE ~ "otherTarget",
    aoi_comp == TRUE ~ "comp",
    aoi_otherComp == TRUE ~ 'otherComp',
    aoi_fillerA == TRUE ~ "fillerA",
    aoi_fillerB == TRUE ~ "fillerB",
    aoi_empty == TRUE ~ 'empty',
    aoi_other == TRUE ~ 'other',
    aoi_goal == TRUE ~ 'goal',
    TRUE ~ 'elsewhere'
  )
  )
  )
fixes = tibble()
for(s in unique(dat$subj)){
  tmp.fix = dat %>% #filter(trial_time >=-1) %>% 
    filter(subj==s)
  for(t in unique(tmp.fix$trial)){
    tmp.trial.fix = tmp.fix %>% 
      filter(trial == t )
    condition = as.character(tmp.trial.fix$condition[1])
    tmp.tar = tmp.trial.fix %>% filter(aoi_target == TRUE)
    tmp.comp = tmp.trial.fix %>% filter(aoi_comp == TRUE)
    object = "target"
    if(nrow(tmp.tar)>0){
      dur_max = max(tmp.tar$duration)
      dur_min = min(tmp.tar$duration)
      dur_median = median(tmp.tar$duration)
      
    } else {
      dur_min = NA
      dur_max = NA
      dur_median = NA
    }
    tmp1 = cbind(s,condition,object,dur_min,dur_max,dur_median)
    object = "comp"
    if(nrow(tmp.comp)>0){
      tmin = min(tmp.comp$trial_time)
      tmax = max(tmp.comp$trial_time)
      tmedian = median(tmp.comp$trial_time)
    } else {
      tmin = NA
      tmax = NA
      tmedian = NA
    }
    tmp2 = cbind(s,condition,object,tmin,tmax,tmedian)
    
    tmp = rbind(tmp1,tmp2)
    fixes = rbind(fixes,tmp)
    
  }
  
}


fixes2 = fixes %>% 
  mutate(tmin= as.numeric(as.character(tmin)),
         tmax = as.numeric(as.character(tmax)),
         tmedian = as.numeric(as.character(tmedian)),
  ) %>%
  as_tibble() %>% 
  drop_na() %>% 
  pivot_longer(cols=names(fixes)[4:ncol(fixes)],
               names_to='value',
               values_to='time')

fixes_mean = fixes %>% 
  drop_na() %>% 
  group_by(s,condition,object) %>% 
  summarize_all(mean)

fixes_grand_mean = fixes_mean %>% 
  drop_na() %>% 
  group_by(condition) %>% select(-s) %>% 
  summarize_all(median)

labels = fixes2 %>%
  drop_na() %>% 
  mutate(value=factor(value,levels=c("tmin","tmax","tmedian"),labels = c("min","max","median")),
         condition = factor(condition,levels=c("confl","no_confl",labels=c("conflict","no-conflict")))) %>% 
  group_by(object,condition,value) %>% 
  summarise_all(mean)
#duration histogram
plot1 <- dat %>% 
  
  select(subj,duration,condition,fix_at) %>%
  filter(fix_at %in% c("comp","target","goal")) %>% 
  mutate(duration = duration/1000) %>% 
  group_by(subj,condition,fix_at) %>%
  mutate(condition = case_when(condition == "conflict" ~ "pair",
                               TRUE ~ "singleton")) %>% 
  drop_na() %>% 
  ggplot(.,aes(x=duration, color=condition, fill=condition)) +
  scale_x_continuous(name="duration (s)")+
  geom_histogram(position="dodge", alpha=0.2,binwidth=0.025)+
  facet_wrap(fix_at~.)+
  theme(
    axis.text=element_text(size=12),
    axis.title=element_text(size=12),
    legend.title=element_blank(), 
    legend.text=element_text(size=12),
    strip.text.x = element_text(size = 12),
    strip.text.y = element_text(size = 12),
    legend.position = "bottom",
    legend.direction="horizontal",
    panel.spacing = unit(2,"lines")
  )+
  scale_color_manual(values=c('red','blue'))+
  scale_fill_manual(values=c('red','blue'))

#trial time histogram
plot2 <-dat %>% 
  select(subj,trial_time,condition,fix_at) %>%
  filter(fix_at %in% c("comp","target","goal")) %>% 
  group_by(subj,condition,fix_at) %>%
  mutate(condition = case_when(condition == "conflict" ~ "pair",
                               TRUE ~ "singleton")) %>% 
  drop_na() %>% 
  ggplot(.,aes(x=trial_time, color=condition, fill=condition)) +
  scale_x_continuous(name="trial time (s)", limits=c(-3.5,3.5), breaks = c(-3.5,0,3.5))+
  geom_histogram(alpha=0.2,position = "dodge",binwidth = 0.25)+
  
  facet_wrap(fix_at~.)+
  theme(
    axis.text=element_text(size=12),
    axis.title=element_text(size=12),
    legend.title=element_blank(), 
    legend.text=element_text(size=12),
    strip.text.x = element_text(size = 12),
    strip.text.y = element_text(size = 12),
    legend.position = "bottom",
    legend.direction="horizontal",
    panel.spacing = unit(2,"lines")
  )+
  scale_color_manual(values=c('red','blue'))+
  scale_fill_manual(values=c('red','blue'))

#by saccadic amplitude
plot3 <- dat %>% 
  select(subj,saccAmpl,condition,fix_at) %>%
  filter(fix_at %in% c("comp","target","goal")) %>% 
  group_by(subj,condition,fix_at) %>%
  mutate(condition = case_when(condition == "conflict" ~ "pair",
                               TRUE ~ "singleton")) %>% 
  group_by(subj,condition,fix_at) %>%
  drop_na() %>% 
  ggplot(.,aes(x=saccAmpl, color=condition, fill=condition)) +
  scale_x_continuous(name="normalized saccadic amplitude")+
  geom_histogram(position="identity", alpha=0.2)+
  facet_wrap(fix_at~.)+
  theme(
    axis.text=element_text(size=12),
    axis.title=element_text(size=12),
    legend.title=element_blank(), 
    legend.text=element_text(size=12),
    strip.text.x = element_text(size = 12),
    strip.text.y = element_text(size = 12),
    legend.position = "bottom",
    legend.direction="horizontal",
    legend.key = element_blank(),
    panel.spacing = unit(2,"lines")
  )+
  scale_color_manual(values=c('red','blue'))+
  scale_fill_manual(values=c('red','blue'))

#compute radians of saccades

dat$dy = NA
dat$dx = NA
# Compute deltas
dat$dx[2:nrow(dat)] <- diff(dat$norm_pos_x)
dat$dy[2:nrow(dat)] <- diff(dat$norm_pos_y)


# Compute angles in radians
dat$angles <- atan2(dat$dy, dat$dx)

# Convert to degrees for easier interpretation
dat$angles_deg <- dat$angles * (360 / pi)



plot4 <- dat %>% 
  select(angles_deg,condition,fix_at) %>% 
  drop_na() %>% 
  group_by(condition,fix_at) %>% 
  mutate(condition = case_when(condition == "conflict" ~ "pair",
                               TRUE ~ "singleton")) %>% 
  filter(fix_at %in% c("goal","elsewhere","target")) %>% 
  ggplot(.,aes(x = angles_deg, fill = condition,color=condition)) +
  geom_histogram(binwidth = 10, position = "dodge",alpha=0.2) +
  coord_polar(start = 300,direction = -1) +
  facet_wrap(~ fix_at) +
  scale_color_manual(values=c('red','blue'))+
  scale_fill_manual(values=c('red','blue'))+
  theme_minimal() +
  labs(
    x = "Angle (degrees)",
    y = "Frequency",
    title = "Angular Distribution of Saccades",
    fill = "condition"
  )+
  scale_x_continuous(limits = c(0,360),
                     breaks = seq(0, 270, by = 90),
                     minor_breaks = seq(0, 360, by = 30)) +
  theme(
    axis.text=element_text(size=12),
    axis.title=element_text(size=12),
    legend.title=element_blank(), 
    legend.text=element_text(size=12),
    strip.text.x = element_text(size = 12),
    strip.text.y = element_text(size = 12),
    legend.position = "bottom",
    legend.direction="horizontal",
    panel.spacing = unit(2,"lines")
  )



combined_plot <- (plot1 | plot2) /
  (plot3 | plot4) +
  plot_annotation(tag_levels = "A")  # tags A, B, C, D
