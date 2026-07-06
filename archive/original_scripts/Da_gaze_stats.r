library(tidyverse)
library(eyetrackingR)

options(contrasts = c("contr.sum","contr.poly"),contr.sum.show.levels = TRUE)

script_path <- rstudioapi::getActiveDocumentContext()$path
dir <- dirname(script_path)
setwd(dir)  #set working directory to script path

gaze_infile = paste(dir,'/gaze_positions/gaze_positions_4analysis.csv',sep="")
gaze_positions_all = read.csv(file = gaze_infile) %>%
  filter(!is.na(condition)) #drop all non-trial (-3.5 to 3.5 seconds after noun onset) data points
names(gaze_positions_all)[29:44] = c('11','12','13','14','21','22','23','24','31','32','33','34','41','42','43','44') #modify names of the columns corresponding to each shelf compartment

gaze2analysis = gaze_positions_all %>% 
  filter(!is.na(trial_time) & trackloss == FALSE) %>% distinct()  %>% 
  mutate(duration = tmax-time) %>%
  mutate(duration = case_when(duration < 0 ~ FALSE,
                              TRUE ~ TRUE)) %>% 
  filter(duration == TRUE) %>% 
  select(-duration) %>% 
  filter(subj %in% c(2,3,5,7,9,10,11,12,14,15,16,17,18,20,22,23,24,25,26,27,28,29,30,31,32,33)) #only those that made it through EEG analysis.

all_words = tibble()
for(s in unique(gaze2analysis$subj)){
  trial = 1
  trial_d = 1
  for(sets in c(1,2)){
    for(pats in c(1,2)){
      if(s < 10){
        word_file_in = paste(dir,'/audio/0',s,'_words2erp_',sets,pats,'.csv',sep="")
        word_file_out = paste(dir,'/audio/0',s,'_words2erp_',sets,pats,'_trialtime.csv',sep="")
      }else{
        word_file_in =  paste(dir,'/audio/',s,'_words2erp_',sets,pats,'.csv',sep="")
        word_file_out =  paste(dir,'/audio/',s,'_words2erp_',sets,pats,'_trialtime.csv',sep="")
        
      }
      word_data = read.csv(word_file_in) %>% 
        mutate(trial = 0)
      
      for(t in 1:nrow(word_data)){
        if(word_data[t,]$pos == 'N'){
          word_data[t,]$trial = trial
          trial = trial+1
        }else if(word_data[t,]$pos == 'D'){
          word_data[t,]$trial = trial_d
          trial_d = trial_d+1  
        }
      }
      subj_data = gaze2analysis %>% filter(subj == s, aoi_target == TRUE,set == sets,pattern==pats) %>% 
        select(trial_time,trial) %>%
        group_by(trial) %>%
        summarize_all(mean)  
      word_data = left_join(word_data,subj_data)
      all_words = rbind(all_words,word_data)
      write.csv(word_data,word_file_out,row.names=FALSE)
    }
  }
  
}

median_d_onset = gaze2analysis %>% 
  select(pos,trial_time) %>%
  filter(pos == "D") %>%
  group_by(pos) %>%
  summarize_all(median) %>% 
  select(trial_time)
median_d_onset = as.numeric(median_d_onset[1,1])*1000

median_n_offset = gaze2analysis %>% filter(condition == "conflict" | condition == "no_conflict") %>% 
  mutate(duration = tmax-time) %>% 
  select(pos,duration) %>%
  filter(pos == "N") %>%
  group_by(pos) %>%
  summarize_all(median) %>% 
  select(duration)
median_n_offset = as.numeric(median_n_offset[1,1])*1000

#transform ET data into format for package
data <- make_eyetrackingr_data(gaze2analysis, 
                               participant_column = "subj",
                               trial_column = "trial",
                               time_column = "trial_time",
                               trackloss_column = "trackloss",
                               aoi_columns = c("aoi_target","aoi_comp","aoi_otherTarget","aoi_otherComp","aoi_fillerB","aoi_fillerA","aoi_goal"),
                               treat_non_aoi_looks_as_missing = F
)

data$trial_time = data$trial_time*1000
data$condition <- as.factor(data$condition)



#transform data for target AOI analysis
response_time <- make_time_sequence_data(data, 
                                           time_bin_size = 100,
                                           predictor_columns = c("condition"),
                                           aois = c('aoi_target','aoi_otherTarget','aoi_comp','aoi_otherComp',"aoi_fillerB","aoi_fillerA",'aoi_goal'),
                                           summarize_by = c("subj")
                                           
)

#statistical threshold for later
num_sub = length(unique((response_time$subj)))

#statistical threshold for later
threshold_t = qt(p = 1 - .05/2, 
                 df = num_sub-1) # pick threshold t based on alpha = .05 two tailed
num_bins = length(unique(response_time$TimeBin))


#data for competitor analysis
response_time_comp <- make_time_sequence_data(data, 
                                              time_bin_size = 100,
                                              predictor_columns = c("condition"),
                                              aois = c('aoi_comp','aoi_otherTarget','aoi_fillerB','aoi_fillerA','aoi_otherComp'),
                                              summarize_by = c("subj")
                                              
)



#mutate data for competitor analysis
response_time_comp = response_time_comp %>% 
  mutate(AOI = case_when(AOI != "aoi_comp" ~ "aoi_AllOther",
                         TRUE ~ AOI))

tmp1 = response_time_comp %>%  filter(AOI == "aoi_comp")
tmp2 = response_time_comp %>%  filter(AOI == "aoi_AllOther") %>% group_by(subj,condition,TimeBin,Time,AOI) %>% summarise_all(mean)

response_time_comp = rbind(tmp1,tmp2)
rm(tmp1,tmp2)
# bootstrap for comp
response_time_comp$aoi_fct = response_time_comp$AOI
response_time_comp = response_time_comp %>% filter(condition == "conflict")
response_time_comp$AOI = "dummy"

response_time_comp$aoi_fct = as.factor(response_time_comp$aoi_fct)

response_time_comp = response_time_comp %>% filter(aoi_fct == "aoi_comp" | aoi_fct == "aoi_AllOther")
#this won't work well with only 1 participant of course....
#test for targets
time_cluster_data_target <- make_time_cluster_data(data = response_time,
                                                   predictor_column = "condition",
                                                   aoi = "aoi_target",
                                                   test = "t.test",
                                                   paired = T,
                                                   threshold = threshold_t
)

clust_analysis_target <- analyze_time_clusters(time_cluster_data_target,samples=4000,within_subj = T,paired=T)

#test for competitors
time_cluster_data_comp <- make_time_cluster_data(data = response_time_comp,
                                                 predictor_column = "aoi_fct",
                                                 aoi = "dummy",
                                                 test = "t.test",
                                                 threshold = threshold_t,
                                                 treatment_level = "aoi_AllOther",
                                                 paired=TRUE
)
clust_analysis_comp <- analyze_time_clusters(time_cluster_data_comp,within_subj = TRUE,paired=TRUE,
                                             samples=4000)



#test for goal

time_cluster_data_goal <- make_time_cluster_data(data = response_time,
                                                 predictor_column = "condition",
                                                 aoi = "aoi_goal",
                                                 test = "t.test",
                                                 threshold = threshold_t,
                                                 treatment_level = "no_conflict",
                                                 paired=TRUE
)


clust_analysis_goal <- analyze_time_clusters(time_cluster_data_goal,within_subj = TRUE,paired=TRUE,
                                             samples=4000)

time_cluster_data_comp <- make_time_cluster_data(data = response_time_comp,
                                                 predictor_column = "timepoint",
                                                 aoi = "aoi_comp",
                                                 test = "t.test",
                                                 paired = T,
                                                 threshold = threshold_t
)


# plot results
plotti1 = response_time %>% 
  #filter(subj %in% c(3,7,12,20,18,23,24,27,31,33)) %>% 
  filter(!(condition == "no_conflict" & AOI == "aoi_comp")) %>% 
  select(Prop,Time,subj,AOI,condition) %>% 
  group_by(Time,AOI,subj,condition) %>% 
  summarize_all(base::mean) %>% 
  ungroup() %>% 
  mutate(condition = case_when(condition == "conflict" ~ 'pair',
                               TRUE ~ 'singleton')) %>% 
  ggplot(aes(x = Time, y = Prop,color=AOI,linetype=condition))+
  geom_smooth(method='loess',span = 0.13,lwd = 1,se = T,level=0.95)+ #je tiefer desto stC$rker smooth+
  theme_light() +
  geom_vline(xintercept = 0,linetype="dotted")+
  geom_vline(xintercept = 471,linetype="dotted")+
  geom_vline(xintercept = median_d_onset,linetype="dotted")+
 #for this one participant, this doesn't make sense, so the next line is commented out (significant time clusters)
  #  geom_rect(data = data.frame(condition = "conflict"), aes(xmin = 1100, xmax = 1900, ymin = 0.51, ymax = 0.52), alpha = 1, fill="black", inherit.aes = FALSE,colour="black")+
   theme(
    axis.text=element_text(size=15),
    axis.title=element_text(size=15),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.text=element_text(size=15),
    strip.text.x = element_text(size = 15),
    strip.text.y = element_text(size = 15),
    legend.background=element_rect(fill = alpha("white", 0))
  )#+

plotti1


# plot results
plotti3 = response_time_comp %>% select(Prop,Time,aoi_fct,subj) %>% 
  group_by(Time,aoi_fct,subj) %>% 
  summarize_all(base::mean) %>% 
  ungroup() %>% 
  dplyr::rename(AOI = aoi_fct) %>% 
  ggplot(aes(x = Time, y = Prop,color=AOI))+
  geom_smooth(method='loess',span = 0.1,lwd = 1,se = T,level=0.95)+ #je tiefer desto stC$rker smooth+
  theme_light() +
  geom_vline(xintercept = 0,linetype="dotted")+
  geom_vline(xintercept = median_n_offset,linetype="dotted")+
  theme(
    axis.text=element_text(size=15),
    axis.title=element_text(size=15),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.text=element_text(size=15),
    strip.text.x = element_text(size = 15),
    strip.text.y = element_text(size = 15),
    legend.background=element_rect(fill = alpha("white", 0))
  )+
  #for this one participant, this doesn't make sense, so the next lines are commented out (significant time clusters), but note how participant 2 looks to the competitor more after noun onset, individual variance is not to underestimate ... ...
  #annotate("rect",xmin = -2900, xmax = -2200, ymin = 0.049, ymax = 0.05,colour="black",fill = "black")+
  #annotate("rect",xmin = 1400, xmax = 3300, ymin = 0.049, ymax = 0.05,colour="black",fill = "black")+
  scale_color_manual(values=c('red','blue'))

plotti3
