
#ADD THE DIRECTOR SURFACE!!!
# this script creates trials from the continuous gaze file, the input is the CORRECTED word file and the gaze file from pupil player (ET software)
library(tidyverse)

script_path <- rstudioapi::getActiveDocumentContext()$path
dir <- dirname(script_path)
setwd(dir)  #set working directory to script path


surface_files = c("gaze_positions_on_surface_11.csv",
                  "gaze_positions_on_surface_12.csv",
                  "gaze_positions_on_surface_13.csv",
                  "gaze_positions_on_surface_14.csv",
                  "gaze_positions_on_surface_21.csv",
                  "gaze_positions_on_surface_22.csv",
                  "gaze_positions_on_surface_23.csv",
                  "gaze_positions_on_surface_24.csv",
                  "gaze_positions_on_surface_31.csv",
                  "gaze_positions_on_surface_32.csv",
                  "gaze_positions_on_surface_33.csv",
                  "gaze_positions_on_surface_34.csv",
                  "gaze_positions_on_surface_41.csv",
                  "gaze_positions_on_surface_42.csv",
                  "gaze_positions_on_surface_43.csv",
                  "gaze_positions_on_surface_44.csv")


subjects = c('02')
done_s = c()

gaze_all_out = paste(dir,'/gaze_positions_all_4analysis.csv',sep="")
gaze_positions_all = tibble()
mean_word_onsets = tibble()
for(s in subjects){
  times_path =  paste(dir,'/helper_files/',sep="")
  gaze_file = paste(dir,'/gaze_positions/gaze_positions.csv',sep="")
  audiopath = paste(dir,'/audio/',sep="")
  surface_folder = paste(dir,'/surfaces/',sep="")
  gaze_before_words_file = paste(dir,'/gaze_positions/gaze_positions_before_words.csv',sep="")
  word_outfile = paste(dir,'/audio/all_words_4analysis.csv',sep="")
  object_positions_file = paste(dir,'/object_positions/object_positions.csv',sep="")
  gaze_subj_out = paste(dir,'/gaze_positions/gaze_positions_4analysis.csv',sep="")
  tmp_gaze_s = paste(dir,'/gaze_positions/tmp_gaze_positions.csv',sep="")
  gaze_positions_s = tibble()
  
  cat(paste("Processing subject",s,"\n"))
  setwd(paste(dir,'/audio/',sep=""))
  file_list = list.files(pattern=paste(s,".*2erp.*.csv",sep=""),path = audiopath)
  
  times_files_list = list.files(pattern = paste("^",s,"_times_.*csv",sep=""),path = times_path)
  timestamps_file_list = list.files(pattern = paste("^",s,"_timestamps_.*csv",sep=""),path = times_path)
  if(s == '08'){file_list = file_list[2:4]
  times_files_list = times_files_list[2:4]
  timestamps_file_list = timestamps_file_list[2:4]}
  
  
  if(!(s %in% done_s)){  
    cat("loading gaze...\n")
    gaze_raw <- as_tibble(read.csv(gaze_file)) %>% 
      select(-gaze_point_3d_x,-gaze_point_3d_y,-gaze_point_3d_z,
             -gaze_normal1_x,-gaze_normal1_y,-gaze_normal1_z,
             -eye_center1_3d_x,-eye_center1_3d_y,-eye_center1_3d_z,
             -gaze_normal0_x,-gaze_normal0_y,-gaze_normal0_z,
             -eye_center0_3d_x,-eye_center0_3d_y,-eye_center0_3d_z)
    words = tibble()
    cat("loading words and combining with gaze...\n")
    for(n_file in 1:length(file_list)){
      file = file_list[n_file]
      
      
      tmp = read.csv(file,sep=",")
      names(tmp)[2] = 'time'
      
      times <- as_tibble(read.csv(paste(dir,'/helper_files/',times_files_list[n_file],sep="")),header=F)
      timestamps = read.csv(paste(dir,'/helper_files/',timestamps_file_list[n_file],sep=""),header=F)
      names(times) = 'time'
      part_gaze = gaze_raw %>% filter(round(gaze_timestamp,7) >= timestamps[1,1]  & round(gaze_timestamp,7) < timestamps[2,1])
      part_gaze = cbind(part_gaze,times)
      #now we have to work out a way to combine the ET data and the words, because they do not exactly match,
      # we do this like this: take the times of the gaze data, add a column id filled with 0 and a column keep filled with 1
      #then we take the times and ids (=word numbers) of the word data frame and create a keep = 1 column.
      #then we combine the two and loop through all times. Whenever id != 0 (i.e. a word) then either the id of the time point before or
      # after (depending on which is closer "if(one_before > one_after & tt < nrow(times_all))") the time point after the current word 
      # (i.e. a gaze time point) gets set to the id of the current word.
      # This means, we do a "closest" join of the two dataframes (words/gaze) with 7 decimal place accuracy.
      times_gaze = part_gaze %>% select(time)
      times_gaze$id = 0
      times_gaze$keep = 1
      times_words=tmp %>% select(time,id)
      times_words$keep = 0
      times_all = rbind(times_gaze,times_words)
      times_all = times_all[order(times_all$time),]
      #  times_all$time[!is.numeric(times_all$time)]
      for(tt in 1:nrow(times_all)){
        if(times_all$id[tt] != 0){
          one_before = abs(times_all$time[tt]-times_all$time[tt-1])
          one_after= abs(times_all$time[tt]-times_all$time[tt+1])
          if(one_before > one_after & tt < nrow(times_all)){
            times_all$id[tt+1] = times_all$id[tt]
          }else{
            times_all$id[tt-1] = times_all$id[tt]
          }  
          
        }
      }
      # finally, filter and join with words by id
      times_all = times_all %>% filter(keep == 1)
      part_gaze = cbind(part_gaze,times_all$id)
      names(part_gaze)[ncol(part_gaze)] = "id"
      tmp_w_to_join = tmp %>% dplyr::rename(audio_time = time)
      part_gaze = left_join(part_gaze,tmp_w_to_join,by="id") %>% mutate(text = tolower(text))
      gaze_positions_s = rbind(gaze_positions_s,part_gaze)
      words = rbind(words,tmp) %>%  as_tibble() %>% mutate(set = as.character(set),pattern = as.character(pattern)) %>% mutate(text = tolower(text))
      #rbind the part_gaze and gaze_positions_s in the end
    } #end loop over word files; output = words joined with gaze of entire experimental session
    cat("loading surfaces...\n")
    for(file in surface_files){
      if(!exists("tmp_surface_positions")){
        tmp_surface_positions = read.csv(paste(surface_folder,file,sep="")) %>% 
          select(-confidence,-world_index,-x_norm,-y_norm,-x_scaled,-y_scaled,-world_timestamp)
        surface = gsub('gaze_positions_on_surface_(.*).csv','\\1',file)
        names(tmp_surface_positions)[ncol(tmp_surface_positions)] = paste(surface,sep='')
      }
      else if(exists("tmp_surface_positions")){
        tmp = read.csv(paste(surface_folder,file,sep="")) %>% 
          select(-confidence,-world_index,-x_norm,-y_norm,-x_scaled,-y_scaled,-world_timestamp)
        surface = gsub('gaze_positions_on_surface_(.*).csv','\\1',file)
        names(tmp)[ncol(tmp)] = paste(surface,sep='')
        tmp_surface_positions = left_join(tmp_surface_positions,tmp,by="gaze_timestamp")
      }
    }
    #join gaze positions and surface positions by timestamp (which are identical of course)
    gaze_positions_s = left_join(gaze_positions_s,tmp_surface_positions,by=c('gaze_timestamp'))
    rm(tmp_surface_positions)
    gaze_positions_s$subj = s
    write.csv(file=gaze_before_words_file,gaze_positions_s,row.names = FALSE)
    write.csv(file=word_outfile,words,row.names=FALSE)
    # now we create the trials in a loop
    gaze_positions_s$trial = NA
    gaze_positions_s$trial_time = NA
    conlist = c('conflict','no_conflict') # could be extended to include open or fillers or middle object
    # trial counter
    trial = 1
    n = nrow(gaze_positions_s)
    #loop through rows
    #sorry for the super slow loop, didn't have time to optimize the code idea.  I tried once but it didn't work... 
    cat(paste("starting_loop"))
    for(pnt in 1:n){
      # if a point is in conlist = a word begins, and if it is a noun and if the point is not na (which is just so the loop works)
      if(gaze_positions_s$condition[pnt] %in% conlist & (!is.na(gaze_positions_s$pos[pnt])&gaze_positions_s$pos[pnt]=='N')){
        #copy trial information
        gaze_positions_s$trial[pnt] = trial
        gaze_positions_s$trial_time[pnt] = 0
        trial_start_time = gaze_positions_s$time[pnt]-3.5
        trial_end_time = gaze_positions_s$time[pnt]+3.5
        #now we copy this info to every time point within +2 and -2 seconds after noun onset
        count_backwards = 1
        count_forwards = 1
        while(gaze_positions_s$time[pnt-count_backwards] > trial_start_time & (gaze_positions_s$time[pnt-count_backwards]-gaze_positions_s$time[pnt] >= -3.5)){
          gaze_positions_s$trial[pnt-count_backwards] = trial
          gaze_positions_s$trial_time[pnt-count_backwards] = gaze_positions_s$time[pnt-count_backwards]-gaze_positions_s$time[pnt]
          gaze_positions_s$condition[pnt-count_backwards] = gaze_positions_s$condition[pnt]
          gaze_positions_s$surface[pnt-count_backwards] = gaze_positions_s$surface[pnt]
          gaze_positions_s$surface_end[pnt-count_backwards] = gaze_positions_s$surface_end[pnt]
          gaze_positions_s$surface_competitor[pnt-count_backwards] = gaze_positions_s$surface_competitor[pnt]
          gaze_positions_s$targetA_surface[pnt-count_backwards] = gaze_positions_s$targetA_surface[pnt]
          gaze_positions_s$targetB_surface[pnt-count_backwards] = gaze_positions_s$targetB_surface[pnt]
          gaze_positions_s$compA_surface[pnt-count_backwards] = gaze_positions_s$compA_surface[pnt]
          gaze_positions_s$compB_surface[pnt-count_backwards] = gaze_positions_s$compB_surface[pnt]
          gaze_positions_s$fillerA_surface[pnt-count_backwards] = gaze_positions_s$fillerA_surface[pnt]
          gaze_positions_s$fillerB_surface[pnt-count_backwards] = gaze_positions_s$fillerB_surface[pnt]
          gaze_positions_s$target_location[pnt-count_backwards] = gaze_positions_s$target_location[pnt]
          
          count_backwards = count_backwards+1
        }
        while(gaze_positions_s$time[pnt+count_forwards] < trial_end_time & (gaze_positions_s$time[pnt+count_forwards]-gaze_positions_s$time[pnt] <= 3.5)){
          gaze_positions_s$trial[pnt+count_forwards] = trial
          gaze_positions_s$trial_time[pnt+count_forwards] = gaze_positions_s$time[pnt+count_forwards]-gaze_positions_s$time[pnt]
          gaze_positions_s$condition[pnt+count_forwards] = gaze_positions_s$condition[pnt]
          gaze_positions_s$surface[pnt+count_forwards] = gaze_positions_s$surface[pnt]
          gaze_positions_s$surface_end[pnt+count_forwards] = gaze_positions_s$surface_end[pnt]
          gaze_positions_s$surface_competitor[pnt+count_forwards] = gaze_positions_s$surface_competitor[pnt]
          gaze_positions_s$targetA_surface[pnt+count_forwards] = gaze_positions_s$targetA_surface[pnt]
          gaze_positions_s$targetB_surface[pnt+count_forwards] = gaze_positions_s$targetB_surface[pnt]
          gaze_positions_s$compA_surface[pnt+count_forwards] = gaze_positions_s$compA_surface[pnt]
          gaze_positions_s$compB_surface[pnt+count_forwards] = gaze_positions_s$compB_surface[pnt]
          gaze_positions_s$fillerA_surface[pnt+count_forwards] = gaze_positions_s$fillerA_surface[pnt]
          gaze_positions_s$fillerB_surface[pnt+count_forwards] = gaze_positions_s$fillerB_surface[pnt]
          gaze_positions_s$target_location[pnt+count_forwards] = gaze_positions_s$target_location[pnt]
          
          count_forwards = count_forwards+1
        }
        percent = round(pnt/n*100,1)
        cat(paste(percent,'% done\n'))
        trial = trial+1
      }
    }
    write.csv(file=tmp_gaze_s,gaze_positions_s,row.names=FALSE)
    
    gaze_positions_s$aoi_target = FALSE
    gaze_positions_s$aoi_otherTarget = FALSE
    gaze_positions_s$aoi_comp = FALSE
    gaze_positions_s$aoi_otherComp = FALSE
    gaze_positions_s$aoi_fillerA = FALSE
    gaze_positions_s$aoi_fillerB = FALSE
    gaze_positions_s$aoi_empty = FALSE
    gaze_positions_s$aoi_other = FALSE
    gaze_positions_s$trackloss = FALSE
    gaze_positions_s$aoi_goal = FALSE
    conlist = c('conflict','no_conflict')
    # the next loop exludes trackloss trials and checks, if at a given time point participants looked at a surface or not  
    tmp = gaze_positions_s %>% filter(!is.na(trial_time) & condition %in% conlist & !is.na(surface) & confidence >= .6)
    tmp2 = gaze_positions_s %>% filter(is.na(trial_time) & !(condition %in% conlist) & is.na(surface) & confidence < .6)
    count = 1
    n = nrow(tmp)
    surface_list=c('11','12','13','14','21','22','23','24','31','32','33','34','41','42','43','44')
    targets_list = unique(tmp$text)
    for(pnt in 1:n){
      if(tmp$condition[pnt] %in% conlist){
        if(!is.na(tmp$surface[pnt])){
          target = paste(tmp$surface[pnt])
          goal = paste(tmp$target_location[pnt])
          if(target == "4"){
            target = "44"
          }
          if(goal == "2"){
            goal = "22"
          }
          if(target == "2"){
            target = "22"
          }
          otherTarget = paste(tmp$targetB_surface[pnt])
          if(otherTarget == "4"){
            otherTarget = "44"
          }
          competitor = paste(tmp$surface_competitor[pnt])
          otherCompetitor = paste(tmp$compB_surface[pnt])
          fillerA = paste(tmp$fillerA_surface[pnt])
          if(fillerA == "4"){
            fillerA = "44"
          }
          fillerB = paste(tmp$fillerB_surface[pnt])
          if(fillerB == "4"){
            fillerB = "44"
          }
          empty = setdiff(surface_list,c(target,competitor,otherTarget,otherCompetitor,fillerA,fillerB,goal))
          #if(!is.na(tmp[pnt,]$director)){  
          if(target != "NA" & target != "Fehler"){
            if(!is.na(tmp[[target]][pnt])){
              if(tmp[[target]][pnt] == 'True'){
                tmp$aoi_target[pnt] = TRUE
              }
            }
          }
          if(goal != "NA" & goal != "Fehler"){
            if(!is.na(tmp[[goal]][pnt])){
              if(tmp[[goal]][pnt] == 'True'){
                tmp$aoi_goal[pnt] = TRUE
              }
            }
          }
          if(otherTarget != "NA" & otherTarget != "Fehler"){
            if(!is.na(tmp[[otherTarget]][pnt])){
              if(tmp[[otherTarget]][pnt] == 'True'){
                tmp$aoi_otherTarget[pnt] = TRUE
              }
            }
          }
          if(competitor!="NA"){
            if(!is.na(tmp[[competitor]][pnt])){
              if(tmp[[competitor]][pnt] == 'True'){
                tmp$aoi_comp[pnt] = TRUE
              }
            }
          }
          if(otherCompetitor!="NA"){
            if(!is.na(tmp[[otherCompetitor]][pnt])){
              if(tmp[[otherCompetitor]][pnt] == 'True'){
                tmp$aoi_otherComp[pnt] = TRUE
              }
            }
          }
          if(fillerA!="NA"){
            if(!is.na(tmp[[fillerA]][pnt])){
              if(tmp[[fillerA]][pnt] == 'True'){
                tmp$aoi_fillerA[pnt] = TRUE
              }
            }
          }
          if(fillerB!="NA"){
            if(!is.na(tmp[[fillerB]][pnt])){
              if(tmp[[fillerB]][pnt] == 'True'){
                tmp$aoi_fillerB[pnt] = TRUE
              }
            }
          }
          if((!is.na(tmp[[empty[1]]][pnt]) & tmp[[empty[1]]][pnt] == 'True' ) | (!is.na(tmp[[empty[2]]][pnt]) & tmp[[empty[2]]][pnt] == 'True' ) | 
             (!is.na(tmp[[empty[3]]][pnt]) & tmp[[empty[3]]][pnt] == 'True' ) | (!is.na(tmp[[empty[4]]][pnt]) & tmp[[empty[4]]][pnt] == 'True' ) |
             (!is.na(tmp[[empty[5]]][pnt]) & tmp[[empty[5]]][pnt] == 'True' ) |(!is.na(tmp[[empty[6]]][pnt]) & tmp[[empty[6]]][pnt] == 'True' ) |
             (!is.na(tmp[[empty[7]]][pnt]) & tmp[[empty[7]]][pnt] == 'True' ) |(!is.na(tmp[[empty[8]]][pnt]) & tmp[[empty[8]]][pnt] == 'True' )){
          tmp$aoi_empty[pnt] = TRUE
          }
          
        }
      }
      percent = round(pnt/n*100,1)
      cat(paste(percent,'% done\n'))
      count = count+1
    }
    tmp2$trackloss = TRUE
    
    
    
    gaze_positions_s = rbind(tmp,tmp2) # add trackloss data again
    gaze_positions_s = gaze_positions_s[order(gaze_positions_s$gaze_timestamp),] #sort by time
    write.csv(file = gaze_subj_out,gaze_positions_s,row.names=FALSE) # save
  }
     gaze_positions_all = rbind(gaze_positions_all,gaze_positions_s) #optional for multiple subjects, or run the code below
  
}

gaze_all_out = paste(dir,'/gaze_positions_all_4analysis.csv',sep="")

write.csv(file = gaze_all_out,gaze_positions_all,row.names=FALSE)

# now the stats script can be run
