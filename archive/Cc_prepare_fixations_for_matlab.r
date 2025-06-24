library(tidyverse)

script_path <- rstudioapi::getActiveDocumentContext()$path
dir <- dirname(script_path)
setwd(dir)  #set working directory to script path
subjects = c('02')

blocks = c('11','12','21','22')
for(s in subjects){
  cat(paste(s,"\n"))

  fixations_file = paste(dir,'/fixations/fixations_4analysis.csv',sep='')

  #read fixations file
  fixations = as_tibble(read.csv(fixations_file,sep=',') )

    fixations = fixations %>% filter(!is.na(trial_time) & !is.na(fixation_id))
  names(fixations)[47:62] = c("AA","AB","AC","AD","BA","BB","BC","BD","CA","CB","CC","CD","EA","EB","EC","ED")
  
  
  pattern = 1
  set = 1
  fixations$block = NA
  fixations$fix_at = NA
  fixations$fix_at[!is.na(fixations$fixation_id)]='elsewhere'
    
  fixations$fix_at[fixations$aoi_target == TRUE] = "target"
  fixations$fix_at[fixations$aoi_target == FALSE & (fixations$aoi_comp == TRUE | fixations$aoi_otherTarget == TRUE | fixations$aoi_otherComp == TRUE | fixations$aoi_fillerA == TRUE | fixations$aoi_fillerB == TRUE)] = 'other'
  fixations$fix_at[fixations$aoi_goal == TRUE] = 'goal'

  for(ii in 1:(nrow(fixations)-1)){
    if(pattern > 2){
      pattern = 1
      set = set+1
    }
    else if((fixations$time[ii] < fixations$time[ii+1]) & pattern < 3){
      fixations$block[ii] = paste(set,pattern,sep='')
    }
    else if((fixations$time[ii] > fixations$time[ii+1])  & pattern < 3){
      fixations$block[ii] = paste(set,pattern,sep='')
      pattern = pattern+1
    }
  }
  for(b in blocks){
    cat(paste(b,"..."))
    fixations = fixations[order(fixations$gaze_timestamp),]
    tmp = fixations %>% filter(block == b )
    fix_out =  paste(dir,'/fixations/fixations_times_',b,'_trials.csv',sep='')
    write.csv(file=fix_out,tmp,row.names=FALSE)
    
  }
  
}

