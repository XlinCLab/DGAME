
#ADD THE DIRECTOR SURFACE!!!
# this script processes the fixation files from pupil player and prepares them for the next script. 
#The name and task suggest similar code to the gaze script Ca, but ithe nature of the fixation data makes an entirely different approach neccessary. See below

library(tidyverse)

script_path <- rstudioapi::getActiveDocumentContext()$path
dir <- dirname(script_path)
setwd(dir)  #set working directory to script path


fixation_files = c("fixations_on_surface_11.csv",
                   "fixations_on_surface_12.csv",
                   "fixations_on_surface_13.csv",
                   "fixations_on_surface_14.csv",
                   "fixations_on_surface_21.csv",
                   "fixations_on_surface_22.csv",
                   "fixations_on_surface_23.csv",
                   "fixations_on_surface_24.csv",
                   "fixations_on_surface_31.csv",
                   "fixations_on_surface_32.csv",
                   "fixations_on_surface_33.csv",
                   "fixations_on_surface_34.csv",
                   "fixations_on_surface_41.csv",
                   "fixations_on_surface_42.csv",
                   "fixations_on_surface_43.csv",
                   "fixations_on_surface_44.csv")

subjects = '02'


studyroot = dir
fix_all_out = paste(studyroot,'/fixations_all_4analysis.csv',sep="")
fixations_all = tibble()
mean_word_onsets = tibble()
for(s in subjects){
  gaze_file = paste(studyroot,'/gaze_positions/gaze_positions_4analysis.csv',sep='')
  surface_folder = paste(studyroot,'/surfaces/',sep="")
  fix_subj_out = paste(studyroot,'/fixations/fixations_4analysis.csv',sep='')

  cat(paste("Processing subject",s,"\n"))


    
    cat("loading fixations...\n")
    gaze <- as_tibble(read.csv(gaze_file,stringsAsFactors = F)) 
    gaze = gaze[,c(1:28,45:ncol(gaze))] %>% select(-norm_pos_x,-norm_pos_y,-base_data)
    for(file in fixation_files){
      if(!exists("tmp_fixation_positions")){
        tmp_fixation_positions = read.csv(paste(surface_folder,file,sep=""),stringsAsFactors = F) %>% 
          group_by(fixation_id) %>%
          filter(row_number()==1) %>% 
          ungroup()
        surface = gsub('fixations_on_surface_(.*).csv','\\1',file)
        names(tmp_fixation_positions)[ncol(tmp_fixation_positions)] = paste(surface,sep='')


      }
      else if(exists("tmp_fixation_positions")){
        tmp = read.csv(paste(surface_folder,file,sep=""),stringsAsFactors = F) %>% 
          group_by(fixation_id) %>%
          filter(row_number()==1) %>% 
          ungroup()
        surface = gsub('fixations_on_surface_(.*).csv','\\1',file)
        names(tmp)[ncol(tmp)] = paste(surface,sep='')
        tmp = tmp %>% select(-world_timestamp,-world_index,-start_timestamp,-duration,-dispersion,-norm_pos_x,-norm_pos_y,-x_scaled,-y_scaled)
        
        tmp_fixation_positions = left_join(tmp_fixation_positions,tmp,by = 'fixation_id')
      }
    }
    tmp_fixation_positions = tmp_fixation_positions[,c(2:ncol(tmp_fixation_positions))] %>% 
      dplyr::rename(gaze_timestamp = start_timestamp)
    
    tmp_fixation_positions$saccAmpl = 0
    
    for(ev in 2:nrow(tmp_fixation_positions)){
      a=c(tmp_fixation_positions$norm_pos_x[ev-1],tmp_fixation_positions$norm_pos_y[ev-1])
      b= c(tmp_fixation_positions$norm_pos_x[ev],tmp_fixation_positions$norm_pos_y[ev])
      tmp_fixation_positions$saccAmpl[ev] = acos( sum(a*b) / ( sqrt(sum(a * a)) * sqrt(sum(b * b)) ) )
      
    }
    #join gaze positions and fixation surface positions by timestamp (which are identical of course)
    gaze_fix = left_join(gaze,tmp_fixation_positions,by="gaze_timestamp") %>% 
      mutate(end_time = gaze_timestamp+(duration/1000)) 
      
    # 1) define the six mappings and your surface list
    surface_list <- c('11','12','13','14',
                      '21','22','23','24',
                      '31','32','33','34',
                      '41','42','43','44')
    
    mappings <- list(
      aoi_target      = "surface",
      aoi_otherTarget = "surface_competitor",
      aoi_comp        = "targetA_surface",
      aoi_otherComp   = "targetB_surface",
      aoi_fillerA     = "fillerA_surface",
      aoi_fillerB     = "fillerB_surface"
    )
    
    # 2) initialize all your AOI flags to FALSE
    aoi_cols <- c(names(mappings), "aoi_empty")
    gaze_fix[aoi_cols] <- FALSE
    
    # 3) for each of the six, vectorised lookup:
    for(aoi_name in names(mappings)) {
      lookup_col <- mappings[[aoi_name]]
      # read off the code (e.g. "11", "23", etc.) for each row
      codes <- as.character(gaze_fix[[lookup_col]])
      # only keep valid codes
      valid <- which(!is.na(gaze_fix$fixation_id) & codes %in% surface_list)
      # set the aoi flag to whatever that column says in the same row
      gaze_fix[[aoi_name]][valid] <- mapply(
        function(r, code) gaze_fix[[code]][r],
        r    = valid,
        code = codes[valid]
      )
    }
    
    # 1) Make sure all your aoi_* are logical
    for(col in c(names(mappings), "aoi_empty")) {
      gaze_fix[[col]] <- as.logical(gaze_fix[[col]])
    }
    
    # 2) Recompute aoi_empty per row (only where there’s a fixation)
    has_fix <- which(!is.na(gaze_fix$fixation_id))
    
    gaze_fix$aoi_empty[has_fix] <- sapply(has_fix, function(i) {
      # check the six other aoi_* flags at row i
      flags <- c(
        gaze_fix$aoi_target[i],
        gaze_fix$aoi_otherTarget[i],
        gaze_fix$aoi_comp[i],
        gaze_fix$aoi_otherComp[i],
        gaze_fix$aoi_fillerA[i],
        gaze_fix$aoi_fillerB[i]
      )
      # empty is TRUE iff *none* of those six are TRUE
      !any(flags, na.rm = TRUE)
    })
    gaze_fix[gaze_fix$aoi_target,]$aoi_fillerA = FALSE
    gaze_fix[gaze_fix$aoi_target,]$aoi_fillerB = FALSE
    write.csv(file = fix_subj_out,gaze_fix,row.names=FALSE)
  
  #   gaze_fix_all = rbind(gaze_fix_all,gaze_fix)
  
}

