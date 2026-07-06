#script to preprocess word lists. Not optimized in any way, a terrible chain of loops, coded as it came to mind, but functional.
library(tidyverse)

# function to retrieve frequency class from leipzig corpus
retrieve_word_data <- function(wordlist, corpus="deu_news_2012_3M", as_dataframe = TRUE){
  if(corpus == "deu_news_2012_3M"){
    url_le <- "http://api.wortschatz-leipzig.de/ws/words/deu_news_2012_3M/word/"
  }
  word_data_list <- vector(mode="list", length = length(wordlist))
  for (ii in 1:length(wordlist)){
    cat(paste(round(ii/length(wordlist)*100,2),'%\n',sep=" "))
    http_le <- paste(url_le, wordlist[[ii]], sep="")
    retrieved_data <- httr::content(httr::GET(http_le, httr::add_headers("Accept: application/json")))
    if (length(retrieved_data)==2){
      word_data_list[[ii]] = list(id = NA,word = wordlist[[ii]],freq = NA, wordRank = NA,frequencyClass = NA)
    }else{word_data_list[[ii]] <- retrieved_data
    }
  }
  if(as_dataframe){
    word_data_df <- do.call(rbind, lapply(word_data_list,unlist))
    word_data_df = as.data.frame(word_data_df)
  }else{
    word_data_list
  }
  
}

script_path <- rstudioapi::getActiveDocumentContext()$path
dir <- dirname(script_path)
setwd(dir)  #set working directory to script path
#studyroot = dir

#subjects and words to count etc.
subjects = c('02')
for(s in subjects){
  audiopath = paste(dir,'audio/',sep="/")
  object_positions_file = paste(dir,'object_positions/object_positions.csv',sep='/')
  setwd(audiopath)
  file_list = list.files(pattern=".*[0-9][0-9]_words_[0-9][0-9].csv",path = audiopath)
  
  pattern = 1
  set = 1
  for(file in file_list){
    cat(paste(file))
    #which object names are possible?
    objects = c('kerze','vase','blume','tasse','rose','Kerze','Vase','Blume','Tasse','Rose')
    filler = c('tube','dose','flasche','pflanze','creme','spritze','paste','Tube','Dose','Flasche','Pflanze','Creme','Spritze','Paste')
    
    
    words = read.csv(file,header=T) %>% distinct() %>% mutate(text = as.character(text)) 
    subject = substr(file,1,2)
    block = substr(file,nchar(file)-5,nchar(file)-4)
    file_out = paste(subject,'_words_',block,'_2analysis.csv',sep="")
    names(words) = c("id","time","text","tmax")
    if(ncol(words)>4){#this is just, because in the original run, we reran the script for participants with "errors" in the word list, e.g. 'äh' before the noun etc. These files all already had a frequency column, and so the script would crash (yes, its a mess)
      words = words %>% select(id,time,text,tmax)
    }
    for(w in 1:nrow(words)){
      if(words$text[w] %in% objects | words$text[w] %in% filler){
        words$text[w] = stringr::str_to_title(words$text[w])
      }
    }
    #retrieve frequency
    if(file != "25_words_22.csv"){
      freqClass = retrieve_word_data(unique(words$text))[,c(2,5)] %>% dplyr::rename(text=word)
      freqClass$text = as.character(freqClass$text)
      words$text = as.character(words$text)
      words = left_join(words,freqClass) %>% rename(freqClass = frequencyClass) %>% mutate(freqClass=as.numeric(as.character(freqClass))) %>% distinct()
    }else{
      words$freqClass =  retrieve_word_data(words$text)[,c(5)]
      words = words %>% mutate(freqClass = as.numeric(freqClass))
      
    }
    words$freqClass[is.na(words$freqClass)]=max(words$freqClass,na.rm = T)+1 #all words without entry receive max frequency class +1
    words = words %>% 
      mutate(text=tolower(text)) %>%
      as_tibble()
    #prepare words, add condition etc.
    objects = c('kerze','vase','blume','tasse','rose')
    filler = c('tube','dose','flasche','pflanze','creme','spritze','paste')
    count = tibble(objects,c(0,0,0,0,0),.name_repair = ~ c("object", "count"))
    count_filler = tibble(filler,c(0,0,0,0,0,0,0),.name_repair = ~ c("object", "count"))
    words$pattern = pattern
    words$set = set
    
    words = words[words$text!="",]
    words$pattern=as.character(words$pattern)
    words$set=as.character(words$set)
    
    words$position = NA
    words$condition = NA
    words$condition_code = NaN
    words$pos = 'word'
    start = 1
    if(file == '11_words_11.csv'){
      start = 300
    }
    for(w in start:nrow(words)){
      if(!((w==268 | w == 264 | w == 279) & file == '02_words_11.csv') & !(file == '33_words_11.csv' & (w == 30 | w == 325 | w == 323))  & !(file == '32_words_21.csv' & (w <= 156)) & !(file == '32_words_12.csv' & (w == 71 | w == 70 )) & !(file == '30_words_21.csv' & (w == 21))  & !(file == '21_words_11.csv' & (w == 42 | w == 43 ))  & !(file == '21_words_22.csv' & (w == 179 | w == 180 )) & !(file == '20_words_21.csv' & (w == 165 | w == 164 )) & !(file == '11_words_22.csv' & (w == 339 | w == 276)) & !(file == '08_words_11.csv' & (w == 127 | w == 129 )) & !(file == '08_words_21.csv' & (w == 175 | w == 183 | w == 188 )) & !(file == '07_words_11.csv' & (w == 258 )) & !(file == '30_words_11.csv' & (w == 367 | w == 389)) & !(w==60 & file == '11_words_21.csv') & !(file == '03_words_11.csv' & (w == 18 | w == 29)) & !(file == '03_words_21.csv' & w == 76) & !(w==148 & file == '31_words_21.csv') & !((w==19 | w == 36) & file == '31_words_11.csv') & !((w==42 | w == 47 | w == 70 | w == 51) & file == '29_words_21.csv') & !(w == 472 & file == '29_words_22.csv') & !(w==19 & file == '29_words_12.csv') & !(w==80 & file == '29_words_11.csv') & !(w==245 & file == '28_words_21.csv') & !(w==156 & file == '27_words_22.csv') & !((w == 9 | w == 14 | w == 12 | w == 92) & file == '27_words_21.csv') & !((w == 424 | w == 426) & file == '25_words_12.csv')  & !(w == 175 & file == '27_words_12.csv')  & !((w == 331 | w == 183 | w == 850) & file == '25_words_21.csv') &!((w == 739) & file == '25_words_22.csv') & !((w == 42 | w == 47 | w == 53) & file == '24_words_12.csv') & !((w == 61 | w == 250) & file == '23_words_11.csv')  & !(( w == 325 | w == 423 | w == 559) & file == '22_words_22.csv') & !(( w == 106 | w == 119 | w == 805 | w == 815) & file == '22_words_11.csv') & !(( w >= 176 & w <= 180) & file == '20_words_22.csv') & !(( w == 273 | w == 317 | w == 318 | w == 460 | w == 461 | w == 688 | w == 689 | w == 482 | w == 483) & file == '07_words_22.csv') & !(( w == 342 | w == 343 | w == 256 | w == 257 | w == 546 | w == 547) & file == '07_words_11.csv') & !(( w == 34 | w == 35 | w == 399 | w == 400) & file == '07_words_12.csv') & !(w==561 & file == '11_words_11.csv')& !(w<=105 & file == '26_words_11.csv') & !(w==36 & file == '18_words_11.csv') & !(w==19 & file == '17_words_21.csv') & !((w==240 | w == 241 | w == 257) & file == '17_words_12.csv') &  !(((w>=47 & w <= 50) | (w>=237 & w <=243)) & file == '10_words_21.csv')  &  !((w>=193 & w <= 199 & file == '12_words_11.csv')) & !(((w == 196 & file == "15_words_22.csv") | ((w>=233 & w <=234) & file == '15_words_12.csv')) | ((w>=346 & w <=347) & file == '15_words_11.csv'))){
        if(words$text[w] %in% objects ){
          if(words$text[w-1]== 'die' | words$text[w-1] == 'der'){
            nback = 1
          }else{
            nback=2
          }
          words$condition[(w-nback):w] = "conflict"
          words$condition_code[w-nback] = 11
          words$condition_code[w] = 12
          words$pos[w+1] = 'next'
          words$pos[w+2] = 'next'
          words$pos[w] = 'N'
          words$pos[w-1] = 'D'
          words$pos[w-2] = 'prev'
          
          words$position[w-nback] = count$count[count$object == words$text[w]]+1
          words$position[w] = count$count[count$object == words$text[w]]+1
          count$count[count$object == words$text[w]] = count$count[count$object == words$text[w]]+1
          
          
        }
        else if(words$text[w] %in% filler & !(file == '05_words_11.csv' & w == 9 )){
          if(words$text[w-1]== 'die' | words$text[w-1] == 'der'){
            nback = 1
          }else{
            nback=2
          }
          words$condition[(w-nback):w] = "no_conflict"
          words$condition_code[w-nback] = 21
          words$condition_code[w] = 22
          words$pos[w+1] = 'next'
          words$pos[w+2] = 'next'
          words$pos[w] = 'N'
          words$pos[w-nback] = 'D'
          words$pos[w-(nback+1)] = 'prev'
          words$pos[w-(nback+2)] = 'prev'
          
          if(nback == 2){
            words$pos[w-1] = 'prev'
          }
          words$position[w-nback] = count_filler$count[count_filler$object == words$text[w]]+1
          words$position[w] = count_filler$count[count_filler$object == words$text[w]]+1
          count_filler$count[count_filler$object == words$text[w]] = count_filler$count[count_filler$object == words$text[w]]+1
        }
      }
    }
    object_positions = as_tibble(read.csv(object_positions_file)) %>% mutate(set = as.character(set),pattern = as.character(pattern),object = as.character(object)) %>% rename(text = object) %>% select(-condition)
    #join words with object position file (created by hand)
    words = left_join(words,object_positions)
    for(w in 1:nrow(words)){
      if(!is.na(words$surface[w])){
        if(words$text[w-1]== 'die' | words$text[w-1] == 'der'){
          nback = 1
        }else{
          nback=2
        }
        words$surface[w-nback] = words$surface[w]
        words$surface_competitor[w-nback] = words$surface_competitor[w]
        words$surface_end[w-nback] = words$surface_end[w]
      }
    }
    #now add other object information to the file
    targets_lc = intersect(unique(object_positions$text[!is.na(object_positions$surface_competitor)]),unique(words$text[words$condition == "conflict" & words$pos=="N"]))
    fillers_lc = intersect(unique(object_positions$text[is.na(object_positions$surface_competitor)]),unique(words$text[words$condition == "no_conflict" & words$pos=="N"]))
    
    where_is_targets = tibble(targets_lc,c(NA,NA),.name_repair = ~c("object","surface"))
    where_is_targets$surface[where_is_targets$object == targets_lc[1]] = object_positions$surface[object_positions$text==targets_lc[1] & object_positions$position == 1 & object_positions$pattern == unique(words$pattern)& object_positions$set == unique(words$set)]
    where_is_targets$surface[where_is_targets$object == targets_lc[2]] = object_positions$surface[object_positions$text==targets_lc[2] & object_positions$position == 1 & object_positions$pattern == unique(words$pattern)& object_positions$set == unique(words$set)]
    
    where_is_comps = tibble(targets_lc,c(NA,NA),.name_repair = ~c("object","surface"))
    where_is_comps$surface[where_is_comps$object == targets_lc[1]] = object_positions$surface_competitor[object_positions$text==targets_lc[1] & object_positions$position == 1 & object_positions$pattern == unique(words$pattern)& object_positions$set == unique(words$set)]
    where_is_comps$surface[where_is_comps$object == targets_lc[2]] = object_positions$surface_competitor[object_positions$text==targets_lc[2] & object_positions$position == 1 & object_positions$pattern == unique(words$pattern)& object_positions$set == unique(words$set)]
    
    where_is_fillers = tibble(fillers_lc,c(NA,NA),.name_repair = ~c("object","surface"))
    where_is_fillers$surface[where_is_fillers$object == fillers_lc[1]] = object_positions$surface[object_positions$text==fillers_lc[1] & object_positions$position == 1 & object_positions$pattern == unique(words$pattern)& object_positions$set == unique(words$set)]
    where_is_fillers$surface[where_is_fillers$object == fillers_lc[2]] = object_positions$surface[object_positions$text==fillers_lc[2] & object_positions$position == 1 & object_positions$pattern == unique(words$pattern)& object_positions$set == unique(words$set)]
    
    #taarget locations
    where_will_fillers = tibble(fillers_lc,c(NA,NA),c(3,3),.name_repair = ~c("object","surface","position"))
    where_will_fillers$surface[where_will_fillers$object == fillers_lc[1]] = object_positions$surface[object_positions$text==fillers_lc[1] & object_positions$position == 2 & object_positions$pattern == unique(words$pattern)& object_positions$set == unique(words$set)]
    where_will_fillers$surface[where_will_fillers$object == fillers_lc[2]] = object_positions$surface[object_positions$text==fillers_lc[2] & object_positions$position == 2 & object_positions$pattern == unique(words$pattern) & object_positions$set == unique(words$set)]
    
    where_will_targets = tibble(targets_lc,c(NA,NA),c(3,3),.name_repair = ~c("object","surface","position"))
    where_will_targets$surface[where_will_targets$object == targets_lc[1]] = object_positions$surface[object_positions$text==targets_lc[1] & object_positions$position == 2 & object_positions$pattern == unique(words$pattern)& object_positions$set == unique(words$set)]
    where_will_targets$surface[where_will_targets$object == targets_lc[2]] = object_positions$surface[object_positions$text==targets_lc[2] & object_positions$position == 2 & object_positions$pattern == unique(words$pattern)& object_positions$set == unique(words$set)]
    
    words$targetA_surface = NA
    words$targetB_surface = NA
    words$fillerA_surface = NA
    words$fillerB_surface = NA
    words$compA_surface = NA
    words$compB_surface = NA
    words$target_location = NA
    
    for(w in 1:nrow(words)){
      if(!is.na(words$position[w])){
        if(words$text[w] %in% targets_lc){
          other_target = setdiff(targets_lc,words$text[w])
          other_comp = setdiff(targets_lc,words$text[w])
          target = intersect(targets_lc,words$text[w])
          comp = intersect(targets_lc,words$text[w])
          if(words$text[w-1] == 'die' | words$text[w-1] == 'der'){
            nback = 1
          }else{
            nback=2
          }
          words$targetA_surface[w-nback] = where_is_targets[where_is_targets$object == target,]$surface
          words$targetB_surface[w-nback] = where_is_targets[where_is_targets$object == other_target,]$surface        
          words$compA_surface[w-nback] = where_is_comps[where_is_comps$object == comp,]$surface 
          words$compB_surface[w-nback] = where_is_comps[where_is_comps$object == other_comp,]$surface        
          words$fillerA_surface[w-nback] = where_is_fillers[1,]$surface
          words$fillerB_surface[w-nback] = where_is_fillers[2,]$surface
          words$targetA_surface[w] = where_is_targets[where_is_targets$object == target,]$surface
          words$targetB_surface[w] = where_is_targets[where_is_targets$object == other_target,]$surface        
          words$compA_surface[w] = where_is_comps[where_is_comps$object == comp,]$surface 
          words$compB_surface[w] = where_is_comps[where_is_comps$object == other_comp,]$surface        
          words$fillerA_surface[w] = where_is_fillers[1,]$surface
          words$fillerB_surface[w] = where_is_fillers[2,]$surface

          where_is_targets$surface[where_is_targets$object == words$text[w]] = words$surface[w]
          where_is_comps$surface[where_is_comps$object == words$text[w]] = words$surface_competitor[w]
        }
        if(words$text[w] %in% fillers_lc){
          other_filler = setdiff(fillers_lc,words$text[w])
          curr_filler = intersect(fillers_lc,words$text[w])
          if(words$text[w-1]== 'die' | words$text[w-1] == 'der'){
            nback = 1
          }else{
            nback=2
          }
          words$targetA_surface[w-nback] = where_is_targets$surface[1]
          words$targetB_surface[w-nback] = where_is_targets$surface[2]        
          words$compA_surface[w-nback] = where_is_comps$surface[1]  
          words$compB_surface[w-nback] = where_is_comps$surface[2]        
          words$fillerA_surface[w-nback] = where_is_fillers[where_is_fillers$object==curr_filler,]$surface
          words$fillerB_surface[w-nback] = where_is_fillers[where_is_fillers$object==other_filler,]$surface
          
          words$targetA_surface[w] = where_is_targets$surface[1]
          words$targetB_surface[w] = where_is_targets$surface[2]        
          words$compA_surface[w] = where_is_comps$surface[1]  
          words$compB_surface[w] = where_is_comps$surface[2]        
          words$fillerA_surface[w] = where_is_fillers[where_is_fillers$object==curr_filler,]$surface
          words$fillerB_surface[w] = where_is_fillers[where_is_fillers$object==other_filler,]$surface
        }
      }
    }
    #this is to enter the goal location
    tar1 = words %>% filter(text == targets_lc[1] & pos == 'N')
    tar2 = words %>% filter(text == targets_lc[2] & pos == 'N')
    fill1 = words %>% filter(text == fillers_lc[1] & pos == 'N')
    fill2 = words %>% filter(text == fillers_lc[2] & pos == 'N')
    rest = words %>% filter(pos != "N")
    
    #enter goal locations
    tar1$target_location= lead(tar1$surface,1)
    tar1$target_location[nrow(tar1)] = tar1$surface_end[1]
    tar2$target_location= lead(tar2$surface,1)
    tar2$target_location[nrow(tar2)] = tar2$surface_end[1]
    fill1$target_location= lead(fill1$surface,1)
    fill1$target_location[nrow(fill1)] = fill1$surface_end[1]
    fill2$target_location= lead(fill2$surface,1)
    fill2$target_location[nrow(fill2)] = fill2$surface_end[1]
    words =rbind(tar1,tar2,fill1,fill2,rest)
    words = words[order(words$id),]
    
    
    
    for(w in 1:nrow(words)){
      if(!is.na(words$target_location[w])){
        if(words$text[w-1]== 'die' | words$text[w-1] == 'der'){
          nback = 1
        }else{
          nback=2
        }
        words$target_location[w-nback] = words$target_location[w]
        
      }
    }
    words$pos[is.na(words$pos)] = 'word'
  #  write_csv(words,file=file_out) #save to file
    if(pattern == 2){
      set = set+1
      pattern = 1
    }else{
      pattern = pattern+1
    }
    
  }
}


####    MANUALLY CORRECT THE OUTPUT BEFORE RUNNING ANY OTHER SCRIPT!
