QC <- function(data, plate, output_dir, na_threshold, rtc_threshold ,normalization_en_ex, report_file){
  
  # Leave a blank line
  cat('', file = report_file, sep = '\n', append = TRUE)
  
  # Samples
  samples_cols <- which(!names(data) %in% c("Well", "ID", "plate"))
  sample_names <- names(data)[which(!names(data) %in% c("Well", "ID", "plate"))]
  
  # Start
  plate_string <- paste("Plate", plate, sep = ' ')
  to_report <- paste('Analysis for plate ', plate, sep ='')
  cat(to_report, file = report_file, sep = '\n', append = TRUE )
  data <- data[c(which(data$plate == plate_string)),]
  
  ## QC
  data.plot <- data[,samples_cols]
  
  # Missing cells plot - plot 1
  myplot <- vis_miss(data.table(data.plot))
  save_as_pdf({print(myplot)}, file.name = paste(output_dir,'/plate',plate,'/vis_missing_cells_plate_',plate,'.pdf', sep = ''), width = 6, height = 4)
  save_image({print(myplot)}, file.name = paste(output_dir,'/plate',plate,'/vis_missing_cells_plate_',plate,'.png', sep = ''))
  
  
  # missing upset - plot 2
  result = tryCatch({
    
    myplot <- gg_miss_upset(data.table(data.plot))
    save_as_pdf({print(myplot)}, file.name = paste(output_dir,'/plate',plate,'/gg_miss_upset_plate_',plate,'.pdf', sep = ''), width = 6, height = 4)
    save_image({print(myplot)}, file.name = paste(output_dir,'/plate',plate,'/gg_miss_upset_plate_',plate,'.png', sep = ''))
    dev.off()
  }, error = function(e) {
    
    dev.off()
    to_report <- paste('Total data gg_miss_upset error for plate ', plate, sep ='')
    cat(to_report, file = report_file, sep = '\n', append = TRUE )
    cat(as.character(e[1]), file = report_file, sep = '\n', append = TRUE )
    
  })
  
  
  # Excluding NA's
  filtered_data.plot <- data.plot
  
  NA_s_per_sample_perc <- colSums(is.na(data.plot))/dim(data)[1]
  samples_to_exclude <- as.numeric(which(NA_s_per_sample_perc > na_threshold))
  
  if (length(samples_to_exclude)>0){
    filtered_data.plot <- filtered_data.plot[, -c(samples_to_exclude)]
    excluded_samples <- names(data.plot)[samples_to_exclude]
    sample_names<- sample_names[which(!sample_names %in% excluded_samples)]
    to_report <- paste('Samples [', paste(unlist(excluded_samples), collapse = ', '), '] were excluded because of their high NAs percentage.', sep = '')
    cat(to_report, file = report_file, sep = '\n', append = TRUE )
  }
  
  if(length(sample_names) == 0){
    print(paste("All samples from plate ", plate, " weere rejected. Check report", sep = ''))
    return(NULL)
  }else if (length(which(startsWith(names(filtered_data.plot ), "Cancer"))) == 0) {
    print(paste("All cancer data from plate ",plate, " were rejected, because of high NA's percentage. Check report.", sep = ''))
  }else if (length(which(startsWith(names(filtered_data.plot ), "Normal"))) == 0) {
    print(paste("All normal data from plate ",plate, " were rejected, because of high NA's percentage. Check report.", sep = ''))
  }
  
  # boxplot quality before normalization
  # TO DO: exclude all the figures on the working directory (ask Nikos for his tool)
  Mnemiopsis_count_data<- filtered_data.plot
  pseudoCount = log2(Mnemiopsis_count_data + 1) # log-transform to make numbers on scale (+1 to avoid zeroes)
  df = melt(pseudoCount, variable.name = "Samples", value.name = "count") # reshape the matrix
  
  #png(file=paste(output_dir,'/plate',plate,'/boxplot_quality_plate_',plate,'.png', sep = ''), width=900, height=600)
  myplot <- ggplot(df, aes(x = Samples, y = count)) + geom_boxplot() + xlab("") +
    ylab(expression(log[2](count + 1)))+
    coord_flip()
  save_as_pdf({print(myplot)}, file.name = paste(output_dir,'/plate',plate,'/boxplot_quality_plate_',plate,'.pdf', sep = ''), width = 6, height = 4)
  save_image({print(myplot)}, file.name = paste(output_dir,'/plate',plate,'/boxplot_quality_plate_',plate,'.png', sep = ''))
  #print(myplot)
  #dev.off()
  
  # TODO
  # Exclude na's
  # Done
  
  # Filtered full data
  if (length(samples_to_exclude)>0){
    filtered_data <- data[,!(names(data) %in% excluded_samples)]
  }else{
    filtered_data <- data
  }
  
  ##normalization/ case
  RTC<- filter(filtered_data, ID == "PPC" | ID == "miRTC")
  how_many_normal <- length(which(startsWith(names(RTC), "Normal")))
  how_many_cancer <- length(which(startsWith(names(RTC), "Cancer")))
  
  RTC.n <- NULL
  RTC.c <- NULL
  
  if (how_many_normal > 0){
    RTC.n<- t(RTC[,c(which(startsWith(names(RTC), "Normal")), length(names(RTC))-1)]%>%  
                group_by(ID) %>%
                summarise_all("mean"))
    RTC.n<- as.data.frame(RTC.n)
    names(RTC.n)<- RTC.n[1,]
    RTC.n<- RTC.n[-1,]
    RTC.n$criterion<- as.numeric(RTC.n$miRTC) - as.numeric(RTC.n$PPC)
    
  }
  
  if (how_many_cancer > 0){
    RTC.c<- t(RTC[,c(which(startsWith(names(RTC), "Cancer")), length(names(RTC))-1)]%>%  
                group_by(ID) %>%
                summarise_all("mean"))
    RTC.c<- as.data.frame(RTC.c)
    names(RTC.c)<- RTC.c[1,]
    RTC.c<- RTC.c[-1,]
    RTC.c$criterion<- as.numeric(RTC.c$miRTC) - as.numeric(RTC.c$PPC)
    
  }
  
  if (!is.null(RTC.c) && !is.null(RTC.n)){
    RTC.all<- rbind(RTC.c, RTC.n)
  }else if (is.null(RTC.c) && !is.null(RTC.n)){
    RTC.all <- RTC.n
  }else if (!is.null(RTC.c) && is.null(RTC.n)){
    RTC.all <- RTC.c
  }else{
    print(paste("All samples from plate ", plate, " weere rejected. Check report", sep = ''))
    return(NULL)
  }
  
  # Pass fail test
  RTC.all$RTC <- paste("Pass")
  failed <- which(RTC.all$criterion >= rtc_threshold)
  RTC.all$RTC[failed] <- paste("Fail")
  
  
  # Drop failed
  if (length(failed) > 0){
    failed_names <- rownames(RTC.all)[failed]
    RTC.all <- RTC.all[-c(failed),]
    sample_names<- sample_names[which(!sample_names %in% failed_names)]
    to_report <- paste('Samples [', paste(unlist(failed_names), collapse = ', '), '] were excluded because they failed RTC criterion.', sep = '')
    cat(to_report, file = report_file, sep = '\n', append = TRUE )
    filtered_data <- filtered_data[,!(names(filtered_data) %in% failed_names)]
    filtered_data.plot <- filtered_data.plot[,!(names(filtered_data.plot) %in% failed_names)]
  }
  
  if(length(sample_names) == 0){
    print(paste("All samples from plate ", plate, " weere rejected. Check report", sep = ''))
    return(NULL)
  }else if (length(sample_names) == 1){
    print(paste("Your dataset has been left only with one sample, the ", sample_names[1] , " sample. All other samples have failed. Check report.", sep = ''))
  }
  
  filtered_data.plot <- as.data.frame(filtered_data.plot)
  names(filtered_data.plot) <- sample_names
  
  filtered_data.plot$ID <- filtered_data$ID
  drops_id <- c("ID") # To drop id column, later
  
  endogenous<- filter(filtered_data.plot, ID == "SNORD61" | ID == "SNORD68" | 
                        ID == "SNORD72" | ID == "SNORD95" | ID=="RNU6-2")
  
  exogenous<- filter(filtered_data.plot, ID == "cel-miR-39")
  
  endogenous.plot <- endogenous[,!(names(endogenous) %in% drops_id)]
  exogenous.plot <- exogenous[,!(names(exogenous) %in% drops_id)]
  
  ## gg miss upset plot - endogenous data
  result = tryCatch({
    #png(file=paste(output_dir,'/plate',plate,'/gg_miss_upset_endogenous_plate_',plate,'.png', sep = ''), width=900, height=600)
    myplot <- gg_miss_upset(data.table(endogenous.plot))
    save_as_pdf({print(myplot)}, file.name = paste(output_dir,'/plate',plate,'/gg_miss_upset_endogenous_plate_',plate,'.pdf', sep = ''), width = 6, height = 4)
    save_image({print(myplot)}, file.name = paste(output_dir,'/plate',plate,'/gg_miss_upset_endogenous_plate_',plate,'.png', sep = ''))
    
    #print(myplot)
    #dev.off()
  }, error = function(e) {
    #dev.off()
    to_report <- paste('Endogenous data gg_miss_upset error for plate ', plate, sep ='')
    cat(to_report, file = report_file, sep = '\n', append = TRUE )
    cat(as.character(e[1]), file = report_file, sep = '\n', append = TRUE )
  })
  
  
  ## gg miss upset plot - exogenous data
  result = tryCatch({
    #png(file=paste(output_dir,'/plate',plate,'/gg_miss_upset_exogenous_plate_',plate,'.png', sep = ''), width=900, height=600)
    myplot <- gg_miss_upset(data.table(exogenous.plot))
    save_as_pdf({print(myplot)}, file.name = paste(output_dir,'/plate',plate,'/gg_miss_upset_exogenous_plate_',plate,'.pdf', sep = ''), width = 6, height = 4)
    save_image({print(myplot)}, file.name = paste(output_dir,'/plate',plate,'/gg_miss_upset_exogenous_plate_',plate,'.png', sep = ''))
    dev.off()
    #print(myplot)
    #dev.off()
  }, error = function(e) {
    #dev.off()
    to_report <- paste('Exogenous data gg_miss_upset error for plate ', plate, sep ='')
    cat(to_report, file = report_file, sep = '\n', append = TRUE )
    cat(as.character(e[1]), file = report_file, sep = '\n', append = TRUE )
  })
  
  
  if (length(sample_names) > 1){
    endogenous<- colMeans(endogenous.plot, na.rm = T)
    exogenous<- colMeans(exogenous.plot, na.rm = T)
  }else{
    endogenous<- mean(endogenous.plot, na.rm = T)
    exogenous<- mean(exogenous.plot, na.rm = T)
  }
  
  
  normalized_data<- filter(filtered_data.plot, ID != "SNORD61" & ID != "SNORD68" & 
                             ID != "SNORD72" & ID != "SNORD95" & 
                             ID!="RNU6-2" &  ID != "cel-miR-39" &
                             ID != "PPC" & ID != "miRTC")
  
  #TO DO: the user could select the endogenous or exogenous normalazitaion
  
  if (normalization_en_ex == 'endogenous'){
    ##normalize with endogenous
    normalized_data[,!(names(normalized_data) %in% drops_id)]<-normalized_data[,!(names(normalized_data) %in% drops_id)]/endogenous
  }else if (normalization_en_ex == 'exogenous'){
    ##normalized with exogenous
    normalized_data[,!(names(normalized_data) %in% drops_id)]<-normalized_data[,!(names(normalized_data) %in% drops_id)]/exogenous
  } else{
    print('Wrong input value in normalization_en_ex variable!')
    return(NULL)
  }
  
  #log transform the data to develop normality in data
  
  colnames(Mnemiopsis_count_data) <- sample_names
  
  ##boxplot with quality
  #corrected data
  if (length(sample_names)>1){
    Mnemiopsis_count_data<- normalized_data[,!(names(normalized_data) %in% drops_id)]
    pseudoCount = log2(Mnemiopsis_count_data + 1) # log-transform to make numbers on scale (+1 to avoid zeroes)
    df = melt(pseudoCount, variable.name = "Samples", value.name = "count") # reshape the matrix
    
    #png(file=paste(output_dir,'/plate',plate,'/boxplot_normalized_',normalization_en_ex,'_plate_',plate,'.png', sep = ''), width=900, height=600)
    
    myplot <- ggplot(df, aes(x = Samples, y = count)) + geom_boxplot() + xlab("") +
      ylab(expression(log[2](count + 1)))+
      coord_flip()
    save_as_pdf({print(myplot)}, file.name = paste(output_dir,'/plate',plate,'/boxplot_normalized_',normalization_en_ex,'_plate_',plate,'.pdf', sep = ''), width = 6, height = 4)
    save_image({print(myplot)}, file.name = paste(output_dir,'/plate',plate,'/boxplot_normalized_',normalization_en_ex,'_plate_',plate,'.png', sep = ''))
    
    #print(myplot)
    #dev.off()
    
  } else {
    Mnemiopsis_count_data<- normalized_data[,!(names(normalized_data) %in% drops_id)]
    pseudoCount = log2(Mnemiopsis_count_data + 1) # log-transform to make numbers on scale (+1 to avoid zeroes)
    df = melt(pseudoCount, variable.name = "Samples", value.name = "count") # reshape the matrix 
    
    #png(file=paste(output_dir,'/plate',plate,'/boxplot_normalized_',normalization_en_ex,'_plate_',plate,'.png', sep = ''), width=900, height=600)
    #myplot <- boxplot(df, ylab = sample_names[1], xlab = "log2(count+1)", horizontal=T)
    save_as_pdf({boxplot(df, ylab = sample_names[1], xlab = "log2(count+1)", horizontal=T)}, file.name = paste(output_dir,'/plate',plate,'/boxplot_normalized_',normalization_en_ex,'_plate_',plate,'.pdf', sep = ''), width = 6, height = 4)
    save_image({boxplot(df, ylab = sample_names[1], xlab = "log2(count+1)", horizontal=T)}, file.name = paste(output_dir,'/plate',plate,'/boxplot_normalized_',normalization_en_ex,'_plate_',plate,'.png', sep = ''))
    
    #dev.off()
  }
  
  return(normalized_data)
}