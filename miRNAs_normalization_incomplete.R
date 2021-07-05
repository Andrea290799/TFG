library(rstatix) 


### FUNCTIONS

data_collection <- function (files_list, output_file){

  ## data_collection ---------------------------------------------------------
  ## Function that gets the data from all the input files and turn some of its 
  ## columns into a list. It also generates a control file. 
  ## files_list: list that contains the input files. 
  ## output_file: control file name. 
  ## It returns a list that contains the data. 
  
  # Number of the column where the specific data is. 
  column_plate_ID <- 3
  column_names <- 4
  column_ct <- 6
  column_A_Rn <- 7
  
  ct <- c() # Cts vector
  names <- c() # miRNA names vector
  plate_IDs <- c() # Plate IDs vector
  A_Rn <- c() # DRn vector
  
  groups <- c() # Groups vector 
  
  output_list <- list()
  
  # Vectors get filled
  for (i in files_list){ 
    
    data <- read.delim(file = i, skip = 4)
    
    for (j in data[,column_names]){
      names <- c(names, j) 
    }
    
    for (j in data[,column_plate_ID]){
      plate_IDs <- c(plate_IDs, j) 
    }

    
    for (j in data[,column_A_Rn]){
      A_Rn <- c(A_Rn, j) 
    }
    
    for (j in data[,column_ct]){
      ct <- c(ct, j) 
    }
  }
  
  # Groups are defined
  for (j in plate_IDs){
    
    if (startsWith(j, "G1")){
      groups <- c(groups, "ICS")
    }
    else if (startsWith(j, "G5")){
      groups <- c(groups, "C")
    }
    else {
      groups <- c(groups, "UC")
    }
  }
  
  output_list[["names"]] <- names
  output_list[["plate_IDs"]] <- plate_IDs
  output_list[["A_Rn"]] <- A_Rn
  output_list[["ct"]] <- ct
  output_list[["groups"]] <- groups
  
  # Control file
  write.table(data.frame("Plate_ID", "miRNA", "Ct", "ARn", "Group"), 
              file = output_file, append = FALSE,  sep = "\t", 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  output <- data.frame(plate_IDs, names, ct, A_Rn, groups)
  
  write.table(output, file = output_file, append = TRUE, 
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  return(output_list)
}



undetermined_IPC_delete <- function (data_list, output_file){

  ## undetermined_IPC_delete -------------------------------------------------
  ## Function that deletes the data of those plates that have 2 or 3 of its IPCs
  ## as undetermined. It also generates a control file.
  ## data list: list that contains the data. 
  ## output_file: control file name. 
  ## It returns a list that contains the data. 
  
  plates_to_remove <- c()
  
  for (j in 1:length(unique(data_list$plate_IDs))){
    # Indexes corresponding plate j that have detector name UniSp3_IPC and
    # undetermined ct value.  
    indexes_to_remove <- which(
      data_list$plate_IDs == unique(data_list$plate_IDs)[j] & 
        data_list$names == "UniSp3_IPC" & data_list$ct == "Undetermined")

    # If there are at least 2 indexes to remove, plate ID j is saved.
    if (length(indexes_to_remove) > 1){
      plates_to_remove <- c(plates_to_remove, 
                           unique(data_list$plate_IDs)[j])
    }
  }
  
  output <- c("Removed plates: ", plates_to_remove, "")

  write.table(output, file = output_file, append = FALSE, 
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  # Plates are deleted
  for (k in plates_to_remove){
    
    # Indexes corresponding to plate k.
    ind <- which(data_list$plate_IDs == k)
    
    data_list$names <- data_list$names[-ind]
    data_list$ct <- data_list$ct[-ind]
    data_list$plate_IDs <- data_list$plate_IDs[-ind]
    data_list$A_Rn <- data_list$A_Rn[-ind]
    data_list$groups <- data_list$groups[-ind]
  }
  
  
  # Control file
  write.table(data.frame("Plate_ID", "miRNA", "Ct", "ARn", "Group"), 
              file = output_file, append = TRUE,  sep = "\t", 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  output <- data.frame(data_list$plate_IDs, data_list$names, data_list$ct,
                      data_list$A_Rn, data_list$groups)
  
  write.table(output, file = output_file, append = TRUE, 
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  return(data_list)
  
}

cts_transformation_NA <- function(data_list, output_file){

  ## cts_transformation_NA ------------------------------------------------
  ## Function that transforms undetermined Ct values to NA. It 
  ## also generates a control file.
  ## data list: list that contains the data. 
  ## output_file: control file name. 
  ## It returns a list that contains the data. 
  
  undetermined = which(data_list$ct == "Undetermined")
  data_list$ct[undetermined] = "NA"
  
  
  # Control file
  write.table(data.frame("Plate_ID", "miRNA", "Ct", "ARn", "Group"), 
              file = output_file, append = TRUE,  sep = "\t", 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  output <- data.frame(data_list$plate_IDs, data_list$names, data_list$ct, 
                      data_list$A_Rn, data_list$groups)
  
  write.table(output, file = output_file, append = TRUE, 
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  return(data_list)
  
}


blank_NA_delete <- function(data_list, output_file){

  ## blank_NA_delete ---------------------------------------------------------
  ## Function that deletes those miRNAs whose cts are 40 or whose miRNA names
  ## are Blank. It also generates a control file.
  ## data list: list that contains the data. 
  ## output_file: control file name. 
  ## It returns a list that contains the data. 
  
  
  to_remove <- which(data_list$names == "Blank_(H2O)" | data_list$ct == "NA")
  
  data_list$plate_IDs <- data_list$plate_IDs[-c(to_remove)]
  data_list$names <- data_list$names[-c(to_remove)]
  data_list$ct <- data_list$ct[-c(to_remove)]
  data_list$groups <- data_list$groups[-c(to_remove)]
  data_list$A_Rn <- data_list$A_Rn[-c(to_remove)]
  
  
  # Control file
  write.table(data.frame("Plate_ID", "miRNA", "Ct", "ARn", "Group"), 
              file = output_file, append = FALSE,  sep = "\t", 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  output <-  data.frame(data_list$plate_IDs, data_list$names, data_list$ct, 
                       data_list$A_Rn, data_list$groups)
  write.table(output, file = output_file, append = TRUE, 
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  
  return(data_list)
  
}


interplate_normalization <- function(data_list, output_file,
                                     reffinder_normfinder_input){
  
  ## interplate_normalization ------------------------------------------------
  ## Function that do interplate_normalization on the data. It also generates a
  ## control file.
  ## data list: list that contains the data. 
  ## output_file: control file name. 
  ## It returns a list that contains the data. 
  
  # This vector will contain relevant information to show in control file
  plate_and_IPC_before_after <- c()
  
  # This vector will contain the means of the 3 IPC of each plate. 
  plate_IPC_means <- c() 
  
  for (i in unique(data_list$plate_IDs)){
    
    # IPCs of each plate
    IPCs <- data_list$ct[which(data_list$names == "UniSp3_IPC" &
                                data_list$plate_IDs == i)]
    
    # IPCs whose Cts differ more than 0.5 Cts are removed
    all_sum_dist <- c()
    
    if (length(IPCs) == 3){
      for (j in 1:length(IPCs)){
        sum_dist <- 0
        dist <- c()
        for (k in 1:length(IPCs)){
          sum_dist <- sum_dist + abs(as.numeric(IPCs[j])-as.numeric(IPCs[k]))
          dist <- c(dist, abs(as.numeric(IPCs[j])-as.numeric(IPCs[k])))
        }
        
        all_sum_dist <- c(all_sum_dist, sum_dist)
      }
      if (length(which(dist > 0.5)) > 0){
        plate_and_IPC_before_after <- c(plate_and_IPC_before_after, i, "Before",
                                       IPCs)
        IPCs <- IPCs[-which.max(all_sum_dist)]
        plate_and_IPC_before_after <- c(plate_and_IPC_before_after, "After",
                                       IPCs, "-----------")
      }
    }

    # Mean of IPCs is calculated
    IPCs_mean <- mean(as.numeric(IPCs)) 
    
    plate_IPC_means <- c(plate_IPC_means, IPCs_mean) 
    
  }
  
  write.table(c(plate_and_IPC_before_after, " "), file = output_file, 
              append = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE,
              quote = FALSE)
  
  plate_IPC_means_mean <- mean(na.omit(plate_IPC_means)) 
  
  # This vector will contain the normalizing factor of each plate
  interplate_norm_factor <- c() 
  
  for (i in 1:length(unique(data_list$plate_IDs))){
    interplate_norm_factor <- c(interplate_norm_factor, plate_IPC_means[i] -
                                 plate_IPC_means_mean)
  }
  
  # Normalizing factor is applied
  for (i in 1:length(interplate_norm_factor)){
    for (j in 1:length(data_list$ct)){
      if (data_list$plate_IDs[j] == unique(data_list$plate_IDs)[i]){
        if (data_list$ct[j] != 40){
          data_list$ct[j] <- as.numeric(data_list$ct[j])-
            interplate_norm_factor[i]
        }
      }
    }
  }
  
  data_list$ct[which(data_list$ct > 40)] <- 40
  
  # Control file
  write.table(data.frame("Plate_ID", "miRNA", "Ct", "ARn", "Group"), 
              file = output_file, append = TRUE,  sep = "\t", 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  output <-  data.frame(data_list$plate_IDs, data_list$names, data_list$ct, 
                       data_list$A_Rn, data_list$groups)
  write.table(output, file = output_file, append = TRUE, 
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  # RefFinder/NormFinder file
  write.table(data.frame("Plate_ID", "miRNA", "Ct", "ARn", "Group"), 
              file = reffinder_normfinder_input, append = FALSE,  sep = "\t", 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  output <-  data.frame(data_list$plate_IDs, data_list$names, data_list$ct, 
                        data_list$A_Rn, data_list$groups)
  write.table(output, file = reffinder_normfinder_input, append = TRUE, 
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  
  
  return(data_list)
}


GME <- function(data_list, output_file){

  ## GME ---------------------------------------------------------------------
  ## Function that applies GME normalization method. It also generates a
  ## control file.
  ## data list: list that contains the data. 
  ## output_file: control file name. 
  ## Returns a list that contains the data. 
  
  all_plates_cts_means <- c() 
  
  # All cts of a plate are listed, doing away with the spikes ones.
  for (i in unique(data_list$plate_IDs)){
    plate_cts <- c()
    for (j in 1:length(data_list$names)){
      
      #print(data_list$plate_IDs[j])
      #print(data_list$names[j])
      #print(data_list$ct[j])
      
      if (data_list$plate_IDs[j] == i & startsWith(data_list$names[j],
                                                   "hsa") == TRUE & 
                                                        data_list$ct[j] != 40 ){
        plate_cts <- c(plate_cts, as.numeric(data_list$ct[j]))
      }
    }
    
    plate_cts_mean <- mean(as.numeric(plate_cts)) 
    all_plates_cts_means <- c(all_plates_cts_means, plate_cts_mean)
    
  }
  
  A_ct <- c()
  # Normalization
  for (i in 1:length(all_plates_cts_means)){
    for (j in 1:length(data_list$ct)){
      if (data_list$plate_IDs[j] == unique(data_list$plate_IDs)[i]){
        A_ct[j] <- as.numeric(data_list$ct[j]) - all_plates_cts_means[i]
      }
    }
  }
  
  # Control file
  
  write.table(data.frame("Plate_ID", "miRNA", "ACt", "ARn", "Group"), 
              file = output_file, append = FALSE,  sep = "\t", 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  output <-  data.frame(data_list$plate_IDs, data_list$names, A_ct, 
                       data_list$A_Rn, data_list$groups)
  write.table(output, file = output_file, append = TRUE, 
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  data_list[["A_ct"]] <- A_ct
  
  return(data_list)
  
}


miRNA_normalization <- function(data_list, output_file, miRNA){

  ## miRNA_normalization -----------------------------------------------------
  ## Function that applies miRNA normalization method and deletes those plates
  ## that don't have ct value of the normalizer miRNA. It also generates a
  ## control file.
  ## data list: list that contains the data. 
  ## output_file: control file name.
  ## miRNA: miRNA chosen to normalize the others.
  ## Returns a list that contains the data. 
  
  # This vector will contain the ct value of the chosen miRNA to normalize
  normalizing_miRNA_cts_all_plates <- c()
  
  # If this value doesn't exist, the plate is lost
  plates_to_remove <- c()
  
  for (i in unique(data_list$plate_IDs)){ 
    
    # Index of the normalizer miRNA
    miRNA_index <- which(data_list$names == miRNA & data_list$plate_IDs == i)

    
    if (length(miRNA_index) == 0){ 
      
      # If there are not value for this miRNA, the plate has to be deleted.
      plates_to_remove <- c(plates_to_remove, i)
    }
    
    else{
      normalizing_miRNA_cts_all_plates <-c(normalizing_miRNA_cts_all_plates,
                                          data_list$ct[miRNA_index])
    }
    
  }
  
  for (i in plates_to_remove){
    indexes_to_remove <- which(data_list$plate_IDs == i)
    
    data_list$plate_IDs <- data_list$plate_IDs[-indexes_to_remove]
    data_list$names <- data_list$names[-indexes_to_remove]
    data_list$ct <- data_list$ct[-indexes_to_remove]
    data_list$A_Rn <- data_list$A_Rn[-indexes_to_remove]
    data_list$groups <- data_list$groups[-indexes_to_remove]
  }
  
  write.table(c("Removed plates:", plates_to_remove, ""), file = output_file,
              append = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE,
              quote = FALSE)
  
  
  # Normalization
  A_ct <- c()
  
  for (i in 1:length(normalizing_miRNA_cts_all_plates)){
    for (j in 1:length(data_list$ct)){
      if (data_list$plate_IDs[j] == unique(data_list$plate_IDs)[i]){
        A_ct[j] <- as.numeric(data_list$ct[j]) -
          as.numeric(normalizing_miRNA_cts_all_plates[i])
      }
    }
  }
  
  # Control file
  
  write.table(data.frame("Plate_ID", "miRNA", "ACt", "ARn", "Group"), 
              file = output_file, append = TRUE,  sep = "\t", 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  output <-  data.frame(data_list$plate_IDs, data_list$names, A_ct, 
                       data_list$A_Rn, data_list$groups)
  write.table(output, file = output_file, append = TRUE, 
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  data_list[["A_ct"]] <- A_ct
  
  return(data_list)
  
}


AACt <- function(data_list, output_file){

  ## AACt --------------------------------------------------------------------
  ## Function that does AACt normalization. It also generates a control file.
  ## data list: list that contains the data. 
  ## output_file: control file name. 
  ## Returns a list that contains the data. 
  
  AAct <- c()

    
    all_Act_group_miRNA_means <- c() 
    
    for (j in 1:length(unique(data_list$names))){
      
      A_ct <-  data_list$A_ct[which(data_list$groups == "C" & data_list$names
                                   == unique(data_list$names)[j])]
      
      Act_group_miRNA_mean <- mean(as.numeric(A_ct)) 
      all_Act_group_miRNA_means <- c(all_Act_group_miRNA_means,
                                    Act_group_miRNA_mean)
      
    }
    
    for (j in 1:length(unique(data_list$names))){
      for (i in 1:length(data_list$groups)){

        if (data_list$names[i] == unique(data_list$names)[j]){
          AAct[i] <- as.numeric(data_list$A_ct[i]) -
            as.numeric(all_Act_group_miRNA_means[j]) 
        }
      }
    }
  
  
  data_list[["AAct"]] <- AAct
  
  # Control file
  output <-  data.frame(data_list$plate_IDs, data_list$names, AAct, 
                       data_list$A_Rn, data_list$groups)
  write.table(output, file = output_file, append = TRUE, 
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  
  return(data_list)
  
  
}

DOS_AAct <- function(data_list, output_file){

  ## DOS_AAct ------------------------------------------------------------
  ## Function that calculates 2^-AACt.It generates the final results file.
  ## data list: list that contains the data. 
  ## output_file: control file name. 

  
  DOS_AAct_value <- c()
  for (i in 1:length(data_list$AAct)){
    DOS_AAct_value<- c(DOS_AAct_value, 2**(-data_list$AAct[i]))
  }
  
  # Final results file
  output <-  data.frame(data_list$plate_IDs, data_list$names, DOS_AAct_value, 
                       data_list$A_Rn, data_list$groups)
  write.table(output, file = output_file, append = TRUE, 
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  

}


main <- function(folder_name, plate_number, normalization_method, miRNA = NULL){

  ## main --------------------------------------------------------------------
  ## Function that executes the main code. 
  ## folder_name: folder that contains the data.
  ## plate_number: number of the analyzed plate. 
  ## normalization_method: method to use to normalize. Options: GME,
  ## normalizing_miRNAs.
  
  if (normalization_method != "GME" & normalization_method != 
      "normalizing_miRNAs"){
  
      stop("The normalization method introduced is not correct. Please, try
         introducing one of these options: GME or normalizing_miRNAs.")
    
  }
  
  if (normalization_method == "normalizing_miRNAs" & is.null(miRNA)){
    stop("You have to introduce the name of the normalizing miRNA.")
  }
  
  path <- getwd()
  
  if (dir.exists(paste(path, folder_name, sep="")) == FALSE){
    stop("Data folder introduced does not exist.")
  }
  
  
  # Files that contain the data are listed
  data_path <- paste(path, folder_name, sep = "")
  list_files <- list.files(data_path, full.names = TRUE)

  
  method_folder <- paste(path, "\\", normalization_method, sep = "")
  
  if (dir.exists(method_folder) == FALSE){dir.create(method_folder)}
  
  plate_folder <- paste(method_folder, "\\Plate_", plate_number, sep = "")
  
  if (dir.exists(plate_folder) == FALSE){dir.create(plate_folder)}
  
  # Output files
  control_file_data_collection <- paste(plate_folder, "\\1_Data_collection_",
                                       plate_number, ".txt", sep = "")
  
  control_file_undetermined_IPC_delete <- paste(plate_folder,
                                               "\\2_Undetermined_IPC_delete_",
                                               normalization_method, "_",
                                               plate_number, ".txt", sep = "")
  
  control_file_cts_transformation_40_NA <- paste(plate_folder,
                                                "\\3_Cts_transformation_NA_", 
                                                normalization_method, "_",
                                                plate_number, ".txt", sep = "")
  
  control_file_blank_NA_delete <- paste(plate_folder, "\\4_Blank_NA_delete_",
                                       normalization_method, "_",
                                       plate_number, ".txt", sep = "")
  
  
  control_file_interplate_normalization <- paste(plate_folder,
                                                "\\5_Interplate_normalization_",
                                                normalization_method, "_",
                                                plate_number, ".txt", sep = "")
  
  file_to_reffinder_normfinder <- paste(getwd(),
                                        "\\to_reffinder_normfinder_file_",
                                        plate_number, ".txt", sep = "")
  
  
  control_file_GME <- paste(plate_folder, "\\6_GME_", normalization_method, "_",
                           plate_number, ".txt", sep = "")
  
  
  control_file_miRNA_normalization <- paste(plate_folder,
                                           "\\6_miRNA_normalization_",
                                           normalization_method, "_",
                                           plate_number, ".txt", sep = "")
  
  
  control_file_AACt <- paste(method_folder, "\\7_AACt_", normalization_method,
                            ".txt", sep = "")
  
  final_results_file <- paste(getwd(), "\\8_RESULTS_", normalization_method,
                             ".txt", sep = "")
  
  if (plate_number == "1"){
    # Headers on common to both plates are written
    write.table(data.frame("Plate_ID", "miRNA", "AACt", "ARn", "Group"), 
                file = control_file_AACt, append = FALSE,  sep = "\t", 
                row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    write.table(data.frame("Plate_ID", "miRNA", "2^-AACt", "ARn", "Group"), 
                file = final_results_file, append = FALSE,  sep = "\t", 
                row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
  
  # All folder files data
  data <- data_collection(list_files, control_file_data_collection)
  
  if (normalization_method == "normalizing_miRNAs" & length(which(data$names ==
                                                                  miRNA)) == 0){
    stop("The introduced miRNA does not exist in the introduced plate.")
  }
  
  cat("\014")
  print("|***_                            | 1/8", quote = FALSE)
  
  # Data without IPC failures
  data <- undetermined_IPC_delete(data, control_file_undetermined_IPC_delete)
  
  cat("\014")
  print("|***_***_                        | 2/8", quote = FALSE)
  
  # NA ct transformation
  data <- cts_transformation_NA(data, control_file_cts_transformation_40_NA)

  cat("\014")
  print("|***_***_***_                    | 3/8", quote = FALSE)
  
  # Blank and NA elimination
  data <- blank_NA_delete(data, control_file_blank_NA_delete)
  
  cat("\014")
  print("|***_***_***_***_                | 4/8", quote = FALSE)
  
  # Interplate normalization
  data <- interplate_normalization(data, control_file_interplate_normalization,
                                   file_to_reffinder_normfinder)
  
  cat("\014")
  print("|***_***_***_***_***_            | 5/8", quote = FALSE)
  
  # ACt
  if (normalization_method == "GME"){
    data <- GME(data, control_file_GME)
    
    cat("\014")
    print("|***_***_***_***_***_***_        | 6/8", quote = FALSE)
  }
  
  else if (normalization_method == "normalizing_miRNAs"){
    
    data <- miRNA_normalization(data, control_file_miRNA_normalization, miRNA)
    
    cat("\014")
    print("|***_***_***_***_***_***_        | 6/8", quote = FALSE)
  }
  
  # AACt
  data <- AACt(data, control_file_AACt)
  
  cat("\014")
  print("|***_***_***_***_***_***_***_    | 7/8", quote = FALSE)
  
  # 2^AACt
  
  DOS_AAct(data, final_results_file)
  
  cat("\014")
  print("|***_***_***_***_***_***_***_***_| 8/8", quote = FALSE)
  
}



### CODE

# GME
main("\\Plate_1_data", "1", "GME")
main("\\Plate_2_data", "2", "GME")

# normalizer miRNAs
#main("\\Plate_1_data", "1", "normalizing_miRNAs", "hsa-miR-198")
#main("\\Plate_2_data", "2", "normalizing_miRNAs", "hsa-miR-200c-5p")
  


