
#FUNCTIONS

data_collection <- function(data_file){
  
  ## data_collection ---------------------------------------------------------
  ## Function that gets the data from all the input files and turns some of its 
  ## columns into a list. 
  ## data_folder: data folder that contains the input files. 
  ## It returns a list that contains the data. 
  
  # Future data list
  output_list = list()
  
  if (file.exists(paste(getwd(),data_file, sep ="//")) == FALSE){
    stop("Introduced file does not exist.")
  }
  
  column_groups <- 5
  column_ct <- 3
  column_names <- 2
  column_plate_ID <- 1
  
  ct = c() 
  names = c() 
  plate_IDs = c() 
  groups = c() 
  


  data <- read.delim(file = data_file)
    
  # Detector names are got
  for (j in data[,column_names]){
    names = c(names, j)
  }
    
  # Repeated UniSp3_IPC names are renamed
  check = 0
  for (x in 1:length(names)){
    if (names[x] == "UniSp3_IPC"){
       if (check == 0){
        names[x] = "UniSp3_IPC_1"
        check = 1
      }
      else if (check == 1){
        names[x] = "UniSp3_IPC_2"
        check = 2
      }
      else{
        names[x] = "UniSp3_IPC_3"
        check=0
          
      }
    }
  }
    
  # Plate Ids are got
  for (j in data[,column_plate_ID]){
    plate_IDs = c(plate_IDs, j)
  }
    
  # Ct values are got. Undetermined values are removed
  for (j in data[,column_ct]){ 
    ct = c(ct,j)
  }
  
  # Groups are got
  for (j in data[,column_groups]){
    groups = c(groups, j)
  }
  
  output_list[["names"]] <- names
  output_list[["plate_IDs"]] <- plate_IDs
  output_list[["ct"]] <- ct
  output_list[["groups"]] <- groups

  return(output_list)
}

main <- function(data_file, method, number){
  
  ## main --------------------------------------------------------------------
  ## Function that executes the main code. 
  ## data_folder: folder that contains the data.
  ## method: method to search for suitable miRNAs to use them for normalize.
  
  data_list = data_collection(data_file)

  
  if (method == "NormFinder"){
    
    column1 = c(" ")
    
    for (i in unique(data_list$names)){
      column1 = c(column1, i)
    }
    
    column1 = c(column1, "group")
    
    matrix = data.frame(column1)

    for (i in unique(data_list$plate_IDs)){
      
      column = c()
      column = c(column, i)
      
      cts_vector = c()
      
      #print(unique(data_list$groups[which(data_list$plate_IDs == i)]))

      for (j in unique(data_list$names)){
        if (length(which(data_list$names == j & data_list$plate == i )) == 1){
          cts_vector = c(cts_vector, data_list$ct[which(data_list$names == j & data_list$plate == i )])
        }
        else{cts_vector = c(cts_vector, NA)}
      }
          
        
      column = c(column, cts_vector) 
        
      column = c(column, unique(data_list$groups[which(data_list$plate_IDs == i)]))
      
      matrix = cbind(matrix, column)
      
    }
    
    # Input file
    output_file = paste(getwd(), "\\NormFinder_input_", number, ".txt", sep = "")
    
    write.table(matrix, file = output_file, append = FALSE, 
                sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    source("r.NormOldStab5.txt")
    
    # Results file
    sink(paste(getwd(), "\\NormFinder_results_", number, ".txt", sep = ""))
    print(Normfinder(output_file))
    sink()
    
  }
  
  else if (method == "RefFinder"){
    
    for (i in unique(data_list$names)){
      
      column = c()
      column = c(column,i)
      
      cts_vector = c()
      
      #print(unique(data_list$groups[which(data_list$plate_IDs == i)]))
      
      for (j in unique(data_list$plate_IDs)){
        if (length(which(data_list$names == i & data_list$plate == j )) == 1){
          cts_vector = c(cts_vector, data_list$ct[which(data_list$names == i & data_list$plate == j )])
        }
        else{cts_vector = c(cts_vector, NA)}
      }
      
      column = c(column, cts_vector)
      
      # First column
      if (i == data_list$names[1]){
        matrix = data.frame(column)
      }
      else{
        matrix = cbind(matrix, column)
      }
      
    }
    
    # Input file
    output_file = paste(getwd(), "\\RefFinder_input_", number, ".txt", sep = "")
    
    write.table(matrix, file = output_file, append = FALSE, 
                sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    
  }
  
  else{stop("Invalid chosen method. Please, introduce one of these methods: 
            'NormFinder' or 'RefFinder'.")}
  
}

# CODE

main("to_reffinder_normfinder_file_1.txt", "NormFinder", "1")
main("to_reffinder_normfinder_file_2.txt", "NormFinder", "2")

main("to_reffinder_normfinder_file_1.txt", "RefFinder", "1")
main("to_reffinder_normfinder_file_2.txt", "RefFinder", "2")
