library("car")
library("ggplot2")
library("rstatix")
library("tidyverse")
library("VennDiagram")
library("WRS2")


### FUNCTIONS

data_collection <- function(file_name){
  
  ## data_collection ----------------------------------------------------
  ## Function that gets the data to analyze. 
  ## file_name: name of file that contains the data to be analyzed. 
  ## Returns a list with the data.
  
  output_list <- list()
  
  path <- getwd()
  data <- read.delim(file = paste(path, file_name, sep = ""), skip = 1)
  
  output_list[["cts"]] <- c(data[,3])
  output_list[["miRNAs"]] <- c(data[,2])
  output_list[["groups"]] <- c(data[,5])
  
  return(output_list)
  
}

statistical_tests_pvalues <- function(cts, groups = NULL, test){

  ## statistical_tests_pvalues -----------------------------------------------
  ## Function that does different statistical tests: 
  ## normality: saphiro.test function
  ## homoscedasticity: leveneTest function
  ## Kruskal-Wallis: kruskal.test function
  ## Pairwise Wilcoxon Rank Sum Tests: pairwise.wilcox.test function
  ## ANOVA: aov function
  ## Pairwise t tests: pairwise.t.test function
  ## Robust one-way ANOVA based on trimmed means: t1way function
  ## Post-hoc t1way: lincon function
  
  ## cts: miRNA cts that will take part of the analysis.
  ## groups: miRNA groups that will take part of the analysis.
  ## It returns a vector with the corresponding p-values.

  # The vector to be returned
  pvalues <- c()

  # NORMALITY
  if (test == "normality"){
  
    # There must be at least 3 values to do this test.
    if (length(cts) > 2){ 
      pvalue_norm <- tryCatch(shapiro.test(unlist(as.numeric(cts)))$p.value, 
                             error = function(e){pvalue_norm = NULL})
        
      if (class(pvalue_norm[1]) != "numeric"){
        pvalue_norm <- NA
      }
    }
        
    else { 
      pvalue_norm <- NA
    }
        
    pvalues <- pvalue_norm

  }
  
  # HOMOCEDASTICITY
  if (test == "homocedasticity"){
    
    # There must be at least 2 groups to do this test.
    if(length(unique(groups))>= 2){ 
      
      pvalue_levene <- leveneTest(y = cts, group = factor(groups), center =
                                   "median")$Pr[1]
    }
    else{
      pvalue_levene <- NA
    }
      
    pvalues <- pvalue_levene
  }
  
  # KRUSCAL-WALLIS
  if (test == "KW"){
    if (length(unique(groups))>= 2){
      pvalue_KW <- tryCatch(kruskal.test(cts, groups)$p.value,
                            error = function(e){pvalue_norm = NA})
    }
    else {pvalue_KW <- NA}
      
     pvalues <-  pvalue_KW
  }
    
  # ANOVA
  if (test == "ANOVA"){
    if (length(unique(groups))>= 2){
      anova_test <- aov(cts ~ groups)
      pvalue_ANOVA <- summary(anova_test)[[1]][["Pr(>F)"]][[1]]
    }
    else {pvalue_ANOVA <- NA}
      
    pvalues <- pvalue_ANOVA
      
  }
  
  # Pairwise Wilcoxon Rank Sum Tests 
  if (test == "post-hoc_KW"){
    if (length(unique(groups)) >= 2){
      pvalue_post_hoc_KW <- pairwise.wilcox.test(cts, groups, 
                                                p.adjust.method = "fdr")$p.value
    }
    pvalues <- pvalue_post_hoc_KW
  }
  
  # Pairwise t tests
  if (test == "post-hoc_ANOVA"){
    if (length(unique(groups)) >= 2){
      pvalue_post_hoc_ANOVA <- pairwise.t.test(cts, groups, 
                                                p.adjust.method = "fdr")$p.value
    }
      
    pvalues <- pvalue_post_hoc_ANOVA
  }
  
  # Robust one-way ANOVA based on trimmed means
  if (test == "robust_methods"){
    
    pvalues <- tryCatch(t1way(formula = cts ~ groups, data =
                  data.frame(groups, cts))$p.value, error = function(e)
                    {pvalues = NULL})

    
  }
  
  if (test == "robust_methods_post_hoc"){
    pvalues <- tryCatch(lincon(formula = cts ~ groups, data =
                       data.frame(groups, cts))$comp[,6], error = function(e)
                       {pvalues = NULL})
    
  }
  

  return(pvalues)
  
}



boxplots <- function(cts, groups, miRNA, boxplot_folder){

  ## boxplots ----------------------------------------------------------------
  ## Function that obtains boxplots for given data.
  ## cts: data ct values.
  ## groups: groups where ct values belong to.
  ## miRNA: miRNA to which ct values belong to.
  
  Grupos <- groups
  
  # We do away with outliers
  to_remove <- which(cts > 40)
    if (length(to_remove) > 0){
      cts <- cts[-to_remove]
      Grupos <- Grupos[-to_remove]
    }

  
  ggplot(data.frame(Grupos, cts), aes(x=as.factor(Grupos), y=cts, fill =
                                        Grupos)) + 
    geom_boxplot(outlier.colour="black", outlier.shape=20,
                 outlier.size=2) +
    labs(x = "Grupos",y = "2^-DDCt", title = miRNA) +
    scale_fill_manual(values=c("#2effd2", "#27c9a3", "#198068")) + 
    theme_classic()
  
  ggsave(paste(boxplot_folder, "/", miRNA, ".pdf", sep=""))
  
}

venn_diagram <- function(file_1, file_2, label_1, label_2, common_method){
  
  ## venn_diagram -------------------------------------------------------------
  ## Function that plots a Venn diagram from file_1, file_2 data. 
  ## file_1, file_2: data files.
  ## label_1, label_2: labels in the diagram.
  ## common method: method that is common to both files
  
  if (file.exists(paste(getwd(),file_1, sep="\\")) == FALSE){
    stop("Introduced data file ", file_1," does not exist.")
  }
  
  if (file.exists(paste(getwd(),file_2, sep="\\")) == FALSE){
    stop("Introduced data file ", file_2," does not exist.")
  }
  
  venn_diagrams_folder <- paste(getwd(), "\\venn_diagrams_folder", sep= "")
  
  if (dir.exists(venn_diagrams_folder) == FALSE){
    dir.create(venn_diagrams_folder)
  }
  
  data <- read.delim(file = file_1)
  method1 <- c(data[,1])
  ICS_C_1 <- c(data[,2])
  UC_C_1 <- c(data[,4])
  UC_ICS_1 <- c(data[,6])
  
  all_1 = list(ICS_C_1, UC_C_1, UC_ICS_1)
  
  data <- read.delim(file = file_2)
  method2 <- c(data[,1])
  ICS_C_2 <- c(data[,2])
  UC_C_2 <- c(data[,4])
  UC_ICS_2 <- c(data[,6])
  
  all_2= list(ICS_C_2, UC_C_2, UC_ICS_2)
  
  titles = c("ICS vs. C", "UC vs. C", "UC vs. ICS")
  names = c("ICS_C", "UC_C", "UC_ICS")
  
  for (i in 1:length(all_1)){
    dif1 = method1[which (all_1[[i]] < 0.05)]
    dif2 = method2[which (all_2[[i]] < 0.05)]

    
    # Common to both groups are counted
    n <- 0
    for (j in dif1){
      n <- n + length(which(dif2 == j))
    }
    
    a1 = length(dif1)
    a2 = length(dif2)
    
    pdf(paste(venn_diagrams_folder, "\\", names[i], "-", label_1, "_", label_2,
              "_venn_diagram_", common_method, ".pdf", sep =""))
    draw.pairwise.venn(area1 = a1, area2 = a2, cross.area = n,
                       category = c(label_1, label_2),
                       lty = rep("blank", 2), fill = c("#2effd2", "#198068"),
                       cat.pos = c(-30,30), main = titles[i])
    dev.off()
    
  }
}



main <- function(file_name, method){

  ## main --------------------------------------------------------------------
  ## Function that executes the main code. 
  ## file_name: name of file that contains the data to be analyzed.
  ## method: "robust" or "not robust" statistical methods. 
  
  if (file.exists(paste(getwd(),file_name, sep="")) == FALSE){
    stop("Introduced data file does not exist.")
  }
  
  data <- data_collection(file_name)
  
  if (method != "not robust" & method != "robust"){

    stop("Method not recognized. Please, introduce 'not robust' or 'robust'.")

  }
  
  else if (method == "not robust"){
    
    boxplot_folder <- paste(getwd(), "\\ANOVA_KW", str_sub(file_name, 11,
                                               length(file_name)-6),
                                                "_BOXPLOTS", sep="")
    
    if (dir.exists(boxplot_folder) == FALSE){dir.create(boxplot_folder)}
    
    # Results file header
    write.table(data.frame("miRNA", "ICS_vs_C_pvalue", "ICS_vs_C_log2(FC)",
                           "UC_vs_C_pvalue", "UC_vs_C_log2(FC)",
                           "UC_vs_ICS_pvalue",
                           "UC_vs_ICS_log2(FC)"), file = 
                  paste("ANOVA_KW_differences", str_sub(file_name, 11, ),
                                                          sep = ""),
                append = FALSE, sep = "\t", row.names = FALSE,
                col.names = FALSE, quote = FALSE)
    
    # Each miRNA population is studied separately.
    for (i in unique(data$miRNAs)){
      # Spikes are not analyzed
      if (startsWith(i, "hsa") == TRUE | startsWith(i, "cel") == TRUE ){
        pvalues_norm <- c()
        pvalue_levene <- ""
        pvalue_ANOVA <- ""
        results_post_hoc_ANOVA <- c()
        pvalue_KW <- ""
        results_post_hoc_KW <- c()
        
        # Indexes corresponding to miRNA i, regardless the group
        all_miRNA_index <- which(data$miRNAs == i) 
        
        all_miRNA_cts <- data$cts[all_miRNA_index]
        all_miRNA_groups <- data$groups[all_miRNA_index]
        
        for (j in c("ICS", "UC", "C")){
          
          # Indexes corresponding to miRNA i, regarding the group j
          group_miRNA_index <- which(data$miRNAs == i & data$groups == j)
    
          group_miRNA_cts <- data$cts[group_miRNA_index]
          
          # If n >= 7, you can aim to do an ANOVA if the data passes the
          # normality and homoscedasticity tests.
          if (length(group_miRNA_cts) >= 7){ 
            
            # Normality test p-values
            pvalues_norm <- c(pvalues_norm, statistical_tests_pvalues(
              group_miRNA_cts, test = "normality"))
    
          }
        }
        
    
        # If the miRNA data passes the normality test in the 3 groups...
        if (is.null(pvalues_norm) == FALSE){
          if (length(pvalues_norm[which(pvalues_norm > 0.05)]) == 3){ 
      
            # Homocedasticity test p-value
            pvalue_levene <- statistical_tests_pvalues(all_miRNA_cts,
                                                      groups = all_miRNA_groups,
                                                      "homocedasticity")
            
          
            # If the miRNA data passes the homocedasticity test...
            if (pvalue_levene > 0.05){
             
              # ANOVA p-value 
              pvalue_ANOVA <- statistical_tests_pvalues(all_miRNA_cts, groups =
                                                         all_miRNA_groups,
                                                        "ANOVA")
              
              if (pvalue_ANOVA <= 0.05 & is.nan(pvalue_ANOVA) == FALSE){
                
                results_post_hoc_ANOVA = statistical_tests_pvalues(
                  all_miRNA_cts, groups = all_miRNA_groups,"post-hoc_ANOVA")
              }
            }
          }
        }
        
        # If the miRNA data does not pass normality or homocedasticity tests...
        if (pvalue_ANOVA == ""){
          
          # KW test p-value
          pvalue_KW <- statistical_tests_pvalues(all_miRNA_cts, groups =
                                                  all_miRNA_groups, "KW")
          if (is.na(pvalue_KW) == FALSE){
            if (pvalue_KW <= 0.05){
              results_post_hoc_KW = statistical_tests_pvalues(all_miRNA_cts,
                                                              groups = 
                                                                all_miRNA_groups,
                                                              "post-hoc_KW")
            }
          }
        }
        
        # Fold Change
        cts_ICS_mean <- mean(all_miRNA_cts[which(all_miRNA_groups == "ICS")])
        cts_UC_mean <- mean(all_miRNA_cts[which(all_miRNA_groups == "UC")])
        cts_C_mean <- mean(all_miRNA_cts[which(all_miRNA_groups == "C")])
        
        
        # Only if one value is <= 0.05, the results are got
        if (length(which(results_post_hoc_ANOVA <= 0.05)) != 0){ 
          
          write.table(data.frame(i, results_post_hoc_ANOVA[1],
                                 log2(cts_ICS_mean/cts_C_mean), 
                                 results_post_hoc_ANOVA[2],
                                 log2(cts_UC_mean/cts_C_mean),
                                 results_post_hoc_ANOVA[4], 
                                 log2(cts_UC_mean/cts_ICS_mean)),
                      file = paste("ANOVA_KW_differences", 
                                   str_sub(file_name, 11, ),
                                   sep = ""), append = TRUE, sep = "\t",
                      row.names = FALSE, col.names = FALSE, quote = FALSE)
          
          boxplots(all_miRNA_cts, all_miRNA_groups, i, boxplot_folder)
          
          
        }
        
        else if (length(which(results_post_hoc_KW <= 0.05)) != 0){ 
          
          write.table(data.frame(i, results_post_hoc_KW[1],
                                 log2(cts_ICS_mean/cts_C_mean),
                                 results_post_hoc_KW[2],
                                 log2(cts_UC_mean/cts_C_mean),
                                 results_post_hoc_KW[4],
                                 log2(cts_UC_mean/cts_ICS_mean)),
                      file =  paste("ANOVA_KW_differences",
                                    str_sub(file_name, 11, ),
                                    sep = ""), append = TRUE, sep="\t",
                      row.names = FALSE, col.names = FALSE, quote = FALSE)
          
          boxplots(all_miRNA_cts, all_miRNA_groups, i, boxplot_folder)
          
        }
      }
    }
  }
  
  else{
    
    boxplot_folder <- paste(getwd(), "\\WRS2", str_sub(file_name, 11, 
                                           length(file_name)-6),
                           "_BOXPLOTS", sep="")
    
    if (dir.exists(boxplot_folder) == FALSE){dir.create(boxplot_folder)}
    
    write.table(data.frame("miRNA", "ICS_vs_C_pvalue", "ICS_vs_C_log2(FC)",
                           "UC_vs_C_pvalue", "UC_vs_C_log2(FC)",
                           "UC_vs_ICS_pvalue", "UC_vs_ICS_log2(FC)"),
                file =  paste("WRS2_differences", str_sub(file_name, 11, ),
                                                          sep = ""),
                append = FALSE, sep = "\t", row.names = FALSE, col.names =
                  FALSE, quote = FALSE)
    
    
    for (i in unique(data$miRNAs)){
      
      if (startsWith(i, "hsa") == TRUE | startsWith(i, "cel") == TRUE){
        
        pvalue_post_hoc <- ""
        
        miRNA_indexes <- which(data$miRNAs == i)
        
        #robust one-way ANOVA based on trimmed means for each miRNA
        pvalue <- statistical_tests_pvalues(data$cts[miRNA_indexes], 
                                           data$groups[miRNA_indexes],
                                           "robust_methods")
        
        
        if (is.null(pvalue) == FALSE){
          if (pvalue <= 0.05){
            
            #post hoc
            pvalue_post_hoc <- statistical_tests_pvalues(
              data$cts[miRNA_indexes], data$groups[miRNA_indexes],
              "robust_methods_post_hoc")
            
            if (is.null(pvalue_post_hoc) == FALSE){
    
              # Fold Change
              cts_ICS_mean <- mean(data$cts[miRNA_indexes][which(
                data$groups[miRNA_indexes] == "ICS")])
              cts_UC_mean <- mean(data$cts[miRNA_indexes][which(
                data$groups[miRNA_indexes] == "UC")])
              cts_C_mean <- mean(data$cts[miRNA_indexes][which(
                data$groups[miRNA_indexes] == "C")])
              
              write.table(data.frame(i, pvalue_post_hoc[1],
                                     log2(cts_ICS_mean/cts_C_mean),
                                     pvalue_post_hoc[2],
                                     log2(cts_UC_mean/cts_C_mean),
                                     pvalue_post_hoc[3],
                                     log2(cts_UC_mean/cts_ICS_mean)),
                          file = paste("WRS2_differences",
                                       str_sub(file_name, 11, ),
                                       sep = ""), append = TRUE, 
                          sep = "\t", row.names = FALSE, col.names = FALSE, 
                          quote = FALSE)
              
              boxplots(data$cts[miRNA_indexes], data$groups[miRNA_indexes], i,
                       boxplot_folder)
            }
          }
        }
      }
    }
  }
}

### CODE

main("\\8_RESULTS_GME.txt", "robust")
main("\\8_RESULTS_GME.txt", "not robust")

venn_diagram("WRS2_differences_GME.txt", "ANOVA_KW_differences_GME.txt",
             "Robustos", "Clásicos", "GME")


#main("\\8_RESULTS_normalizing_miRNAs.txt", "robust")
#main("\\8_RESULTS_normalizing_miRNAs.txt", "not robust")
