# This script reproduces the single cell RNAseq analysis of human pancreatic 
# samples and PBMCs in Figures 2A, 2B, 2C, 2D, 2E, 2F, 3A, 3B, 3C, 3D, 3E, 3F,
# 4A, 4B, 4C, 4D, 4E, 4F, 4G, 5A, 5B, 5C, 5D, 5E, 5F, 6A, 6B, 6C, 6D, 6E, 6F, 6G 
# and Supplementary Figures S2B, S2C, S2D, S3A, S3B, S3C, S3D, S3E, S4A, S4B, S4C, 
# S5A, S5B, S5C, S5D, S5E of the paper "Multimodal Mapping of the Tumor 
# and Peripheral Blood Immune Landscape in Human Pancreatic Cancer" 
# The raw data was processed in line with the 
# Seurat workflow outlined on the Satija Lab website 
# (https://satijalab.org/seurat/v3.2/pbmc3k_tutorial.html) as well as in the 
# following reference:
#
# Stuart T, Butler A, Hoffman P, Hafemeister C, Papalexi E, III WMM, Hao Y, 
# Stoeckius M, Smibert P, Satija R (2019). "Comprehensive Integration of 
# Single-Cell Data." Cell, 177, 1888-1902. doi: 10.1016/j.cell.2019.05.031, 
# https://doi.org/10.1016/j.cell.2019.05.031.

#Load required packages

library(Seurat)
library(dplyr)
library(magrittr)
library(data.table)
library(Matrix)
library(devtools)
library(RcppArmadillo)
library(Rcpp)
library(scales)
library(pheatmap)
library(gplots)
library(ggplot2)
library(cowplot)
library(tibble)
library(xlsx)
library(data.table)


#Install batch correction package harmony
install_github("immunogenomics/harmony")
devtools::install_github("immunogenomics/harmony")
library(harmony)

#If using processed data from GEO/dbGaP, skip ahead to section Figures/Data Analysis

# Table of Contents

# 1. Functions
# 2. Data Pre-processing
# 3. Interactome Analysis
# 4. Figures/Data Analysis

#____________________________________________________________________________________________________________________________________________________________________#

# 1. Functions

#Interactome Code

check_genes <- function(genes, database, object_genes) {
  'Check to make sure that wanted genes are in reference you provide and in object
  
   Args:
   genes (chr vector): list of potential wanted genes
   database (data.frame): table with reference ligand/receptor pairs
   object_genes (chr vector): list of genes present in Seurat object
   
   Returns:
   not_found_list (chr list): list of genes not found in database and/or object
  '
  
  database_genes <- as.vector(unlist(database))
  not_found_database <- c()
  not_found_object <- c()
  
  for (i in 1:length(genes)) {
    if (!(genes[i] %in% database_genes)) {
      not_found_database <- c(not_found_database, genes[i])
    }
    
    if (!(genes[i] %in% object_genes)) {
      not_found_object <- c(not_found_object, genes[i])
    }
  }
  
  not_found_list <- list(database=not_found_database, object=not_found_object)
  return(not_found_list)
}

make_LR_pairs <- function(ligands, receptors, database) {
  'Make all LR pairs based on wanted ligands and/or receptors and database provided
  
   Args:
   ligands (chr vector): list of wanted ligands
   receptors (chr vector): list of wanted receptors
   database (data.frame): table with reference ligand/receptor pairs
   
   Returns:
   wanted_LR (data.frame): data.frame with ligand/receptor pairs from wanted ligands and receptors
  '
  
  wanted_LR <- data.frame(Ligand = character(0), Receptor = character(0))
  
  for (ligand in ligands){
    # list of corresponding receptors
    corresponding_receptors <- unique(as.character(database[,2][grep(ligand, database[,1])]))
    
    for (receptor in corresponding_receptors) {
      LR_row <- data.frame(Ligand = ligand, Receptor = receptor)
      wanted_LR <- rbind(wanted_LR, LR_row)
    }
  }
  
  # filter out unwanted receptors
  wanted_LR <- wanted_LR[which(wanted_LR$Receptor %in% wanted_receptors),]
  
  return(wanted_LR)
}

create_LR_table <- function(ligands, receptors, cell_types, LRs, avg0, avg1) {
  'Create table w/ ligand, receptor, source and target cells, average expressions, and IDs
   
   Args:
   ligands (chr vector): list of wanted ligands
   receptors (chr vector): list of wanted receptors
   cell_types (chr vector): list of common cell types between two objects
   LRs (data.frame): table with with wanted ligand/receptor pairs (from make_LR_pairs function)
   avg0 (list of num data.frame): average expression table of object 0
   avg1 (list of num data.frame): average expression table of object 1
   
   Returns:
   LR_table (data.frame): table of potential ligand/receptor pairs
    Contains ligand, receptor, source, target, average expression of ligand and receptors from 
    source and target cells. Also contains IDs (important for Cytoscape). 
  '
  
  LR_table <- data.frame()
  count <- 0
  
  for (ligand in ligands) {
    known_receptors <- LRs$Receptor[which(LRs$Ligand == ligand)]
    
    for (receptor in known_receptors) {
      for (i in c(1:length(cell_types))) {
        for (j in c(1:length(cell_types))) {
          LR_table <- rbind(LR_table,data.frame(ligands = ligand, receptors = receptor, 
                                                source = cell_types[i], target = cell_types[j], 
                                                avg_lig_0 = avg0[ligand, cell_types[i]],
                                                avg_rec_0 = avg0[receptor, cell_types[j]],
                                                avg_lig_1 = avg1[ligand, cell_types[i]],
                                                avg_rec_1 = avg1[receptor, cell_types[j]]))
        }
      }
    }
    
    cat(count, ' ')
    count <- count+1
  }
  
  # Create IDs
  source_ID <- c()
  target_ID <- c()
  
  # Find optimal cell type IDs
  ids_key <- c()
  
  for (i in 1:length(cell_types)) {
    for (j in 1:min(nchar(cell_types))) {
      name_1 <- substring(cell_types[i], 1, j)
      name_list <- c()
      name_list <- c(sapply(cell_types[-i], function(x) name_list <- c(name_list, substring(x, 1, j))))
      
      if(name_1 %in% name_list) {
        next
      }
      else {
        ids_key <- c(ids_key, name_1) 
        break
      }
    }
  }
  
  names(ids_key) <- cell_types
  
  for (i in c(1:length(LR_table[,1]))) {
    letter1 <- as.character(ids_key[LR_table$source[i]])
    letter2 <- as.character(ids_key[LR_table$target[i]])
    n1 <- which(ligands == LR_table$ligands[i])
    n2 <- which(receptors == LR_table$receptors[i])
    
    source_ID <- c(source_ID, paste(letter1, 'L', n1, sep = ''))
    target_ID <- c(target_ID, paste(letter2, 'R', n2, sep = ''))
  }
  
  LR_table <- cbind(LR_table, data.frame(source_ID = source_ID, target_ID = target_ID))
  
  return(LR_table)
}

avg_LR_filt <- function(table, threshold) {
  'Calculates states (ON/OFF --> 1/0) and filter if there is no expression in both groups
  
   Args:
   table (data.frame): table with potential ligand/receptors (from create_LR_table)
   threshold (num): average expression threshold 
   
   Returns:
   table (data.frame): filtered table based on average expression
  '
  
  # Find states of pairs in each group
  LR_states <- data.frame(lig_0 = character(0), rec_0 = character(0),
                          lig_1 = character(0), rec_1 = character(0),
                          is_on = logical(0))
  
  for (i in c(1:length(LR_table$avg_lig_0))) {
    row_states <- data.frame(lig_0 = table$avg_lig_0[i] > threshold, 
                             rec_0 = table$avg_rec_0[i] > threshold,
                             lig_1 = table$avg_lig_1[i] > threshold, 
                             rec_1 = table$avg_rec_1[i] > threshold)
    row_states$is_on <- (row_states$lig_0 & row_states$rec_0) | (row_states$lig_1 & row_states$rec_1)
    
    LR_states <- rbind(LR_states, row_states)
  }
  
  table <- cbind(table, ifelse((LR_states$lig_0 & LR_states$rec_0), 1, 0))
  table <- cbind(table, ifelse((LR_states$lig_1 & LR_states$rec_1), 1, 0))
  
  colnames(table)[11] <- 'is_on_0'
  colnames(table)[12] <- 'is_on_1'
  
  # Filter out pairs if pairs in both group are less than threshold 
  table <- table[which(LR_states$is_on == TRUE),]
  
  return(table)
}

LR_diff <- function(table, data_0, data_1, genes, label, alpha = 0.05) {
  'Calculate Wilcoxon-rank test on ligands and receptors (separately) between groups
    Order of test: data_0 --> data_1
  
    Args:
    table (data.frame): table with potential ligand/receptors (from create_LR_table)
    data_0 (data.frame): table of gene expression from object 0
    data_1 (data.frame): table of gene expression from object 1
    genes (chr vector): list of wanted genes
    label (chr): name of meta.data slot in Seurat object
    alpha (num): alpha level for Wilcoxon-rank test
    
    Returns:
    table (data.frame): filtered table based on Wilcoxon-rank test
  '
  
  table$lig_diff <- rep(0, length(table[,1]))
  table$rec_diff <- rep(0, length(table[,1]))
  
  table$lig_diff_p <- rep(0, length(table[,1]))
  table$rec_diff_p <- rep(0, length(table[,1]))
  
  for (i in 1:length(table[,1])) {
    ligand <- table$ligands[i]
    receptor <- table$receptors[i]
    source <- table$source[i]
    target <- table$target[i]
    
    lig_0_data <- data_0[which(data_0[,label] == source), ligand]
    rec_0_data <- data_0[which(data_0[,label] == target), receptor]
    
    lig_1_data <- data_1[which(data_1[,label] == source), ligand]
    rec_1_data <- data_1[which(data_1[,label] == target), receptor]
    
    lig_wilcox <- wilcox.test(lig_0_data, lig_1_data, exact = F, paired = F)
    rec_wilcox <- wilcox.test(rec_0_data, rec_1_data, exact = F, paired = F)
    
    table$lig_diff_p[i] <- lig_wilcox$p.value
    table$rec_diff_p[i] <- rec_wilcox$p.value
  }
  
  table$lig_diff_p_adj <- p.adjust(table$lig_diff_p, method = 'bonferroni')
  table$rec_diff_p_adj <- p.adjust(table$rec_diff_p, method = 'bonferroni')
  
  # If not significant, then 0
  # If significant, then 1
  for (i in 1:length(table[,1])) {
    table$lig_diff[i] <- ifelse(table$lig_diff_p_adj[i] < alpha, 1, 0)
    table$rec_diff[i] <- ifelse(table$rec_diff_p_adj[i] < alpha, 1, 0)
  }
  
  # # If there is difference, then find if increase/decrease
  # # for ligands
  # for (i in 1:length(table[,1])) {
  #   if (table$lig_diff[i] == 1) {
  #     table$lig_diff[i] <- ifelse(table$avg_lig_0[i] > table$avg_lig_1[i], yes = -1, no = 1)
  #   }
  #   
  #   if (table$rec_diff[i] == 1) {
  #     table$rec_diff[i] <- ifelse(table$avg_rec_0[i] > table$avg_rec_1[i], yes = -1, no = 1)
  #   }
  # }
  
  return(table)
}

generate_supplement <- function(LR_table, cytoscape_nodes) {
  'Generate supplemental information for nodes to input back into Cytoscape
    1) ligand/receptor 2) ligand/receptor name 
  
   Args:
    LR_table (data.frame): table with ligand/receptor pairs
    cytoscape_ids (data.frame): exported cytoscape node table
    
   Returns:
    node_names_ids (data.frame): table with gene names matching IDs and if ligand/receptor
  '
  
  # List of node names
  ids <- c(as.character(LR_table$source_ID), as.character(LR_table$target_ID))
  gene_names <- c(as.character(LR_table$ligands), as.character(LR_table$receptors))
  node_names_ids <- data.frame(ID = ids, names = gene_names, LR = rep(c('L', 'R'),each=length(LR_table$source_ID)),stringsAsFactors = F)
  node_names_ids <- unique(node_names_ids)
  
  # sort supplemental names in order of what is in cytoscape
  cytoscape_ids <- gsub('""', '', cytoscape_nodes$`shared name`)
  order <- match(cytoscape_ids, node_names_ids$ID)
  node_names_ids <- node_names_ids[order,]
  node_names_ids$ID <- paste('"', node_names_ids$ID, '"', sep = '')
  
  return(node_names_ids)
}

# AutomatedClusterMarkerTable returns FindAllMarkers table with extra bits of useful information
# and an educated guess about cluster identity

AutomatedClusterMarkerTable <- function(Seurat_Object){
  library(dplyr)
  library(tibble)
  library(Seurat)
  ClusterList <- list()
  Idents(object = Seurat_Object) <- "seurat_clusters"
  current.cluster.ids <- sort(as.numeric(levels(Seurat_Object@active.ident)))
  new.cluster.ids <- c()
  
  
  for(i in current.cluster.ids){
    List_Position <- i + 1
    ClusterList[[List_Position]] <- FindMarkers(object = Seurat_Object, ident.1 = i, min.pct = 0.25, only.pos = TRUE)
    Positive_Genes <- rownames(ClusterList[[List_Position]])
    Num_Positive_Genes <- length(Positive_Genes)
    
    RPS_Num <- length(grep(pattern = "^RPS", x = Positive_Genes))
    RPL_Num <- length(grep(pattern = "^RPL", x = Positive_Genes))
    RP_Percent <- sum(RPS_Num, RPL_Num)/length(Positive_Genes)*100
    RP_Label <- paste("RP%:", RP_Percent, sep = " ")
    
    Mito_Num <- length(grep(pattern = "^MT-", x = Positive_Genes))
    Mito_Percent <- Mito_Num/length(Positive_Genes)*100
    Mito_Label <- paste("Mito%:", RP_Percent, sep = " ")
    
    ClusterCells <- WhichCells(object = Seurat_Object, idents = i)
    Cell_Barcodes <- unlist(Seurat_Object@assays$RNA@counts@Dimnames[2])
    Cell_Number <- c()
    for(k in 1:length(ClusterCells)){
      Cell_Position <- grep(pattern = ClusterCells[k],x = Cell_Barcodes, value = FALSE)
      Cell_Number <- c(Cell_Number,Cell_Position)
    }
    
    
    S_Score <- Seurat_Object@meta.data$S.Score
    G2M_Score <- Seurat_Object@meta.data$G2M.Score
    Cluster_S_Score <- S_Score[Cell_Number]
    Cluster_G2M_Score <- G2M_Score[Cell_Number]
    Avg_Cluster_S_Score <- mean(Cluster_S_Score)
    Avg_Cluster_G2M_Score <- mean(Cluster_G2M_Score)
    Cluster_S_Score_Range <- range(Cluster_S_Score)
    Cluster_G2M_Score_Range <- range(Cluster_G2M_Score)
    
    nFeature <- Seurat_Object@meta.data$nFeature_RNA
    nCount <- Seurat_Object@meta.data$nCount_RNA
    Mito <- Seurat_Object@meta.data$percent.mt
    
    Cluster_nFeature <- nFeature[Cell_Number]
    Cluster_nCount <- nCount[Cell_Number]
    Cluster_Mito <- Mito[Cell_Number]
    
    Avg_Cluster_nFeature <- as.integer(mean(Cluster_nFeature))
    Avg_Cluster_nCount <- as.integer(mean(Cluster_nCount))
    Max_Cluster_Mito <- max(Cluster_Mito)
    
    Cell_Types <- c("Epi","T Cell","Myeloid","B Cell","Fibroblast","RBC","NK", "Endo","Acinar")
    
    Epi_Markers <- c("KRT7","KRT8","KRT18","KRT19","EPCAM","CDH1")
    T_Cell_Markers <- c("CD3E","CD3G","CD3D","CD4","IL7R","CD8A","LEF1")
    Myeloid_Markers <- c("CD14","ITGAM","MNDA","MPEG1","ITGAX")
    B_Cell_Markers <- c("CD79A","MS4A1","CD19")
    Fibroblast_Markers <- c("CDH11","PDGFRA","PDGFRB","ACTA2")
    RBC_Markers <- c("HBA1","HBB","HBA2")
    NK_Markers <- c("NCR3","FCGR3A","NCAM1","KLRF1","KLRC1","CD38","KLRC1")
    Endo_Markers <- c("CDH5","PECAM1")
    Acinar_Markers <- c("TRY4","SPINK1","AMY2A")
    All_Markers <- list(Epi_Markers,T_Cell_Markers,Myeloid_Markers,B_Cell_Markers,Fibroblast_Markers,RBC_Markers,NK_Markers,Endo_Markers,Acinar_Markers)
    
    Epi_Score <- 0
    T_Cell_Score <- 0
    Myeloid_Score <- 0
    B_Cell_Score <- 0
    Fibroblast_Score <- 0
    RBC_Score <- 0
    NK_Score <- 0
    Endo_Score <- 0
    Acinar_Score <- 0 
    All_Scores <- list(Epi_Score,T_Cell_Score,Myeloid_Score,B_Cell_Score,Fibroblast_Score,RBC_Score,NK_Score,Endo_Score,Acinar_Score)
    Weighted_Scores <- c()
    Score_Weights <- c(1.85,1.85,2.22,3.7,2.78,3.7,1.85,5.56,3.7) 
    
    for(h in 1:length(All_Markers)){
      Markers_to_Test<- All_Markers[[h]]
      Marker_Row <- h
      for(j in 1:length(Markers_to_Test)){
        Gene_Found <- 0
        Gene_Found <- length(grep(pattern = Markers_to_Test[j], x = Positive_Genes))
        if(Gene_Found > 0 ){
          All_Scores[[Marker_Row]] <- All_Scores[[Marker_Row]]+1
        }
      }
      Weighted_Scores[Marker_Row] <- All_Scores[[Marker_Row]]*Score_Weights[Marker_Row]
    }
    
    ClusterID <- which(Weighted_Scores >= 5.5)
    if(length(ClusterID) > 0){
      if(length(ClusterID) > 1){
        ID <- "Multiple"
      }else{
        ID <- Cell_Types[ClusterID]
      }
    }else{
      ID <- i
    } 
    if(RP_Percent > 30){
      ID <- paste("RP_",ID,sep = "")
    }
    if(Avg_Cluster_S_Score > 0.01 | Avg_Cluster_G2M_Score > 0.01){
      CellCycleID <- "Cycling"
      ID <- paste("Cycling_",ID,sep = "")
    }else{
      CellCycleID <- "N/A"
    }
    if(Avg_Cluster_nCount < 700){
      ID <- paste("G_",ID,sep = "")
    }
    new.cluster.ids <- c(new.cluster.ids,ID)
    
    Label_Row <- length(Positive_Genes) + 1
    Label_Row2 <- length(Positive_Genes) + 2
    Label_Row3 <- length(Positive_Genes) + 3
    Label_Row4 <- length(Positive_Genes) + 4
    Label_Row5 <- length(Positive_Genes) + 5
    
    Label1 <- c("Summary:",paste("Cluster",i, sep = " "), paste("ID:",ID, sep = " "),paste("Mito%:",Mito_Percent, sep = " "),paste("RP%:",RP_Percent, sep = " "))
    Label2 <- c("Immune Summary",paste("T Cell Score:",All_Scores[[2]], sep = " "),paste("Myeloid Score:",All_Scores[[3]], sep = " "),paste("B Cell Score:",All_Scores[[4]], sep = " "),
                paste("NK Score:",All_Scores[[7]], sep = " "))
    Label3 <- c(paste("Epi Score:",All_Scores[[1]], sep = " "),paste("Fib Score:",All_Scores[[5]], sep = " "),paste("Acinar Score:",All_Scores[[9]], sep = " "),
                paste("Endo Score:",All_Scores[[8]], sep = " "),paste("RBC Score:",All_Scores[[6]], sep = " "))
    Label4 <- c("Avg S Score:", Avg_Cluster_S_Score, "Avg G2M Score:", Avg_Cluster_G2M_Score, CellCycleID)
    Label5 <- c("Filter Info",paste("Avg. nGene:", Avg_Cluster_nFeature, sep = " "),paste("Avg. nCounts:", Avg_Cluster_nCount, sep = " ")
                , paste("Highest Mito:",Max_Cluster_Mito, sep = " "), paste("# Cells:",length(Cell_Number), sep = " "))
    
    ClusterList[[List_Position]][Label_Row,] <- Label1
    ClusterList[[List_Position]][Label_Row2,] <- Label2
    ClusterList[[List_Position]][Label_Row3,] <- Label3
    ClusterList[[List_Position]][Label_Row4,] <- Label4
    ClusterList[[List_Position]][Label_Row5,] <- Label5
    
    ClusterList[[List_Position]] <- rownames_to_column(.data = ClusterList[[List_Position]],var = "Gene")
    ClusterList[[List_Position]][Label_Row,"Gene"] <- "Summary1"
    ClusterList[[List_Position]][Label_Row2,"Gene"] <- "Summary2"
    ClusterList[[List_Position]][Label_Row3,"Gene"] <- "Summary3"
    ClusterList[[List_Position]][Label_Row4,"Gene"] <- "Summary4"
    ClusterList[[List_Position]][Label_Row5,"Gene"] <- "Summary5"
    
    
  }
  
  ClusterDataFrame <- bind_rows(ClusterList, .id = "column_label")
  ClusterDataFrame <- ClusterDataFrame[,-1]
  ClusterPackage <- list(ClusterDataFrame, new.cluster.ids)
  return(ClusterPackage)
}

#Circos Code

# CircosFunctions will output the text files required to run Circos
#SetIdents to Circos Labels Prior to Starting
CircosFunctions <- function(InteractomeData, SeuratObject, CellTypeVector, POI, Lig_or_Rec, LR_List, Species, Cutoff_Lig, Cutoff_Rec){
  #InteractomeData: Filtered Interactome Data
  #SeuratObject: Object used to make interactome 
  #CellTypeVector: Vector of cell types to be included in interactome
  #POI: A number representing the cell population of interest to single out as the Ligand in a receptor plot or the receptor in a ligand plot
  #Lig_or_Rec: T = POI as Ligand in receptor plot, F = POI as receptor in ligand plot
  #LR_List: A dataframe of ligands and paired receptors, ex. The Ramilowski List (organized human, mouse, suborganized ligand-receptor)
  #Species: T = Human, F = Mouse
  #Cutoff_Lig: Cutoff ligand expression value for the interactome
  #Cutoff_Rec: Cutoff receptor expression value for the interactome
  #Cutoff Values:
  #Positive = Cutoff determined by summary function ex Min, 1st Q, Median, Mean, 3rd Q, Max (1-6)
  #Negative = Arbitrary cutoff (set -0.05 for a 0.05 cutoff)
  
  #POI Numbers:
  #1 - CD8 T Cells
  #2 - CD4 T Cells
  #3 - T Cells
  #4 - myCAF Fibroblasts
  #5 - iCAF Fibroblasts
  #6 - Fibroblasts
  #7 - Epithelial
  #8 - Acinar
  #9 - Endothelial
  #10 - Mast Cells
  #11 - Granulocytes
  #12 - Macrophages
  #13 - MDSCs 
  #14 - Myeloid
  #15 - NK Cells
  #16 - Dendritic Cells
  #17 - Endocrine
  #18 - Perivascular
  #19 - B Cells
  
  CircosFiles <- list()
  
  CellTypes <- c("CD8TCells","CD4TCells","TCells","myCAFFibroblast","iCAFFibroblast","Fibroblasts","Epithelial","Acinar","Endothelial","MastCells",
                 "Granulocytes","Macrophages","MDSCs","Myeloid","NKCells","DendriticCells","Endocrine","Perivascular","BCells")
  
  CellTypeColors <- c("darkgreen", "limegreen","forestgreen","darkslategray3","darkcyan","darkcyan","red3","deeppink",
                      "mediumvioletred","gold","orange1","darkorange2", "tan2","darkorange2","darkorchid","chocolate4","darkred","lightskyblue","gold1")
  
  #Karyotype Function
  
  Karyotype_File <- data.frame(
    "Chr" = "chr -",
    "Name" = "Epithelial",
    "Label" = "Epithelial",
    "Start" = 0,
    "End" = 0,
    "Color" = "chr0"
  )
  
  LigList_File <- data.frame(
    "Chr" = "A",
    "Start" = 0,
    "End" = 0,
    "Name" = "A"
  )
  
  RecList_File <- data.frame(
    "Chr" = "A",
    "Start" = 0,
    "End" = 0,
    "Name" = "A"
  )
  
  DataVec_Source <- c()
  
  for(i in 1:length(CellTypeVector)){
    CellType <- CellTypeVector[i]
    DataVec_Source <- c(DataVec_Source, which(x = InteractomeData$source == CellType))
  }
  
  Data_Source <- InteractomeData[DataVec_Source,]
  
  DataVec_Target <- c()
  
  for(i in 1:length(CellTypeVector)){
    CellType <- CellTypeVector[i]
    DataVec_Target <- c(DataVec_Target, which(x = Data_Source$target == CellType))
  }
  
  Data_Target <- Data_Source[DataVec_Target,]
  
  Data <- Data_Target[which(Data_Target$source != Data_Target$target),]
  
  if(Cutoff_Lig > 0 & Cutoff_Rec > 0){
    Data <- Data[which(Data$avg_lig_0 > summary(Data$avg_lig_0)[1] | Data$avg_rec_0 > summary(Data$avg_rec_0)[Cutoff_Rec]),]
  }else{
    if(Cutoff_Lig < 0 & Cutoff_Lig < 0){
      Data <- Data[which(Data$avg_lig_0 > -1*-.01),]
      Data <- Data[which(Data$avg_rec_0 > -1*Cutoff_Rec),]
    }else{
      if(Cutoff_Lig > 0){
        Data <- Data[which(Data$avg_lig_0 > summary(Data$avg_lig_0)[Cutoff_Lig]),]
        Data <- Data[which(Data$avg_rec_0 > -1*Cutoff_Rec),]
      }else{
        Data <- Data[which(Data$avg_rec_0 > summary(Data$avg_rec_0)[Cutoff_Rec]),]
        Data <- Data[which(Data$avg_lig_0 > -1*Cutoff_Lig),]
      }
    }
  }
  
  for (i in 1:length(table(Data$source))) {
    #Determine if population is present after filtering by significance
    if( length(which(Data$source == rownames(table(Data$source))[i])) > 0 ){
      #Subset Population of Interest
      POI_Lig_Data <- Data[which(Data$source == rownames(table(Data$source))[i]),]
      POI_Rec_Data <- Data[which(Data$target == rownames(table(Data$source))[i]),]
      POI_Lig <- unique(as.character(POI_Lig_Data$ligands))
      POI_Rec <-  unique(as.character(POI_Rec_Data$receptors))
      POI_Lig_Num <- length(POI_Lig)
      POI_Rec_Num <- length(POI_Rec)
      if(POI_Lig_Num > POI_Rec_Num){
        BandSize <- POI_Lig_Num*10-1
      }else{
        BandSize <- POI_Rec_Num*10-1
      }
      InsertDF <- data.frame(
        "Chr" = "chr -",
        "Name" = as.character(rownames(table(Data$source))[i]),
        "Label" = as.character(rownames(table(Data$source))[i]),
        "Start" = 0,
        "End" = BandSize,
        "Color" = paste("chr",i,sep = "")
      ) 
      Karyotype_File <- merge(x = Karyotype_File, y = InsertDF, by = c("Chr","Name","Label","Start","End","Color"), all.y = T, all.x = T, sort = F)  
      
      if(POI_Lig_Num > POI_Rec_Num){
        
        InsertDF_LigList <- data.frame(
          "Chr" = rep(x = as.character(rownames(table(Data$source))[i]), times = POI_Lig_Num),
          "Start" = seq(from = 0, to = POI_Lig_Num*10-10, by = 10),
          "End" = seq(from = 9, to = POI_Lig_Num*10-1, by = 10),
          "Name" = POI_Lig
        ) 
        
        InsertDF_RecList <- data.frame(
          "Chr" = rep(x = as.character(rownames(table(Data$target))[i]), times = POI_Rec_Num),
          "Start" = ceiling(seq(from = 0, to = POI_Lig_Num*10-(POI_Lig_Num*10-10-0)/POI_Rec_Num, by = (POI_Lig_Num*10-1)/POI_Rec_Num)),
          "End" = floor(seq(from = (POI_Lig_Num*10-1)/POI_Rec_Num, to = POI_Lig_Num*10-1, by = (POI_Lig_Num*10-1)/POI_Rec_Num)),
          "Name" = POI_Rec
        ) 
        
      }else{
        
        InsertDF_LigList <- data.frame(
          "Chr" = rep(x = as.character(rownames(table(Data$source))[i]), times = POI_Lig_Num),
          "Start" = ceiling(seq(from = 0, to = POI_Rec_Num*10-(POI_Rec_Num*10-10-0)/POI_Lig_Num, by = (POI_Rec_Num*10-1)/POI_Lig_Num)),
          "End" = floor(seq(from = (POI_Rec_Num*10-1)/POI_Lig_Num, to = POI_Rec_Num*10-1, by = (POI_Rec_Num*10-1)/POI_Lig_Num)),
          "Name" = POI_Lig
        ) 
        
        InsertDF_RecList <- data.frame(
          "Chr" = rep(x = as.character(rownames(table(Data$target))[i]), times = POI_Rec_Num),
          "Start" = seq(from = 0, to = POI_Rec_Num*10-10, by = 10),
          "End" = seq(from = 9, to = POI_Rec_Num*10-1, by = 10),
          "Name" = POI_Rec
        ) 
      }
      
      Karyotype_File <- merge(x = Karyotype_File, y = InsertDF, by = c("Chr","Name","Label","Start","End","Color"), all.y = T, all.x = T)
      LigList_File <- merge(x = LigList_File, y = InsertDF_LigList, by = c("Chr","Start","End","Name"), all.y = T, all.x = T, sort = F)
      RecList_File <- merge(x = RecList_File, y = InsertDF_RecList, by = c("Chr","Start","End","Name"), all.y = T, all.x = T, sort = F)
      
    } 
  }
  
  Karyotype_File <- Karyotype_File[-1,]
  LigList_File <- LigList_File[-1,]
  RecList_File <- RecList_File[-1,]
  
  sentenceString <- Karyotype_File$Name
  searchString <- ' '
  replacementString <- ''
  sentenceString = sub(searchString,replacementString,sentenceString)
  sentenceString = sub(searchString,replacementString,sentenceString)
  Karyotype_File$Name <- sentenceString
  Karyotype_File$Label <- sentenceString
  
  sentenceString <- LigList_File$Chr
  searchString <- ' '
  replacementString <- ''
  sentenceString = sub(searchString,replacementString,sentenceString)
  sentenceString = sub(searchString,replacementString,sentenceString)
  LigList_File$Chr <- sentenceString
  
  
  sentenceString <- RecList_File$Chr
  searchString <- ' '
  replacementString <- ''
  sentenceString = sub(searchString,replacementString,sentenceString)
  sentenceString = sub(searchString,replacementString,sentenceString)
  RecList_File$Chr <- sentenceString
  
  Band_Color <- c()
  
  for (i in 1:dim(Karyotype_File)[1]) {
    CellType <- as.character(Karyotype_File[i,"Name"])
    Band_Color <- c(Band_Color, CellTypeColors[which(CellTypes == CellType)])
  }
  Karyotype_File[,"Color"] <- Band_Color
  
  LigList <- LigList_File
  RecList <- RecList_File
  
  if(Species == T){
    LigRecList <- LR_List[,1:2]
  }else{
    LigRecList <- LR_List[,3:4]
  }
  
  
  
  
  #Expression Function
  
  LigExpression <- data.frame(
    "Chr"   = "A",
    "Name" = "A",
    "Expression" = 0
  )
  
  RecExpression <- data.frame(
    "Chr"   = "A",
    "Name" = "A",
    "Expression" = 0
  )
  
  Avg_Expression <- AverageExpression(object = SeuratObject, features = c(LigList$Name,RecList$Name))
  Avg_Expression <- Avg_Expression[[1]]
  
  for(i in 1:length(CellTypeVector)){
    
    CellSearch <- CellTypeVector[i]
    sentenceString <- CellSearch
    searchString <- ' '
    replacementString <- ''
    sentenceString = sub(searchString,replacementString,sentenceString)
    sentenceString = sub(searchString,replacementString,sentenceString)
    CellSearch <- sentenceString
    
    ident1 <- CellTypeVector[i]
    CellExpression <- subset.data.frame(x = Avg_Expression, select = ident1)
    
    LigFeatures <- as.character(LigList[which(LigList$Chr == CellSearch),"Name"])
    RecFeatures <- as.character(RecList[which(RecList$Chr == CellSearch),"Name"])
    
    Lig_DE <- data.frame(
      "Exp" = 0,
      "Name" = "a"
    )
    Rec_DE <- data.frame(
      "Exp" = 0,
      "Name" = "a"
    )
    for(k in 1:length(LigFeatures)){
      
      Lig_Exp <- subset.data.frame(x = CellExpression, subset = rownames(Avg_Expression) == LigFeatures[k])
      Lig_Exp[1,2] <- LigFeatures[k]
      colnames(Lig_Exp)[1] <- "Exp"
      colnames(Lig_Exp)[2] <- "Name"
      
      Lig_DE <- merge(x = Lig_DE, y = Lig_Exp, by = c("Exp","Name"), all.y = T, all.x = T, sort = F)
      
    }
    
    Lig_DE <- Lig_DE[-1,]
    
    for(t in 1:length(RecFeatures)){
      
      Rec_Exp <- subset.data.frame(x = CellExpression, subset = rownames(Avg_Expression) == RecFeatures[t])
      Rec_Exp[1,2] <- RecFeatures[t]
      colnames(Rec_Exp)[1] <- "Exp"
      colnames(Rec_Exp)[2] <- "Name"
      
      Rec_DE <- merge(x = Rec_DE, y = Rec_Exp, by = c("Exp","Name"), all.y = T, all.x = T, sort = F)
    }
    
    Rec_DE <- Rec_DE[-1,] 
    
    Lig_DF <- data.frame(
      "Chr" = rep(x = CellTypeVector[i], times = dim(Lig_DE)[1]),
      "Name" = Lig_DE$Name,
      "Expression" = Lig_DE$Exp
    )
    Rec_DF <- data.frame(
      "Chr" = rep(x = CellTypeVector[i], times = dim(Rec_DE)[1]),
      "Name" = Rec_DE$Name,
      "Expression" = Rec_DE$Exp
    )
    
    LigExpression <- merge(x = LigExpression, y = Lig_DF, by = c("Chr","Name","Expression"), all.y = T, all.x = T, sort = F)
    RecExpression <- merge(x = RecExpression, y = Rec_DF, by = c("Chr","Name","Expression"), all.y = T, all.x = T, sort = F)
    
    
  }
  
  LigExpression <- LigExpression[-1,]
  RecExpression <- RecExpression[-1,]
  
  sentenceString <- LigExpression$Chr
  searchString <- ' '
  replacementString <- ''
  sentenceString = sub(searchString,replacementString,sentenceString)
  sentenceString = sub(searchString,replacementString,sentenceString)
  LigExpression$Chr <- sentenceString
  
  
  sentenceString <- RecExpression$Chr
  searchString <- ' '
  replacementString <- ''
  sentenceString = sub(searchString,replacementString,sentenceString)
  sentenceString = sub(searchString,replacementString,sentenceString)
  RecExpression$Chr <- sentenceString
  
  ExpressionOutput <- merge(x = LigExpression, y = RecExpression, by = c("Chr","Name","Expression"), all.y = T, all.x = T, sort = F)
  
  #Squeeze expression values
  ExpressionInput <- ExpressionOutput
  Normalized_Expression <- log(ExpressionOutput$Expression*1000)
  ScaledExpression <- ((Normalized_Expression - range(Normalized_Expression)[1])/
                         (range(Normalized_Expression)[2]- range(Normalized_Expression)[1]))*4
  ExpressionInput$Expression <- ScaledExpression
  
  #Text and Link Function
  
  CellNumCounter <- 1
  
  if(Lig_or_Rec == T){
    
    if(POI == 1){
      CellType <- CellTypes[1]
      CellTypeColor <- CellTypeColors[1]
    }else{
      while (POI != CellNumCounter) {
        CellNumCounter <- CellNumCounter + 1
      }
      CellType <- CellTypes[CellNumCounter]
      CellTypeColor <- CellTypeColors[CellNumCounter]
    }
    
    POI_Lig_Data <- LigList[which(LigList$Chr == CellType),]
    POI_Rec_Data <- RecList[which(RecList$Chr != CellType),]
    Text_File <- merge(x = POI_Lig_Data, y = POI_Rec_Data, by = c("Chr","Start","End","Name"), all.y = T, all.x = T, sort = F)
    
    Fil_LigRecList_Vec <- c()
    
    for(i in 1:dim(LigRecList)[1]){
      LigRecList_Lig <- as.character(LigRecList[i,1])
      LigRecList_Rec <- as.character(LigRecList[i,2])
      Lig_Present <- grep(pattern = paste("^", LigRecList_Lig,"$", sep = ""), x = POI_Lig_Data$Name)
      Rec_Present <- grep(pattern = paste("^", LigRecList_Rec,"$", sep = ""), x = POI_Rec_Data$Name)
      LigList_Row <- which(LigRecList[,1] == LigRecList_Lig)
      RecList_Row <- which(LigRecList[,2] == LigRecList_Rec)
      if(length(Lig_Present) > 0 & length(Rec_Present) > 0){
        Fil_LigRecList_Vec <- c(Fil_LigRecList_Vec, intersect(LigList_Row,RecList_Row))
      }
    }
    
    LigRecList_Fil <- LigRecList[Fil_LigRecList_Vec,]
    
    Link_File <- data.frame(
      "Chr" = "A",
      "Start" = 0,
      "End" =  0,
      "Chr1" = "A",
      "Start1" = 0,
      "End1" = 0
    )
    
    for (i in 1:length(POI_Lig_Data$Name)) {
      Rec_Data <- data.frame(
        "Chr" = "A",
        "Start" = 0,
        "End" = 0
      )
      
      Lig <- as.character(POI_Lig_Data$Name[i])
      Receptors <- as.character(LigRecList_Fil[,2][grep(x = LigRecList_Fil[,1], pattern = paste("^", Lig,"$", sep = ""),value = F)])
      
      if(length(Receptors) == 1){
        Rec_Data_Insert <- POI_Rec_Data[grep(x = POI_Rec_Data$Name, pattern = paste("^", Receptors,"$", sep = ""),value = F),1:3]
        Rec_Data <- merge(x = Rec_Data, y = Rec_Data_Insert, by = c("Chr","Start","End"), all.y = T, all.x = T, sort = F)
      }else{
        for (l in 1:length(Receptors)) {
          Receptor <- Receptors[l]
          Rec_Data_Insert <- POI_Rec_Data[grep(x = POI_Rec_Data$Name, pattern = paste("^", Receptor,"$", sep = ""),value = F),1:3]
          Rec_Data <- merge(x = Rec_Data, y = Rec_Data_Insert, by = c("Chr","Start","End"), all.y = T, all.x = T, sort = F)
        }
      }
      Rec_Data <- Rec_Data[-1,]
      
      Lig_Data <- POI_Lig_Data[i,1:3]
      
      if(length(Rec_Data$Chr) > 1){
        Lig_Data <- rbind(Lig_Data, Lig_Data[rep(1, length(Rec_Data$Chr)-1), ])
      }
      
      Lig_Rec_Merge <- bind_cols(x = Lig_Data, y = Rec_Data)
      
      Link_File <- rbind(x = Link_File, y = Lig_Rec_Merge)
    }
    
    
    
    Link_File <- Link_File[-1,]
    Link_File[,(dim(Link_File)[2]+1)] <-rep(x = paste("color=",CellTypeColor,sep = ""), times = dim(Link_File)[1])
    
    Expression_File <- Text_File
    Expression_File[,5] <- rep(0, times = dim(Expression_File)[1])
    
    for(k in 1:dim(Expression_File)[1]){
      GeneCell <- Expression_File[k,"Chr"]
      Gene <- as.character(Expression_File[k,"Name"])
      ExpressionPull <- ExpressionInput[which(ExpressionInput$Chr == GeneCell & ExpressionInput$Name == Gene),]
      Expression_File[k,5] <- ExpressionPull$Expression
    }
    
    Expression_File <- Expression_File[,-4]
    
  }else{
    
    if(POI == 1){
      CellType <- CellTypes[1]
    }else{
      while (POI != CellNumCounter) {
        CellNumCounter <- CellNumCounter + 1
      }
      CellType <- CellTypes[CellNumCounter]
      CellTypeColor <- CellTypeColors[CellNumCounter]
    }
    
    POI_Lig_Data <- LigList[which(LigList$Chr != CellType),]
    POI_Rec_Data <- RecList[which(RecList$Chr == CellType),]
    Text_File <- merge(x = POI_Rec_Data, y = POI_Lig_Data, by = c("Chr","Start","End","Name"), all.y = T, all.x = T, sort = F)
    
    Fil_LigRecList_Vec <- c()
    
    for(i in 1:dim(LigRecList)[1]){
      LigRecList_Lig <- as.character(LigRecList[i,1])
      LigRecList_Rec <- as.character(LigRecList[i,2])
      Lig_Present <- grep(pattern = paste("^", LigRecList_Lig,"$", sep = ""), x = POI_Lig_Data$Name)
      Rec_Present <- grep(pattern = paste("^", LigRecList_Rec,"$", sep = ""), x = POI_Rec_Data$Name)
      LigList_Row <- which(LigRecList[,1] == LigRecList_Lig)
      RecList_Row <- which(LigRecList[,2] == LigRecList_Rec)
      if(length(Lig_Present) > 0 & length(Rec_Present) > 0){
        Fil_LigRecList_Vec <- c(Fil_LigRecList_Vec, intersect(LigList_Row,RecList_Row))
      }
    }
    
    LigRecList_Fil <- LigRecList[Fil_LigRecList_Vec,]    
    
    Link_File <- data.frame(
      "Chr" = "A",
      "Start" = 0,
      "End" =  0,
      "Chr1" = "A",
      "Start1" = 0,
      "End1" = 0
    )
    
    
    
    for (i in 1:length(POI_Rec_Data$Name)) {
      Lig_Data <- data.frame(
        "Chr" = "A",
        "Start" = 0,
        "End" = 0
      )
      
      Rec <- as.character(POI_Rec_Data$Name[i])
      Ligands <- as.character(LigRecList_Fil[,1][grep(x = LigRecList_Fil[,2], pattern = paste("^", Rec,"$", sep = ""),value = F)])
      
      if(length(Ligands == 1)){
        Lig_Data_Insert <- POI_Lig_Data[grep(x = POI_Lig_Data$Name, pattern = paste("^", Ligands,"$", sep = ""),value = F),1:3]
        Lig_Data <- merge(x = Lig_Data, y = Lig_Data_Insert, by = c("Chr","Start","End"), all.y = T, all.x = T, sort = F)
      }else{
        for (l in 1:length(Ligands)) {
          Ligand <- Ligands[l]
          Lig_Data_Insert <- POI_Lig_Data[grep(x = POI_Lig_Data$Name, pattern = paste("^", Ligand,"$", sep = ""),value = F),1:3]
          Lig_Data <- merge(x = Lig_Data, y = Lig_Data_Insert, by = c("Chr","Start","End"), all.y = T, all.x = T, sort = F)
        }
      }
      Lig_Data <- Lig_Data[-1,]
      
      Rec_Data <- POI_Rec_Data[i,1:3]
      
      if(length(Lig_Data$Chr) > 1){
        Rec_Data <- rbind(Rec_Data, Rec_Data[rep(1, length(Lig_Data$Chr)-1), ])
      }
      Lig_Rec_Merge <- bind_cols(x = Rec_Data, y = Lig_Data)
      
      Link_File <- rbind(x = Link_File, y = Lig_Rec_Merge)
    }
    
    
    
    Link_File <- Link_File[-1,]
    
    Link_Color <- c()
    
    for (i in 1:dim(Link_File)[1]) {
      CellType <- as.character(Link_File[i,"Chr1"])
      Link_Color <- c(Link_Color, CellTypeColors[which(CellTypes == CellType)])
    }
    Link_Color <- paste("color=",Link_Color,"_a4",sep = "")
    Link_File[,(dim(Link_File)[2]+1)] <- Link_Color
    
    Expression_File <- Text_File
    Expression_File[,5] <- rep(0, times = dim(Expression_File)[1])
    
    for(k in 1:dim(Expression_File)[1]){
      GeneCell <- Expression_File[k,"Chr"]
      Gene <- as.character(Expression_File[k,"Name"])
      ExpressionPull <- ExpressionInput[which(ExpressionInput$Chr == GeneCell & ExpressionInput$Name == Gene),]
      Expression_File[k,5] <- ExpressionPull$Expression
    }
    
    Expression_File <- Expression_File[,-4]
    
  }
  
  CircosFiles[[1]] <- Karyotype_File
  CircosFiles[[2]] <- Text_File
  CircosFiles[[3]] <- Expression_File
  CircosFiles[[4]] <- Link_File 
  
  return(CircosFiles)
}


#____________________________________________________________________________________________________________________________________________________________________#

# 2. Data Pre-processing

#Tissue Object

#Load in all filtered runs using Read10X or Read10X_h5 functions
#   Data.data <- Read10X("filepath/")/Read10X_h5("filepath/file.ext")

PDAC_TISSUE_1.data <- Read10X("~/Desktop/PDAC_TISSUE_1/filtered_feature_bc_matrix/")
PDAC_TISSUE_14.data <- Read10X_h5("~/Desktop/PDAC_TISSUE_14/filtered_feature_bc_matrix.h5")

#Generate Seurat objects using CreateSeuratObject function, min.cells=3, min.features=100
# Data <- CreateSeuratObject(counts = Data.data, project = "PDA", min.cells=3, min.features=100)

PDAC_TISSUE_1 <- CreateSeuratObject(counts = PDAC_TISSUE_1.data, project = 'PDAC_TISSUE_1', min.cells = 3, min.features = 100)

#Add Metadata

#ID Metadata
PDAC_TISSUE_16$ID <- "PDAC_TISSUE_16"
PDAC_TISSUE_15$ID <- "PDAC_TISSUE_15"
PDAC_TISSUE_13$ID <- "PDAC_TISSUE_13"
PDAC_TISSUE_14$ID <- "PDAC_TISSUE_14"
PDAC_TISSUE_12$ID <- "PDAC_TISSUE_12"
PDAC_TISSUE_11B$ID <- "PDAC_TISSUE_11B" 
PDAC_TISSUE_11A$ID <- "PDAC_TISSUE_11A" 
PDAC_TISSUE_10$ID <- "PDAC_TISSUE_10"
PDAC_TISSUE_9$ID <- "PDAC_TISSUE_9"
PDAC_TISSUE_8$ID <- "PDAC_TISSUE_8"
PDAC_TISSUE_7$ID <- "PDAC_TISSUE_7"
PDAC_TISSUE_6$ID <- "PDAC_TISSUE_6"
PDAC_TISSUE_5$ID <- "PDAC_TISSUE_5"
PDAC_TISSUE_4$ID <- "PDAC_TISSUE_4"
PDAC_TISSUE_3$ID <- "PDAC_TISSUE_3"
PDAC_TISSUE_2$ID <- "PDAC_TISSUE_2"
PDAC_TISSUE_1$ID <- "PDAC_TISSUE_1"
AdjNorm_TISSUE_1$ID <-"AdjNorm_TISSUE_1"
AdjNorm_TISSUE_3$ID <- "AdjNorm_TISSUE_3"
AdjNorm_TISSUE_2$ID <- "AdjNorm_TISSUE_2"

#Disease State Metadata
PDAC_TISSUE_16$DiseaseState <-"PDAC"
PDAC_TISSUE_15$DiseaseState <-"PDAC"
PDAC_TISSUE_13$DiseaseState <-"PDAC"
PDAC_TISSUE_14$DiseaseState <-"PDAC"
PDAC_TISSUE_12$DiseaseState <-"PDAC"
PDAC_TISSUE_11B$DiseaseState <-"PDAC"  
PDAC_TISSUE_11A$DiseaseState <-"PDAC"
PDAC_TISSUE_10$DiseaseState <-"PDAC"
PDAC_TISSUE_9$DiseaseState <-"PDAC"
PDAC_TISSUE_8$DiseaseState <-"PDAC"
PDAC_TISSUE_7$DiseaseState <-"PDAC"
PDAC_TISSUE_6$DiseaseState <-"PDAC"
PDAC_TISSUE_5$DiseaseState <-"PDAC"
PDAC_TISSUE_4$DiseaseState <-"PDAC"
PDAC_TISSUE_3$DiseaseState <-"PDAC"
PDAC_TISSUE_2$DiseaseState <-"PDAC"
PDAC_TISSUE_1$DiseaseState <-"PDAC"
AdjNorm_TISSUE_1$DiseaseState <- "AdjNorm"
AdjNorm_TISSUE_3$DiseaseState <- "AdjNorm"
AdjNorm_TISSUE_2$DiseaseState <- "AdjNorm"

#Merge all objects
TotalTissue.combined <- merge(PDAC_TISSUE_1, y = c(PDAC_TISSUE_2, PDAC_TISSUE_3, PDAC_TISSUE_4, 
                                                   PDAC_TISSUE_5, PDAC_TISSUE_6, PDAC_TISSUE_7, 
                                                   PDAC_TISSUE_8, PDAC_TISSUE_9, PDAC_TISSUE_10, 
                                                   PDAC_TISSUE_11A, PDAC_TISSUE_11B, PDAC_TISSUE_12,
                                                   PDAC_TISSUE_13, PDAC_TISSUE_14, PDAC_TISSUE_15, 
                                                   PDAC_TISSUE_16))
#Check all objects present
table(TotalTissue.combined$orig.ident)

# Changing between meta.data for identities- you can change this by altering what you input into your metadata
Idents(object = TotalTissue.combined) <- 'ID'
# Check active identity
levels(TotalTissue.combined)

#NORMALIZE DATA
TotalTissue.combined <- NormalizeData(object = TotalTissue.combined, normalization.method = "LogNormalize", 
                                      scale.factor = 10000)
#Percent Mitochondrial Genes
#QC Metric used to remove cells with overabundant Mitochondrial genes, typically associated with nuclear wash out during sequencing
TotalTissue.combined[["percent.mt"]] <- PercentageFeatureSet(TotalTissue.combined, pattern = "^MT-")

#Plot the nFeatures/counts/% Mito to get general idea about the quality of your data
VlnPlot(TotalTissue.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = .0)

#FIND VARIABLE GENES
TotalTissue.combined<- FindVariableFeatures(object = TotalTissue.combined, mean.function = ExpMean, dispersion.function = LogVMR, 
                                            x.low.cutoff = 0.0125, y.cutoff = 0.5)


#Calculate Cell Cycle Score (S-G2M Difference)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
TotalTissue.combined<- CellCycleScoring(TotalTissue.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
TotalTissue.combined$CC.Difference <- TotalTissue.combined$S.Score - TotalTissue.combined$G2M.Score

#THIS STEP MAY TAKE A VERY LONG TIME
#Scale Data
TotalTissue.combined<- ScaleData(object = TotalTissue.combined, vars.to.regress = "nCount_RNA", features = rownames(TotalTissue.combined))

#It is recommended that you save this object so you do not have to re-run the scaling step more than necessary
save(TotalTissue.combined,file="TotalTissue.combined.RData")
#____________________________________________________________________________________________________________________________________________________________________#

#Data Visualization

#Run PCA and Determine Dimensions for 90% Variance
TotalTissue.combined <- RunPCA(object = TotalTissue.combined, features = VariableFeatures(object = TotalTissue.combined))
stdev <- TotalTissue.combined@reductions$pca@stdev
var <- stdev^2

EndVar = 0

for(i in 1:length(var)){
  total <- sum(var)
  numerator <- sum(var[1:i])
  expvar <- numerator/total
  if(End == 0){
    if(expvar > 0.9){
      EndVar <- EndVar + 1
      PCNum <- i
    }
  }
}
#Confirm #PC's determined explain > 90% of variance
sum(var[1:PCNum])/ sum(var)


#Find Neighbors + Find CLusters (without harmony batch correction)
TotalTissue.combined <- FindNeighbors(object = TotalTissue.combined, dims = 1:PCNum)
TotalTissue.combined <- FindClusters(object = TotalTissue.combined, resolution = 1.2)

#Run UMAP and get unlabelled cluster UMAP and violin plot (without harmony batch correction)
TotalTissue.combined <- RunUMAP(object = TotalTissue.combined, dims = 1:PCNum)
DimPlot(object = TotalTissue.combined, reduction = "umap", label = TRUE, pt.size = 0.5)
TotalTissue.combined[["UMAP_Clusters"]] <- Idents(object = TotalTissue.combined)

#UMAP split by group
DimPlot(object = TotalTissue.combined, reduction = "umap", label = TRUE, pt.size = 0.5, split.by = "ID", ncol=5)

# Changing between meta.data for identities- you can change this by altering what you input into your metadata - will need to do this to make the right plots
Idents(object = TotalTissue.combined) <- 'UMAP_Clusters'
# Check active identity
levels(TotalTissue.combined)

#Umap with overlaid group
Idents(object = TotalTissue.combined) <- 'ID'
DimPlot(object = TotalTissue.combined, reduction = "umap", label = TRUE, pt.size = 0.5)
VlnPlot(TotalTissue.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Find neighbors and clusters WITH harmony batch correction
options(repr.plot.height = 2.5, repr.plot.width = 6)
TotalTissue.combined <- TotalTissue.combined %>% 
  RunHarmony("ID", plot_convergence = TRUE)

TotalTissue.combined.harmony <- FindNeighbors(object = TotalTissue.combined, dims = 1:PCNum, reduction ="harmony")
TotalTissue.combined.harmony <- FindClusters(object = TotalTissue.combined, resolution = 1.2, reduction ="harmony")

#Run UMAP and get unlabelled cluster UMAP and violin plot
TotalTissue.combined.harmony <- RunUMAP(object = TotalTissue.combined.harmony, dims = 1:PCNum, reduction = "harmony")
DimPlot(object = TotalTissue.combined.harmony, reduction = "umap", label = TRUE, pt.size = 0.5)

#UMAP split by group
DimPlot(object = TotalTissue.combined.harmony, reduction = "umap", label = TRUE, pt.size = 0.5, split.by = "ID", ncol=5)

# Changing between meta.data for identities- you can change this by altering what you input into your metadata - will need to do this to make the right plots
Idents(object = TotalTissue.combined.harmony) <- 'group'
# Check active identity
levels(TotalTissue.combined.harmony)

#Umap with overlaid group
Idents(object = TotalTissue.combined.harmony) <- 'ID'
DimPlot(object = TotalTissue.combined.harmony, reduction = "umap", label = TRUE, pt.size = 0.5)
VlnPlot(TotalTissue.combined.harmony, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Make the active identity UMAP clusters again
Idents(object = TotalTissue.combined.harmony) <- 'UMAP_Clusters'
# Check active identity
levels(TotalTissue.combined.harmony)

#MAKE CLUSTER MARKER TABLE

# This function automates the FindMarkers function and uses the list of markers to broadly
# identify cell types based on a preselected list of markers. Markers chosen for human samples.
# The output is a dataframe containing the FindMarkers output.

ClusterMarkerTable <- function(ClustNum,Data,NumMarkers,PCT) {
  ClusterMarks <- FindMarkers(object = Data, ident.1 = ClustNum[1], min.pct = PCT)
  ClusterTable <- as.data.frame(head(x = ClusterMarks, n = NumMarkers))
  MarkerLabelRow <- length(ClusterTable$p_val)+1
  MarkerLabel <- c("Cluster", ClustNum[1], "Markers", "Are", "Above")
  ClusterTable[MarkerLabelRow,] <- MarkerLabel
  
  if(length(ClustNum) > 2) {
    for(i in ClustNum[2:length(ClustNum)]){
      ClusterMarks <- FindMarkers(object = Data, ident.1 = i, min.pct = PCT)
      ClusterTable <- rbind(ClusterTable,head(x = ClusterMarks, n = NumMarkers))
      MarkerLabelRow <- length(ClusterTable$p_val)+1
      MarkerLabel <- c("Cluster", i, "Markers", "Are", "Above")
      ClusterTable[MarkerLabelRow,] <- MarkerLabel
    }
  }
  View(ClusterTable)
  return(ClusterTable)
}

#change cluster marker range to match desired cluster numbers from UMAP
TotalTissue.combined.harmony_clusters <- ClusterMarkerTable(c(44:50),TotalTissue.combined.harmony,100,0.25) 
write.csv(TotalTissue.combined.harmony_clusters, file = "~/Desktop/TotalTissue.combined.harmony_clusters101719clusters44-50.csv")

#Save the batch corrected object
save(TotalTissue.combined.harmony, file = 'TotalTissue.combined.harmony_scaled.RData')

#Define Clusters Based on Marker Expression

#Epithelial
FeaturePlot(object = TotalTissue.combined.harmony, features = c('PRSS1', 'CTRB2',
                                                                'REG1A','CLU','MKI67','KRT8','SPINK1','KRT19','KRT18', 'KRT18','TFF1','MUC1'), cols = c("grey", "deeppink"), reduction = "umap", pt.size = .5)
#Fibroblasts
FeaturePlot(object = TotalTissue.combined.harmony, features = c('ACTA2', 'CDH11', 'PDGFRB', 'COL1A1', 'COL3A1', 
                                                                'RGS5', 'IGFBP7', 'PDPN','DCN','MCAM','IL6','APOE','GLI1','GLI2','GLI3',
                                                                'PDGFA'), cols = c("grey", "deeppink"), reduction = "umap", pt.size = .5)
#ENDOTHELIAL
FeaturePlot(object = TotalTissue.combined.harmony, features = c('VWF','CDH5'), cols = c("grey", "deeppink"), reduction = "umap", pt.size = .5)

#MAST CELLS
FeaturePlot(object = TotalTissue.combined.harmony, features = c('TPSAB1','CPA3'), cols = c("grey", "deeppink"), reduction = "umap", pt.size = .5)

#MYELOID
FeaturePlot(object = TotalTissue.combined.harmony, features = c('CD14','ITGAM','FCGR3A','FCGR3B','APOE',
                                                                'C1QA','MARCO','LYZ','HLA-DRA'), cols = c("grey", "deeppink"), reduction = "umap", pt.size = .5)
#DCS
FeaturePlot(object = TotalTissue.combined.harmony_fil, features = c('ITGAE','LYZ','CLEC9A','BATF3','IRF8','IDO1','CD207','CD1A'
                                                                    ,'CD1C', 'HLA-DRA','CCL22','LAMP3','IL22RA2','CD101'), cols = c("grey", "deeppink"), reduction = "umap", pt.size = .5)
#TCELLS&NK
FeaturePlot(object = TotalTissue.combined.harmony, features = c('CD2', 'CD3D','CD3E','NCAM1','NKG7','CD4','CD8A','PRF1'
                                                                ,'IFNG', 'GZMB','CD69','FOXP3','TIGIT','TOP2A','FCGR3A'), cols = c("grey", "deeppink"), reduction = "umap", pt.size = .5)
#BCELLS
FeaturePlot(object = TotalTissue.combined.harmony, features = c('CD79A', 'MS4A1','CD20','CD138','IGJ','IGLL5','CXCR4','CD117'
                                                                ,'CD27','HLA-DRA'), cols = c("grey", "deeppink"), reduction = "umap", pt.size = .5)
#LABEL THE CLUSTERS round 1 for collapsed populations
Idents(object = TotalTissue.combined.harmony) <- 'UMAP_Clusters'
levels(TotalTissue.combined.harmony)
current.cluster.ids <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50)
new.cluster.ids <- c("CD4 T Cells","CD8 T Cells","Myeloid","Epithelial","Epithelial","Epithelial", "Myeloid","Myeloid","Myeloid","Mast Cells","Myeloid","Acinar",
                     "Myeloid","NK Cells", "Plasma Cells","Fibroblast","Epithelial", "Junk", "Acinar", "Epithelial", "B Cells","Myeloid","Fibroblast","Epithelial","CD4 T Cells",
                     "Fibroblast",  "Acinar", "Epithelial", "Epithelial", "Cycling", "Endothelial", "CD8 T Cells", 
                     "Epithelial", "Epithelial", "Dendritic Cell", "Fibroblast", "CD8 T Cells", "Myeloid", "Cycling", "Myeloid", "Junk", "Junk", "RBC", "Dendritic Cell", "Endocrine", 
                     "Epithelial", "Endothelial", "Epithelial", "Epithelial", "Epithelial", "Plasma Cells")
names(x = new.cluster.ids) <- levels(x = TotalTissue.combined.harmony)
TotalTissue.combined.harmony <- RenameIdents(object = TotalTissue.combined.harmony, new.cluster.ids)
DimPlot(object =TotalTissue.combined.harmony, reduction = "umap", pt.size = 1, label = T)

#Put these new cluster labels as metadata
TotalTissue.combined.harmony[["Collapsed_labels"]] <- Idents(object = TotalTissue.combined.harmony)

# Changing between meta.data for identities- you can change this by altering what you input into your metadata - will need to do this to make the right plots
Idents(object = TotalTissue.combined.harmony) <- 'Collapsed_labels'
# Check active identity
levels(TotalTissue.combined.harmony)

#LABEL THE CLUSTERS round 2 for expanded populations
Idents(object = TotalTissue.combined.harmony) <- 'UMAP_Clusters'
levels(TotalTissue.combined.harmony)

current.cluster.ids <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50)
new.cluster.ids <- c("CD4 T Cells","CD8 T Cells","Macrophage","Epithelial","Epithelial","Epithelial", "Granulocyte","Macrophage","Granulocyte","Mast Cells","Macrophage","Acinar",
                     "Macrophage","NK Cells", "Plasma Cells","iCAF Fibroblast","Epithelial","Junk","Acinar", "Epithelial", "B Cells","Macrophage","iCAF Fibroblast","Epithelial","T Reg",
                     "myCAF Fibroblast", "Acinar", "Epithelial", "Epithelial", "Cycling", "Endothelial", "CD8 T Cells", 
                     "Epithelial", "Epithelial", "Dendritic Cell", "myCAF Fibroblast", "CD8 T Cells", "Granulocyte", "Cycling", "Macrophage", "Junk", "Junk", "RBC", "Dendritic Cell", "Endocrine", 
                     "Epithelial", "Endothelial", "Epithelial", "Epithelial", "Epithelial", "Plasma Cells")
names(x = new.cluster.ids) <- levels(x = TotalTissue.combined.harmony)
TotalTissue.combined.harmony <- RenameIdents(object = TotalTissue.combined.harmony, new.cluster.ids)

#Put these new cluster labels as metadata
TotalTissue.combined.harmony[["Expanded_labels"]] <- Idents(object = TotalTissue.combined.harmony)

DimPlot(object =TotalTissue.combined.harmony, reduction = "umap", pt.size = 1, label = T)

# Changing between meta.data for identities- you can change this by altering what you input into your metadata - will need to do this to make the right plots
Idents(object = TotalTissue.combined.harmony) <- 'Expanded_labels'
# Check active identity
levels(TotalTissue.combined.harmony)

#remove the junk and rbcs and save as a filtered object - '_fil'
TotalTissue.combined.harmony_fil <- SubsetData(TotalTissue.combined.harmony,ident.remove = c("RBC", "Cycling", "Junk"))

# save final object
save(TotalTissue.combined.harmony_fil, file = 'TotalTissue.combined.harmony_fil.RData')
#____________________________________________________________________________________________________________________________________________________________________#

#Subset Individual Tissue Populations

#CD8 T Cells 

Idents(object = TotalTissue.combined.harmony_fil) <- 'Expanded_labels'
levels(TotalTissue.combined.harmony_fil)
CD8_T_cell_subset <- subset(TotalTissue.combined.harmony_fil, idents = "CD8 T Cells")

CD8_T_cell_subset <- NormalizeData(object = CD8_T_cell_subset, normalization.method = "LogNormalize", scale.factor = 10000)
CD8_T_cell_subset <- FindVariableFeatures(object = CD8_T_cell_subset, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, y.cutoff = 0.5)
CD8_T_cell_subset <- ScaleData(object = CD8_T_cell_subset, vars.to.regress = "nCount_RNA", features = rownames(CD8_T_cell_subset))
CD8_T_cell_subset <- RunPCA(object = CD8_T_cell_subset, pc.genes = CD8_T_cell_subset@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)

#Find number of PCs that gives 90% variance
st_dev <- CD8_T_cell_subset@reductions$pca@stdev
var <- st_dev^2
sum(var[1:38])/sum(var)
#Harmony Batch
#Find neighbors and clusters WITH harmony batch correction
options(repr.plot.height = 2.5, repr.plot.width = 6)
CD8_T_cell_subset <- CD8_T_cell_subset %>% 
  RunHarmony("ID", plot_convergence = TRUE)
#Find clusters
CD8_T_cell_subset <- FindNeighbors(object = CD8_T_cell_subset, dims = 1:38, save.SNN = TRUE, force.recalc = T, reduction = "harmony")
CD8_T_cell_subset <- FindClusters(object = CD8_T_cell_subset, resolution = 1.2, verbose = F, reduction = "harmony")
#Run UMAP
CD8_T_cell_subset <- RunUMAP(object = CD8_T_cell_subset, dims = 1:38, reduction = 'harmony')
DimPlot(object = CD8_T_cell_subset, reduction = "umap", pt.size = 2, label = T)

#Assess Cluster Quality with FindAllMarkers Function
T_cell.markers <- FindAllMarkers(CD8_T_cell_subset)

#Label the clusters
current.cluster.ids <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19)
new.cluster.ids <- c("G","CD8","CD8","G","G","G","G","CD8","G","CD8","CD8","CD8","CD8","CD8", "G","CD8","G","G","CD8","CD8")
names(x = new.cluster.ids) <- levels(x = CD8_T_cell_subset)
CD8_T_cell_subset <- RenameIdents(object = CD8_T_cell_subset, new.cluster.ids)
CD8_T_cell_subset[["CD8_T_cell_subset_labels"]]<- Idents(object = CD8_T_cell_subset)

#Filter CD8 T cells
CD8_T_cell_subset_fil <- SubsetData(CD8_T_cell_subset, ident.use = "CD8")
levels(CD8_T_cell_subset_fil)
DimPlot(CD8_T_cell_subset_fil, reduction = "umap", label = F)

#Renormalize
CD8_T_cell_subset_fil <- NormalizeData(object = CD8_T_cell_subset_fil, normalization.method = "LogNormalize", scale.factor = 10000)
CD8_T_cell_subset_fil <- FindVariableFeatures(object = CD8_T_cell_subset_fil, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, y.cutoff = 0.5)
CD8_T_cell_subset_fil <- ScaleData(object = CD8_T_cell_subset_fil, vars.to.regress = "nCount_RNA", features = rownames(CD8_T_cell_subset_fil))
CD8_T_cell_subset_fil <- RunPCA(object = CD8_T_cell_subset_fil, pc.genes = CD8_T_cell_subset_fil@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)

save(CD8_T_cell_subset_fil, file = 'CD8_T_cell_subset_fil.RData')

#Find number of PCs that gives 90% variance
st_dev <- CD8_T_cell_subset_fil@reductions$pca@stdev
var <- st_dev^2
sum(var[1:38])/sum(var)
#Harmony Batch
#Find neighbors and clusters WITH harmony batch correction
options(repr.plot.height = 2.5, repr.plot.width = 6)
CD8_T_cell_subset_fil <- CD8_T_cell_subset_fil %>% 
  RunHarmony("ID", plot_convergence = TRUE)
#Find clusters
CD8_T_cell_subset_fil <- FindNeighbors(object = CD8_T_cell_subset_fil, dims = 1:38, save.SNN = TRUE, force.recalc = T, reduction = "harmony")
CD8_T_cell_subset_fil <- FindClusters(object = CD8_T_cell_subset_fil, resolution = 1.2, verbose = F, reduction = "harmony")
#Run UMAP
CD8_T_cell_subset_fil <- RunUMAP(object = CD8_T_cell_subset_fil, dims = 1:38, reduction = 'harmony')

Idents(object = TotalTissue.combined) <- 'seurat_clusters'
current.cluster.ids <- c(0:5)
new.cluster.ids <- c("CD8 T cell 1", "CD8 T cell 2","CD8 T cell 3","CD8 T cell 4","CD8 T cell 5","CD8 T cell 6")
names(x = new.cluster.ids) <- levels(x = CD8_T_cell_subset_fil)
CD8_T_cell_subset_fil <- RenameIdents(object = CD8_T_cell_subset_fil, new.cluster.ids)
CD8_T_cell_subset_fil[["CD8_expanded"]]<- Idents(object = CD8_T_cell_subset_fil)

#Find all markers to categorize CD8 T cells
T_cell.markers_eff_mem <- FindAllMarkers(CD8_T_cell_subset_fil)
write.csv(T_cell.markers_eff_mem, file = "~/Desktop/T.cell_markers_eff_mem.csv")

#Label T cell populations
Idents(object = TotalTissue.combined) <- 'CD8_expanded'
current.cluster.ids <- c("CD8 T cell 1", "CD8 T cell 2","CD8 T cell 3","CD8 T cell 4","CD8 T cell 5","CD8 T cell 6")
new.cluster.ids <- c("Effector","Exhausted","Effector","Memory","Exhausted","Exhausted")
names(x = new.cluster.ids) <- levels(x = CD8_T_cell_subset_fil)
CD8_T_cell_subset_fil <- RenameIdents(object = CD8_T_cell_subset_fil, new.cluster.ids)
CD8_T_cell_subset_fil[["CD8_collapsed"]]<- Idents(object = CD8_T_cell_subset_fil)

#Save file
save(CD8_T_cell_subset_fil, file = 'CD8_Cells_fil.RData')

#____________________________________________________________________________________________________________________________________________________________________#

#NK Cells
Idents(object = TotalTissue.combined.harmony_fil) <- 'Expanded_labels'
levels(TotalTissue.combined.harmony_fil)
NK_Cells <- subset(TotalTissue.combined.harmony_fil, idents = "NK Cells")
NK_Cells <- NormalizeData(object = NK_Cells, normalization.method = "LogNormalize", scale.factor = 10000)
NK_Cells <- FindVariableFeatures(object = NK_Cells, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, y.cutoff = 0.5)
NK_Cells <- ScaleData(object = NK_Cells, vars.to.regress = "nCount_RNA", features = rownames(NK_Cells))
NK_Cells <- RunPCA(object = NK_Cells, pc.genes = NK_Cells@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)

#Find number of PCs that gives 90% variance
st_dev <- NK_Cells@reductions$pca@stdev
var <- st_dev^2
sum(var[1:42])/sum(var)
#Harmony Batch
#Find neighbors and clusters WITH harmony batch correction
options(repr.plot.height = 2.5, repr.plot.width = 6)
NK_Cells <- NK_Cells %>% 
  RunHarmony("ID", plot_convergence = TRUE)
#Find clusters
NK_Cells <- FindNeighbors(object = NK_Cells, dims = 1:42, save.SNN = TRUE, force.recalc = T, reduction = "harmony")
NK_Cells <- FindClusters(object = NK_Cells, resolution = 1.2, verbose = F, reduction = "harmony")
#Run UMAP
NK_Cells <- RunUMAP(object = NK_Cells, dims = 1:42, reduction = 'harmony')

#Label clusters
current.cluster.ids <- c(0,1,2,3,4,5,6,7,8)
new.cluster.ids <- c("Cluster 1","NK cell 1","NK cell 1","NK cell 1","NK cell 2","Cluster 1","Cluster 2","g","NK cell 3")
names(x = new.cluster.ids) <- levels(x = NK_Cells)
NK_Cells <- RenameIdents(object = NK_Cells, new.cluster.ids)
NK_Cells[["NK_clusters_labels"]]<- Idents(object = NK_Cells)

#Filter NK cells
NK_cell_filter <- SubsetData(NK_Cells, ident.remove = c("g"))
NK_cell_filter[["NK_filter_clusters_labels"]]<- Idents(object = NK_cell_filter)

#Save file
save(NK_cell_filter, file = 'NK_cell_filter.RData')

#____________________________________________________________________________________________________________________________________________________________________#

#CD4 T Cells
Idents(object = TotalTissue.combined.harmony_fil) <- 'Expanded_labels'
levels(TotalTissue.combined.harmony_fil)
CD4_Cells <- subset(TotalTissue.combined.harmony_fil, idents = "CD4 T Cells")
CD4_Cells <- NormalizeData(object = CD4_Cells, normalization.method = "LogNormalize", scale.factor = 10000)
CD4_Cells <- FindVariableFeatures(object = CD4_Cells, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, y.cutoff = 0.5)
CD4_Cells <- ScaleData(object = CD4_Cells, vars.to.regress = "nCount_RNA", features = rownames(CD4_Cells))
CD4_Cells <- RunPCA(object = CD4_Cells, pc.genes = CD4_Cells@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
CD4_Cells_save <- CD4_Cells

#Find number of PCs that gives 90% variance
st_dev <- CD4_Cells@reductions$pca@stdev
var <- st_dev^2
sum(var[1:40])/sum(var)
#Harmony Batch
#Find neighbors and clusters WITH harmony batch correction
options(repr.plot.height = 2.5, repr.plot.width = 6)
CD4_Cells <- CD4_Cells %>% 
  RunHarmony("ID", plot_convergence = TRUE)
#Find clusters
CD4_Cells <- FindNeighbors(object = CD4_Cells, dims = 1:40, save.SNN = TRUE, force.recalc = T, reduction = "harmony")
CD4_Cells <- FindClusters(object = CD4_Cells, resolution = 1.2, verbose = F, reduction = "harmony")
#Run UMAP
CD4_Cells <- RunUMAP(object = CD4_Cells, dims = 1:40, reduction = 'harmony')

#Assess Cluster Quality with ClusterMarkerTable Function
CD4_Cells_clusters <- ClusterMarkerTable(c(0:17),CD4_Cells,100,0.25) #change cluster marker range to match number of clusters shown in UMAP
write.csv(CD4_Cells_clusters, file = "~/Desktop/CD4clusters.csv")

#Label the clusters
Idents(CD4_Cells) <- 'seurat_clusters'
current.cluster.ids <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17)
new.cluster.ids <- c("CD4","CD4","CD4","CD4","CD4","CD4","G","CD4","CD4","CD4","G","CD4","G","CD4", "CD4","CD4","G","G")
names(x = new.cluster.ids) <- levels(x = CD4_Cells)
CD4_Cells <- RenameIdents(object = CD4_Cells, new.cluster.ids)

#Filtered CD4 T cells
CD4_Cells_fil <- SubsetData(CD4_Cells, ident.use = "CD4")

#Renormalize
CD4_Cells_fil <- NormalizeData(object = CD4_Cells_fil, normalization.method = "LogNormalize", scale.factor = 10000)
CD4_Cells_fil <- FindVariableFeatures(object = CD4_Cells_fil, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, y.cutoff = 0.5)
CD4_Cells_fil <- ScaleData(object = CD4_Cells_fil, vars.to.regress = "nCount_RNA", features = rownames(CD4_Cells_fil))
CD4_Cells_fil <- RunPCA(object = CD4_Cells_fil, pc.genes = CD4_Cells_fil@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)

save(CD4_Cells_fil, file = 'CD4_Cells_fil.RData')

#Find number of PCs that gives 90% variance
st_dev <- CD4_Cells_fil@reductions$pca@stdev
var <- st_dev^2
sum(var[1:39])/sum(var)
#Harmony Batch
#Find neighbors and clusters WITH harmony batch correction
options(repr.plot.height = 2.5, repr.plot.width = 6)
CD4_Cells_fil <- CD4_Cells_fil %>% 
  RunHarmony("ID", plot_convergence = TRUE)
#Find clusters
CD4_Cells_fil <- FindNeighbors(object = CD4_Cells_fil, dims = 1:39, save.SNN = TRUE, force.recalc = T, reduction = "harmony")
CD4_Cells_fil <- FindClusters(object = CD4_Cells_fil, resolution = 1.2, verbose = F, reduction = "harmony")
#Run UMAP
CD4_Cells_fil <- RunUMAP(object = CD4_Cells_fil, dims = 1:39, reduction = 'harmony')

#Save  file
save(CD4_Cells_fil, file = 'CD4_Cells_fil.RData')

#____________________________________________________________________________________________________________________________________________________________________#

#Myeloid
Idents(object = TotalTissue.combined.harmony_fil) <- 'Collapsed_labels'
levels(TotalTissue.combined.harmony_fil)
Myeloid_subset <- subset(TotalTissue.combined.harmony_fil, idents = "Myeloid")
Myeloid_subset <- NormalizeData(object = Myeloid_subset, normalization.method = "LogNormalize", scale.factor = 10000)
Myeloid_subset <- FindVariableFeatures(object = Myeloid_subset, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, y.cutoff = 0.5)
Myeloid_subset <- ScaleData(object = Myeloid_subset, vars.to.regress = "nCount_RNA", features = rownames(Myeloid_subset))
Myeloid_subset <- RunPCA(object = Myeloid_subset, pc.genes = Myeloid_subset@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)

#Find number of PCs that gives 90% variance
st_dev <- Myeloid_subset@reductions$pca@stdev
var <- st_dev^2
sum(var[1:33])/sum(var)
#Harmony Batch
#Find neighbors and clusters WITH harmony batch correction
options(repr.plot.height = 2.5, repr.plot.width = 6)
Myeloid_subset <- Myeloid_subset %>% 
  RunHarmony("ID", plot_convergence = TRUE)
#Find clusters
Myeloid_subset <- FindNeighbors(object = Myeloid_subset, dims = 1:33, save.SNN = TRUE, force.recalc = T, reduction = "harmony")
Myeloid_subset <- FindClusters(object = Myeloid_subset, resolution = 1.2, verbose = F, reduction = "harmony")
#Run UMAP
Myeloid_subset <- RunUMAP(object = Myeloid_subset, dims = 1:33, reduction = 'harmony')

#Assess Cluster Quality with FindAllMarkers Function
Myeloid.markers <- FindAllMarkers(Myeloid_subset)
write.csv(Myeloid.markers, file = "~/Desktop/Myeloidcellmarkers.csv")

#Label the clusters 
Idents(Myeloid_subset) <- "Myeloid_subset_clusters_numbered"
levels(Myeloid_subset)
current.cluster.ids <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)
new.cluster.ids <- c("Granulocyte 1", "Classical Macrophages (CCR2+)", "Resident Macrophages", "Alternatively Activated Macrophages","Alternatively Activated Macrophages", 
                     "Granulocyte 1", "Classical Macrophages (CCR2+)", "Granulocyte 1", "Alternatively Activated Macrophages", "Classical Macrophages (CCR2+)", 
                     "Classical Macrophages (CCR2+)", "Granulocyte 1", "Junk", "Junk", "Granulocyte 1", "Alternatively Activated Macrophages",
                     "Macrophage", "Alternatively Activated Macrophages","Resident Macrophages",
                     "Granulocyte 1","Alternatively Activated Macrophages","Junk","Granulocyte 2","Junk", "Mast Cells")
names(x = new.cluster.ids) <- levels(x = Myeloid_subset)
Myeloid_subset <- RenameIdents(object = Myeloid_subset, new.cluster.ids)

#Filter out Junk and Mast cells and save as filtered object
Idents(Myeloid_subset_fil) <- "Myeloid_subset_labels_expanded"
levels(Myeloid_subset_fil)
Myeloid_subset_fil <- SubsetData(Myeloid_subset,ident.remove = c("Junk", "Mast Cells"))

#Save file
save(Myeloid_subset_fil, file = 'Myeloid_subset_fil.RData')

#____________________________________________________________________________________________________________________________________________________________________#

#Dendritic Cells

Idents(object = TotalTissue.combined.harmony_fil) <- 'Expanded_labels'
levels(TotalTissue.combined.harmony_fil)
Dendritic_Cells <- subset(TotalTissue.combined.harmony_fil, idents = "Dendritic Cell")
Dendritic_Cells <- NormalizeData(object = Dendritic_Cells, normalization.method = "LogNormalize", scale.factor = 10000)
Dendritic_Cells <- FindVariableFeatures(object = Dendritic_Cells, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, y.cutoff = 0.5)
Dendritic_Cells <- ScaleData(object = Dendritic_Cells, vars.to.regress = "nCount_RNA", features = rownames(Dendritic_Cells))
Dendritic_Cells <- RunPCA(object = Dendritic_Cells, pc.genes = Dendritic_Cells@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)

#Find number of PCs that gives 90% variance
st_dev <- Dendritic_Cells@reductions$pca@stdev
var <- st_dev^2
sum(var[1:35])/sum(var)
#Harmony Batch
#Find neighbors and clusters WITH harmony batch correction
options(repr.plot.height = 2.5, repr.plot.width = 6)
Dendritic_Cells <- Dendritic_Cells %>% 
  RunHarmony("ID", plot_convergence = TRUE)
#Find clusters
Dendritic_Cells <- FindNeighbors(object = Dendritic_Cells, dims = 1:35, save.SNN = TRUE, force.recalc = T, reduction = "harmony")
Dendritic_Cells <- FindClusters(object = Dendritic_Cells, resolution = 1.2, verbose = F, reduction = "harmony")
#Run UMAP
Dendritic_Cells <- RunUMAP(object = Dendritic_Cells, dims = 1:35, reduction = 'harmony')

#Assess Cluster Quality with ClusterMarkerTable Function
Dendritic_Cells_clusters <- ClusterMarkerTable(c(0:9),Dendritic_Cells,100,0.25) #change cluster marker range to match number of clusters shown in UMAP
write.csv(Dendritic_Cells_clusters, file = "~/Desktop/DCclusters.csv")

#Label the clusters
Idents(Dendritic_Cells) <- 'seurat_clusters'
current.cluster.ids <- c(0,1,2,3,4,5,6,7,8,9)
new.cluster.ids <- c("Dendritic_Cells","Dendritic_Cells","Dendritic_Cells","Dendritic_Cells","Dendritic_Cells","Dendritic_Cells","Dendritic_Cells","Dendritic_Cells","Dendritic_Cells","J")
names(x = new.cluster.ids) <- levels(x = Dendritic_Cells)
Dendritic_Cells <- RenameIdents(object = Dendritic_Cells, new.cluster.ids)

#Filter Dendrtiic cells
Dendritic_Cells_fil <- SubsetData(Dendritic_Cells, ident.use = "Dendritic_Cells")
levels(Dendritic_Cells_fil)
DimPlot(Dendritic_Cells_fil, reduction = "umap", label = T)

#Renormalize
Dendritic_Cells_fil <- NormalizeData(object = Dendritic_Cells_fil, normalization.method = "LogNormalize", scale.factor = 10000)
Dendritic_Cells_fil <- FindVariableFeatures(object = Dendritic_Cells_fil, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, y.cutoff = 0.5)
Dendritic_Cells_fil <- ScaleData(object = Dendritic_Cells_fil, vars.to.regress = "nCount_RNA", features = rownames(Dendritic_Cells_fil))
Dendritic_Cells_fil <- RunPCA(object = Dendritic_Cells_fil, pc.genes = Dendritic_Cells_fil@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)

save(Dendritic_Cells_fil, file = 'Dendritic_Cells_fil.RData')

#Find number of PCs that gives 90% variance
st_dev <- Dendritic_Cells_fil@reductions$pca@stdev
var <- st_dev^2
sum(var[1:38])/sum(var)
#Harmony Batch
#Find neighbors and clusters WITH harmony batch correction
options(repr.plot.height = 2.5, repr.plot.width = 6)
Dendritic_Cells_fil <- Dendritic_Cells_fil %>% 
  RunHarmony("ID", plot_convergence = TRUE)
#Find clusters
Dendritic_Cells_fil <- FindNeighbors(object = Dendritic_Cells_fil, dims = 1:38, save.SNN = TRUE, force.recalc = T, reduction = "harmony")
Dendritic_Cells_fil <- FindClusters(object = Dendritic_Cells_fil, resolution = 1.2, verbose = F, reduction = "harmony")

#Run UMAP
Dendritic_Cells_fil <- RunUMAP(object = Dendritic_Cells_fil, dims = 1:38, reduction = 'harmony')

#Find all markers to categorize Dendritic cells
Dendritic_Cells_top_fil <- FindAllMarkers(Dendritic_Cells_fil, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#Label clusters
Idents(Dendritic_Cells_fil) <- 'seurat_clusters'
current.cluster.ids <- c(0,1,2,3,4,5,6,7,8,9)
new.cluster.ids <- c("pDC1",
                     "cDC2_A",
                     "Langerhan_like_DC1",
                     "cDC1",
                     "Langerhan_like_DC2",
                     "cDC2_B",
                     "cDC2_B",
                     "pDC2",
                     "Activated_DC1",
                     "Activated_DC2")
names(x = new.cluster.ids) <- levels(x = Dendritic_Cells_fil)
Dendritic_Cells_fil[["Dendritic_official"]]<- Idents(object = Dendritic_Cells_fil)

#Save file
save(Dendritic_Cells_fil, file = 'Dendritic_Cells_fil.RData')

#____________________________________________________________________________________________________________________________________________________________________#

#PBMC Object 

#Load in all filtered runs using Read10X or Read10X_h5 functions
#   Data.data <- Read10X("filepath/")

PDAC_PBMC_1.data <- Read10X("~/Desktop/PDAC_PBMC_1/filtered_feature_bc_matrix/")

#Generate Seurat objects using CreateSeuratObject function, min.cells=3, min.features=100
# Data <- CreateSeuratObject(counts = Data.data, project = "PDA", min.cells=3, min.features=100)

PDAC_PBMC_1 <- CreateSeuratObject(counts = PDAC_PBMC_1.data, project = 'PDAC_PBMC_1', min.cells = 3, min.features = 100)

#Add Metadata

#ID Metadata
PDAC_PBMC_16$ID <- "PDAC_PBMC_16"
PDAC_PBMC_15$ID <- "PDAC_PBMC_15"
PDAC_PBMC_13$ID <- "PDAC_PBMC_13"
PDAC_PBMC_14$ID <- "PDAC_PBMC_14"
PDAC_PBMC_12$ID <- "PDAC_PBMC_12"
PDAC_PBMC_10B$ID <- "PDAC_PBMC_10B" 
PDAC_PBMC_10A$ID <- "PDAC_PBMC_10A" 
PDAC_PBMC_11$ID <- "PDAC_PBMC_11"
PDAC_PBMC_9$ID <- "PDAC_PBMC_9"
PDAC_PBMC_8$ID <- "PDAC_PBMC_8"
PDAC_PBMC_7$ID <- "PDAC_PBMC_7"
PDAC_PBMC_6$ID <- "PDAC_PBMC_6"
PDAC_PBMC_5$ID <- "PDAC_PBMC_5"
PDAC_PBMC_4$ID <- "PDAC_PBMC_4"
PDAC_PBMC_3$ID <- "PDAC_PBMC_3"
PDAC_PBMC_2$ID <- "PDAC_PBMC_2"
PDAC_PBMC_1$ID <- "PDAC_PBMC_1"
Healthy_PBMC_1$ID <-"Healthy_PBMC_1"
Healthy_PBMC_2$ID <- "Healthy_PBMC_2"
Healthy_PBMC_3$ID <- "Healthy_PBMC_3"
Healthy_PBMC_4$ID <- "Healthy_PBMC_4"

#Disease State Metadata
PDAC_PBMC_16$DiseaseState <- "PDAC"
PDAC_PBMC_15$DiseaseState <- "PDAC"
PDAC_PBMC_13$DiseaseState <- "PDAC"
PDAC_PBMC_14$DiseaseState <- "PDAC"
PDAC_PBMC_12$DiseaseState <- "PDAC"
PDAC_PBMC_10B$DiseaseState <- "PDAC" 
PDAC_PBMC_10A$DiseaseState <- "PDAC" 
PDAC_PBMC_11$DiseaseState <- "PDAC"
PDAC_PBMC_9$DiseaseState <- "PDAC"
PDAC_PBMC_8$DiseaseState <- "PDAC"
PDAC_PBMC_7$DiseaseState <- "PDAC"
PDAC_PBMC_6$DiseaseState <- "PDAC"
PDAC_PBMC_5$DiseaseState <- "PDAC"
PDAC_PBMC_4$DiseaseState <- "PDAC"
PDAC_PBMC_3$DiseaseState <- "PDAC"
PDAC_PBMC_2$DiseaseState <- "PDAC"
PDAC_PBMC_1$DiseaseState <- "PDAC"
Healthy_PBMC_1$DiseaseState <-"Healthy"
Healthy_PBMC_2$DiseaseState <- "Healthy"
Healthy_PBMC_3$DiseaseState <- "Healthy"
Healthy_PBMC_4$DiseaseState <- "Healthy"

#Merge all objects
PBMC_Merge <- merge(PDAC_PBMC_1, y = c(PDAC_PBMC_2, PDAC_PBMC_3, PDAC_PBMC_4, 
                                                   PDAC_PBMC_5, PDAC_PBMC_6, PDAC_PBMC_7, 
                                                   PDAC_PBMC_8, PDAC_PBMC_9, PDAC_PBMC_11, 
                                                   PDAC_PBMC_10A, PDAC_PBMC_10B, PDAC_PBMC_12,
                                                   PDAC_PBMC_13, PDAC_PBMC_14, PDAC_PBMC_15, 
                                                   PDAC_PBMC_16))

#Calculate MT %
PBMC_Merge[["percent.mt"]] <- PercentageFeatureSet(object = PBMC_Merge, pattern = "^MT-")

#Normalize Data
PBMC_Merge <- NormalizeData(object = PBMC_Merge, normalization.method = "LogNormalize", scale.factor = 10000)

#Find Variable Genes
PBMC_Merge <- FindVariableFeatures(object = PBMC_Merge, selection.method = "vst", nfeatures = 2000)

#Calculate Cell Cycle Score (S-G2M Difference)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
PBMC_Merge <- CellCycleScoring(PBMC_Merge, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
PBMC_Merge$CC.Difference <- PBMC_Merge$S.Score - PBMC_Merge$G2M.Score

#THIS STEP MAY TAKE A VERY LONG TIME
#Scale Data
all.genes <- rownames(x = PBMC_Merge)
PBMC_Merge <- ScaleData(PBMC_Merge, vars.to.regress = c("nCount_RNA","CC.Difference"), features = all.genes)

#It is recommended that you save this object so you do not have to re-run the scaling step more than necessary
# save Seurat object
save(PBMC_Merge, file = 'PBMC_Merge.RData')

#Run PCA and Determine Dimensions for 90% Variance
PBMC_Merge <- RunPCA(object = PBMC_Merge, features = VariableFeatures(object = PBMC_Merge))
stdev <- PBMC_Merge@reductions$pca@stdev
var <- stdev^2
sum(var[1:30])/ sum(var)


#Batch Correction
PBMC_Merge <- RunHarmony(PBMC_Merge, group.by.vars = c("Patient_ID"), dims.use = 1:30, verbose = F)

#Find Neighbors + Find Clusters
PBMC_Merge <- FindNeighbors(object = PBMC_Merge, dims = 1:30, reduction = "harmony")
PBMC_Merge <- FindClusters(object = PBMC_Merge, resolution = 1.2)

#Run UMAP and get unlabeled cluster UMAP and violin plot
PBMC_Merge <- RunUMAP(object = PBMC_Merge, dims = 1:30, reduction = "harmony")
DimPlot(object = PBMC_Merge_Fil, reduction = "umap", label = TRUE, pt.size = 0.5)

#Determine n Cutoffs
VlnPlot(PBMC_Merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
PBMC_Merge_Fil <- subset(x = PBMC_Merge, subset = nCount_RNA < 25000 )
PBMC_Merge_Fil <- subset(x = PBMC_Merge, subset = nCount_RNA < 25000 & percent.mt < 15 )
PBMC_Merge_Fil <- subset(x = PBMC_Merge, subset = nCount_RNA > 1000 & nCount_RNA < 25000 & percent.mt < 15)
VlnPlot(PBMC_Merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
table(PBMC_Merge_Fil@active.ident)

#Subset Seurat Object
PBMC_Merge_Fil <- subset(x = PBMC_Merge, idents =  c(0:4,6:12,14,15,17:20,23:27,29:32))

#Normalize Data
PBMC_Merge_Fil <- NormalizeData(object = PBMC_Merge_Fil, normalization.method = "LogNormalize", scale.factor = 10000)

#Find Variable Genes
PBMC_Merge_Fil <- FindVariableFeatures(object = PBMC_Merge_Fil, selection.method = "vst", nfeatures = 2000)

#Calculate Cell Cycle Score (S-G2M Difference)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
PBMC_Merge_Fil <- CellCycleScoring(PBMC_Merge_Fil, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
PBMC_Merge_Fil$CC.Difference <- PBMC_Merge_Fil$S.Score - PBMC_Merge_Fil$G2M.Score

#Scale Data
all.genes <- rownames(x = PBMC_Merge_Fil)
PBMC_Merge_Fil <- ScaleData(PBMC_Merge_Fil, vars.to.regress = c("nCount_RNA","CC.Difference"), features = all.genes)

#Run PCA and Determine Dimensions for 90% Variance
PBMC_Merge_Fil <- RunPCA(object = PBMC_Merge_Fil, features = VariableFeatures(object = PBMC_Merge_Fil))
stdev <- PBMC_Merge_Fil@reductions$pca@stdev
var <- stdev^2
sum(var[1:26])/ sum(var)

KillSwitch = 0

for(i in 1:length(var)){
  total <- sum(var)
  numerator <- sum(var[1:i])
  expvar <- numerator/total
  if(KillSwitch == 0){
    if(expvar > 0.9){
      KillSwitch <- KillSwitch + 1
      PCNum <- i
    }
  }
}

#Batch Correction
PBMC_Merge_Fil <- RunHarmony(PBMC_Merge_Fil, group.by.vars = c("Patient_ID"), dims.use = 1:PCNum, verbose = F)

#Find Neighbors + Find CLusters
PBMC_Merge_Fil <- FindNeighbors(object = PBMC_Merge_Fil, dims = 1:PCNum, reduction = "harmony")
PBMC_Merge_Fil <- FindClusters(object = PBMC_Merge_Fil, resolution = 1.2)

#Run UMAP and get unlabelled cluster UMAP and violin plot
PBMC_Merge_Fil <- RunUMAP(object = PBMC_Merge_Fil, dims = 1:31, reduction = "harmony")

#Run ACMT Function for Annotating
CMT_PBMC_Merge_Fil <- AutomatedClusterMarkerTable(Seurat_Object = PBMC_Merge_Fil)
write.csv(CMT_PBMC_Merge_Fil[[1]], file = "~/Desktop/CMT_PBMC_Merge_Fil.csv") 
View(CMT_PBMC_Merge_Fil[[1]])

#Final Labelling
Idents(object = PBMC_Merge_Fil) <- "seurat_clusters"
levels(PBMC_Merge_Fil)
PBMC_Merge_Fil <- RenameIdents(PBMC_Merge_Fil, 
                               "0" = "CD4 T Cells",
                               "1" = "CD4 T Cells",
                               "2" = "Granulocyte",
                               "3" = "Granulocyte",
                               "4" = "NK Cells",
                               "5" = "CD8 T Cells",
                               "6" = "CD8 T Cells",
                               "7" = "Monocyte",
                               "8" = "Granulocyte",
                               "9" = "Monocyte",
                               "10" = "Monocyte",
                               "11" = "Granulocyte",
                               "12" = "Granulocyte",
                               "13" = "B Cells",
                               "14" = "Monocyte",
                               "15" = "CD8 T Cells",
                               "16" = "CD4 T Cells",
                               "17" = "CD4 T Cells",
                               "18" = "NK Cells",
                               "19" = "Monocyte",
                               "20" = "Plasma Cells",
                               "21" = "Junk",
                               "22" = "CD8 T Cells",
                               "23" = "Junk",
                               "24" = "Junk",
                               "25" = "CD4 T Cells",
                               "26" = "RBC",
                               "27" = "pDCs",
                               "28" = "CD4 T Cells")
PBMC_Merge_Fil <- SubsetData(PBMC_Merge_Fil, ident.remove = c("Junk", "RBC"))

#Save object
save(PBMC_Merge_Fil, file = 'PBMC_Merge_Fil.RData')

#____________________________________________________________________________________________________________________________________________________________________#

# 3. Interactome Analysis - Supplementary Figure S5D

#Subset normal panc and PDA tissue
Idents(object = TotalTissue.combined.harmony_fil) <- 'DiseaseState'
object_0 <- subset(TotalTissue.combined.harmony_fil,subset = "AdjNorm")  
object_1 <- subset(TotalTissue.combined.harmony_fil,subset = "PDAC")

# make sure that it has 2 columns (ligand-receptor)
LR_Pairs_database <- as.data.frame(fread("Literature Supported Receptor Ligand Pairs Ramilowski.csv", header = T, stringsAsFactors = F))
LR_Pairs_database <- LR_Pairs_database[,1:2]

# list of wanted genes
ligand_genes <- as.character(LR_Pairs_database$Ligand.ApprovedSymbol)
receptor_genes <- as.character(LR_Pairs_database$Receptor.ApprovedSymbol)

# Create merged Seurat object
meta_label <- 'Expanded_labels'
# Switch ids of objects to meta_label
Idents(object_0) <- meta_label
Idents(object_1) <- meta_label
merge <- merge(object_0, object_1)

# Find common cell type populations and rename if need
object0_cell_types <- as.character(levels(object_0))
object1_cell_types <- as.character(levels(object_1))
cell_types <- object0_cell_types[which(object0_cell_types %in% object1_cell_types)]
View(cell_types)
cell_types <- cell_types[c(1,3,6,7,10)] # where you choose cell types want to look at

# Get info for only common cell types
Idents(merge) <- meta_label
merge <- subset(merge, idents = cell_types)

# Get gene list in object 
merge_genes <- rownames(merge)

# Create group Seurat objects
Idents(merge) <- 'DiseaseState'
normal <- subset(merge, idents = 'AdjNorm')
pda <- subset(merge, idents = 'PDAC')
Idents(normal) <- meta_label
Idents(pda) <- meta_label

'-------------------------------------------------------------------------------'
not_found_ligands <- check_genes(ligand_genes, LR_Pairs_database, merge_genes) 
not_found_receptors <- check_genes(receptor_genes, LR_Pairs_database, merge_genes)

# separate genes into ligands and receptors and get rid of genes not found
not_in_lig <- unique(not_found_ligands$object)
wanted_ligands <- unique(ligand_genes[!ligand_genes %in% not_in_lig])

not_in_rec <- unique(not_found_receptors$object)
wanted_receptors <- unique(receptor_genes[!receptor_genes %in% not_in_rec])

genes <- unique(c(wanted_ligands, wanted_receptors))

# get need info from object into table
# change second vars to what meta.data refers to labels
normal_data <- FetchData(object = normal, vars = c(genes, meta_label))
pda_data <- FetchData(object = pda, vars = c(genes, meta_label))

# Make all possible LR pairs
LRs <- make_LR_pairs(wanted_ligands, wanted_receptors, LR_Pairs_database)

# Calculate average expression for all cell types for each group
avg_0 <- AverageExpression(object = normal, features = genes, verbose = F, alpha = 0.05)
avg_1 <- AverageExpression(object = pda, features = genes, verbose = F, alpha = 0.05)

# Create table w/ average expressions for every ligand/receptor and source/target cell combination
LR_table <- create_LR_table(wanted_ligands, wanted_receptors, cell_types, LRs, avg_0$RNA, avg_1$RNA)
LR_table[,1:4] <- data.frame(apply(LR_table[,1:4], 2, as.character), stringsAsFactors = FALSE)

# Choose threshold --> here I'm using median
summary(c(LR_table$avg_lig_0, LR_table$avg_lig_1, LR_table$avg_rec_0, LR_table$avg_rec_1))
threshold <- median(c(LR_table$avg_lig_0, LR_table$avg_lig_1, LR_table$avg_rec_0, LR_table$avg_rec_1))

# Filter LR pairs based on average expresions
LR_table <- avg_LR_filt(LR_table, threshold)

# Calculate p-value between groups
LR_table <- LR_diff(LR_table, normal_data, pda_data, genes, meta_label)

# Filter LR pairs where expression is constant between groups
LR_table <- LR_table[which(LR_table$lig_diff != 0 | LR_table$rec_diff != 0),]

# Sort by ligand
LR_table <- arrange(LR_table, desc(lig_diff != 0 & rec_diff != 0), lig_diff_p_adj)

write.csv(LR_table, file = 'human_pda_tissue_LR_lig_sort.csv', row.names = F)

# list of genes (Does not go into cytoscape, can be used for functional annotation)
genes_list <- c(LR_table$ligands, LR_table$receptors)
genes_list <- unique(genes_list)
write.csv(genes_list, file = 'human_pda_tissue_LR_genes.csv')

# get supplemental data from Cytoscape (need to have at least shared name column)
cytoscape_nodes <- fread('human_pda_tissue_LR_lig_sort.csv default node.csv')
node_supplement <- generate_supplement(LR_table, cytoscape_nodes)
write.csv(node_supplement, file = 'human_pda_tissue_LR_names.csv', row.names = F)


#____________________________________________________________________________________________________________________________________________________________________#


# 4. Figures/Data Analysis

# Figure 2A

Idents(TotalTissue.combined.harmony_fil) <- "Collapsed_labels"
DimPlot(object = TotalTissue.combined.harmony_fil, reduction = "umap", split.by = "DiseaseState", pt.size = 0.5, label = F, ncol=5, 
        cols = c("CD8 T Cells" = "darkgreen", 
                 "CD4 T Cells" = "limegreen", 
                 "T Reg" = "darkseagreen4",
                 "Macrophage" = "darkorange2", 
                 "Granulocyte" = "orange1",
                 "Dendritic Cell" = "chocolate4",
                 "Mast Cells" = "gold",
                 "NK Cells" = "darkorchid",
                 "B Cells" = "dodgerblue1",
                 "Plasma Cells" = "navy",
                 "Fibroblasts" = "darkcyan",
                 "Acinar" = "deeppink",
                 "Endothelial" = "mediumvioletred",
                 'Epithelial' = "red3",
                 "Endocrine" = "darkred",
                 "Plasma Cells" = "lemonchiffon"))

# Figure 2B

Idents(TotalTissue.combined.harmony_fil) <- "Collapsed_labels"
DotPlot(TotalTissue.combined.harmony_fil, features = c('MS4A1', 'CD79A','IGJ','LYZ', 'APOE', 'ITGAM','ITGAX', 'GZMB', 'HLA-DRA', 'CD14', 'CD3D',
                                                       'NKG7','NCAM1', 'CD3E','PDGFRB', 'RGS5', 'IRF7','PLVAP', 'VWF', 'CDH5','TPSAB1', 'CPA3', 'LUM', 'DCN','PRSS1', 'CTRB2',
                                                       'REG1A','SPP1','MMP7', 'CLU','KRT8','KRT19', 'KRT18','TFF1','CCL22','LAMP3')
        , cols = c('red', 'blue'), dot.scale = 10) + RotatedAxis() + FontSize(x.text = 17, y.text = 17)

# Figure 2C



# Figure 2D



# Figure 2E



# Figure 2F



# Figure 3A

Idents(CD8_T_cell_subset_fil) <- "CD8_collapsed"
DimPlot(object = CD8_T_cell_subset_fil, reduction = "umap", pt.size = 2, label = F, split.by = "DiseaseState", 
        ncol = 4,cols = c("CD8 T cell 1"="deeppink","CD8 T cell 3"="springgreen4",
                          "CD8 T cell 4"="deepskyblue3","CD8 T cell 5"="olivedrab2",
                          "CD8 T cell 6"="darkturquoise","CD8 T cell 2"="plum"))

# Figure 3B

Idents(CD8_T_cell_subset_fil) <- 'seurat_clusters'
T_cell.markers_fil <- FindAllMarkers(CD8_T_cell_subset_fil)
write.csv(T_cell.markers, file = "~/Desktop/Tcellmarkers_fil.csv")
top_10_fil <- T_cell.markers_fil%>%group_by(cluster)%>%top_n(n=10, wt = avg_logFC)
DoHeatmap(CD8_T_cell_subset_fil, features = top_10_fil$gene, size = 5.5)

# Figure 3C



# Figure 3D/Supplementary Figure S4B

table(CD8_T_cell_subset_fil$ID, CD8_T_cell_subset_fil$CD8_collapsed)
table(CD8_T_cell_subset_fil$ID)

# Figure 3E 

Idents(CD8_T_cell_subset_fil) <- "CD8_collapsed"
Eff_CD8 <- subset(CD8_T_cell_subset_fil, idents = "Effector")
Idents(Eff_CD8) <- "DiseaseState"
VlnPlot(Eff_CD8, features = c("GZMA","KLF2","RORA","KLRG1","CD74","LTB","SLC4A10","TNFRSF25","HLA-C","HIF1A","GZMB","ITGA1"), cols = c("red","blue"), pt.size = 0, ncol = 2)

# Figure 3F 

Idents(CD8_T_cell_subset_fil) <- "CD8_collapsed"
Exh_CD8 <- subset(CD8_T_cell_subset_fil, idents = "Exhausted")
Idents(Exh_CD8) <- "DiseaseState"
VlnPlot(Exh_CD8, features = c("EOMES","KLF2","GIMAP7","TCF7","ARHGAP25","SELPLG","STAT1","HLA-DPB1","IFITM2","SYNE2","CD84","MBP"), cols = c("red","blue"), pt.size = 0, ncol = 2)

# Figure 4A

Idents(NK_cell_filter) <- "NK_filter_clusters_labels"
DimPlot(object = NK_cell_filter, reduction = "umap", pt.size = 1.3, label = F, 
        cols = c("NK cell 1" = "darkorchid", "Cluster 1" = "darkcyan", "NK cell 2" = "lightsalmon2", 
                 "Cluster 2" = "navy", "NK cell 3" = "lightpink1"), 
        split.by = "DiseaseState")

# Figure 4B

Idents(NK_cell_filter) <- "NK_filter_clusters_labels"
NK_cell_filter_markers <- FindAllMarkers(NK_cell_filter, only.pos = T)
write.csv(NK_cell_filter_markers, file = "NK_cell_filter_markers.CSV")
NK_cell_filter_top10 <- NK_cell_filter_markers%>%group_by(cluster)%>%top_n(n=10, wt = avg_logFC)
DoHeatmap(NK_cell_filter, features = NK_cell_filter_top10$gene, size = 5.5)

# Figure 4C

Idents(NK_cell_filter) <- "NK_filter_clusters_labels"
VlnPlot(NK_cell_filter, features = c("NKG7","TIGIT","FCGR3A","LAG3","NCAM1","PDCD1",
                                     "NCR1","HAVCR2","CD3E","PRF1","CD8A","GZMB"),
        cols = c("Cluster 1" = "darkcyan", "NK cell 1" = "darkorchid", "NK cell 2" = "lightsalmon2", 
                 "Cluster 2" = "navy", "NK cell 3" = "lightpink1"), pt.size = 0.5, ncol = 2, y.max = 2)

# Figure 4D



# Figure 4E

Idents(CD4_Cells_fil) <- "seurat_clusters"
DimPlot(object = CD4_Cells_fil, reduction = "umap", pt.size = 2, label = F, 
        cols = c("0"="gray70","1"="darkseagreen4","2"="aquamarine","3"="darkblue","4"="darkolivegreen",
                 "5"="mediumseagreen","6"="palegreen","7"="seagreen","8"="palegreen4","9"="lightgreen","10"="green4",
                 "11"="forestgreen","12"="olivedrab1"))

# Figure 4F

Idents(CD4_Cells_fil) <- "seurat_clusters"
CD4_topgenes_fil <- FindAllMarkers(CD4_Cells_fil, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
rp_mt_genes <- CD4_topgenes_fil$gene[grep("^RP|^MT-", CD4_topgenes_fil$gene)]
CD4_topgenes_fil_1 <- CD4_topgenes_fil %>% filter(!gene %in% rp_mt_genes)
top10_CD4 <- CD4_topgenes_fil_1 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(CD4_Cells_fil, features = top10_CD4$gene)

# Figure 4G

FeaturePlot(CD4_Cells_fil, c("CD4","FOXP3","CTLA4","TIGIT"), cols = c('grey77','red'), pt.size = 1.5)

# Figure 5A



# Figure 5B



# Figure 5C

FeaturePlot(Myeloid_subset_fil, c("LGALS9","CD274","PVR","CSF1R","SIRPA","HLA-DQ1"), cols = c('grey77','red'), pt.size = 1.5)

# Figure 5D

FeaturePlot(Myeloid_subset_fil, c("MARCO","FCGR3B","CD68"), cols = c('grey77','red'), pt.size = 1.5)
Idents(Myeloid_subset_fil) <- "DiseaseState"
VlnPlot(Myeloid_subset_fil, features = c("LGALS9","SIRPA","PVR"),
        cols = c("red","blue"), pt.size = 0.5, ncol = 2, y.max = 2)

# Figure 5E

Idents(Dendritic_Cells_fil) <- "Dendritic_official"
DimPlot(object = Dendritic_Cells_fil, reduction = "umap", pt.size = 2, label = F)

# Figure 5F

Idents(Dendritic_Cells_fil) <- "Dendritic_official"
top10_Dendritic_fil <- Dendritic_Cells_top_fil %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(Dendritic_Cells_fil, features = top10_Dendritic_fil$gene)

# Figure 6A

Idents(Dendritic_Cells_fil) <- "Dendritic_official"
VlnPlot(Myeloid_subset_fil, features = c("CLEC9A","BATF3","CD1A","CD207","CCL22",
                                         "LAMP3","GZMB","IRF8","SIRPA","ITGAE",
                                         "ITGAM","HLA-DQ1"),pt.size = 0.5, ncol = 2)


# Figure 6B



# Figure 6C

Idents(TotalTissue.combined.harmony_fil) <- "Collapsed_labels"
Mac_Lig_Gran_DC_CD4_CD8_NK_Rec <- CircosFunctions(InteractomeData = Full_list_LR_012720_lig_sort_collapsed_labels, SeuratObject = TotalTissue.combined.harmony_fil, 
                                             CellTypeVector = c("Macrophage","Granulocyte","Dendritic Cells","CD4 T Cells","CD8 T Cells", "NK Cells"), POI = 12, Lig_or_Rec = T,
                                             LR_List = Literature.Supported.Receptor.Ligand.Pairs.Ramilowski, Species = T, 
                                             Cutoff_Lig = -0.5, Cutoff_Rec = -0.5)

write.table(x = Mac_Lig_Gran_DC_CD4_CD8_NK_Rec[[1]], file = "Mac_Lig_Gran_DC_CD4_CD8_NK_Rec_Karyotype.txt", row.names = F, quote = F, col.names = F)
write.table(x = Mac_Lig_Gran_DC_CD4_CD8_NK_Rec[[2]], file = "Mac_Lig_Gran_DC_CD4_CD8_NK_Rec_Text.txt", row.names = F, quote = F, col.names = F)
write.table(x = Mac_Lig_Gran_DC_CD4_CD8_NK_Rec[[3]], file = "Mac_Lig_Gran_DC_CD4_CD8_NK_Rec_Expression.txt", row.names = F, quote = F, col.names = F)
write.table(x = Mac_Lig_Gran_DC_CD4_CD8_NK_Rec[[4]], file = "Mac_Lig_Gran_DC_CD4_CD8_NK_Rec_Links.txt", row.names = F, quote = F, col.names = F)

# Figure 6D

Idents(TotalTissue.combined.harmony_fil) <- "Collapsed_labels"
Gran_Lig_Mac_DC_CD4_CD8_NK_Rec <- CircosFunctions(InteractomeData = Full_list_LR_012720_lig_sort_collapsed_labels, SeuratObject = TotalTissue.combined.harmony_fil, 
                                             CellTypeVector = c("Macrophage","Granulocyte","Dendritic Cells","CD4 T Cells","CD8 T Cells", "NK Cells"), POI = 11, Lig_or_Rec = T,
                                             LR_List = Literature.Supported.Receptor.Ligand.Pairs.Ramilowski, Species = T, 
                                             Cutoff_Lig = -0.5, Cutoff_Rec = -0.5)

write.table(x = Gran_Lig_Mac_DC_CD4_CD8_NK_Rec[[1]], file = "Gran_Lig_Mac_DC_CD4_CD8_NK_Rec_Karyotype.txt", row.names = F, quote = F, col.names = F)
write.table(x = Gran_Lig_Mac_DC_CD4_CD8_NK_Rec[[2]], file = "Gran_Lig_Mac_DC_CD4_CD8_NK_Rec_Text.txt", row.names = F, quote = F, col.names = F)
write.table(x = Gran_Lig_Mac_DC_CD4_CD8_NK_Rec[[3]], file = "Gran_Lig_Mac_DC_CD4_CD8_NK_Rec_Expression.txt", row.names = F, quote = F, col.names = F)
write.table(x = Gran_Lig_Mac_DC_CD4_CD8_NK_Rec[[4]], file = "Gran_Lig_Mac_DC_CD4_CD8_NK_Rec_Links.txt", row.names = F, quote = F, col.names = F)

# Figure 6E

Idents(TotalTissue.combined.harmony_fil) <- "Collapsed_labels"
DC_Lig_Mac_Gran_CD4_CD8_NK_Rec <- CircosFunctions(InteractomeData = Full_list_LR_012720_lig_sort_collapsed_labels, SeuratObject = TotalTissue.combined.harmony_fil, 
                                             CellTypeVector = c("Macrophage","Granulocyte","Dendritic Cells","CD4 T Cells","CD8 T Cells", "NK Cells"), POI = 16, Lig_or_Rec = T,
                                             LR_List = Literature.Supported.Receptor.Ligand.Pairs.Ramilowski, Species = T, 
                                             Cutoff_Lig = -0.5, Cutoff_Rec = -0.5)

write.table(x = DC_Lig_Mac_Gran_CD4_CD8_NK_Rec[[1]], file = "DC_Lig_Mac_Gran_CD4_CD8_NK_Rec_Karyotype.txt", row.names = F, quote = F, col.names = F)
write.table(x = DC_Lig_Mac_Gran_CD4_CD8_NK_Rec[[2]], file = "DC_Lig_Mac_Gran_CD4_CD8_NK_Rec_Text.txt", row.names = F, quote = F, col.names = F)
write.table(x = DC_Lig_Mac_Gran_CD4_CD8_NK_Rec[[3]], file = "DC_Lig_Mac_Gran_CD4_CD8_NK_Rec_Expression.txt", row.names = F, quote = F, col.names = F)
write.table(x = DC_Lig_Mac_Gran_CD4_CD8_NK_Rec[[4]], file = "DC_Lig_Mac_Gran_CD4_CD8_NK_Rec_Links.txt", row.names = F, quote = F, col.names = F)

# Figure 6F

Idents(TotalTissue.combined.harmony_fil) <- "Collapsed_labels"
Endo_Lig_CD4_CD8_NK_Rec <- CircosFunctions(InteractomeData = Full_list_LR_012720_lig_sort_collapsed_labels, SeuratObject = TotalTissue.combined.harmony_fil, 
                                             CellTypeVector = c("Endothelial","CD4 T Cells","CD8 T Cells", "NK Cells"), POI = 9, Lig_or_Rec = T,
                                             LR_List = Literature.Supported.Receptor.Ligand.Pairs.Ramilowski, Species = T, 
                                             Cutoff_Lig = -0.5, Cutoff_Rec = -0.5)

write.table(x = Endo_Lig_CD4_CD8_NK_Rec[[1]], file = "Endo_Lig_CD4_CD8_NK_Rec_Karyotype.txt", row.names = F, quote = F, col.names = F)
write.table(x = Endo_Lig_CD4_CD8_NK_Rec[[2]], file = "Endo_Lig_CD4_CD8_NK_Rec_Text.txt", row.names = F, quote = F, col.names = F)
write.table(x = Endo_Lig_CD4_CD8_NK_Rec[[3]], file = "Endo_Lig_CD4_CD8_NK_Rec_Expression.txt", row.names = F, quote = F, col.names = F)
write.table(x = Endo_Lig_CD4_CD8_NK_Rec[[4]], file = "Endo_Lig_CD4_CD8_NK_Rec_Links.txt", row.names = F, quote = F, col.names = F)

# Figure 6G 

Idents(TotalTissue.combined.harmony_fil) <- "Collapsed_labels"
Epi_Lig_CD4_CD8_NK_Rec <- CircosFunctions(InteractomeData = Full_list_LR_012720_lig_sort_collapsed_labels, SeuratObject = TotalTissue.combined.harmony_fil, 
                                             CellTypeVector = c("Epithelial","CD4 T Cells","CD8 T Cells", "NK Cells"), POI = 7, Lig_or_Rec = T,
                                             LR_List = Literature.Supported.Receptor.Ligand.Pairs.Ramilowski, Species = T, 
                                             Cutoff_Lig = -0.5, Cutoff_Rec = -0.5)

write.table(x = Epi_Lig_CD4_CD8_NK_Rec[[1]], file = "Epi_Lig_CD4_CD8_NK_Rec_Karyotype.txt", row.names = F, quote = F, col.names = F)
write.table(x = Epi_Lig_CD4_CD8_NK_Rec[[2]], file = "Epi_Lig_CD4_CD8_NK_Rec_Text.txt", row.names = F, quote = F, col.names = F)
write.table(x = Epi_Lig_CD4_CD8_NK_Rec[[3]], file = "Epi_Lig_CD4_CD8_NK_Rec_Expression.txt", row.names = F, quote = F, col.names = F)
write.table(x = Epi_Lig_CD4_CD8_NK_Rec[[4]], file = "Epi_Lig_CD4_CD8_NK_Rec_Links.txt", row.names = F, quote = F, col.names = F)

# Figure 6H 

Idents(CD8_T_cell_subset_fil) <- "DiseaseState"
VlnPlot(CD8_T_cell_subset_fil, features = c("TIGIT","CD96","PVRIG","CD226"), 
        pt.size = 0, ncol = 2)

# Figure 6I 

Idents(TotalTissue.combined.harmony_fil) <- "Collapsed_labels"
DotPlot(TotalTissue.combined.harmony_fil, 
        features = c('PVRIG', 'CD226','CD96','TIGIT', 'PVRL2', 'PVR'), 
        cols = c('red', 'blue'), dot.scale = 10) + RotatedAxis() + FontSize(x.text = 17, y.text = 17)

# Supplementary Figure S2B

Idents(TotalTissue.combined) <- "ID"
DimPlot(object = TotalTissue.combined, reduction = "umap", label = F, pt.size = 0.5)
Idents(TotalTissue.combined.harmony) <- "ID"
DimPlot(object = TotalTissue.combined.harmony, reduction = "umap", pt.size = 0.5, label = F)

# Supplementary Figure S2C/D

Idents(TotalTissue.combined.harmony_fil) <- "Collapsed_labels"
DimPlot(object = TotalTissue.combined.harmony_fil, reduction = "umap", split.by = "ID", pt.size = 0.5, label = F)

# Supplementary Figure S3A

DimPlot(object = PBMC_Merge_Fil, reduction = "umap", label = F, pt.size = 0.5, 
        cols = c("CD4 T Cells" = "limegreen", "CD8 T Cells" = "darkgreen", "Monocyte" = "darkorange2", "Granulocyte" = "orange1", 
                 "B Cells" = "gold1", "Plasma Cells" = "lightgoldenrod1", "pDCs" = "darkturquoise", "NK Cells" = "darkorchid"), split.by = "DiseaseState")


# Supplementary Figure S3B

Idents(PBMC_Merge_Fil) <- "Collapsed_labels"
DotPlot(PBMC_Merge_Fil, features = c( "CD4","FCGR3B","NKG7", "NCAM1", "CD8A", "CD3E",
                                      "ITGAM", "CD14" ,"IGLL5",'CD19', "CD79A","IRF7", "GZMB","PTPRC"), 
        cols = c("blue", "red"), dot.scale = 10) + RotatedAxis() + FontSize(x.text = 17, y.text = 17)

# Supplementary Figure S3C



# Supplementary Figure S3D



# Supplementary Figure S3E



# Supplementary Figure S4A

FeaturePlot(CD8_T_cell_subset_fil, c("CD3D","PDCD1","LAG3","TIGIT","HAVCR2",
                                     "IFNG","GZMB","GZMK","EOMES"), 
            cols = c('grey77','red'), pt.size = 1.5)

# Supplementary Figure S4C



# Supplementary Figure S5A

Idents(Myeloid_subset_fil) <- "Myeloid_subset_labels_expanded"
Myeloid_subset_fil <- RenameIdents(Myeloid_subset_fil,"Resident Macrophages_PDAC" = "1",
                                   "Resident Macrophages_AdjNorm" = "2", 
                                   "Classical Macrophages (CCR2+)_PDAC" = "3", 
                                   "Classical Macrophages (CCR2+)_AdjNorm" = "4", 
                                   "Alternatively Activated Macrophages_PDAC" = "5", 
                                   "Alternatively Activated Macrophages_AdjNorm" = "6", "Granulocyte 2_PDAC" = "7",                          
                                   "Granulocyte 2_AdjNorm" = "8",                       
                                   "Macrophage_PDAC" = "9",                            
                                   "Macrophage_AdjNorm" = "10", 
                                   "Granulocyte 1_PDAC" = "11", 
                                   "Granulocyte 1_AdjNorm" = "12")
Myeloid_subset_fil[["Myeloid_VlnPlot_Labels"]]<- Idents(object = Myeloid_subset_fil)

VlnPlot(Myeloid_subset_fil, features = c("CD274", "LGALS9","PVR","TNFSF4","SIRPA","HLA-DQA1",
                                         "ICOSLG","TNFSF18","CD80","CSF1R","CD40","CD70"), 
        cols = c("11" = "orange1", "12" = "orange", "3" = "steelblue1","4" = "steelblue", 
                 "1" = "darkgreen","2" = "green4", "5" = "deeppink2", "6" = "deeppink", 
                 "9"= "red3", "10" = "red", "7" = "gold2", "8" = "goldenrod"), pt.size = 0, y.max = 5, ncols = 2) 

# Supplementary Figure S5B



# Supplementary Figure S5C



# Supplementary Figure S5E

Idents(TotalTissue.combined.harmony_fil) <- "Collapsed_labels"
DotPlot(TotalTissue.combined.harmony_fil, 
        features = c('ADORA3', 'ADORA2B','ADORA2A','ADORA1'), 
        cols = c('red', 'blue'), dot.scale = 10) + RotatedAxis() + FontSize(x.text = 17, y.text = 17)



