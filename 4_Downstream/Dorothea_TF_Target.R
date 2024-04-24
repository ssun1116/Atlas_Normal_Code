rm(list = ls())
library(dplyr)

## TF collection
TF_collection = read.csv("Tables/CollecTRI_TF_collection.csv")
length(unique(TF_collection$source)) # 1186

TF_collection2 <- TF_collection %>%
  group_by(source) %>%
  summarise(target = paste(unique(target), collapse = ","))
colnames(TF_collection2) = c("names", "target")

## TF activity
TF_activity = read.csv("Tables/Dorothea_TF_activity_by_Leiden.csv")
TF_activity = TF_activity[,2:8]

TF_activity2 = merge(TF_activity, TF_collection2, by = "names")
TF_activity2 = TF_activity2[TF_activity2$pvals_adj < 0.05,]

## marker gene
marker = read.csv("Tables/DEG_rank_genes_groups_wilcoxon_pct.csv")

cols_per_group <- 4

# Split the dataframe into a list of dataframes with 4 columns each
df_list <- lapply(seq(1, ncol(marker), cols_per_group), function(i) {
  df_subset <- marker[, i:min(i + cols_per_group - 1, ncol(marker))]
  df_subset$cluster <- paste0("C", (i-1)/4)  # Add a new column named "cluster" with the index of the list
  return(df_subset)
})

# Ensure all data frames have the same column names
df_list <- lapply(df_list, function(df_subset) {
  colnames(df_subset) <- c("gene", "logfc", "FDR", "pct", "cluster")
  return(df_subset)
})

# Bind the dataframes row-wise
marker_df <- do.call(rbind, df_list)
marker_df = marker_df[marker_df$FDR < 0.05,]
marker_df = marker_df[abs(marker_df$logfc) > 0.1 & marker_df$pct > 0.01,]


### Test for C0
TF_activity2_test = TF_activity2[TF_activity2$group == "C0",]
TF_activity2_test = TF_activity2_test %>% top_n(., n = 10, meanchange)

marker_df_test = marker_df[marker_df$cluster == "C0",]

# Desired set of genes to retain
desired_genes <- marker_df_test$gene

# Function to filter genes based on the desired set
filter_genes <- function(gene_string, desired_genes) {
  genes <- unlist(strsplit(gene_string, ",", fixed = TRUE))
  filtered_genes <- genes[genes %in% desired_genes]
  return(paste(filtered_genes, collapse = ","))
}

# Apply the function to the "target" column
TF_activity2_test$target_marker <- sapply(TF_activity2_test$target, filter_genes, desired_genes = desired_genes)


### for all clusters
TF_activity_list = list()
for (cluster in unique(marker_df$cluster)){
  TF_activity2_tmp = TF_activity2[TF_activity2$group == cluster,]
  TF_activity2_tmp = TF_activity2_tmp %>% top_n(., n = 10, statistic)
  
  marker_df_tmp = marker_df[marker_df$cluster == cluster,]
  
  # Desired set of genes to retain
  desired_genes <- marker_df_tmp$gene
  
  # Function to filter genes based on the desired set
  filter_genes <- function(gene_string, desired_genes) {
    genes <- unlist(strsplit(gene_string, ",", fixed = TRUE))
    filtered_genes <- genes[genes %in% desired_genes]
    return(paste(filtered_genes, collapse = ","))
  }
  
  # Apply the function to the "target" column
  TF_activity2_tmp$target_marker <- sapply(TF_activity2_tmp$target, filter_genes, desired_genes = desired_genes)
  TF_activity_list[[cluster]] = TF_activity2_tmp
}

TF_activity_final = do.call(rbind, TF_activity_list)
TF_activity_final$group_names = paste(TF_activity_final$group, TF_activity_final$names, sep = "_")

TF_activity_final = TF_activity_final[,c(10, 1:9)]
TF_activity_final$group <- factor(TF_activity_final$group)
TF_activity_final = TF_activity_final %>% group_by(group) %>% arrange(-statistic, .by_group=TRUE)

write.csv(TF_activity_final, "Tables/Dorothea_TF_Leiden_Target_Marker.csv")
