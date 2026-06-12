.libPaths('/scratch/midway2/qhauck/conda_env/rstudio-server/lib/R/library/')
library(tidyverse)
library(readr)
library(fastTopics)
library(nnls)
library(GenomicRanges)

args <- commandArgs(TRUE)

arg_names <- c(
  "input_file", "meta_file", "threads", "min_samples_per_junc",
  "outdir", "outprefix", "elbow_cutoff", "max_factors",
  "avg_clust_read_count_min", "gencode_rds", "avg_sample_reads_per_cluster",
  "gene_range_bound"
)
if (length(args) != length(arg_names)) {
  stop(sprintf(
    "tcga_LF_saving.R expected %d arguments (%s); got %d.",
    length(arg_names), paste(arg_names, collapse = ", "), length(args)
  ))
}

input_file <- args[1]
meta_file <- args[2]
threads <- as.numeric(args[3])
min_samples_per_junc <- as.numeric(args[4])
outdir <- args[5]
outprefix <- args[6]
elbow_cutoff <- as.numeric(args[7])
max_factors <- as.numeric(args[8])
avg_clust_read_count_min <- as.numeric(args[9])

#to use to remove junctions that are far away from annotated gene model
gencode <- readRDS(args[10])
avg_sample_reads_per_cluster <- as.numeric(args[11])
# maximum allowed junction coordinate away from gencode established gene boundaries
gene_range_bound <- as.numeric(args[12])

numeric_args <- c(
  threads = threads,
  min_samples_per_junc = min_samples_per_junc,
  elbow_cutoff = elbow_cutoff,
  max_factors = max_factors,
  avg_clust_read_count_min = avg_clust_read_count_min,
  avg_sample_reads_per_cluster = avg_sample_reads_per_cluster,
  gene_range_bound = gene_range_bound
)
if (any(is.na(numeric_args))) {
  stop("Failed to parse numeric arguments: ",
       paste(names(numeric_args)[is.na(numeric_args)], collapse = ", "))
}


# input_file <- 'all_genes/chrX/snaptron_output/SPIN2B_snaptron.tsv'
# meta_file <- 'tcga_v2_samples.tsv'
# threads <- 22
# min_samples_per_junc <- 5
# outdir <- 'testing/'
# outprefix <- 'SPIN2B'
# elbow_cutoff <- 0.0001
# max_factors <- 32
# avg_clust_read_count_min <- 10
# gene_range_bound <- 2000
# #to use to remove junctions that are far away from annotated gene model
# gencode <- readRDS('hg38_granges.RDS')
# avg_sample_reads_per_cluster <- 5


gene_name <- str_split_i(basename(input_file), '_', 1)
gene_data <- gencode[gencode$gene_name == gene_name & gencode$type %in% c("exon")]
t_models <- rtracklayer::split(gene_data, gene_data$transcript_id)

lower_bound <- min(t_models@unlistData@ranges@start)
upper_bound <- max(as.numeric(t_models@unlistData@ranges@start) + as.numeric(t_models@unlistData@ranges@width))

if(!dir.exists(outdir)){
  dir.create(outdir)
}

if (!dir.exists(paste0(outdir, outprefix))){
  dir.create(paste0(outdir, outprefix))
}

if (file.size(input_file) < 100L){
  #do nothing
  # L_output_path <- paste0(outdir, '/', outprefix, '/', 'L.tsv')
  # write_tsv(data.frame(), L_output_path)
  #
  # F_output_path <- paste0(outdir, '/', outprefix, '/', 'F.tsv')
  # write_tsv(data.frame(), F_output_path)

  message("Skipping ", outprefix, ": input file < 100 bytes")
  saveRDS(list(), paste0(outdir, outprefix, '/', 'res.RDS'))

}else{
  
  tcga_metadata <- read_tsv(meta_file)
  
  min_junctions_for_pca <- 3
  
  gene_table <- read_tsv(input_file)
  
  #at least one gene has duplicates for some reason
  #Get rid of any junctions that probably don't actually correspond to the gene
  #get rid of junctions on opposite strand
  gene_table <- gene_table %>%  
    filter(samples_count > min_samples_per_junc) %>%
    unique() %>%
    filter(start > (lower_bound - gene_range_bound)) %>%
    filter(end < (upper_bound + gene_range_bound)) %>% 
    filter(strand %in% unique(as.character(strand(gene_data))))
  
  normal_rail_ids <- filter(tcga_metadata, grepl('Normal', cgc_sample_sample_type))$rail_id
  
  if(nrow(gene_table) < min_junctions_for_pca){
    message("Skipping ", outprefix, ": only ", nrow(gene_table),
            " junctions left after initial sample/strand/range filter (< ",
            min_junctions_for_pca, ")")
    saveRDS(list(), paste0(outdir, outprefix, '/', 'res.RDS'))
  }else{

    #functions to parse the snaptron samples column
    split_to_ids <- function(tem){
      return(str_split(tem, ':')[[1]][1])
    }
    split_to_counts <- function(tem){
      return(str_split(tem, ':')[[1]][2])
    }
    
    #setup matrix of samples by junction counts
    just_counts <- matrix(0, nrow = nrow(gene_table), ncol = length(tcga_metadata$rail_id))
    colnames(just_counts) <- tcga_metadata$rail_id
    rownames(just_counts) <- as.character(1:nrow(just_counts))
    
    #updating just_counts with counts corresponding to junctions of interest
    for (i in 1:nrow(gene_table)){
      pre_split <- str_split(gene_table$samples[i], ',')
      pre_split <- pre_split[[1]][2:length(pre_split[[1]])]
      
      rail_ids <- as.vector(sapply(pre_split, split_to_ids))
      counts <- as.vector(sapply(pre_split, split_to_counts))
      
      just_counts[cbind(i, rail_ids)] <- counts
    }
    
    #ensure proper behavior of matrix as numeric
    colnames_to_keep <- colnames(just_counts)
    just_counts <- matrix(as.numeric(just_counts), nrow = nrow(just_counts), ncol = ncol(just_counts))
    colnames(just_counts) <- colnames_to_keep
    
    #as.matrix needed to avoid 1 row issue where matrix is just numeric
    #deprecated as we don't care about 1 row matrices here anyway
    just_counts <- just_counts[,!colnames(just_counts) %in% normal_rail_ids]
    
    #filter out junctions that previously passed due to many normal samples having counts of them
    gene_table <- gene_table[!rowSums(just_counts > 0) < min_samples_per_junc, ]
    just_counts <- just_counts[!rowSums(just_counts > 0) < min_samples_per_junc, ]
    
    #now chance of having too few junctions
    if(nrow(gene_table) < min_junctions_for_pca){
      message("Skipping ", outprefix, ": only ", nrow(gene_table),
              " junctions left after dropping normal samples + low-support junctions (< ",
              min_junctions_for_pca, ")")
      saveRDS(list(), paste0(outdir, outprefix, '/', 'res.RDS'))
    }else{

      junction_by_count <- cbind(gene_table, just_counts)
      
      
      # leafcutter-esque clustering ------------------------------------------------------
      
      clustering <- junction_by_count %>%
        group_by(start) %>%
        mutate(ID = cur_group_id())
      
      
      clusters <- clustering$ID
      
      for (i in 2:max(clusters)){
        ends <- clustering$end[clusters==i]
        if (sum(ends %in% clustering$end[clusters < i]) > 0){
          index_of_matching_end <- c(1:length(clusters))[(clustering$end %in% ends & clustering$ID < i)][1]
          clusters[clusters == i] <- clusters[index_of_matching_end]
        }
      }
      
      junction_by_count$cluster <- clusters
      
      intron_cluster <- junction_by_count %>%
        group_by(cluster)%>%
        mutate(across(colnames(just_counts), ~sum(.x), .names = "{.col}_clust_sum"))
      
      #end of leafcutter-esque clustering--------------------------------------------
      
      #remove clusters that do not have an average of at least x reads 
      #retain only junction info columns, counts and cluster column
      clust_sums <- intron_cluster[,((ncol(intron_cluster)-1)-ncol(just_counts)+1):ncol(intron_cluster)] %>%
        unique()
      #kept cluster to avoid throwing away clusters w/same dist of junctions (unlikely), now need to remove it
      clust_sums <- subset(clust_sums, select = -cluster)
      clust_to_remove <- (rowSums(clust_sums)/ncol(clust_sums) < avg_clust_read_count_min)
      filtered_intron_cluster <- intron_cluster %>% 
        filter(!cluster %in% unique(clusters)[clust_to_remove]) %>% 
        .[,(1:ncol(junction_by_count))]
      
      #remove any very low coverage clusters prior to removing low coverage samples
      clust_sums <- clust_sums[!clust_to_remove, ]
      
      #need to filter out genes that all clusters are thrown away
      if(nrow(clust_sums) == 0){
        message("Skipping ", outprefix,
                ": no clusters passed avg_clust_read_count_min = ",
                avg_clust_read_count_min)
        saveRDS(list(), paste0(outdir, outprefix, '/', 'res.RDS'))
      }else{
        
        junction_names <- paste0(filtered_intron_cluster$start, '-', filtered_intron_cluster$end)
        
        #remove samples that do not have an average of at least 5 reads per cluster
        samples_to_remove <- (colSums(clust_sums)/nrow(clust_sums) < avg_sample_reads_per_cluster)
        filtered_intron_cluster <- filtered_intron_cluster[,c(rep(F, 18),!as.vector(samples_to_remove), T)]
        
        
        #We may have no samples that fit this criterion
        if(ncol(filtered_intron_cluster) == 1){
          message("Skipping ", outprefix,
                  ": no samples passed avg_sample_reads_per_cluster = ",
                  avg_sample_reads_per_cluster)
          saveRDS(list(), paste0(outdir, outprefix, '/', 'res.RDS'))
        }else{
          
          #remove any junctions that in final dat have no support, not counting cluster column in RowSums
          
          junction_names <- junction_names[(rowSums(filtered_intron_cluster[1:(ncol(filtered_intron_cluster)-1)]) >  0)]
          filtered_intron_cluster <- filtered_intron_cluster[(rowSums(filtered_intron_cluster[1:(ncol(filtered_intron_cluster)-1)]) > 0),]
          
          #want something actually relevant to perform PCA on
          if (nrow(filtered_intron_cluster) < min_junctions_for_pca || ncol(filtered_intron_cluster)
              < min_junctions_for_pca){

            message("Skipping ", outprefix, ": after final rowSum filter, ",
                    nrow(filtered_intron_cluster), " junctions x ",
                    ncol(filtered_intron_cluster), " cols (< ",
                    min_junctions_for_pca, ")")
            saveRDS(list(), paste0(outdir, outprefix, '/', 'res.RDS'))

          }else{
            
            #want samples to at least have a read, highly unlikely this happens but possible if
            #a sample only has reads in an unsupported junction
            filtered_intron_cluster <- filtered_intron_cluster[, colSums(filtered_intron_cluster) > 0]
            
            #normalize junctions by total counts to 'equally weight' samples in PCA
            normalized_junctions <- t(t(filtered_intron_cluster[1:(ncol(filtered_intron_cluster)-1)])/colSums(filtered_intron_cluster[1:(ncol(filtered_intron_cluster)-1)]))
            
            #use number of principal components to reach a certain explained variance to bound how large of a k to use
            pc_for_elbow <- prcomp(t(normalized_junctions))
            proportion_of_var <- summary(pc_for_elbow)$importance[2,]
            proportion_of_var <- as.numeric(proportion_of_var)
            cumul_var <- cumsum(proportion_of_var)
            
            ideal_factors <- length(proportion_of_var) - sum(cumul_var > (1-elbow_cutoff))
            
            #CHOOSING NUMBER OF FACTORS
            num_factors <- max(2, min(max_factors, ideal_factors))
            
            dat <- t(as.matrix(filtered_intron_cluster[,1:(ncol(filtered_intron_cluster)-1)]))
            res <- fit_poisson_nmf(dat, num_factors, control = list('nc' = threads))
            
            rownames(res$F) <- junction_names
            
            
            dat <- as.data.frame(dat)
            colnames(dat) <- junction_names
            saveRDS(res, paste0(outdir, outprefix, '/', 'res.RDS'))
            write_tsv(dat, paste0(outdir, outprefix, '/', 'dat.tsv'))
          }
        }
      }
    }
  }
}