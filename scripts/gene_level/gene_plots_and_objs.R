.libPaths('/scratch/midway2/qhauck/conda_env/rstudio-server/lib/R/library/')

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(readr)
  library(fastTopics)
  library(ggbio)
  library(patchwork)
  library(ggforce)
  library(htmlwidgets)
  library(plotly)
  library(nnls)
  library(data.table)
  library(patchwork)
  library(bit64)
  library(rtracklayer)
})

args <- commandArgs(TRUE)

#arguments for generating dat
min_junctions_for_pca <- as.numeric(args[1])
gene_range_bound <- as.numeric(args[2])
avg_clust_read_count_min <- as.numeric(args[3])
input_file <- args[4]
min_samples_per_junc <- as.numeric(args[5])
avg_sample_reads_per_cluster <- as.numeric(args[6])

#arguments for generating res from dat
gene_name <- args[7]
chr <- args[8]
gencode <- readRDS(args[9])
tcga_metadata <- fread(args[10], data.table = F)
threads <- as.numeric(args[11])
elbow_cutoff <- as.numeric(args[12])
max_factors <- as.numeric(args[13])
lfc_thresh <- as.numeric(args[14])
abs_thresh <- as.numeric(args[15])
sval_thresh <- as.numeric(args[16])


if(!dir.exists(paste0( chr, '/', gene_name))){
  dir.create(paste0( chr, '/', gene_name))
}

file.create(paste0( chr, '/', gene_name, '/de_res.RDS'))
file.create(paste0( chr, '/', gene_name, '/whole_factor.html'))
file.create(paste0( chr, '/', gene_name, '/res.RDS'))
file.create(paste0( chr, '/', gene_name, '/dat.tsv'))


gene_data <- gencode[gencode$gene_name == gene_name & gencode$type %in% c("exon")]
t_models <- rtracklayer::split(gene_data, gene_data$transcript_id)

lower_bound <- min(t_models@unlistData@ranges@start)
upper_bound <- max(as.numeric(t_models@unlistData@ranges@start) + as.numeric(t_models@unlistData@ranges@width))


if (file.size(input_file) < 100L){
  
  quit(save = 'no', status = 0)
  
}else{
  
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
    quit(save = 'no', status = 0)
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
      quit(save = 'no', status = 0)
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
          clusters[clusters == i] = clusters[index_of_matching_end]
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
        quit(save = 'no', status = 0)
      }else{
        
        junction_names <- paste0(filtered_intron_cluster$start, '-', filtered_intron_cluster$end)
        
        #remove samples that do not have an average of at least 5 reads per cluster
        samples_to_remove <- (colSums(clust_sums)/nrow(clust_sums) < avg_sample_reads_per_cluster)
        filtered_intron_cluster <- filtered_intron_cluster[,c(rep(F, 18),!as.vector(samples_to_remove), T)]
        
        
        #We may have no samples that fit this criterion
        if(ncol(filtered_intron_cluster) == 1){
          quit(save = 'no', status = 0)
        }else{
          
          #remove any junctions that in final dat have no support, not counting cluster column in RowSums
          
          junction_names <- junction_names[(rowSums(filtered_intron_cluster[1:(ncol(filtered_intron_cluster)-1)]) >  0)]
          filtered_intron_cluster <- filtered_intron_cluster[(rowSums(filtered_intron_cluster[1:(ncol(filtered_intron_cluster)-1)]) > 0),]
          
          #want something actually relevant to perform PCA on
          if (nrow(filtered_intron_cluster) < min_junctions_for_pca || ncol(filtered_intron_cluster)
              < min_junctions_for_pca){
            
            quit(save = 'no', status = 0)
            
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
            dat <- as.data.frame(dat)
            
            colnames(dat) <- junction_names
          }
        }
      }
    }
  }
}
write.table(dat, paste0(chr, '/', gene_name, '/dat.tsv'), row.names = T)

dat <- as.matrix(dat)

if(ncol(dat) < 3 || nrow(dat) < 1000){
  quit(save = 'no', status = 0)
}

palet <- c(
  'red', 'blue', 'green', 'purple', 'orange', 'cyan2', 'goldenrod4', 'pink', 
       'magenta', 'grey', 'brown', 'darkolivegreen',
       'skyblue2', 'yellow', 'chartreuse4', 'black'
)

palet_labels <- c(
  '1'='red', '2'='blue', '3'='green', '4'='purple', '5'='orange', 
  '6'='cyan2', '7'='goldenrod4', '8'='pink', '9'='magenta', '10'='grey', 
  '11'='brown', '12'='darkolivegreen', 
  '13'='skyblue2', '14'='yellow', '15'='chartreuse4', '16'='black'
)


# helper functions --------------------------------------------------------
round_down_to_thousand <- function(value) {
  # Divide the value by 1000, apply floor, and multiply back by 1000
  rounded_value <- floor(value / 1000) * 1000
  return(rounded_value)
}


round_up_to_thousand <- function(value) {
  # Divide the value by 1000, apply floor, and multiply back by 1000
  rounded_value <- ceiling(value / 1000) * 1000
  return(rounded_value)
}

best_fit_L <- function(results, data_to_fit){
  F_mat <- results$F
  
  L <- matrix(0, nrow = nrow(data_to_fit), ncol = ncol(F_mat))
  
  # Solve NNLS for each sample
  for (i in 1:nrow(data_to_fit)) {
    fit <- nnls(F_mat, data_to_fit[i, ])
    L[i, ] <- fit$x
  }
  return(L)
}



# data processing ---------------------------------------------------------


#choosing number of factors to use for visualization
normalized_junctions <- dat/rowSums(dat)
pc_for_elbow <- prcomp(normalized_junctions)
proportion_of_var <- summary(pc_for_elbow)$importance[2,]
proportion_of_var <- as.numeric(proportion_of_var)
cumul_var <- cumsum(proportion_of_var)

ideal_factors <- length(proportion_of_var) - sum(cumul_var > elbow_cutoff)

#CHOOSING NUMBER OF FACTORS
num_factors <- max(2, min(max_factors, ideal_factors))

#setting up res objects and saving them
new_res <- poisson2multinom(fit_poisson_nmf(dat, num_factors, control = list('nc' = threads)))

de_res <- de_analysis(new_res, dat, control = list('nc' = threads))
s_vals <- de_res$svalue
s_vals[is.na(s_vals)] <- 1
lfc <- replace_na(de_res$postmean, 0)
lfc[is.na(lfc)] <- 0

saveRDS(new_res, paste0( chr, '/', gene_name, '/res.RDS'))
saveRDS(de_res, paste0( chr, '/', gene_name, '/de_res.RDS'))


gene_data <- gencode[gencode$gene_name == gene_name & gencode$type %in% c("exon")]
t_models <- rtracklayer::split(gene_data, gene_data$transcript_id)

#setting up plotting for curves to be added to
min_coord <- round_down_to_thousand(min(c(start(gene_data@ranges), end(gene_data@ranges))))
max_coord <- round_up_to_thousand(max(c(start(gene_data@ranges), end(gene_data@ranges))))
p <- ggplot() + ggbio::geom_alignment(t_models, gap.geom = 'arrow', label = F, exon.rect.h = 0.4) + theme_minimal()+
  scale_x_continuous(breaks = seq(min_coord, max_coord, by = 10000))

build <- ggplot_build(p)
max_y <- max(build$layout$panel_params[[1]]$y.range)

whole_line_df <- data.frame(x = numeric(), y=numeric(), xend = numeric(), yend = numeric(), width = numeric(), me = character())
line_df <- data.frame(x = numeric(), y=numeric(), xend = numeric(), yend = numeric(), width = numeric(), factor = character())

cancer_types <- tcga_metadata[tcga_metadata$rail_id %in% rownames(dat),]$gdc_cases.project.project_id

# num_clusters <- 5
# avg_clust_sample_min <- 5
# #filtering out samples with low expression
# low_expression_rail_ids <- rownames(dat[rowSums(dat) < num_clusters*avg_clust_sample_min, ])
# cancer_types[rownames(dat) %in% low_expression_rail_ids] = 'Z-LOW EXPRESSION'

struc_plot <- structure_plot(new_res, grouping = cancer_types, gap = 100, colors = palet[1:num_factors], n = 5000)+
  theme(aspect.ratio = 1/5)+
  labs(title = paste0('TCGA Factorization of GENE: ', gene_name))


factor_info <- data.frame(new_res$F, rownames(new_res$F))
colnames(factor_info) <- c(1:num_factors, 'junction')

factor_plot <- factor_info %>%
  mutate(start = as.numeric(str_split_i(junction, '-', 1))) %>%
  mutate(end = as.numeric(str_split_i(junction, '-', 2)))

k=1
#iterate through factors and add junctions from them as curves to line_df, only if they pass DE thresholds
for (j in 1:num_factors){
  this_whole_line_df <- data.frame(x = factor_plot$start, y = max_y+5*j, xend = factor_plot$end, yend = max_y+5*j, width = factor_plot[,j], factor = as.character(j)) %>%
    filter(width > abs_thresh/10)
  
  whole_line_df <- rbind(whole_line_df, this_whole_line_df)
  
  #isolate only the DE junctions that are significant and positively upregulated relative to all other factors
  max_lfc <- apply(lfc, 1, max)
  
  this_line_df <- data.frame(x = factor_plot$start, y = max_y+5*k, xend = factor_plot$end, yend = max_y+5*k, width = factor_plot[,j], factor = as.character(j)) %>%
    cbind(s_vals) %>%
    cbind(lfc) %>%
    .[.[,6+j+num_factors] > lfc_thresh & .[,6+j+num_factors]/max_lfc == 1,] %>%
    .[.[, 'width'] > abs_thresh,] %>%
    .[.[,6+j] < sval_thresh,] %>%
    select(c(x,y,xend,yend,width,factor))
  
  if(nrow(this_line_df) != 0){
    k = k+1
  }
  
  line_df <- rbind(line_df, this_line_df)
}

#add curves to plot and plot factors along w/structure plot
whole_p <- p +
  geom_curve(whole_line_df, mapping = aes(x=x, y=y, xend=xend, yend=yend, linewidth = width, colour = factor), curvature = -0.3)


whole_factor_plot <- whole_p+ylim(0, 25+5*length(palet)+5) + guides(linewidth = 'none') + scale_color_manual(values = palet_labels[1:num_factors])

whole_combined_plot <- struc_plot / whole_factor_plot

#print(whole_combined_plot)


de_p <- p +
  geom_curve(line_df, mapping = aes(x=x, y=y, xend=xend, yend=yend, linewidth = width, colour = factor), curvature = -0.3)

factors_included <- unique(line_df$factor)

de_factor_plot <- de_p+ylim(0, max_y+5*length(factors_included)+5) + guides(linewidth = 'none') + scale_color_manual(values = palet_labels[factors_included])

de_combined_plot <- struc_plot / de_factor_plot


##WHOLE FACTOR PLOT
running_whole_plot <- p
height_scaling=0.7
for (i in 1:nrow(whole_line_df)){
  y_base <- max_y + as.numeric(whole_line_df$factor[i])*4
  x_vals <- seq(whole_line_df$x[i], whole_line_df$xend[i], length.out = 100)
  y_top <- sin((x_vals - whole_line_df$x[i]) / (whole_line_df$xend[i] - whole_line_df$x[i]) * pi)*height_scaling + y_base + whole_line_df$width[i]/max(whole_line_df$width)*3.5
  
  # Generate the bottom arc (flat or another function)
  y_bottom <- sin((x_vals - whole_line_df$x[i]) / (whole_line_df$xend[i] - whole_line_df$x[i]) * pi)*height_scaling + y_base
  
  # Combine into a polygon (top arc + reversed bottom arc to close the shape)
  polygon_df <- data.frame(
    x = c(x_vals, rev(x_vals)),
    y = c(y_top, rev(y_bottom))
  )
  
  running_whole_plot <- running_whole_plot + geom_polygon(data = polygon_df, aes(x = x, y = y), fill = palet[as.numeric(whole_line_df$factor[i])], alpha = 1)
}

running_whole_plot <- running_whole_plot + scale_y_continuous(breaks = c(0,10), labels = c('',''))
ggp_whole_factor <- ggplotly(running_whole_plot, dynamicTicks = T)
polygon_traces <- 4:(3 + nrow(whole_line_df))
whole_line_df$label <- paste0(whole_line_df$x, "-", whole_line_df$xend)

# Apply custom text to each trace
ggp_whole_factor <- purrr::reduce(
  .x = seq_along(polygon_traces),
  .f = function(plotly_obj, i) {
    style(
      plotly_obj,
      hoverinfo = "text",
      text = whole_line_df$label[i],  # Assign ONE label per trace
      traces = polygon_traces[i]     # Target one trace at a time
    )
  },
  .init = ggp_whole_factor   # Your original plotly object
)



running_de_plot <- p

#establishing arrow for graph before any graphing could occur
strand = as.character(t_models@unlistData@strand[1])
if (strand == '+'){
  arrow_start <- min_coord
  arrow_end <- min_coord + (max_coord-min_coord)/10
}else{
  arrow_start <- max_coord
  arrow_end <- max_coord - (max_coord-min_coord)/10
}

ggp_struc <- ggplotly(struc_plot)

#DE FACTOR PLOT
if(nrow(line_df) != 0){
  for (i in 1:nrow(line_df)){
    y_base <- max_y + as.numeric(line_df$factor[i])*4
    x_vals <- seq(line_df$x[i], line_df$xend[i], length.out = 100)
    y_top <- sin((x_vals - line_df$x[i]) / (line_df$xend[i] - line_df$x[i]) * pi)*height_scaling + y_base + line_df$width[i]/max(whole_line_df$width)*3.5
    
    # Generate the bottom arc (flat or another function)
    y_bottom <- sin((x_vals - line_df$x[i]) / (line_df$xend[i] - line_df$x[i]) * pi)*height_scaling + y_base
    
    # Combine into a polygon (top arc + reversed bottom arc to close the shape)
    polygon_df <- data.frame(
      x = c(x_vals, rev(x_vals)),
      y = c(y_top, rev(y_bottom))
    )
    
    label = paste0(line_df$x[i], '-', line_df$xend[i])
    running_de_plot <- running_de_plot + geom_polygon(data = polygon_df, aes(x = x, y = y), fill = palet[as.numeric(line_df$factor[i])], alpha = 1)
  }
  
  running_de_plot <- running_de_plot + scale_y_continuous(breaks = c(0,10), labels = c('',''))
  ggp_de_factor <- ggplotly(running_de_plot, dynamicTicks = T)
  polygon_traces <- 4:(3 + nrow(line_df))
  line_df$label <- paste0(line_df$x, "-", line_df$xend)
  
  # Apply custom text to each trace
  ggp_de_factor <- purrr::reduce(
    .x = seq_along(polygon_traces),
    .f = function(plotly_obj, i) {
      style(
        plotly_obj,
        hoverinfo = "text",
        text = line_df$label[i],  # Assign ONE label per trace
        traces = polygon_traces[i]     # Target one trace at a time
      )
    },
    .init = ggp_de_factor  # Your original plotly object
  )
  
  
  #Adding Gene Direction arrow to plot
  ggp_de_factor <- ggp_de_factor  %>%
    add_annotations(
      x = arrow_end,        # Arrow end point (xend)
      y = -1.4,              # Arrow end point (yend)
      ax = arrow_start,      # Arrow start point (x)
      ay = -1.4,
      xref = "x2",
      yref = "y2",
      axref = "x2",
      ayref = "y2",
      text = "",           # No label on arrow
      showarrow = TRUE,
      arrowhead = 10,       # Style (1-7, 2=default arrow)
      arrowsize = 3,     # Smaller arrow size
      arrowwidth = 0.5,    # Thinner arrow line
      arrowcolor = "black" # Color
    ) %>%
    # Add text label above arrow
    add_annotations(
      x = (arrow_start + arrow_end)/2,  # Midpoint of arrow
      y = -0.3,                       # Position above arrow
      text = "gene direction",
      showarrow = FALSE,
      font = list(size = 9),         # Adjust font size
      xanchor = "center"              # Center-align text
    )
  
  htmlwidgets::saveWidget(subplot(ggp_struc, ggp_de_factor, nrows = 2, margin = 0.07), paste0( chr, '/', gene_name, '/de_factor.html'))
  
}
ggp_whole_factor <- ggp_whole_factor  %>%
  add_annotations(
    x = arrow_end,        # Arrow end point (xend)
    y = -1.4,              # Arrow end point (yend)
    ax = arrow_start,      # Arrow start point (x)
    ay = -1.4,
    xref = "x2",
    yref = "y2",
    axref = "x2",
    ayref = "y2",
    text = "",           # No label on arrow
    showarrow = TRUE,
    arrowhead = 10,       # Style (1-7, 2=default arrow)
    arrowsize = 3,     # Smaller arrow size
    arrowwidth = 0.5,    # Thinner arrow line
    arrowcolor = "black" # Color
  ) %>%
  # Add text label above arrow
  add_annotations(
    x = (arrow_start + arrow_end)/2,  # Midpoint of arrow
    y = -0.3,                       # Position above arrow
    text = "gene direction",
    showarrow = FALSE,
    font = list(size = 9),         # Adjust font size
    xanchor = "center"              # Center-align text
  )



###{r just TCGA}
htmlwidgets::saveWidget(subplot(ggp_struc, ggp_whole_factor, nrows = 2, margin = 0.07), paste0( chr, '/', gene_name, '/whole_factor.html'))
unlink(paste0( chr, '/', gene_name, '/whole_factor_files/'), recursive = TRUE)
unlink(paste0( chr, '/', gene_name, '/de_factor_files/'), recursive = TRUE)
