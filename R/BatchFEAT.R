require(phyloseq)
require(rhdf5)
require(biom)
require(shiny)
require(dplyr)
require(DT)
require(ggplot2)
require(vegan)
require(lazyeval)
ggsave <- ggplot2::ggsave; body(ggsave) <- body(ggplot2::ggsave)[-2]
options(digits = 4)


#' Command-line FEAT function to perform analysis in batch for a given experiment.
#' @param biom_table BIOM table
#' @param mapping_file Mapping File
#' @param FMT_pairs tsv indicating FMT pairs (see example)
#' @param input_params input parameters in tsv
#' @param output_dir output directory to leave output files
#' @param ... taxonomy file (optional)
#' @return many files containing QC and summary statistics, as well as FMT visualizations
#' @export
#'
batchFEAT <- function(biom_file, mapping_file, FMT_pairs, input_params, output_dir, taxonomy=NULL) {
  start_time <- proc.time()
  cat("----Batch FMT Efficacy Analysis Toolkit-----\n")
  cat("***Starting***\n")
  cat("Loading OTU Table & Mapping file...\n")
  # Load table + mapping
  mapping <- read.delim(file = mapping_file)
  mapping$X.SampleID <- as.character(mapping$X.SampleID)
  biom_table <- biom::read_biom(biom_file)
  table_data <- as.matrix(biom_data(biom_table))
  biom_only <- as.data.frame(t(table_data))
  raw_table_otu_count <- ncol(biom_only)
  colnames(biom_only) <- paste("OTU", colnames(biom_only), sep = "_")
  biom_only$X.SampleID <- as.character(row.names(biom_only))
  row.names(biom_only) <- NULL

  cat("OTU Table & Mapping file loaded...\n")
  
  if(!is.null(taxonomy)){
	  cat("Taxonomy provided, loading Taxonomy map...\n")
	tax_map <- read.delim(taxonomy, header = FALSE)
	colnames(tax_map) <- c("OTU", "Taxon", "Quality_Score")
	tax_map$OTU <- paste("OTU", tax_map$OTU, sep = "_")
	
	tax_check <- TRUE
	tax_file <- taxonomy

	cat("Taxonomy map loaded successfully...\n")
	}
	else if (is.null(taxonomy)){
		tax_check <- FALSE
		tax_file <- 'N/A'
	}
  
  cat("Adding Metadata to OTU table...\n")

  metadata_table <- merge(biom_only, mapping, by = 'X.SampleID')

  cat("Metadata added to OTU table... \n")

  # Get total number of OTUs to start
  n_otus_starting <- raw_table_otu_count
  cat("# of OTUs to start:", n_otus_starting, '\n')

  # Load Parameters
  parameters <- read.delim(file = input_params)
  metadata_category <- as.character(parameters[parameters$Parameter == 'metadata', 'Value'])
  min_fraction <- as.numeric(as.character(parameters[parameters$Parameter == "min_fraction", "Value"]))
  min_abundance <- as.numeric(as.character(parameters[parameters$Parameter == "min_abundance", "Value"]))
  comparison <- as.character(parameters[parameters$Parameter == "comparison", "Value"])


  transplant_pairs <- read.delim(file = FMT_pairs, header = TRUE)

  table_out_items <- c("Comment", "TransplantID", "Taxonomy_Added", "Taxonomy_Filename", "OTU_Table_Filename", "N_OTUs_starting",
                       "Mapping_Filename", "Metadata_Category_FMT",
                       "DonorID", "N_Donor_Samples", "RecipientID", "N_Recipient_Samples", "Post_FMT_ID", "N_Post_FMT", "Rel_Abundance_Filter", "N_OTUs_Relative_Abundance_Filtered",
                       "Fleeting_OTUs_Filter", "N_OTUs_Fleeting_Filtered", "Comparison_Metric",
                       "N_don", "N_rec", "N_post_fmt", "N_post_fmt_excl",
                       "N_post_fmt_unique", "N_shared",
                       "N_P_Donor", "D_Engraft", "P_Donor",
                       "N_P_Recipient", "R_Persist", "P_Recipient", "P_Shared", "P_Unique",
                       "D_Engraft_excl", "P_Donor_excl",
                       "R_Persist_excl", "P_Recipient_excl")
					   
  slim_table_out_items <- c("Comment", "TransplantID", "D_Engraft", "R_Persist", "P_Donor", "P_Recipient", "P_Shared", "P_Unique")
  slim_table_out <- data.frame(t(slim_table_out_items))
  colnames(slim_table_out) <- make.names(slim_table_out_items)
  slim_table_out <- slim_table_out[-1,]
  
  table_out <- data.frame(t(table_out_items))
  colnames(table_out) <- make.names(table_out_items)
  table_out <- table_out[-1, ]
  dir.create(file.path(getwd(), output_dir), showWarnings = FALSE)
  setwd(file.path(getwd(), output_dir))

  for (i in 1:nrow(transplant_pairs)) {
    donor <- as.character(transplant_pairs[i, "Donor"])
    recipient <- as.character(transplant_pairs[i, "Recipient"])
    post_fmt <- as.character(transplant_pairs[i, "Post_FMT"])
    comment <- as.character(transplant_pairs[i, "Comment"])
	output_name <- paste(post_fmt, comment, sep = "_")
    dir.create(file.path(getwd(), output_name), showWarnings = FALSE)
    setwd(file.path(getwd(), output_name))
    cat("-----Start Iteration (", i, "of", nrow(transplant_pairs),") -----\n")
    transplant_id <- paste(donor, " into ", recipient, " (", post_fmt,")",sep = "")
	
    cat("Performing analysis for the following transplant:", transplant_id, "--", comment, "--\n")
    cat("Generating transplant-specific table\n")
    experiment_specific <- split_into_experiment_batch(metadata_table, mapping, metadata_category, donor, recipient, post_fmt)

    n_donor_samples <- nrow(experiment_specific[experiment_specific[,metadata_category] == donor, ])
    n_recipient_samples <- nrow(experiment_specific[experiment_specific[,metadata_category] == recipient, ])
    n_post_fmt_samples <- nrow(experiment_specific[experiment_specific[,metadata_category] == post_fmt, ])
    cat("# Donor Samples:", n_donor_samples, "\n")
    cat("# Recipient Samples:", n_recipient_samples, "\n")
    cat("# Post-FMT Samples:", n_post_fmt_samples, "\n")

    cat("Normalizing transplant-specific table\n")
    cat("Filtering OTUs that do not have a relative abundance of at least", min_abundance, "...\n")

    normalized_filtered_list <- normalize_and_filter_batch(experiment_specific, min_abundance)

    unnormalized_experiment_specific <- normalized_filtered_list[[2]] # For QC
    normalized_unfiltered <- normalized_filtered_list[[3]] # For QC
    normalized_filtered <- normalized_filtered_list[[1]]


    dir.create(file.path(getwd(), 'Tables'), showWarnings = FALSE)
    setwd(file.path(getwd(), 'Tables'))

    full_table <- inner_join(normalized_filtered, mapping, by = "X.SampleID")
    write.csv(full_table, paste(output_name, "relative_abundance_filtered_experiment_specific_table.csv", sep = "_"), row.names = FALSE)
    n_otus_after_relative_filter <- ncol(full_table) - ncol(mapping[,c(-1,-2,-3)])

    cat("# of OTUs that pass this filter:", n_otus_after_relative_filter, '\n')

    donor_only_table <- full_table[full_table[,metadata_category] == donor,]
    write.csv(donor_only_table, paste(output_name, "donor_only_table.csv", sep = "_"), row.names = FALSE)
    recipient_only_table <- full_table[full_table[,metadata_category] == recipient,]
    write.csv(recipient_only_table, paste(output_name, "recipient_only_table.csv", sep = "_"), row.names = FALSE)
    post_fmt_only_table <- full_table[full_table[,metadata_category] == post_fmt,]
    write.csv(post_fmt_only_table, paste(output_name, "post_fmt_only_table.csv", sep = "_"), row.names = FALSE)
	
	
    cat("Filtering OTUs that are fleeting (are present in less than", percent(min_fraction), "of samples of each group)\n")
    donor_nonzero_table <- filter_function(nonzero_fraction_getter(donor_only_table), min_fraction, donor)
    recipient_nonzero_table <- filter_function(nonzero_fraction_getter(recipient_only_table), min_fraction, recipient)
    post_fmt_nonzero_table <- filter_function(nonzero_fraction_getter(post_fmt_only_table), min_fraction, post_fmt)
	
	if(!is.null(taxonomy)) {
		cat("Collapsing OTUs into Taxa...\n")
		donor_nonzero_table <- collapse_taxonomy(inner_join(donor_nonzero_table, tax_map, by = "OTU"))
		recipient_nonzero_table <- collapse_taxonomy(inner_join(recipient_nonzero_table, tax_map, by = "OTU"))
		post_fmt_nonzero_table <- collapse_taxonomy(inner_join(post_fmt_nonzero_table, tax_map, by = "OTU"))
		
		full_nonzero <- bind_rows(bind_rows(donor_nonzero_table, recipient_nonzero_table), post_fmt_nonzero_table)
		
		distinct_otus_nonzero <- distinct(as.data.frame(full_nonzero$Taxon))
		
	} else if (is.null(taxonomy)) {
		full_nonzero <- bind_rows(bind_rows(donor_nonzero_table, recipient_nonzero_table), post_fmt_nonzero_table)
		
		distinct_otus_nonzero <- distinct(as.data.frame(full_nonzero$OTU))
	}


    write.csv(full_nonzero, paste(output_name, "full_nonzero_table.csv", sep = "_"), row.names = FALSE)
    n_otus_after_fleeting_filter <- nrow(distinct_otus_nonzero)
    cat("# of OTUs/Taxa that remain after this filter:", n_otus_after_fleeting_filter, "\n")

    cat("Writing tables...\n")

	if(!is.null(taxonomy)){
		joining_category <- "Taxon"
	} else if (is.null(taxonomy)) {
		joining_category <- "OTU"
	}

    # Make Tables of OTUs unique to each condition
    donor_unique <- anti_join(donor_nonzero_table, recipient_nonzero_table, by = joining_category)
    write.csv(donor_unique, paste(output_name, "donor_unique_table.csv", sep = "_"), row.names = FALSE)
    cat("donor_unique.csv written\n")
    recipient_unique <- anti_join(recipient_nonzero_table, donor_nonzero_table, by = joining_category)
    write.csv(recipient_unique, paste(output_name, "recipient_unique_table.csv", sep = "_"), row.names = FALSE)
    cat("recipient_unique.csv written\n")
    post_fmt_unique <- anti_join(anti_join(post_fmt_nonzero_table, donor_nonzero_table, by = joining_category), recipient_nonzero_table, by = joining_category)
    write.csv(post_fmt_unique, paste(output_name, "post_fmt_unique_table.csv", sep = "_"), row.names = FALSE)
    cat("post_fmt_unique.csv written\n")
    post_fmt_full <- post_fmt_nonzero_table
    write.csv(post_fmt_full, paste(output_name, "post_fmt_full_table.csv", sep = "_"), row.names = FALSE)
    cat("post_fmt_full.csv written\n")

    post_fmt_donor <- semi_join(post_fmt_nonzero_table, donor_unique, by = joining_category)
    post_fmt_rec <- semi_join(post_fmt_nonzero_table, recipient_unique, by = joining_category)
    post_fmt_excl <- bind_rows(post_fmt_donor, post_fmt_rec)
    write.csv(post_fmt_excl, paste(output_name, "pos_fmt_excl_table.csv", sep = "_"), row.names = FALSE)
    cat("post_fmt_excl.csv written\n")
    shared_pre <- semi_join(donor_nonzero_table, recipient_nonzero_table, by = joining_category)
    write.csv(shared_pre, paste(output_name, "shared_pre_table.csv", sep = "_"), row.names = FALSE)
    cat('shared_pre.csv written\n')
    shared_throughout <- semi_join(shared_pre, post_fmt_nonzero_table, by = joining_category)
    write.csv(shared_throughout, paste(output_name, "shared_throughout_table.csv", sep = "_"), row.names = FALSE)
    cat("shared_throughout.csv written\n")
	
    # Get numbers of OTUs unique to each condition
    N_otus_unique_donor <- nrow(donor_unique)
    N_otus_unique_recipient <- nrow(recipient_unique)
    N_otus_unique_post_fmt <- nrow(post_fmt_unique)
    N_otus_shared_pre <- nrow(shared_pre)
    N_otus_shared_throughout <- nrow(shared_throughout)
    N_otus_post_fmt_full <- nrow(post_fmt_full)
    N_otus_post_fmt_excl <- nrow(post_fmt_excl)

    #Create a combined table of OTUs unique for the different conditions, including nonzero fraction
    OTUs_unique_combined <- bind_rows(bind_rows(donor_unique, recipient_unique), post_fmt_unique) %>% arrange(Specificity)
    OTUs_unique_combined_excluded <- bind_rows(donor_unique, recipient_unique) %>% arrange(Specificity)
	

	donor_relative_abundance <- metric_function(donor_only_table, paste("Donor", donor, sep = "_"), comparison)
	recipient_relative_abundance <- metric_function(recipient_only_table, paste("Recipient", recipient, sep = "_"), comparison)
	post_fmt_relative_abundance <- metric_function(post_fmt_only_table, paste("Post_FMT", post_fmt, sep = "_"), comparison)
	
	if (!is.null(taxonomy)){
		relative_abundance_unfiltered <- inner_join(full_join(full_join(donor_relative_abundance[,-3], recipient_relative_abundance[,-3], by = "OTU"), post_fmt_relative_abundance[,-3], by = "OTU"), tax_map, by = "OTU")
		full_relative_abundance <- collapse_taxonomy(inner_join(OTUs_unique_combined, relative_abundance_unfiltered, by = joining_category))
		full_relative_abundance_excl <- collapse_taxonomy(inner_join(OTUs_unique_combined_excluded, relative_abundance_unfiltered, by = joining_category))
	
	} else if (is.null(taxonomy)) {
		relative_abundance_unfiltered <- full_join(full_join(donor_relative_abundance[,-3], recipient_relative_abundance[,-3], by = "OTU"), post_fmt_relative_abundance[,-3], by = "OTU")
		full_relative_abundance <- inner_join(OTUs_unique_combined, relative_abundance_unfiltered, by = joining_category)
		full_relative_abundance_excl <- inner_join(OTUs_unique_combined_excluded, relative_abundance_unfiltered, by = joining_category)
	
	}
	
	write.csv(full_relative_abundance, paste(output_name, 'full_relative_abundance.csv', sep = "_"), row.names = FALSE)
	write.csv(full_relative_abundance_excl, paste(output_name, 'excl_relative_abundance.csv', sep = "_"), row.names = FALSE)
	
	
    setwd('../')

   
    cat("Computing metrics...\n")

    dir.create(file.path(getwd(), 'Individual_Metrics'), showWarnings = FALSE)
    setwd(file.path(getwd(), 'Individual_Metrics'))

    # P_Donor_table, the table of OTUs in the post-transplant samples that came from the donor.
    P_Donor_table <- semi_join(post_fmt_full, donor_unique, by = joining_category)
    write.csv(P_Donor_table, paste(output_name, "P_Donor_table.csv", sep = "_"), row.names = FALSE)
    P_Donor_table_excluded <- semi_join(post_fmt_excl, donor_unique, by = joining_category)
    write.csv(P_Donor_table_excluded, paste(output_name, "P_Donor_table_excl.csv", sep = "_"), row.names = FALSE) ## may be unnecessary

    # N_P_Donor, the number of OTUs in the post-transplant samples that came from the donor
    N_P_Donor <- nrow(P_Donor_table)
    N_P_Donor_excluded <- nrow(P_Donor_table_excluded)

    # D_Engraft the proportion of donor OTUs that made it into the post-transplant samples
    D_Engraft <- N_P_Donor/N_otus_unique_donor
    D_Engraft_excluded <- N_P_Donor_excluded/N_otus_unique_donor

    # P_Donor, the proportion of OTUs in post-transplant samples that came from the donor
    P_Donor <- N_P_Donor/N_otus_post_fmt_full
    P_Donor_excluded <- N_P_Donor_excluded/N_otus_post_fmt_excl
	
    # P_Recipient_table, the table of OTUs in the post-transplant samples that came from the recipient
    P_Recipient_table <- semi_join(post_fmt_full, recipient_unique, by = joining_category)
    write.csv(P_Recipient_table, paste(output_name, "P_Recipient_table.csv", sep = "_"), row.names = FALSE)
    P_Recipient_table_excluded <- semi_join(post_fmt_excl, recipient_unique, by = joining_category)
    write.csv(P_Recipient_table_excluded, paste(output_name, "P_Recipient_table_excl.csv", sep = "_"), row.names = FALSE)

    setwd('../')
	
    # P_Recipient, the number of OTUs in the post-transplant samples that came from the recipient
    N_P_Recipient <- nrow(P_Recipient_table)
    N_P_Recipient_excluded <- nrow(P_Recipient_table_excluded)

    # R_Persist the proportion of recipient OTUs that remained in the post-transplant samples
    R_Persist <- N_P_Recipient/N_otus_unique_recipient
    R_Persist_excluded <- N_P_Recipient_excluded/N_otus_unique_recipient

    # P_Recipient, the proportion of OTUs in post-transplant samples that came from the recipient
    P_Recipient <- N_P_Recipient/N_otus_post_fmt_full
    P_Recipient_excluded <- N_P_Recipient_excluded/N_otus_post_fmt_excl
	
	# P_Shared and P_Unique
	P_Shared <- N_otus_shared_throughout/N_otus_post_fmt_full
	P_Unique <- N_otus_unique_post_fmt/N_otus_post_fmt_full

    cat('Compiling metrics...\n')

    Item <- c("Comment", "TransplantID", "Taxonomy_Added", "Taxonomy_Filename", "OTU_Table_Filename", "N_OTUs_starting",
                       "Mapping_Filename", "Metadata_Category_FMT",
                       "DonorID", "N_Donor_Samples", "RecipientID", "N_Recipient_Samples", "Post_FMT_ID", "N_Post_FMT", "Rel_Abundance_Filter", "N_OTUs_Relative_Abundance_Filtered",
                       "Fleeting_OTUs_Filter", "N_OTUs_Fleeting_Filtered", "Comparison_Metric",
                       "N_don", "N_rec", "N_post_fmt", "N_post_fmt_excl",
                       "N_post_fmt_unique", "N_shared",
                       "N_P_Donor", "D_Engraft", "P_Donor",
                       "N_P_Recipient", "R_Persist", "P_Recipient", "P_Shared", "P_Unique",
                       "D_Engraft_excl", "P_Donor_excl",
                       "R_Persist_excl", "P_Recipient_excl")
    Value <- c(comment, transplant_id, tax_check, tax_file, biom_file, n_otus_starting,
               mapping_file, metadata_category,
               donor, n_donor_samples, recipient, n_recipient_samples, post_fmt, n_post_fmt_samples, min_abundance, n_otus_after_relative_filter,
               min_fraction, n_otus_after_fleeting_filter, comparison,
               N_otus_unique_donor, N_otus_unique_recipient, N_otus_post_fmt_full, N_otus_post_fmt_excl,
               N_otus_unique_post_fmt, N_otus_shared_throughout,
               N_P_Donor, D_Engraft, P_Donor,
               N_P_Recipient, R_Persist, P_Recipient, P_Shared, P_Unique,
               D_Engraft_excluded, P_Donor_excluded,
               R_Persist_excluded, P_Recipient_excluded)
    output_table <- data.frame(Item, Value)
    pre_out <- as.data.frame(t(output_table))
    colnames(pre_out) <- make.names(Item)
    row.names(pre_out) <- NULL
	
	slim_items <- c("Comment", "TransplantID", "D_Engraft", "R_Persist", "P_Donor", "P_Recipient", "P_Shared", "P_Unique")
	slim_values <- c(comment, transplant_id, D_Engraft, R_Persist, P_Donor, P_Recipient, P_Shared, P_Unique)
	
	slim_output <- data.frame(slim_items, slim_values)
	slim_pre_out <- as.data.frame(t(slim_output))
	colnames(slim_pre_out) <- make.names(slim_items)
	row.names(slim_pre_out) <- NULL
	
	slim_table_out <- rbind(slim_table_out, slim_pre_out[-1,])
	

    table_out <- rbind(table_out, pre_out[-1, ])

    cat("Writing output metric table, plot, & QC...\n")
	
    write.csv(output_table, file = paste(output_name,"metric_summary_table.csv", sep = "_"), row.names = FALSE)
    metric_plot <- visualize_metrics_batch(N_otus_unique_donor, N_otus_unique_recipient, N_otus_post_fmt_full, N_P_Donor, N_P_Recipient, N_otus_unique_post_fmt, N_otus_shared_throughout, output_name, paste('Full', joining_category, output_name, sep = "_"), comment)
    metric_plot_excl <- visualize_metrics_batch(N_otus_unique_donor, N_otus_unique_recipient, N_otus_post_fmt_excl, N_P_Donor, N_P_Recipient, 0, 0, output_name, paste('Excl', joining_category, output_name, sep = "_"), comment)

    ## QC Tables & Metrics

    dir.create(file.path(getwd(), 'QC'), showWarnings = FALSE)
    setwd(file.path(getwd(), 'QC'))

    QC_absolute_filtered_OTU_distribution <- colSums(biom_only[,-length(biom_only)])
    QC_absolute_filtered_OTU_distribution_plot <- qplot(QC_absolute_filtered_OTU_distribution, xlab = 'Total OTU Representation (Count)', geom='density', main = "Initial OTU Table Count Distribution")
    pdf(paste("QC_absolute_filtered_OTU_distribution_plot", output_name, '.pdf', sep = '_'), width = 7, height = 7)
    print(QC_absolute_filtered_OTU_distribution_plot)
    dev.off()

    QC_experiment_specific_OTU_distribution <- colSums(unnormalized_experiment_specific)
    QC_experiment_specific_OTU_distribution_plot <- qplot(QC_experiment_specific_OTU_distribution, xlab = 'Total OTU Representation (Count)', geom='density', main = "Experiment-Specific OTU Table Count Distribution")
    pdf(paste("QC_experiment_specific_OTU_distribution_plot", output_name, '.pdf', sep = '_'), width = 7, height = 7)
    print(QC_experiment_specific_OTU_distribution_plot)
    dev.off()

    QC_relative_filtered_OTU_distribution <- colSums(normalized_filtered[, -length(normalized_filtered)])
    QC_relative_filtered_OTU_distribution_plot <- qplot(QC_relative_filtered_OTU_distribution, xlab = 'Total OTU Representation (Sum of Relative Abundances by OTU)', geom='density', main = "Normalized, Relative Abundance Filtered OTU Table Count Distribution", xlim = c(0,1))
    pdf(paste("QC_relative_filtered_OTU_distribution_plot", output_name, '.pdf', sep = '_'), width = 7, height = 7)
    print(QC_relative_filtered_OTU_distribution_plot)
    dev.off()

    QC_absolute_filtered_Sample_Depth_distribution <- rowSums(biom_only[,-length(biom_only)])
    QC_absolute_filtered_Sample_Depth_distribution_plot <- qplot(QC_absolute_filtered_Sample_Depth_distribution, geom="density", xlab = "Sequences Per Sample", ylab = "Count", main = "Initial OTU Table Sample Representation")
    pdf(paste("QC_absolute_filtered_Sample_Depth_distribution_plot", output_name, '.pdf', sep = '_'), width = 7, height = 7)
    print(QC_absolute_filtered_Sample_Depth_distribution_plot)
    dev.off()

    QC_experiment_specific_Sample_Depth_distribution <- rowSums(unnormalized_experiment_specific)
    QC_experiment_specific_Sample_Depth_distribution_plot <- qplot(QC_experiment_specific_Sample_Depth_distribution, geom="density", xlab = "Sequences Per Sample", ylab = "Count", main = "Experiment-Specific OTU Table Sample Representation")
    pdf(paste("QC_experiment_specific_Sample_Depth_distribution_plot", output_name, '.pdf', sep = '_'), width = 7, height = 7)
    print(QC_experiment_specific_Sample_Depth_distribution_plot)
    dev.off()

    QC_relative_filtered_Sample_Depth_distribution <- rowSums(normalized_filtered[,-length(normalized_filtered)])
    QC_relative_filtered_Sample_Depth_distribution_plot <- qplot(QC_relative_filtered_Sample_Depth_distribution, geom="density", xlab = "Proportion of Sequences Remaining Per Sample", ylab = "Count", main = "Normalized, Relative Abundance Filtered OTU Table Sample Representation",xlim =  c(0,1))
    pdf(paste("QC_relative_filtered_Sample_Depth_distribution_plot", output_name, '.pdf', sep = '_'), width = 7, height = 7)
    print(QC_relative_filtered_Sample_Depth_distribution_plot)
    dev.off()

    setwd('../')
	

    setwd('../')
    cat("Finished writing outputs.\n")

    cat("-----End Iteration-----\n")
  }

  write.csv(table_out, file = "All_Transplant_Summary.csv", row.names = FALSE)
  write.csv(slim_table_out, file = "Slim_All_Transplant_Summary.csv", row.names = FALSE)
  
  setwd('../')
  end_time <- proc.time()
  elapsed <- (end_time - start_time)[['elapsed']]
  cat("***Finished***\n")
  cat(paste("Total Time: ", round(elapsed, digits = 2), " seconds\n", sep = ""))

}
