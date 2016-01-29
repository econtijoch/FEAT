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
#' @return many files containing QC and summary statistics, as well as FMT visualizations
#' @export
#'
batchFEAT <- function(biom_file, mapping_file, FMT_pairs, input_params, output_dir) {

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
  cat("Adding Metadata to OTU table...\n")

  metadata_table <- merge(biom_only, mapping, by = 'X.SampleID')

  write.csv(metadata_table, "testing.csv", row.names = F)

  cat("Metadata added to OTU table... \n")

  # Get total number of OTUs to start
  n_otus_starting <- raw_table_otu_count
  cat("# of OTUs to start:", n_otus_starting, '\n')

  # Load Parameters
  parameters <- read.delim(file = input_params)
  min_count <- as.numeric(as.character(parameters[parameters$Parameter == 'min_count', 'Value']))
  metadata_category <- as.character(parameters[parameters$Parameter == 'metadata', 'Value'])
  min_fraction <- as.numeric(as.character(parameters[parameters$Parameter == "min_fraction", "Value"]))
  min_abundance <- as.numeric(as.character(parameters[parameters$Parameter == "min_abundance", "Value"]))
  comparison <- as.character(parameters[parameters$Parameter == "comparison", "Value"])



  transplant_pairs <- read.delim(file = FMT_pairs, header = TRUE)

  table_out_items <- c("Comment", "TransplantID", "Input OTU Table", "# OTUs in raw table",
                       "Mapping File", "Metadata Category to Select FMT Details",
                       "Donor", "# Donor Samples", "Recipient", "# Recipient Samples", "Post-FMT Recipient", "# Post-FMT Recipient Samples", "Minimum Relative Abundance Filter", "# OTUs After Minimum Relative Abundance Filter",
                       "Fleeting OTU Filter Threshold", "# OTUs After Fleeting OTU Filter (Final # OTUs)", "Comparison Metric",
                       "N_don: # specific & unique to donor", "N_rec: # specific & unique to recipient", "N_post_fmt: # specific to post-FMT", "N_post_fmt_excl: # specific to post-FMT, excluding shared and unique",
                       "N_post_fmt_unique: # Unique to post-fmt (should be low/zero)", "N_shared: # shared across donor and recipient throughout",
                       "FMT_don: # from donor in post-FMT", "D_Frac_FMT: proportion of donor in post-FMT", "FMT_FracD: proportion of post-FMT from donor",
                       "FMT_rec: # from recipient in post-FMT", "R_Frac_FMT: proportion of recipient in post-FMT", "FMT_FracR: propotion of post-FMT from recipient",
                       "FMT_don_excl: # from donor in post-FMT (excluding shared & unique)", "D_Frac_FMT_excl: proportion of donor in post-FMT (excluding shared & unique)", "FMT_FracD_excl: proportion of post-FMT from donor (excluding shared & unique)",
                       "FMT_rec_excl: # from recipient in post-FMT (excluding shared & unique)", "R_Frac_FMT_excl : proportion of recipient in post-FMT (excluding shared & unique)", "FMT_FracR_excl: propotion of post-FMT from recipient (excluding shared & unique)")

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
    output_name <- post_fmt
    dir.create(file.path(getwd(), output_name), showWarnings = FALSE)
    setwd(file.path(getwd(), output_name))
    cat("-----Start Iteration (", i, "of", nrow(transplant_pairs),") -----\n")
    transplant_id <- paste(donor, "into", recipient, sep = " ")
    cat("Performing analysis for the following transplant:", transplant_id, "(", comment, ")\n")
    cat("Generating transplant-specific table\n")
    experiment_specific <- split_into_experiment_batch(metadata_table, mapping, metadata_category, donor, recipient, post_fmt)
    write.csv(experiment_specific, paste(output_name, "experiment_specific_table.csv", sep = "_"), row.names = FALSE)


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

    full_table <- inner_join(normalized_filtered, mapping, by = "X.SampleID")
    write.csv(full_table, paste(output_name, "full_table_plus_metadata.csv", sep = "_"), row.names = FALSE)
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

    full_nonzero <- bind_rows(bind_rows(donor_nonzero_table, recipient_nonzero_table, post_fmt_nonzero_table))
    write.csv(full_nonzero, paste(output_name, "full_nonzero_table.csv", sep = "_"), row.names = FALSE)

    distinct_otus_nonzero <- distinct(as.data.frame(full_nonzero$OTU))
    n_otus_after_fleeting_filter <- nrow(distinct_otus_nonzero)
    cat("# of OTUs that remain after this filter:", n_otus_after_fleeting_filter, "\n")

    cat("Writing tables...\n")
    # Make Tables of OTUs unique to each condition
    donor_unique <- anti_join(donor_nonzero_table, recipient_nonzero_table, by = "OTU")
    write.csv(donor_unique, paste(output_name, "donor_unique_table.csv", sep = "_"), row.names = FALSE)
    cat("donor_unique.csv written\n")
    recipient_unique <- anti_join(recipient_nonzero_table, donor_nonzero_table, by = "OTU")
    write.csv(recipient_unique, paste(output_name, "recipient_unique_table.csv", sep = "_"), row.names = FALSE)
    cat("recipient_unique.csv written\n")
    post_fmt_unique <- anti_join(anti_join(post_fmt_nonzero_table, donor_nonzero_table, by = "OTU"), recipient_nonzero_table, by = "OTU")
    write.csv(post_fmt_unique, paste(output_name, "post_fmt_unique_table.csv", sep = "_"), row.names = FALSE)
    cat("post_fmt_unique.csv written\n")
    post_fmt_full <- post_fmt_nonzero_table
    write.csv(post_fmt_full, paste(output_name, "post_fmt_full_table.csv", sep = "_"), row.names = FALSE)
    cat("post_fmt_full.csv written\n")
    post_fmt_excl <- semi_join(post_fmt_nonzero_table, union(donor_unique, recipient_unique), by = "OTU")
    write.csv(post_fmt_excl, paste(output_name, "pos_fmt_excl_table.csv", sep = "_"), row.names = FALSE)
    cat("post_fmt_excl.csv written\n")
    shared_pre <- semi_join(donor_nonzero_table, recipient_nonzero_table, by = "OTU")
    write.csv(shared_pre, paste(output_name, "shared_pre_table.csv", sep = "_"), row.names = FALSE)
    cat('shared_pre.csv written\n')
    shared_throughout <- semi_join(shared_pre, post_fmt_nonzero_table, by = "OTU")
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


    cat("Computing metrics...\n")

    # FMT_don_table, the table of OTUs in the post-transplant samples that came from the donor.
    FMT_don_table <- semi_join(post_fmt_full, donor_unique, by = "OTU")
    write.csv(FMT_don_table, paste(output_name, "FMT_don_table.csv", sep = "_"), row.names = FALSE)
    FMT_don_table_excluded <- semi_join(post_fmt_excl, donor_unique, by = "OTU")
    write.csv(FMT_don_table_excluded, paste(output_name, "FMT_don_table_excl.csv", sep = "_"), row.names = FALSE) ## may be unnecessary

    # FMT_don, the number of OTUs in the post-transplant samples that came from the donor
    FMT_don <- nrow(FMT_don_table)
    FMT_don_excluded <- nrow(FMT_don_table_excluded)

    # D_Frac_FMT the proportion of donor OTUs that made it into the post-transplant samples
    D_Frac_FMT <- FMT_don/nrow(donor_unique)
    D_Frac_FMT_excluded <- FMT_don_excluded/nrow(donor_unique)

    # FMT_FracD, the proportion of OTUs in post-transplant samples that came from the donor
    FMT_FracD <- FMT_don/nrow(post_fmt_full)
    FMT_FracD_excluded <- FMT_don_excluded/nrow(post_fmt_excl)

    # FMT_rec_table, the table of OTUs in the post-transplant samples that came from the recipient
    FMT_rec_table <- semi_join(post_fmt_full, recipient_unique, by = "OTU")
    write.csv(FMT_rec_table, paste(output_name, "FMT_rec_table.csv", sep = "_"), row.names = FALSE)
    FMT_rec_table_excluded <- semi_join(post_fmt_excl, recipient_unique, by = "OTU")
    write.csv(FMT_rec_table_excluded, paste(output_name, "FMT_rec_table_excl.csv", sep = "_"), row.names = FALSE)


    # FMT_rec, the number of OTUs in the post-transplant samples that came from the recipient
    FMT_rec <- nrow(FMT_rec_table)
    FMT_rec_excluded <- nrow(FMT_rec_table_excluded)

    # R_Frac_FMT the proportion of recipient OTUs that remained in the post-transplant samples
    R_Frac_FMT <- FMT_rec/nrow(recipient_unique)
    R_Frac_FMT_excluded <- FMT_rec_excluded/nrow(recipient_unique)

    # FMT_FracR, the proportion of OTUs in post-transplant samples that came from the recipient
    FMT_FracR <- FMT_rec/nrow(post_fmt_full)
    FMT_FracR_excluded <- FMT_rec_excluded/nrow(post_fmt_excl)

    cat('Compiling metrics...\n')

    Item <- c("Comment", "TransplantID", "Input OTU Table", "# OTUs in raw table",
              "Mapping File", "Metadata Category to Select FMT Details",
              "Donor", "# Donor Samples", "Recipient", "# Recipient Samples", "Post-FMT Recipient", "# Post-FMT Recipient Samples", "Minimum Relative Abundance Filter", "# OTUs After Minimum Relative Abundance Filter",
              "Fleeting OTU Filter Threshold", "# OTUs After Fleeting OTU Filter (Final # OTUs)", "Comparison Metric",
              "N_don: # specific & unique to donor", "N_rec: # specific & unique to recipient", "N_post_fmt: # specific to post-FMT", "N_post_fmt_excl: # specific to post-FMT, excluding shared and unique",
              "N_post_fmt_unique: # Unique to post-fmt (should be low/zero)", "N_shared: # shared across donor and recipient throughout",
              "FMT_don: # from donor in post-FMT", "D_Frac_FMT: proportion of donor in post-FMT", "FMT_FracD: proportion of post-FMT from donor",
              "FMT_rec: # from recipient in post-FMT", "R_Frac_FMT: proportion of recipient in post-FMT", "FMT_FracR: propotion of post-FMT from recipient",
              "FMT_don_excl: # from donor in post-FMT (excluding shared & unique)", "D_Frac_FMT_excl: proportion of donor in post-FMT (excluding shared & unique)", "FMT_FracD_excl: proportion of post-FMT from donor (excluding shared & unique)",
              "FMT_rec_excl: # from recipient in post-FMT (excluding shared & unique)", "R_Frac_FMT_excl : proportion of recipient in post-FMT (excluding shared & unique)", "FMT_FracR_excl: propotion of post-FMT from recipient (excluding shared & unique)")
    Value <- c(comment, transplant_id, biom_file, n_otus_starting,
               mapping_file, metadata_category,
               donor, n_donor_samples, recipient, n_recipient_samples, post_fmt, n_post_fmt_samples, min_abundance, n_otus_after_relative_filter,
               min_fraction, n_otus_after_fleeting_filter, comparison,
               N_otus_unique_donor, N_otus_unique_recipient, N_otus_post_fmt_full, N_otus_post_fmt_excl,
               N_otus_unique_post_fmt, N_otus_shared_throughout,
               FMT_don, D_Frac_FMT, FMT_FracD,
               FMT_rec, R_Frac_FMT, FMT_FracR,
               FMT_don_excluded, D_Frac_FMT_excluded, FMT_FracD_excluded,
               FMT_rec_excluded, R_Frac_FMT_excluded, FMT_FracR_excluded)
    output_table <- data.frame(Item, Value)
    pre_out <- as.data.frame(t(output_table))
    colnames(pre_out) <- make.names(Item)
    row.names(pre_out) <- NULL

    table_out <- rbind(table_out, pre_out[-1, ])

    cat("Writing output metric table, plot, & QC...\n")

    write.csv(output_table, file = paste(output_name,"metric_table.csv", sep = "_"), row.names = FALSE)
    metric_plot <- visualize_metrics_batch(N_otus_unique_donor, N_otus_unique_recipient, N_otus_post_fmt_full, FMT_don, FMT_rec, N_otus_unique_post_fmt, N_otus_shared_throughout, output_name, paste('Full', output_name, sep = "_"))
    metric_plot_excl <- visualize_metrics_batch(N_otus_unique_donor, N_otus_unique_recipient, N_otus_post_fmt_excl, FMT_don, FMT_rec, 0, 0, output_name, paste('Excl', output_name, sep = "_"))

    ## QC Tables & Metrics


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
    cat("Table written to working directory.\n")

    cat("-----End Iteration-----\n")
  }

  write.csv(table_out, file = "All_Transplant_Summary.csv")

  setwd('../')


}
