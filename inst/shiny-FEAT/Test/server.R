require(phyloseq)
require(rhdf5)
require(biom)
require(shiny)
require(DT)
require(ggplot2)
require(vegan)
require(lazyeval)
require(FEAT)
require(shinyjs)
require(stringr)
require(reshape2)
require(dplyr)

# source('helper_funcitons.R')
ggsave <- ggplot2::ggsave
body(ggsave) <- body(ggplot2::ggsave)[-2]
options(digits = 4)

# Define server logic
shinyServer(function(input, output, session) {
  # Read OTU Table
  biom_table <- reactive({
    validate(need(
      !is.null(input$input_otu_table) &&
        !is.null(input$metadata_file),
      "Please load an OTU table and metadata file."
    ))
    raw_table <- input$input_otu_table
    withProgress(message = "1) Loading OTU table", value = 0.75, {
      biom_table <- read_biom(raw_table$datapath)
      data <- as(biom_data(biom_table), "matrix")
      biom_only <- as.data.frame(t(data))
      raw_table_otu_count <<- ncol(biom_only)
      colnames(biom_only) <-
        paste("OTU", colnames(biom_only), sep = "_")
      biom_only$X.SampleID <-
        as.character(row.names(biom_only))
      row.names(biom_only) <- NULL
    })
    return(biom_only)
  })
  
  # Read metadata
  metadata <- reactive({
    validate(need(
      !is.null(input$input_otu_table) &&
        !is.null(input$metadata_file),
      "Please load an OTU table and metadata file."
    ))
    map_file <- input$metadata_file
    if (is.null(map_file)) {
      return()
    }
    map <- read.delim(file = map_file$datapath)
    map$X.SampleID <- as.character(map$X.SampleID)
    if ("biomass_ratio" %in% colnames(map)) {
      map <- map[!is.na(map$biomass_ratio), ]
    } else {
      stop("Biomass data not in metadata file")
    }
    
    return(map)
  })
  
  # Create an output to toggle FMT detail selection UI
  output$files_loaded <- reactive({
    return(!(is.null(metadata()) | is.null(taxonomy())))
  })
  outputOptions(output, 'files_loaded', suspendWhenHidden = FALSE)
  
  output$depth <- renderUI({
    if (input$toggle_taxonomy) {
      selectInput(
        'selected_depth',
        label = "Select Taxonomic Depth",
        choices = c(
          "Kingdom",
          "Phylum",
          "Class",
          "Order",
          "Family",
          "Genus",
          "Species"
        ),
        selected = 'Species'
      )
    } else{
      return()
    }
  })
  
  # Get number of OTUs to start
  num_otus_raw <- reactive({
    if (is.null(biom_table())) {
      return()
    }
    return(raw_table_otu_count)
  })
  # Return number of OTUs to start
  output$raw_otu_count <- renderText({
    return(paste("# of OTUs in raw table:", num_otus_raw()))
  })
  
  # Read in taxonomy table, allow for using either greengenes or user-created
  taxonomy <- reactive({
    if (is.null(input$id_tax_file)) {
      return()
    }
    withProgress(message = "Loading Taxonomy", value = 0, {
      raw <- input$id_tax_file
      table <- read.delim(raw$datapath, header = FALSE)
      incProgress(0.5, detail = "Cleaning Up")
      if (ncol(table) == 2) {
        colnames(table) <- c("OTU_ID", "Taxon")
      } else if (ncol(table) == 3) {
        colnames(table) <- c("OTU_ID", "Taxon", "QualityScore")
      }
      table$OTU_ID <- paste("OTU", table$OTU, sep = "_")
      output <-
        as.data.frame(apply(table, 2, function(x)
          gsub("\\s+", "", x)))
      output$Taxon <- as.character(output$Taxon)
      output$OTU_ID <- as.character(output$OTU_ID)
    })
    return(output)
  })
  
  # Create a data table that combines the metadata and OTU data
  biom_merged <- reactive({
    return(dplyr::inner_join(metadata(), biom_table(), by = 'X.SampleID'))
  })
  
  ##########################################################################################
  # Define relevant FMT conditions                                                         #
  ##########################################################################################
  
  # Define metadata conditions to compare to populate drop-down menus
  conditions_of_compare <- reactive({
    if (!is.null(metadata()) && !is.null(input$comparison)) {
      category <- metadata()[[input$comparison]]
      return(levels(category))
    } else{
      return()
    }
  })
  
  ## Help select FMT details
  # Select metadata category to pick FMT details from
  output$compare_options <- renderUI({
    selectInput("comparison",
                "Select Metadata Category",
                names(metadata()),
                selected = "Compare")
  })
  # Select Donor
  output$donor <- renderUI({
    selectInput("donor",
                "Select donor condition",
                conditions_of_compare(),
                selected = conditions_of_compare()[1])
  })
  # Select Recipient
  output$recipient <- renderUI({
    selectInput(
      "recipient",
      "Select pre-transplant recipient condition",
      conditions_of_compare(),
      selected = conditions_of_compare()[2]
    )
  })
  # Select Post-FMT
  output$post_fmt <- renderUI({
    selectInput(
      "post_fmt",
      "Select post-transplant recipient condition",
      conditions_of_compare(),
      selected = conditions_of_compare()[3]
    )
  })
  
  # Create labels to keep things straight
  donor <-
    eventReactive(input$go, {
      return(input$donor)
    })
  recipient <-
    eventReactive(input$go, {
      return(input$recipient)
    })
  post_fmt <-
    eventReactive(input$go, {
      return(input$post_fmt)
    })
  output$donor_id <- renderText({
    donor()
  })
  output$recipient_id <- renderText({
    recipient()
  })
  output$post_fmt_id <- renderText({
    post_fmt()
  })
  output$test_metric <- renderText({
    input$comparison_test
  })
  
  # Create Experiment-specific table - only smaples from selected conditions
  biom_experiment_specific <-
    eventReactive(input$go, {
      validate(need(
        !is.null(input$input_otu_table) &&
          !is.null(input$metadata_file),
        "Please load an OTU table and metadata file."
      ))
      withProgress(message = "2) Creating Experiment-Specific Table", value = 0.5, {
        table_out <-
          biom_merged()[biom_merged()[, input$comparison] %in% c(input$donor, input$recipient, input$post_fmt),]
      })
      return(table_out)
    })
  
  # Pull number of donor samples, create output
  N_donor_samples <-
    eventReactive(input$go, {
      if (is.null(biom_experiment_specific())) {
        return()
      }
      table <- biom_experiment_specific()
      donor <- table[table[, input$comparison] == input$donor,]
      output <- nrow(donor)
      return(output)
    })
  text_donor_samples <- reactive({
    output <- paste("# Donor Samples Used:", N_donor_samples())
    return(output)
  })
  output$num_donor_samples <- renderText({
    text_donor_samples()
  })
  
  # Pull number of recipient samples, create output
  N_recipient_samples <-
    eventReactive(input$go, {
      if (is.null(biom_experiment_specific())) {
        return()
      }
      table <- biom_experiment_specific()
      recipient <-
        table[table[, input$comparison] == input$recipient,]
      output <- nrow(recipient)
      return(output)
    })
  text_recipient_samples <- reactive({
    output <- paste("# Recipient Samples Used:", N_recipient_samples())
    return(output)
  })
  output$num_recipient_samples <-
    renderText({
      text_recipient_samples()
    })
  
  # Pull number of post-FMT samples, create output
  N_post_fmt_samples <-
    eventReactive(input$go, {
      if (is.null(biom_experiment_specific())) {
        return()
      }
      table <- biom_experiment_specific()
      post_fmt <-
        table[table[, input$comparison] == input$post_fmt,]
      output <- nrow(post_fmt)
      return(output)
    })
  text_post_fmt_samples <- reactive({
    output <- paste("# Post-FMT Samples Used:", N_post_fmt_samples())
    return(output)
  })
  output$num_post_fmt_samples <-
    renderText({
      text_post_fmt_samples()
    })
  
  # Create an output to toggle quick FMT summary
  output$data_analyzed <- reactive({
    return(!is.null(biom_experiment_specific()))
  })
  outputOptions(output, 'data_analyzed', suspendWhenHidden = FALSE)
  
  # Create table of OTU data only (specific for selected samples)
  otus_only <- eventReactive(input$go, {
    merged <- biom_experiment_specific()
    row.names(merged) <- merged$X.SampleID
    table <- merged[, grepl("OTU_", names(merged))]
    return(table)
  })
  
  # Create table of metadata only (specific for selected samples)
  metadata_only <- eventReactive(input$go, {
    table <-
      biom_experiment_specific()[, !grepl("OTU_", names(biom_experiment_specific()))]
    table$X.SampleID <- as.character(table$X.SampleID)
    return(table)
  })
  
  ## Normalize and Filter experiment-specific table
  relative <- eventReactive(input$go, {
    withProgress(message = '3) Normalizing and Filtering', value = 0.2, {
      normalized <- otus_only() / rowSums(otus_only())
      incProgress(0.4)
      max_otu_fraction <- apply(normalized, 2, max)
      filtered <-
        normalized[, max_otu_fraction > input$min_OTU_fraction]
    })
    filtered_table_otu_count <<- ncol(filtered)
    return(filtered)
  })
  
  
  # Create Output for filtered OTU table OTU count
  num_otus_relative_filtered <- eventReactive(input$go, {
    if (is.null(biom_table())) {
      return()
    }
    return(filtered_table_otu_count)
  })
  output$num_otus_after_relative_filter <- renderText({
    return(paste(
      "# of OTUs after rel. abundance filter:",
      num_otus_relative_filtered()
    ))
  })
  
  ## Scale Normalized table
  scaled <- eventReactive(input$go, {
    withProgress(message = "4) Scaling Abundances", value = 0.8, {
      output <- relative() * metadata_only()$biomass_ratio
    })
    return(output)
  })
  
  ### Taxonomy Addition/Compression
  # Melt data to add taxonomy
  abundance_selected_melted <- eventReactive(input$go, {
    if (input$abundance_type == 'absolute') {
      scaled_only <- scaled()
      scaled_only$X.SampleID <- row.names(scaled_only)
      melted <-
        melt(
          data = scaled_only,
          id.vars = 'X.SampleID',
          measure.vars = colnames(scaled_only)[-length(colnames(scaled_only))]
        )
      # QC to elimnate errors and spurious negative abundance (spot-checked, this occurs only once from all samples, and is minimally negative anyways...)
      if (length(melted[melted$value < 0, ]$value) > 0) {
        melted[melted$value < 0, ]$value <- 0
      }
    } else {
      relative_only <- relative()
      relative_only$X.SampleID <- row.names(relative_only)
      melted <-
        melt(
          data = relative_only,
          id.vars = 'X.SampleID',
          measure.vars = colnames(relative_only)[-length(colnames(relative_only))]
        )
    }
    names(melted)[names(melted) == 'variable'] <- 'OTU_ID'
    melted$OTU_ID <- as.character(melted$OTU_ID)
    return(melted)
  })
  
  
  # Add Taxonomy and collapse upon common taxa
  
  phylogeny <-
    c("Kingdom",
      "Phylum",
      "Class",
      "Order",
      "Family",
      "Genus",
      "Species")
  
  abundance_selected_taxonomy <- eventReactive(input$go, {
    withProgress(message = 'Collapsing by Taxonomy', value = 0.5, {
      tax_added <-
        dplyr::inner_join(abundance_selected_melted(), taxonomy(), by = "OTU_ID")
      incProgress(0.3)
      grouped <- dplyr::group_by(tax_added, X.SampleID, Taxon)
      collapsed <- dplyr::summarize(grouped, abundance = sum(value))
      split_tax <-
        tidyr::separate(
          collapsed,
          Taxon,
          phylogeny,
          sep = ";.__",
          remove = F,
          fill = 'right'
        )
      split_tax$Kingdom <- "Bacteria"
      split_tax$Abundance_Type <- input$abundance_type
    })
    return(split_tax)
  })
  
  # Create a semi-static variable to copture current abundance type
  abundance_type <- eventReactive(input$go, {
    return(input$abundance_type)
  })
  
  # Create a table of per-sample-per-taxa abundances at the selected depth
  sample_abundance_by_depth <- reactive({
    i <- grep(input$selected_depth, phylogeny)
    label <-
      paste(phylogeny[i - 1], input$selected_depth, sep = ".")
    dots <- lapply(phylogeny[1:i], as.symbol)
    grouped1 <-
      dplyr::group_by(abundance_selected_taxonomy(), X.SampleID)
    grouped2 <- dplyr::group_by_(grouped1, .dots = dots, add = T)
    table <-
      dplyr::summarize(grouped2, abundance = sum(abundance))
    if (i > 1) {
      table$short_label <-
        paste(table[[phylogeny[i - 1]]], table[[phylogeny[i]]], sep = ".")
      table <-
        tidyr::unite_(table,
                      'long_label',
                      phylogeny[1:i],
                      sep = ".",
                      remove = F)
    } else {
      table$short_label <- table[[phylogeny[i]]]
      table$long_label <- table[[phylogeny[i]]]
    }
    table$Abundance_Type <- abundance_type()
    table$Depth <- phylogeny[i]
    return(table)
  })
  
  # Create a compact version of this table that contains as many rows as variables, and columns with metadata + taxa abundances
  compact_sample_abundance_by_depth <- reactive({
    compact_table <-
      dcast(
        sample_abundance_by_depth(),
        X.SampleID + Depth + Abundance_Type ~ long_label,
        value.var = 'abundance'
      )
    output <-
      dplyr::left_join(metadata_only(), compact_table, by = "X.SampleID")
    return(output)
  })
  
  # Create label switch for OTUs vs Taxa
  OTU_taxa_label <- reactive({
    if (input$toggle_taxonomy) {
      output <- "Taxa"
    } else {
      output <- "OTUs"
    }
    return(output)
  })
  
  abundance_selected_metadata_added <- eventReactive(input$go, {
    if (input$abundance_type == 'absolute') {
      scaled_only <- scaled()
      scaled_only$X.SampleID <- as.character(row.names(scaled_only))
      output <-
        dplyr::left_join(metadata_only(), scaled_only, by = 'X.SampleID')
    } else {
      relative_only <- relative()
      relative_only$X.SampleID <-
        as.character(row.names(relative_only))
      output <-
        dplyr::left_join(metadata_only(), relative_only, by = "X.SampleID")
    }
    if (input$toggle_taxonomy) {
      return(compact_sample_abundance_by_depth())
    } else {
      return(output)
    }
  })
  
  detail_text <- eventReactive(input$go, {
    if (input$toggle_taxonomy) {
      output <- paste(input$selected_depth,
                      "-Level Analysis (Taxonomy Added)",
                      sep = "")
    } else {
      output <- "OTU-based Analysis (Taxonomy NOT Added)"
    }
    return(output)
  })
  output$other_details <- renderText({
    detail_text()
  })
  
  ## FEAT metrics
  # OTUs specific to Donor, in minimum fraction of samples
  abundance_selected_donor_specific <- eventReactive(input$go, {
    donor_only <-
      abundance_selected_metadata_added()[abundance_selected_metadata_added()[, input$comparison] == donor(),]
    if (input$toggle_taxonomy) {
      table <-
        as.data.frame(t(summarise_each(donor_only[, grepl("Bacteria", names(donor_only))], funs(
          sum(. > 0) / length(.)
        ))))
      table$Taxon <- row.names(table)
      colnames(table) <- c('Fraction', 'Taxon')
      output <-
        dplyr::filter(table_reorder_first(table, 'Taxon'),
                      Fraction >= input$min_fraction)
      output$Specificity <- donor()
    } else {
      table <-
        as.data.frame(t(summarise_each(donor_only[, grepl("OTU_", names(donor_only))], funs(
          sum(. > 0) / length(.)
        ))))
      table$OTU <- row.names(table)
      colnames(table) <- c('Fraction', 'OTU_ID')
      output <-
        dplyr::filter(table_reorder_first(table, 'OTU_ID'),
                      Fraction >= input$min_fraction)
      output$Specificity <- donor()
    }
    return(output)
  })
  
  # OTUs specific to Recipient, in minimum fraction of samples
  abundance_selected_recipient_specific <- eventReactive(input$go, {
    recipient_only <-
      abundance_selected_metadata_added()[abundance_selected_metadata_added()[, input$comparison] == recipient(),]
    if (input$toggle_taxonomy) {
      table <-
        as.data.frame(t(summarise_each(
          recipient_only[, grepl("Bacteria", names(recipient_only))], funs(sum(. > 0) / length(.))
        )))
      table$Taxon <- row.names(table)
      colnames(table) <- c('Fraction', 'Taxon')
      output <-
        dplyr::filter(table_reorder_first(table, 'Taxon'),
                      Fraction >= input$min_fraction)
      output$Specificity <- recipient()
    } else {
      table <-
        as.data.frame(t(summarise_each(
          recipient_only[, grepl("OTU_", names(recipient_only))], funs(sum(. > 0) / length(.))
        )))
      table$OTU <- row.names(table)
      colnames(table) <- c('Fraction', 'OTU_ID')
      output <-
        dplyr::filter(table_reorder_first(table, 'OTU_ID'),
                      Fraction >= input$min_fraction)
      output$Specificity <- recipient()
    }
    return(output)
  })
  
  # OTUs specific to Post-FMT, in minimum fraction of samples
  abundance_selected_post_fmt_specific <- eventReactive(input$go, {
    post_fmt_only <-
      abundance_selected_metadata_added()[abundance_selected_metadata_added()[, input$comparison] == post_fmt(),]
    if (input$toggle_taxonomy) {
      table <-
        as.data.frame(t(summarise_each(
          post_fmt_only[, grepl("Bacteria", names(post_fmt_only))], funs(sum(. > 0) / length(.))
        )))
      table$Taxon <- row.names(table)
      colnames(table) <- c('Fraction', 'Taxon')
      output <-
        dplyr::filter(table_reorder_first(table, 'Taxon'),
                      Fraction >= input$min_fraction)
      output$Specificity <- post_fmt()
    } else {
      table <-
        as.data.frame(t(summarise_each(
          post_fmt_only[, grepl("OTU_", names(post_fmt_only))], funs(sum(. > 0) / length(.))
        )))
      table$OTU <- row.names(table)
      colnames(table) <- c('Fraction', 'OTU_ID')
      output <-
        dplyr::filter(table_reorder_first(table, 'OTU_ID'),
                      Fraction >= input$min_fraction)
      output$Specificity <- post_fmt()
    }
    return(output)
  })
  
  # Full table of OTUs that are specific and non-fleeting
  full_nonfleeting <- eventReactive(input$go, {
    if (is.null(abundance_selected_donor_specific()) |
        is.null(abundance_selected_recipient_specific()) |
        is.null(abundance_selected_post_fmt_specific())) {
      return()
    }
    withProgress(message = paste('Filtering Fleeting', OTU_taxa_label()),
                 value = 0.5,
                 {
                   out <-
                     bind_rows(
                       bind_rows(
                         abundance_selected_donor_specific(),
                         abundance_selected_recipient_specific()
                       ),
                       abundance_selected_post_fmt_specific()
                     )
                 })
    return(out)
  })
  
  # List of distinct OTUs that are specific and non-fleeting (i.e. relevant OTUs)
  relevant_OTUs <- eventReactive(input$go, {
    if (is.null(abundance_selected_donor_specific()) |
        is.null(abundance_selected_recipient_specific()) |
        is.null(abundance_selected_post_fmt_specific())) {
      return()
    }
    if (input$toggle_taxonomy) {
      otus <- full_nonfleeting()$Taxon
    } else {
      otus <- full_nonfleeting()$OTU
    }
    return(unique(otus))
  })
  # Create an output to track how many OTUs survive this filter
  N_otus_after_fleeting_filter <- eventReactive(input$go, {
    return(length(relevant_OTUs()))
  })
  n_otus_fleeting_filter_text <- eventReactive(input$go, {
    output <- paste(
      '# of',
      OTU_taxa_label(),
      'after fleeting filter:',
      N_otus_after_fleeting_filter()
    )
    return(output)
  })
  output$num_otus_after_fleeting_filter <- renderText({
    return(n_otus_fleeting_filter_text())
  })
  
  ## Create Table of OTUs unique  to Donor
  donor_unique <- eventReactive(input$go, {
    if (input$toggle_taxonomy) {
      out <-
        dplyr::anti_join(
          abundance_selected_donor_specific(),
          abundance_selected_recipient_specific(),
          by = "Taxon"
        )
    } else {
      out <-
        dplyr::anti_join(
          abundance_selected_donor_specific(),
          abundance_selected_recipient_specific(),
          by = "OTU_ID"
        )
    }
    return(out)
  })
  ## Create Table of OTUs unique  to Recipient
  recipient_unique <- eventReactive(input$go, {
    if (input$toggle_taxonomy) {
      out <-
        dplyr::anti_join(
          abundance_selected_recipient_specific(),
          abundance_selected_donor_specific(),
          by = "Taxon"
        )
    } else {
      out <-
        dplyr::anti_join(
          abundance_selected_recipient_specific(),
          abundance_selected_donor_specific(),
          by = "OTU_ID"
        )
    }
    return(out)
  })
  ## Create Table of OTUs unique to Post-FMT (should be small)
  post_fmt_unique <- eventReactive(input$go, {
    if (input$toggle_taxonomy) {
      out <-
        dplyr::anti_join(
          dplyr::anti_join(
            abundance_selected_post_fmt_specific(),
            abundance_selected_donor_specific(),
            by = "Taxon"
          ),
          abundance_selected_recipient_specific(),
          by = "Taxon"
        )
    } else {
      out <-
        dplyr::anti_join(
          dplyr::anti_join(
            abundance_selected_post_fmt_specific(),
            abundance_selected_donor_specific(),
            by = "OTU_ID"
          ),
          abundance_selected_recipient_specific(),
          by = "OTU_ID"
        )
    }
    return(out)
  })
  ## Create Table of OTUs shared before transplant
  shared_pre <- eventReactive(input$go, {
    if (input$toggle_taxonomy) {
      out <-
        dplyr::semi_join(
          abundance_selected_donor_specific(),
          abundance_selected_recipient_specific(),
          by = "Taxon"
        )
    } else {
      out <-
        dplyr::semi_join(
          abundance_selected_donor_specific(),
          abundance_selected_recipient_specific(),
          by = "OTU_ID"
        )
    }
    return(out)
  })
  ## Create Table of OTUs Shared throughout
  shared_throughout <- eventReactive(input$go, {
    if (input$toggle_taxonomy) {
      output <-
        dplyr::semi_join(shared_pre(),
                         abundance_selected_post_fmt_specific(),
                         by = "Taxon")
      output$Specificity <- "Shared"
    } else {
      output <-
        dplyr::semi_join(shared_pre(),
                         abundance_selected_post_fmt_specific(),
                         by = "OTU_ID")
      output$Specificity <- "Shared"
    }
    return(output)
  })
  
  # Create Table of OTUs shared but lost
  shared_lost <- eventReactive(input$go, {
    if (input$toggle_taxonomy) {
      output <-
        dplyr::anti_join(shared_pre(),
                         abundance_selected_post_fmt_specific(),
                         by = "Taxon")
    } else {
      output <-
        dplyr::anti_join(shared_pre(),
                         abundance_selected_post_fmt_specific(),
                         by = "OTU_ID")
    }
    if (nrow(output) > 0) {
      output$Specificity <- "Shared_lost"
    }
    else {
      output <- NULL
    }
    return(output)
  })
  
  ## Return numbers of OTUs unique to each condition
  # Number donor unique
  N_Donor <- eventReactive(input$go, {
    if (is.null(donor_unique())) {
      return()
    }
    return(nrow(donor_unique()))
  })
  
  output$N_Donor <- renderUI({
    HTML(paste(
      "<strong>N",
      tags$sub("Donor"),
      ":</strong>",
      " ",
      N_Donor(),
      sep = ""
    ))
  })
  
  # Number recipient unique
  N_Recipient <- eventReactive(input$go, {
    if (is.null(recipient_unique())) {
      return()
    }
    return(nrow(recipient_unique()))
  })
  output$N_Recipient <- renderUI({
    HTML(paste(
      "<strong>N",
      tags$sub("Recipient"),
      ":</strong>",
      " ",
      N_Recipient(),
      sep = ""
    ))
  })
  
  # Number post-fmt unique
  N_P_Unique <- eventReactive(input$go, {
    if (is.null(post_fmt_unique())) {
      return()
    }
    return(nrow(post_fmt_unique()))
  })
  output$N_P_Unique <- renderUI({
    HTML(paste(
      "<strong>N",
      tags$sub("Post-FMT (Unique)"),
      ":</strong>",
      " ",
      N_P_Unique(),
      sep = ""
    ))
  })
  
  # Number shared pre-transplant
  N_otus_shared_pre <- eventReactive(input$go, {
    if (is.null(shared_pre())) {
      return()
    }
    return(nrow(shared_pre()))
  })
  output$num_otus_shared_pre <- renderText({
    N_otus_shared_pre()
  })
  
  # Number shared throughout
  N_P_Shared <- eventReactive(input$go, {
    if (is.null(shared_throughout())) {
      return()
    }
    return(nrow(shared_throughout()))
  })
  output$N_P_Shared <- renderUI({
    HTML(paste(
      "<strong>N",
      tags$sub("Post-FMT (Shared)"),
      ":</strong>",
      " ",
      N_P_Shared(),
      sep = ""
    ))
  })
  
  # Number shared but lost after transplant
  N_otus_shared_lost <- eventReactive(input$go, {
    if (is.null(shared_lost())) {
      return(0)
    }
    return(nrow(shared_lost()))
  })
  output$N_shared_lost <- renderUI(({
    HTML(paste(
      "<strong>N",
      tags$sub("Lost | Shared"),
      ":</strong>",
      " ",
      N_otus_shared_lost(),
      sep = ""
    ))
  }))
  
  # Total post-FMT
  N_P_Total <- eventReactive(input$go, {
    if (is.null(abundance_selected_post_fmt_specific())) {
      return()
    }
    return(nrow(abundance_selected_post_fmt_specific()))
  })
  output$N_P_Total <- renderUI({
    HTML(paste(
      "<strong>N",
      tags$sub("Post-FMT (Total)"),
      ":</strong>",
      " ",
      N_P_Total(),
      sep = ""
    ))
  })
  
  #### FMT Metrics
  
  
  P_Donor_table <- eventReactive(input$go, {
    if (input$toggle_taxonomy) {
      output <-
        semi_join(abundance_selected_post_fmt_specific(),
                  donor_unique(),
                  by = "Taxon")
    } else {
      output <-
        semi_join(abundance_selected_post_fmt_specific(),
                  donor_unique(),
                  by = "OTU_ID")
    }
    return(output)
  })
  output$P_Donor_table <- renderDataTable({
    P_Donor_table()
  })
  
  # N_P_Donor, the number of OTUs in the post-transplant samples that came from the donor
  N_P_Donor <- eventReactive(input$go, {
    if (is.null(P_Donor_table())) {
      return()
    }
    return(nrow(P_Donor_table()))
  })
  output$N_P_Donor <- renderUI({
    HTML(paste(
      "<strong>N",
      tags$sub("Donor| Post-FMT"),
      ":</strong>",
      " ",
      N_P_Donor(),
      sep = ""
    ))
  })
  
  # D_Engraft the proportion of donor OTUs that made it into the post-transplant samples
  D_Engraft <- eventReactive(input$go, {
    if (is.null(donor_unique())) {
      return()
    }
    fraction <- N_P_Donor() / N_Donor()
    return(fraction)
  })
  output$D_Engraft <- renderUI({
    HTML(paste(
      "<strong>D",
      tags$sub("Engraft"),
      ":</strong>",
      " ",
      round(D_Engraft(), 3),
      sep = ""
    ))
  })
  
  # P_Donor, the proportion of OTUs in post-transplant samples that came from the donor
  P_Donor <- eventReactive(input$go, {
    if (is.null(abundance_selected_post_fmt_specific())) {
      return()
    }
    fraction <- N_P_Donor() / N_P_Total()
    return(fraction)
  })
  output$P_Donor <- renderUI({
    HTML(paste(
      "<strong>P",
      tags$sub("Donor"),
      ":</strong>",
      " ",
      round(P_Donor(), 3),
      sep = ""
    ))
  })
  
  # P_Recipient_table, the table of OTUs in the post-transplant samples that came from the recipient
  P_Recipient_table <- eventReactive(input$go, {
    if (input$toggle_taxonomy) {
      output <-
        semi_join(abundance_selected_post_fmt_specific(),
                  recipient_unique(),
                  by = "Taxon")
    } else {
      output <-
        semi_join(abundance_selected_post_fmt_specific(),
                  recipient_unique(),
                  by = "OTU_ID")
    }
    return(output)
  })
  output$P_Recipient_table <- renderDataTable({
    P_Recipient_table()
  })
  
  # P_Recipient, the number of OTUs in the post-transplant samples that came from the recipient
  N_P_Recipient <- eventReactive(input$go, {
    if (is.null(P_Recipient_table())) {
      return()
    }
    return(nrow(P_Recipient_table()))
  })
  output$N_P_Recipient <- renderUI({
    HTML(paste(
      "<strong>N",
      tags$sub("Recipient | Post-FMT"),
      ":</strong>",
      " ",
      N_P_Recipient(),
      sep = ""
    ))
  })
  
  # R_Persist the proportion of recipient OTUs that remained in the post-transplant samples
  R_Persist <- eventReactive(input$go, {
    if (is.null(recipient_unique())) {
      return()
    }
    fraction <- N_P_Recipient() / N_Recipient()
    return(fraction)
  })
  output$R_Persist <- renderUI({
    HTML(paste(
      "<strong>R",
      tags$sub("Persist"),
      ":</strong>",
      " ",
      round(R_Persist(), 3),
      sep = ""
    ))
  })
  
  # P_Recipient, the proportion of OTUs in post-transplant samples that came from the recipient
  P_Recipient <- eventReactive(input$go, {
    if (is.null(abundance_selected_post_fmt_specific())) {
      return()
    }
    fraction <- N_P_Recipient() / N_P_Total()
    return(fraction)
  })
  output$P_Recipient <- renderUI({
    HTML(paste(
      "<strong>P",
      tags$sub("Recipient"),
      ":</strong>",
      " ",
      round(P_Recipient(), 3),
      sep = ""
    ))
  })
  
  # P_Shared , proporiton of OTUs in post-transplant samples that are shared
  P_Shared <- eventReactive(input$go, {
    if (is.null(abundance_selected_post_fmt_specific())) {
      return()
    }
    fraction <- N_P_Shared() / N_P_Total()
    return(fraction)
  })
  
  output$P_Shared <- renderUI({
    HTML(paste(
      "<strong>P",
      tags$sub("Shared"),
      ":</strong>",
      " ",
      round(P_Shared(), 3),
      sep = ""
    ))
  })
  
  # P_Unique, proporiton of OTUs in post-transplant samples that are unique
  P_Unique <- eventReactive(input$go, {
    if (is.null(abundance_selected_post_fmt_specific())) {
      return()
    }
    fraction <- N_P_Unique() / N_P_Total()
    return(fraction)
  })
  output$P_Unique <- renderUI({
    HTML(paste(
      "<strong>P",
      tags$sub("Unique"),
      ":</strong>",
      " ",
      round(P_Unique(), 3),
      sep = ""
    ))
  })
  
  # Metric Vizualization
  metric_vis <- eventReactive(input$go, {
    p <-
      visualize_metrics(
        N_Donor(),
        N_Recipient(),
        N_P_Total(),
        N_P_Donor(),
        N_P_Recipient(),
        N_P_Unique(),
        N_P_Shared(),
        post_fmt()
      )
    return(p)
    
  })
  output$metric_visualization <- renderPlot({
    metric_vis()
  })
  
  
  final_otu_matrix <- eventReactive(input$go, {
    if (input$toggle_taxonomy) {
      table <-
        abundance_selected_metadata_added()[, c(T, grepl(
          'Bacteria',
          colnames(abundance_selected_metadata_added())
        )[-1])]
    } else {
      table <-
        abundance_selected_metadata_added()[, c(T, grepl('OTU_', colnames(
          abundance_selected_metadata_added()
        ))[-1])]
    }
    matrix_otus <-  as.matrix(table[, -1])
    rownames(matrix_otus) <- table[, 1]
    return(matrix_otus)
  })
  
  pc_table <- reactive({
    distance_table <-
      as.matrix(vegan::vegdist(final_otu_matrix(), method = input$distance_method))
    pca <-
      cmdscale(
        distance_table,
        k = nrow(final_otu_matrix()) - 1,
        eig = TRUE,
        add = TRUE
      )
    return(pca)
  })
  
  data_viz_table <- reactive({
    shannon <- vegan::diversity(final_otu_matrix(), index = 'shannon')
    simpson <-
      vegan::diversity(final_otu_matrix(), index = 'simpson')
    pc_only <-
      bind_cols(
        as.data.frame(pc_table()$points[, 1]),
        as.data.frame(pc_table()$points[, 2]),
        as.data.frame(pc_table()$points[, 3]),
        as.data.frame(pc_table()$points[, 4]),
        as.data.frame(pc_table()$points[, 5])
      )
    colnames(pc_only) <- c("PC1", "PC2", "PC3", "PC4", "PC5")
    pc_only$X.SampleID <-
      as.character(row.names(final_otu_matrix()))
    pc_only$Shannon <- shannon
    pc_only$Simpson <- simpson
    if (input$toggle_taxonomy) {
      pre_output <-
        dplyr::left_join(pc_only, abundance_selected_metadata_added(), by = 'X.SampleID')
    } else {
      pre_output <-
        dplyr::left_join(pc_only, abundance_selected_metadata_added(), by = 'X.SampleID')
    }
    output_table <- droplevels.data.frame(pre_output)
    return(output_table)
  })
  
  
  output$pc1 <- renderText({
    eig <- pc_table()$eig
    return(percent(eig[1] / sum(eig)))
  })
  output$pc2 <- renderText({
    eig <- pc_table()$eig
    return(percent(eig[2] / sum(eig)))
  })
  output$pc3 <- renderText({
    eig <- pc_table()$eig
    return(percent(eig[3] / sum(eig)))
  })
  output$pc4 <- renderText({
    eig <- pc_table()$eig
    return(percent(eig[4] / sum(eig)))
  })
  output$pc5 <- renderText({
    eig <- pc_table()$eig
    return(percent(eig[5] / sum(eig)))
  })
  
  output$plot_x <- renderUI({
    selectInput('plot_x', 'X', names(data_viz_table()))
  })
  output$plot_y <- renderUI({
    selectInput('plot_y',
                'Y',
                names(data_viz_table()),
                names(data_viz_table())[[2]])
  })
  output$plot_color_by <- renderUI({
    selectInput('plot_color_by',
                'Color By',
                c('None', names(data_viz_table())),
                selected = input$comparison)
  })
  output$plot_facet_row <- renderUI({
    selectInput('facet_row', 'Facet Row', c(None = '.', names(data_viz_table())))
  })
  output$plot_facet_col <- renderUI({
    selectInput('facet_col', 'Facet Column', c(None = '.', names(data_viz_table())))
  })
  
  plot_inputs <- reactive({
    if (!is.null(input$plot_x) && !is.null(input$plot_y)) {
      output <-
        paste(
          input$distance_method,
          input$plot_x,
          input$plot_y,
          input$plot_color_by,
          input$facet_row,
          input$facet_col,
          input$point_size
        )
    }
    else {
      output <- NULL
    }
    return(output)
  })
  
  spree_plot <- eventReactive(plot_inputs() , {
    x_spree <-
      paste("PC", as.character(seq(1, length(pc_table(
        
      )$eig))), sep = " ")
    y_spree <- (pc_table()$eig / sum(pc_table()$eig)) * 100
    
    spree_data <- data.frame(x_spree, y_spree)
    
    p <-
      ggplot(data = spree_data[1:5,], aes(x = x_spree, y = y_spree)) + geom_bar(stat = 'identity') + labs(x = "", y = "Percent Explained") + theme_classic() + coord_cartesian(ylim = c(0, 100)) + theme(
        axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.line.x = element_line(size = 1, color = 'black'),
        axis.line.y = element_line(size = 1, color = 'black'),
        axis.ticks = element_line(size = 1.5)
      ) + scale_y_continuous(expand = c(0, 0))
    return(p)
  })
  output$spree_plot <- renderPlot({
    spree_plot()
  })
  
  data_viz_plot <- eventReactive(plot_inputs(), {
    xval <- input$plot_x
    yval <- input$plot_y
    p <-
      ggplot2::ggplot(data = data_viz_table(), aes_string(x = xval, y = yval)) + geom_point(size = input$point_size) +
      theme_classic() +
      theme(
        axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.line.x = element_line(size = 1, color = 'black'),
        axis.line.y = element_line(size = 1, color = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_line(size = 1.5)
      )
    
    
    if (!is.null(input$plot_color_by)) {
      if (input$plot_color_by != 'None') {
        p <- p + aes_string(color = input$plot_color_by)
        if (input$plot_color_by == input$comparison) {
          colors <- c("blue3", "firebrick3", "darkgreen")
          names(colors) <- c(donor(), recipient(), post_fmt())
          p <- p + scale_colour_manual(values = colors)
        }
      }
    }
    
    if (!is.null(input$facet_row) & !is.null(input$facet_col)) {
      facets <- paste(input$facet_row, '~', input$facet_col, sep = " ")
      if (facets != '. ~ .') {
        p <- p + facet_grid(facets)
      }
    }
    return(p)
  })
  output$data_plot <- renderPlot({
    data_viz_plot()
  })
  
  
  
  
  # Dummy Output
  #output$test <- renderDataTable({
  #  sample_abundance_by_depth()
  #})
  # output$test2 <- renderPrint({
  #   data_viz_table()
  # })
  
})