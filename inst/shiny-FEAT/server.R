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
require(BiomassWorkflow)

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
      map$biomass_ratio <- 1
    }
    
    return(map)
  })
  
  # Create an output to toggle FMT detail selection UI
  output$files_loaded <- reactive({
    return(!(is.null(metadata()) | is.null(taxonomy())))
  })
  outputOptions(output, 'files_loaded', suspendWhenHidden = FALSE)
  
  output$depth <- renderUI({
    if (!is.null(input$id_tax_file)) {
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
    } else {
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
      return(data.frame())
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
  
  output$check_taxonomy <- renderUI({
    selectInput(
      'toggle_taxonomy',
      label = 'Incorporate Taxonomy?',
      choices = c(
        "Work with Taxonomy" = TRUE,
        "Work with OTUs" = FALSE
      ),
      selected = !is.null(input$id_tax_file)
    )
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
    if (i > 2) {
      table$short_label <-
        paste(table[[phylogeny[i - 2]]], table[[phylogeny[i - 1]]], table[[phylogeny[i]]], sep = ".")
      table <-
        tidyr::unite_(table,
                      'long_label',
                      phylogeny[1:i],
                      sep = ".",
                      remove = F)
    } else if (i > 1) {
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
  output$abundance_selected_taxonomy <-
    renderDataTable({
      abundance_selected_donor_specific()
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
  
  # Engraftment metric - 2 log of P_Donor/P_Recipient
  Engraftment <- eventReactive(input$go, {
    if (is.null(abundance_selected_post_fmt_specific())) {
      return()
    }
    engraftment_metric <- log(P_Donor() / P_Recipient(), 2)
    return(engraftment_metric)
  })
  
  output$Engraftment <- renderUI({
    HTML(paste("<strong>E:</strong>",
               round(Engraftment(), 3),
               sep = " "))
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
  
  
  donor_unique_table <- eventReactive(input$go, {
    metadata_added <-
      dplyr::left_join(sample_abundance_by_depth(), metadata_only(), by = "X.SampleID")
    
    donor_full <-
      metadata_added[metadata_added[, input$comparison] ==  donor(), ]
    donor_grouped <-
      dplyr::group_by(donor_full, long_label, short_label)
    donor_summarized <-
      dplyr::summarize(donor_grouped, MeanAbundance = mean(abundance))
    donor_output <-
      dplyr::filter(donor_summarized, long_label %in% donor_unique()[, 1])
    if (input$abundance_type == 'absolute') {
      missing_abundance <-
        mean(donor_full$biomass_ratio) - sum(donor_output$MeanAbundance)
    } else {
      missing_abundance <- 1 - sum(donor_output$MeanAbundance)
    }
    missing_table <-
      data.frame(
        long_label = "Shared/Other",
        short_label = "Shared/Other",
        MeanAbundance = missing_abundance
      )
    output <- dplyr::bind_rows(donor_output, missing_table)
    output$Condition <- "Donor"
    colnames(output) <-
      c("Taxon", "ShortName", "MeanAbundance", "Condition")
    output$MeanAbundance <- round(output$MeanAbundance, 4)
    return(output)
  })
  
  recipient_unique_table <- eventReactive(input$go, {
    metadata_added <-
      dplyr::left_join(sample_abundance_by_depth(), metadata_only(), by = "X.SampleID")
    
    recipient_full <-
      metadata_added[metadata_added[, input$comparison] ==  recipient(), ]
    recipient_grouped <-
      dplyr::group_by(recipient_full, long_label, short_label)
    recipient_summarized <-
      dplyr::summarize(recipient_grouped, MeanAbundance = mean(abundance))
    recipient_output <-
      dplyr::filter(recipient_summarized, long_label %in% recipient_unique()[, 1])
    if (input$abundance_type == 'absolute') {
      missing_abundance <-
        mean(recipient_full$biomass_ratio) - sum(recipient_output$MeanAbundance)
    } else {
      missing_abundance <- 1 - sum(recipient_output$MeanAbundance)
    }
    missing_table <-
      data.frame(
        long_label = "Shared/Other",
        short_label = "Shared/Other",
        MeanAbundance = missing_abundance
      )
    output <- dplyr::bind_rows(recipient_output, missing_table)
    output$Condition <- "Recipient"
    colnames(output) <-
      c("Taxon", "ShortName", "MeanAbundance", "Condition")
    output$MeanAbundance <- round(output$MeanAbundance, 4)
    return(output)
  })
  
  post_fmt_table <- eventReactive(input$go, {
    metadata_added <-
      dplyr::left_join(sample_abundance_by_depth(), metadata_only(), by = "X.SampleID")
    
    post_fmt_full <-
      metadata_added[metadata_added[, input$comparison] ==  post_fmt(), ]
    post_fmt_grouped <-
      dplyr::group_by(post_fmt_full, long_label, short_label)
    post_fmt_summarized <-
      dplyr::summarize(post_fmt_grouped, MeanAbundance = mean(abundance))
    post_fmt_output <-
      dplyr::filter(post_fmt_summarized,
                    long_label %in% abundance_selected_post_fmt_specific()[, 1])
    if (input$abundance_type == 'absolute') {
      missing_abundance <-
        mean(post_fmt_full$biomass_ratio) - sum(post_fmt_output$MeanAbundance)
    } else {
      missing_abundance <- 1 - sum(post_fmt_output$MeanAbundance)
    }
    missing_table <-
      data.frame(
        long_label = "Lost/Other",
        short_label = "Lost/Other",
        MeanAbundance = missing_abundance
      )
    
    output <- dplyr::bind_rows(post_fmt_output, missing_table)
    output$Condition <- "Post-FMT"
    colnames(output) <-
      c("Taxon", "ShortName", "MeanAbundance", "Condition")
    output$MeanAbundance <- round(output$MeanAbundance, 4)
    return(output)
  })
  
  donor_engrafted_table <- eventReactive(input$go, {
    output <-
      donor_unique_table()[donor_unique_table()[["Taxon"]] %in% P_Donor_table()[, 1],]
    return(output)
  })
  
  recipient_persisted_table <- eventReactive(input$go, {
    output <-
      recipient_unique_table()[recipient_unique_table()[["Taxon"]] %in% P_Recipient_table()[, 1],]
    return(output)
  })
  
  post_fmt_donor_table <- eventReactive(input$go, {
    output <-
      post_fmt_table()[post_fmt_table()[["Taxon"]] %in% P_Donor_table()[, 1],]
    return(output)
  })
  post_fmt_recipient_table <- eventReactive(input$go, {
    output <-
      post_fmt_table()[post_fmt_table()[["Taxon"]] %in% P_Recipient_table()[, 1],]
    return(output)
  })
  post_fmt_unique_table <- eventReactive(input$go, {
    output <-
      post_fmt_table()[post_fmt_table()[["Taxon"]] %in% post_fmt_unique()[, 1],]
    return(output)
  })
  shared_throughout_table <- eventReactive(input$go, {
    output <-
      post_fmt_table()[post_fmt_table()[["Taxon"]] %in% shared_throughout()[, 1],]
    return(output)
  })
  
  output$click_info <- renderDataTable({
    if (!is.null(input$id_tax_file)) {
      if (input$toggle_taxonomy) {
        click <- input$plot_click
        maximum <- max(c(N_Donor(), N_Recipient(), N_P_Total()))
        buffer <- max(c(1, maximum / 5))
        
        if (!is.null(click)) {
          if (click$x > buffer &
              click$x < (N_Donor() + buffer) &
              click$y < 45 & click$y > 30) {
            # Donor Unique
            output <- donor_unique_table()
          } else if (click$x > 2 * buffer + N_Donor() &
                     click$x < (N_Recipient() + N_Donor() + 2 * buffer) &
                     click$y < 45 & click$y > 30) {
            # Recipient Unique
            output <- recipient_unique_table()
          } else if (click$x > 3 * buffer + N_Donor() + N_Recipient() &
                     click$x < (N_P_Total() + N_Recipient() + N_Donor() + 3 * buffer) &
                     click$y < 45 & click$y > 30) {
            # Post-FMT
            output <- post_fmt_table()
          } else if (click$x > buffer &
                     click$x < (N_P_Donor() + buffer) &
                     click$y < 25 & click$y > 10) {
            # Donor Engrafted
            output <- donor_engrafted_table()
          } else if (click$x > (N_Recipient() + N_Donor() + 2 * buffer - N_P_Recipient()) &
                     click$x < (N_Recipient() + N_Donor() + 2 * buffer) &
                     click$y < 25 & click$y > 10) {
            # Recipient Engrafted
            output <- recipient_persisted_table()
          } else if (click$x > 3 * buffer + N_Donor() + N_Recipient() &
                     click$x < (N_P_Donor() + N_Recipient() + N_Donor() + 3 * buffer) &
                     click$y < 25 & click$y > 10) {
            # Post-FMT Donor
            output <- post_fmt_donor_table()
          } else if (click$x > 3 * buffer + N_Donor() + N_Recipient() + N_P_Donor() &
                     click$x < (N_P_Donor() + N_P_Recipient() + N_Recipient() + N_Donor() + 3 *
                                buffer) &
                     click$y < 25 & click$y > 10) {
            # Post-FMT Recipient
            output <- post_fmt_recipient_table()
          } else if (click$x > 3 * buffer + N_Donor() + N_Recipient() + N_P_Donor() + N_P_Recipient() &
                     click$x < (
                       N_P_Donor() + N_P_Recipient() + N_P_Unique() + N_Recipient() + N_Donor() + 3 *
                       buffer
                     ) & click$y < 25 & click$y > 10) {
            # Post-FMT Unique
            output <- post_fmt_unique_table()
          } else if (click$x > 3 * buffer + N_Donor() + N_Recipient() + N_P_Donor() + N_P_Recipient() + N_P_Unique() &
                     click$x < (N_P_Total() + N_Recipient() + N_Donor() + 3 * buffer) &
                     click$y < 25 & click$y > 10) {
            # Post-FMT Shared
            output <- shared_throughout_table()
          } else {
            return()
          }
          
          return(output)
        } else {
          return()
        }
      } else {
        return()
      }
    } else {
      return()
    }
    
  },
  options = list(
    lengthMenu = c(50, 100, 200),
    pageLength = 50,
    orderClasses = TRUE
  ))
  
  
  output$click_plot <- renderPlot({
    if (!is.null(input$id_tax_file)) {
      if (input$toggle_taxonomy) {
        click <- input$plot_click
        maximum <- max(c(N_Donor(), N_Recipient(), N_P_Total()))
        buffer <- max(c(1, maximum / 5))
        
        if (!is.null(click)) {
          if (click$x > buffer &
              click$x < (N_Donor() + buffer) &
              click$y < 45 & click$y > 30) {
            # Donor Unique
            plot_table <- donor_unique_table()
          } else if (click$x > 2 * buffer + N_Donor() &
                     click$x < (N_Recipient() + N_Donor() + 2 * buffer) &
                     click$y < 45 & click$y > 30) {
            # Recipient Unique
            plot_table <- recipient_unique_table()
          } else if (click$x > 3 * buffer + N_Donor() + N_Recipient() &
                     click$x < (N_P_Total() + N_Recipient() + N_Donor() + 3 * buffer) &
                     click$y < 45 & click$y > 30) {
            # Post-FMT
            plot_table <- post_fmt_table()
          } else if (click$x > buffer &
                     click$x < (N_P_Donor() + buffer) &
                     click$y < 25 & click$y > 10) {
            # Donor Engrafted
            plot_table <- donor_engrafted_table()
          } else if (click$x > (N_Recipient() + N_Donor() + 2 * buffer - N_P_Recipient()) &
                     click$x < (N_Recipient() + N_Donor() + 2 * buffer) &
                     click$y < 25 & click$y > 10) {
            # Recipient Engrafted
            plot_table <- recipient_persisted_table()
          } else if (click$x > 3 * buffer + N_Donor() + N_Recipient() &
                     click$x < (N_P_Donor() + N_Recipient() + N_Donor() + 3 * buffer) &
                     click$y < 25 & click$y > 10) {
            # Post-FMT Donor
            plot_table <- post_fmt_donor_table()
          } else if (click$x > 3 * buffer + N_Donor() + N_Recipient() + N_P_Donor() &
                     click$x < (N_P_Donor() + N_P_Recipient() + N_Recipient() + N_Donor() + 3 *
                                buffer) &
                     click$y < 25 & click$y > 10) {
            # Post-FMT Recipient
            plot_table <- post_fmt_recipient_table()
          } else if (click$x > 3 * buffer + N_Donor() + N_Recipient() + N_P_Donor() + N_P_Recipient() &
                     click$x < (
                       N_P_Donor() + N_P_Recipient() + N_P_Unique() + N_Recipient() + N_Donor() + 3 *
                       buffer
                     ) & click$y < 25 & click$y > 10) {
            # Post-FMT Unique
            plot_table <- post_fmt_unique_table()
          } else if (click$x > 3 * buffer + N_Donor() + N_Recipient() + N_P_Donor() + N_P_Recipient() + N_P_Unique() &
                     click$x < (N_P_Total() + N_Recipient() + N_Donor() + 3 * buffer) &
                     click$y < 25 & click$y > 10) {
            # Post-FMT Unique
            plot_table <- shared_throughout_table()
          } else {
            plot_table <- NULL
          }
          
          if (!is.null(plot_table)) {
            if (input$abundance_type == 'relative') {
              plot_table$MeanAbundance <- plot_table$MeanAbundance * 100
              ylabel <- "Relative Abundance (% of OTUs)\n"
              ylimit <- 100
            } else {
              ylabel <-
                "Absolute Abundance \n(ug DNA for each Taxa per mg Feces)\n"
              ylimit <- max(metadata_only()$biomass_ratio)
            }
            
            i <- grep(input$selected_depth, phylogeny)
            legend_label <-
              paste(
                "Taxon (",
                phylogeny[i - 2],
                ".",
                phylogeny[i - 1],
                ".",
                input$selected_depth,
                ")",
                sep = ""
              )
            plot <-
              ggplot(data = plot_table,
                     aes(
                       x = Condition,
                       y = MeanAbundance,
                       fill = ShortName
                     )) + geom_bar(color =
                                     "black", stat = 'identity') + EJC_theme() + scale_fill_manual(name = legend_label, values = EJC_colors) + labs(x = "",
                                                                                                                                                    y = ylabel,
                                                                                                                                                    title = "") + theme(
                                                                                                                                                      legend.position = 'right',
                                                                                                                                                      legend.direction = 'vertical',
                                                                                                                                                      aspect.ratio = 2
                                                                                                                                                    ) + guides(fill = guide_legend(reverse = TRUE)) + coord_cartesian(ylim = c(0, ylimit))
            final_output <- plot
            
            return(final_output)
          }
        } else {
          return()
        }
      } else {
        return()
      }
    } else {
      return()
    }
  })
  
  ########
  
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
  
  percent_explained_table <- reactive({
    eig <- pc_table()$eig
    PC <- c("PC1", "PC2", "PC3", "PC4", "PC5")
    PercentExplained <- percent(eig / sum(eig))[1:5]
    df <- data.frame(PC, PercentExplained)
    return(df)
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
          input$point_size,
          input$plot_type
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
    plot <-
      ggplot2::ggplot(data = data_viz_table(), aes_string(x = xval, y = yval))
    if ('Boxplot' %in% input$plot_type) {
      plot <- plot + geom_boxplot()
    }
    
    if ('Scatter' %in% input$plot_type) {
      plot <- plot + geom_point(size = input$point_size)
    }
    p <- plot +
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
          p <-
            p + scale_colour_manual(values = colors)
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
  
  # Parameter Summary Table
  
  param_summary_table <- reactive({
    validate(need(
      !is.null(input$input_otu_table) &&
        !is.null(input$metadata_file),
      "Please load an OTU table and mapping file."
    ))
    Item <- c(
      "Input OTU Table",
      "Metadata File",
      'Taxonomy File',
      "Metadata Category",
      "Donor Condition",
      "Recipient Condition",
      "Post-FMT Condition",
      "Minimum Relative Abundance Filter",
      "Fleeting Filter",
      "Abundance Type",
      "Incorporate Taxonomy?",
      "Taxonomic Depth",
      "Number Donor Samples",
      "Number Recipient Samples",
      "Number Post-FMT Samples",
      "Starting OTUs",
      "OTUs After Relative Abundance Filer",
      "OTUs (or Taxa) after Fleeting OTU Filter",
      "Distance Metric for PCA"
    )
    
    Value <-
      c(
        input$input_otu_table$name,
        input$metadata_file$name,
        input$id_tax_file$name,
        input$comparison,
        donor(),
        recipient(),
        post_fmt(),
        input$min_OTU_fraction,
        input$min_fraction,
        input$abundance_type,
        input$toggle_taxonomy,
        input$selected_depth,
        N_donor_samples(),
        N_recipient_samples(),
        N_post_fmt_samples(),
        num_otus_raw(),
        num_otus_relative_filtered(),
        N_otus_after_fleeting_filter(),
        input$distance_method
      )
    
    df <- data.frame(Item, Value)
    return(df)
  })
  
  # Metric Summary Table
  metric_summary_table <- reactive({
    validate(need(
      !is.null(input$input_otu_table) &&
        !is.null(input$metadata_file),
      "Please load an OTU table and mapping file."
    ))
    Item <- c(
      "Donor",
      "Recipient",
      "Post-FMT",
      "Taxonomy Added?",
      "N_Donor",
      "N_Recipient",
      "N_P_Total",
      "N_P_Unique",
      "N_P_Shared",
      "N_Shared_Lost",
      "N_P_Donor",
      "N_P_Recipient",
      "D_Engraft",
      "R_Persist",
      "P_Donor",
      "P_Recipient",
      "P_Unique",
      "P_Shared",
      "Engraft"
    )
    
    
    Value <- c(
      donor(),
      recipient(),
      post_fmt(),
      input$toggle_taxonomy,
      N_Donor(),
      N_Recipient(),
      N_P_Total(),
      N_P_Unique(),
      N_P_Shared(),
      N_otus_shared_lost(),
      N_P_Donor(),
      N_P_Recipient(),
      D_Engraft(),
      R_Persist(),
      P_Donor(),
      P_Recipient(),
      P_Unique(),
      P_Shared(),
      Engraftment()
    )
    
    df <- data.frame(Item, Value)
    return(df)
    
  })
  
  ## Download Handler
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(donor(),
            "_into_",
            recipient(),
            "_output_files",
            ".zip",
            sep = "")
    },
    content = function(fname) {
      fs <- c()
      tmpdir <- tempdir()
      setwd(tempdir())
      
      # Parameter Summary
      param_summary_path <- 'Parameter_Summary.csv'
      write.csv(x = param_summary_table(), param_summary_path, row.names = FALSE)
      
      # Metric Summary
      metric_summary_path <- 'Metric_Summary.csv'
      write.csv(x = metric_summary_table(),
                metric_summary_path,
                row.names = FALSE)
      
      #Metric Vizualization
      metric_vis_path <- 'Transplant_Visualization.pdf'
      save_plot(metric_vis_path,
                metric_vis(),
                base_height = 8,
                base_width = 15)
      
      # Donor Unique Table
      donor_unique_path <- "Donor_Unique_Taxa.csv"
      write.csv(x = donor_unique_table() , donor_unique_path, row.names = FALSE)
      
      # Recipient Unique Table
      recipient_unique_path <- "Recipient_Unique_Taxa.csv"
      write.csv(x = recipient_unique_table() ,
                recipient_unique_path,
                row.names = FALSE)
      
      # Post-FMT Table
      post_fmt_full_path <- "Post-FMT_Full_Taxa.csv"
      write.csv(x = post_fmt_table() , post_fmt_full_path, row.names = FALSE)
      
      # Donor Engrafted Table
      donor_engrafted_path <- "Donor_Engrafted_Taxa.csv"
      write.csv(x = donor_engrafted_table() ,
                donor_engrafted_path,
                row.names = FALSE)
      
      # Recipient Persisted Table
      recipient_persisted_path <- "Recipient_Persisted_Taxa.csv"
      write.csv(x = recipient_persisted_table() ,
                recipient_persisted_path,
                row.names = FALSE)
      
      # Post-FMT Donor Table
      post_fmt_donor_path <- "Post-FMT_Donor_Taxa.csv"
      write.csv(x = post_fmt_donor_table() ,
                post_fmt_donor_path,
                row.names = FALSE)
      
      # Post-FMT Recipient Table
      post_fmt_recipient_path <- "Post-FMT_Recipient_Taxa.csv"
      write.csv(x = post_fmt_recipient_table() ,
                post_fmt_recipient_path,
                row.names = FALSE)
      
      # Post-FMT Unique Table
      post_fmt_unique_path <- "Post-FMT_Unique_Taxa.csv"
      write.csv(x = post_fmt_unique_table() ,
                post_fmt_unique_path,
                row.names = FALSE)
      
      # Shared Throughout Table
      shared_throughout_path <- "Shared_Taxa.csv"
      write.csv(x = shared_throughout_table() ,
                shared_throughout_path,
                row.names = FALSE)
      
      # PCA Plot
      pca_plot_path <- "PCA_Plot.pdf"
      pdf(pca_plot_path, height = 6, width = 8)
      colors <- c("blue3", "firebrick3", "darkgreen")
      names(colors) <- c(donor(), recipient(), post_fmt())
      p <-
        ggplot2::ggplot(data = data_viz_table(), aes(x = PC1, y = PC2, colour = Compare)) + geom_point(size = 4) + theme_classic() +
        theme(
          axis.text = element_text(size = 14, face = "bold"),
          axis.title = element_text(size = 16, face = "bold"),
          axis.line.x = element_line(size = 1, color = 'black'),
          axis.line.y = element_line(size = 1, color = 'black'),
          axis.text.x = element_text(angle = 0, hjust = 0.5),
          axis.ticks = element_line(size = 1.5)
        ) + scale_colour_manual(values = colors)
      print(p)
      dev.off()
      
      # Percent Explained
      percent_explained_path <- "PCA_Percent_Explained.csv"
      write.csv(percent_explained_table(),
                percent_explained_path,
                row.names = FALSE)
      
      # Plotting_Data
      plot_data_path <- "Plot_Data_Table.csv"
      write.csv(data_viz_table(), plot_data_path, row.names = FALSE)
      
      fs <-
        c(
          param_summary_path,
          metric_summary_path,
          metric_vis_path,
          donor_unique_path,
          recipient_unique_path,
          post_fmt_full_path,
          donor_engrafted_path,
          recipient_persisted_path,
          post_fmt_donor_path,
          post_fmt_recipient_path,
          post_fmt_unique_path,
          shared_throughout_path,
          pca_plot_path,
          percent_explained_path,
          plot_data_path
        )
      
      zip(zipfile = fname, files = fs)
    },
    contentType = "application/zip"
  )
  
  
  
  
  # Dummy Output
  #output$test <- renderDataTable({
  #  sample_abundance_by_depth()
  #})
  # output$test2 <- renderPrint({
  #   data_viz_table()
  # })
  
})