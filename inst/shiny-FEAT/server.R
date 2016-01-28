require(rhdf5)
require(biom)
require(shiny)
require(dplyr)
require(DT)
require(ggplot2)
require(vegan)
# source('helper_funcitons.R')
ggsave <- ggplot2::ggsave; body(ggsave) <- body(ggplot2::ggsave)[-2]
options(digits = 4)

# Define server logic
shinyServer(function(input, output, session) {


  ##########################################################################################
  # Initial table import & pre-processing                                                  #
  ##########################################################################################
  # Import mapping file
  mapping <- reactive({
    map_file <- input$mapping_file
    if(is.null(map_file)){return()}
    map <- read.delim(file = map_file$datapath)
    map$X.SampleID <- as.character(map$X.SampleID)
    return(map)
  })

  # Import biom table
  biom_first_filtered <- reactive({
    validate(need(!is.null(input$input_otu_table), "Please load an OTU table."))
    raw_table <- input$input_otu_table
    withProgress(message = 'Loading OTU table', value = 0.35, {
    biom_table <- read_biom(raw_table$datapath)
    data <- as(biom_data(biom_table), "matrix")
    biom_only <- as.data.frame(t(data))
    raw_table_otu_count <<- ncol(biom_only)
    colnames(biom_only) <- paste("OTU", colnames(biom_only), sep = "_")
    biom_only$X.SampleID <- as.character(row.names(biom_only))
    row.names(biom_only) <- NULL
    })
    withProgress(message = 'Adding Metadata to OTU table', value = 0.75, {
    output <- inner_join(biom_only, mapping(), by = 'X.SampleID')
    })
    return(output)
  })
  # Get number of OTUs
  num_otus_raw <- reactive({
    if (is.null(biom_first_filtered())) {return()}
    return(raw_table_otu_count)})
  output$raw_otu_count <- renderText({num_otus_raw()})
  # Report Success
  output$zeroeth_check <- renderText({
    if (!is.null(num_otus_raw())) {
    output <- "Table loaded successfully..."
    }
  })

  # Create an experiment-specific OTU table
  biom_experiment_specific <- eventReactive(input$split_into_experiment, {
    validate(need(!is.null(input$input_otu_table) && !is.null(input$mapping_file), "Please load an OTU table and mapping file."))
    tables <- split_into_experiment(biom_first_filtered(), mapping(), input$comparison, input$donor, input$recipient, input$post_fmt)
    return(tables)
  })

  N_donor_samples <- eventReactive(input$split_into_experiment, {
    if (is.null(biom_experiment_specific())) {return()}
    table <- biom_experiment_specific()
    donor <- table[table[,input$comparison] == input$donor, ]
    output <- nrow(donor)
    return(output)
  })
  text_donor_samples <- reactive({
    output <- paste("# Donor Samples Used:", N_donor_samples())
    return(output)
  })
  N_recipient_samples <- eventReactive(input$split_into_experiment, {
    if (is.null(biom_experiment_specific())) {return()}
    table <- biom_experiment_specific()
    recipient <- table[table[,input$comparison] == input$recipient, ]
    output <- nrow(recipient)
    return(output)
  })
  text_recipient_samples <- reactive({
    output <- paste("# Recipient Samples Used:", N_recipient_samples())
    return(output)
  })
  N_post_fmt_samples <- eventReactive(input$split_into_experiment, {
    if (is.null(biom_experiment_specific())) {return()}
    table <- biom_experiment_specific()
    post_fmt <- table[table[,input$comparison] == input$post_fmt, ]
    output <- nrow(post_fmt)
    return(output)
  })
  text_post_fmt_samples <- reactive({
    output <- paste("# Post-FMT Samples Used:", N_post_fmt_samples())
    return(output)
  })

  output$num_donor_samples <- renderText({text_donor_samples()})
  output$num_recipient_samples <- renderText({text_recipient_samples()})
  output$num_post_fmt_samples <- renderText({text_post_fmt_samples()})

  # Report success
  output$second_check <- renderText({
    if (!is.null(biom_experiment_specific())) {
    output <- "Experiment-specific table created..."
    }
    else {
      output <- "..."
    }
  })

  # Normalize and filter the experiment-specific table
  biom_normalized_filtered <- eventReactive(input$normalize_filter, {
    table <- normalize_and_filter(biom_experiment_specific(), input$min_OTU_fraction)
  })
  #Report success
  output$third_check <- renderText({
    if (!is.null(biom_normalized_filtered())) {
    output <- "Experiment-specific table has been normalized and filtered..."
    }
  })

  # Add mapping information
  # Table now has experiment-specific samples from the original biom file, all OTUs that pass preprocessing filters, and all metadata
  biom_table <- eventReactive(input$normalize_filter, {
    validate(need(!is.null(input$input_otu_table) && !is.null(input$mapping_file), "Please load an OTU table and mapping file."))
    output <- merge(biom_normalized_filtered(), mapping()[,c(-1,-2,-3)])
    return(output)
  })

  # After this processing step, report back how many OTUs survive the preprocessing filters
  num_otus_preprocess <- eventReactive(input$normalize_filter, {
    validate(need(!is.null(input$input_otu_table) && !is.null(input$mapping_file), "Please load an OTU table and mapping file."))
    n_col_biom_table <- ncol(biom_table())
    n_col_mapping <- ncol(mapping()[, c(-1,-2,-3)])
    n_otus <- n_col_biom_table - n_col_mapping
    return(ncol(biom_table()) - ncol(mapping()[,c(-1,-2,-3)]))
  })
  output$n_otus_preprocess <- renderText({num_otus_preprocess()})


  ##########################################################################################
  # Define relevant FMT conditions                                                         #
  ##########################################################################################

  # Define mapping conditions to compare to populate drop-down menus
  conditions_of_compare <- reactive({
    if(!is.null(mapping()) && !is.null(input$comparison)){
      category <- mapping()[[input$comparison]]
      return(levels(category))
    } else{return()}
  })

  # Help select FMT details
  output$compare_options <- renderUI({
  selectInput("comparison", "Select comparison category (from mapping file)", names(mapping()), selected = "Compare")
    })
  output$donor <- renderUI({
    selectInput("donor", "Select donor condition", conditions_of_compare(), selected = conditions_of_compare()[1])
  })
  output$recipient <- renderUI({
    selectInput("recipient", "Select pre-transplant recipient condition", conditions_of_compare(), selected = conditions_of_compare()[2])
  })
  output$post_fmt <- renderUI({
    selectInput("post_fmt", "Select post-transplant recipient condition", conditions_of_compare(),selected = conditions_of_compare()[3])
  })


  # Create labels to keep things straight
  donor <- eventReactive(input$split_into_experiment, {return(input$donor)})
  recipient <- eventReactive(input$split_into_experiment, {return(input$recipient)})
  post_fmt <- eventReactive(input$split_into_experiment, {return(input$post_fmt)})
  output$donor_id <- renderText({donor()})
  output$recipient_id <- renderText({recipient()})
  output$post_fmt_id <- renderText({post_fmt()})
  output$donor_id_taxa <- renderText({donor()})
  output$recipient_id_taxa <- renderText({recipient()})
  output$post_fmt_id_taxa <- renderText({post_fmt()})
  output$test_metric <- renderText({input$comparison_test})

  ##########################################################################################
  # Produce table of settings to report back what is being studied                         #
  ##########################################################################################

  # Populate settings table, render it to UI, and allow for download of it
  settings <- reactive({
    validate(need(!is.null(input$input_otu_table) && !is.null(input$mapping_file), "Please load an OTU table and mapping file."))
    Setting <- c("Input OTU Table", "# OTUs in input table",  "Mapping File", "Metadata Category to Select FMT Details", "Donor", "# Donor Samples", "Recipient", "# Recipient Samples", "Post-FMT Recipient", "# Post-FMT Recipient Samples", "Minimum Relative Abundance Filter", "# OTUs After Minimum Relative Abundance Filter", "Fleeting OTU Filter Threshold", "# OTUs After Fleeting OTU Filter (Final # OTUs)", "Comparison Metric", "Exclude OTUs Found Only in Post-FMT Samples or Shared")
    Value <- c(input$input_otu_table$name, num_otus_raw(), input$mapping_file$name, input$comparison, donor(), N_donor_samples(), recipient(), N_recipient_samples(), post_fmt(), N_post_fmt_samples(), input$min_OTU_fraction, num_otus_preprocess(), input$min_fraction, N_otus_after_nonzero_filter(), input$comparison_test, input$remove_OTUs_test_specific)
    df <- data.frame(Setting, Value)
    return(df)
  })
  output$settings <- renderTable({settings()})
  output$x0 <- downloadHandler(filename = function() {'settings.csv'}, content = function(con) {
    write.csv(settings(), con, row.names = FALSE)
  })

  ##########################################################################################
  # Start creating tables specific for the selected conditions                             #
  ##########################################################################################

  # Create tables that are specific for each given condition
  donor_only_table <- eventReactive(input$normalize_filter, {
    validate(need(!is.null(input$input_otu_table) && !is.null(input$mapping_file), "Please load an OTU table and mapping file."))
    table <- get_condition_only_table_function(biom_table(), input$comparison, donor())
    return(table)
  })
  recipient_only_table <- eventReactive(input$normalize_filter, {
    validate(need(!is.null(input$input_otu_table) && !is.null(input$mapping_file), "Please load an OTU table and mapping file."))
    table <- get_condition_only_table_function(biom_table(), input$comparison, recipient())
    return(table)
  })
  post_fmt_only_table <- eventReactive(input$normalize_filter, {
    validate(need(!is.null(input$input_otu_table) && !is.null(input$mapping_file), "Please load an OTU table and mapping file."))
    table <- get_condition_only_table_function(biom_table(), input$comparison, post_fmt())
    return(table)
  })

  # Create a merged table of the given conditions only, and place the column indicating the condition first, send to UI and allow for download
  full_table <- eventReactive(input$normalize_filter, {
    if(is.null(donor_only_table()) | is.null(recipient_only_table()) | is.null(post_fmt_only_table()) | is.null(input$comparison)){return()}
    in_table <- bind_rows(bind_rows(donor_only_table(),recipient_only_table()), post_fmt_only_table())
    out <- table_reorder_first(in_table, input$comparison)
    output <- droplevels(out)
    return(output)
  })
  output$experiment_specific_full_table <- renderDataTable({full_table()}, options = list(lengthMenu = c(50,100,200), pageLength = 50, orderClasses = TRUE))
  output$x1 <- downloadHandler(filename = function() {
       'raw_table.csv'
      }, content = function(con) {write.csv(full_table(), con, row.names = FALSE)})

  ##########################################################################################
  # Query whether OTUs are reliably present in a sufficient number of samples              #
  ##########################################################################################

  # Generate tables for each condition with the fraction of samples in which each OTU has non-zero relative abundances, filtered by the slider input fraction
  donor_nonzero <- reactive({
    in_table <- nonzero_fraction_getter(donor_only_table())
    out <- filter_function(in_table, input$min_fraction, donor())
    return(out)
  })
  recipient_nonzero <- reactive({
    in_table <- nonzero_fraction_getter(recipient_only_table())
    out <- filter_function(in_table, input$min_fraction, recipient())
    return(out)
  })
  post_fmt_nonzero <- reactive({
    in_table <- nonzero_fraction_getter(post_fmt_only_table())
    out <- filter_function(in_table, input$min_fraction, post_fmt())
    return(out)
  })

  # Generate merged table with the fraction of samples in which each OTU has non-zero values filtered by the slider input fraction
  full_nonzero <- reactive({
    if(is.null(donor_nonzero()) | is.null(recipient_nonzero()) | is.null(post_fmt_nonzero())){return()}
    out <- bind_rows(bind_rows(donor_nonzero(), recipient_nonzero()), post_fmt_nonzero())
    return(out)
  })
  output$fraction_samples_with_nonzero_values_for_OTUs <- renderDataTable({full_nonzero()}, options = list(lengthMenu = c(50,100,200), pageLength = 50, orderClasses = TRUE))
  output$x2 <- downloadHandler(filename = function() {
    'fraction_nonzero.csv'
  }, content = function(con) {write.csv(full_nonzero(), con, row.names = FALSE)})

  # Get unique OTUs that survive this filter
  non_overlapping_nonzero <- reactive({
    if(is.null(donor_nonzero()) | is.null(recipient_nonzero()) | is.null(post_fmt_nonzero())){return()}
    otus <- as.data.frame(full_nonzero()$OTU)
    distinct <- distinct(otus)
    return(distinct)
  })

  # Create an output to track how many OTUs survive this filter
  N_otus_after_nonzero_filter <- reactive({
    validate(need(!is.null(input$input_otu_table) && !is.null(input$mapping_file), "Please load an OTU table and mapping file."))
    return(nrow(non_overlapping_nonzero()))
  })
  output$num_otus_after_nonzero_filter <- renderText({N_otus_after_nonzero_filter()})

  ##########################################################################################
  # Use logic joining of these OTUs to define unique and shared OTUs                       #
  ##########################################################################################

  # Produce the tables of OTUs that are unique for each condition
  donor_unique <- reactive({
    out <- anti_join(donor_nonzero(), recipient_nonzero(), by = "OTU")
    return(out)
  })
  recipient_unique <- reactive({
    out <- anti_join(recipient_nonzero(), donor_nonzero(), by = "OTU")
    return(out)
  })
  post_fmt_unique <- reactive({
    out <- anti_join(anti_join(post_fmt_nonzero(), donor_nonzero(), by = "OTU"), recipient_nonzero(), by = "OTU")
    return(out)
  })
  post_fmt_nonunique <- reactive({
    a <- semi_join(post_fmt_nonzero(), donor_unique(), by = "OTU")
    b <- semi_join(post_fmt_nonzero(), recipient_unique(), by = "OTU")
    out <- bind_rows(a, b)
    return(out)
  })
  shared_pre <- reactive({
    output <- semi_join(donor_nonzero(), recipient_nonzero(), by = "OTU")
    return(output)
  })
  shared_throughout <- reactive({
    output <- semi_join(shared_pre(), post_fmt_nonzero(), by = "OTU")
    return(output)
  })

  post_fmt_selected <- reactive({
    if(input$remove_OTUs_test_specific) {
      return(post_fmt_nonunique())
    }
    else {
      return(post_fmt_nonzero())
    }
  })

  # Return numbers of OTUs unique to each condition
  N_otus_unique_donor <- reactive({
    if (is.null(donor_unique())) {return()}
    return(nrow(donor_unique()))})
  output$num_otus_unique_donor <- renderText({N_otus_unique_donor()})

  N_otus_unique_recipient <- reactive({
    if (is.null(recipient_unique())) {return()}
    return(nrow(recipient_unique()))})
  output$num_otus_unique_recipient <- renderText({N_otus_unique_recipient()})

  N_otus_unique_post_fmt <- reactive({
    if (is.null(post_fmt_unique())) {return()}
    return(nrow(post_fmt_unique()))})
  output$num_otus_unique_post_fmt <- renderText({N_otus_unique_post_fmt()})

  N_otus_nonunique_post_fmt <- reactive({
    if (is.null(post_fmt_nonqunique())) {return()}
    return(nrow(post_fmt_nonunique()))})
  output$num_otus_nonunique_post_fmt <- renderText({N_otus_nonunique_post_fmt()})

  N_otus_shared_pre <- reactive({
    if (is.null(shared_pre())) {return()}
    return(nrow(shared_pre()))})
  output$num_otus_shared_pre <- renderText({N_otus_shared_pre()})

  N_otus_shared_throughout <- reactive({
    if (is.null(shared_throughout())) {return()}
    return(nrow(shared_throughout()))})
  output$num_otus_shared_throughout <- renderText({N_otus_shared_throughout()})

  N_otus_post_fmt_nonzero <- reactive({
    if (is.null(post_fmt_nonzero())) {return()}
    return(nrow(post_fmt_nonzero()))
  })
  output$num_otus_post_fmt_nonzero <- renderText({N_otus_post_fmt_nonzero()})

  N_otus_post_fmt_selected <- reactive({
    if (is.null(post_fmt_selected())) {return()}
    return(nrow(post_fmt_selected()))
  })
  output$num_otus_post_fmt_selected <- renderText({N_otus_post_fmt_selected()})

  output$remove_OTUs <- renderText({
    if (input$remove_OTUs_test_specific){
      "Yes"
    }
    else {
      "No"
    }
  })

  #Create a combined table of OTUs unique for the different conditions, including nonzero fraction. Send to UI and allow download
  OTUs_unique_combined <- reactive({
  # Merge these tables together
    if (input$remove_OTUs_test_specific) {
      merged <- bind_rows(donor_unique(), recipient_unique())
    } else {
    merged <- bind_rows(bind_rows(donor_unique(), recipient_unique()), post_fmt_unique())
    }
    output <- merged %>% arrange(Specificity)
    return(output)
  })
  output$otus_unique_to_condition <- renderDataTable({OTUs_unique_combined()}, options = list(lengthMenu = c(50,100,200), pageLength = 50, orderClasses = TRUE))
  output$x3 <- downloadHandler(filename = function() {
    'otus_unique_to_each_condition_combined.csv'
  }, content = function(con) {write.csv(OTUs_unique_combined(), con, row.names = FALSE)})

  ##########################################################################################
  # Compute mean or median OTU relative abundance in all samples of a condition            #
  ##########################################################################################

  # Generate tables with the mean or median OTU relative abundance samples in a given condition. This will be calculated for ALL OTUs, not just those specific for each condition
  donor_relative_abundance <- reactive({
    output <- metric_function(donor_only_table(), "Donor", donor(), input$comparison_test)
    return(output)
  })
  recipient_relative_abundance <- reactive({
    output <- metric_function(recipient_only_table(), "Recipient", recipient(), input$comparison_test)
    return(output)
  })
  post_fmt_relative_abundance <- reactive({
    output <- metric_function(post_fmt_only_table(), "Post_FMT", post_fmt(), input$comparison_test)
    return(output)
  })

  # Combine the tables of mean or median OTU relative abundances, and filter to only display the combined unique OTUs (+/- the unique to post-fmt)
  full_relative_abundance_unique <- reactive({
    mid <- full_join(full_join(donor_relative_abundance()[,-3], recipient_relative_abundance()[,-3], by = "OTU"),  post_fmt_relative_abundance()[,-3], by = "OTU")
    output <- inner_join(OTUs_unique_combined(), mid, by = "OTU")
    return(output)
  })

  ##########################################################################################
  # Compute differences mean or median OTU relative abundances                             #
  ##########################################################################################

  # Compute pairwise differences of the relative abundances
  full_relative_abundance_unique_plus_differences <- reactive({
    output <- full_relative_abundance_unique() %>% mutate(Donor_VS_Recipient = abs(Donor-Recipient)) %>% mutate(Donor_VS_Post_FMT = abs(Donor-Post_FMT)) %>% mutate(Recipient_VS_Post_FMT = abs(Recipient-Post_FMT))
    return(output)
  })

  # Filter this table to only display OTUs that have a minimum differences of relative abundance as given by slider input. Push table to UI and allow download
  filtered_relative_abundance_unique_plus_differences <- reactive({
    output <- filter_by_metric(full_relative_abundance_unique_plus_differences(), input$min_diff, list("Donor_VS_Recipient", "Donor_VS_Post_FMT", "Recipient_VS_Post_FMT"))
    return(output)
  })
  output$filtered_relative_abundance_of_unique_OTUs_plus_differences <- renderDataTable({filtered_relative_abundance_unique_plus_differences()}, options = list(lengthMenu = c(50,100,200), pageLength = 50, orderClasses = TRUE))
  output$x5 <- downloadHandler(filename = function() {
    'filtered_relative_abundance_unique_plus_differences.csv'
  }, content = function(con) {
    write.csv(filtered_relative_abundance_unique_plus_differences(), con, row.names = FALSE)
  })

  # Pull number of OTUs that remain after these steps
  N_otus_after_abundance_difference_filter <- reactive({
    validate(need(!is.null(input$input_otu_table) && !is.null(input$mapping_file), "Please load an OTU table and mapping file."))
    return(nrow(filtered_relative_abundance_unique_plus_differences()))
  })
  output$num_otus_after_abundance_difference_filter <- renderText({N_otus_after_abundance_difference_filter()})


  ##########################################################################################
  # Compute Relevant FMT metrics by OTU                                                    #
  ##########################################################################################

  # FMT_don_table, the table of OTUs in the post-transplant samples that came from the donor.
  FMT_don_table <- reactive({
    validate(need(!is.null(input$input_otu_table) && !is.null(input$mapping_file), "Please load an OTU table and mapping file."))
      output <- semi_join(post_fmt_selected(), donor_unique(), by = "OTU")
    return(output)
  })
  output$FMT_don_table <- renderDataTable({FMT_don_table()})

  # FMT_don, the number of OTUs in the post-transplant samples that came from the donor
  FMT_don <- reactive({
    if (is.null(FMT_don_table())) {return()}
    return(nrow(FMT_don_table()))
  })
  output$FMT_don <- renderText({FMT_don()})

  # D_Frac_FMT the proportion of donor OTUs that made it into the post-transplant samples
  D_Frac_FMT <- reactive({
    if (is.null(donor_unique())) {return()}
    num_donor <- nrow(donor_unique())
    fraction <- FMT_don()/num_donor
    return(fraction)
  })
  output$D_Frac_FMT <- renderText({D_Frac_FMT()})

  # FMT_FracD, the proportion of OTUs in post-transplant samples that came from the donor
  FMT_FracD <- reactive({
    if (is.null(post_fmt_selected())) {return()}
    num_post_fmt <- nrow(post_fmt_selected())
    fraction <- FMT_don()/num_post_fmt
    return(fraction)
  })
  output$FMT_FracD <- renderText({FMT_FracD()})

  # FMT_rec_table, the table of OTUs in the post-transplant samples that came from the recipient
  FMT_rec_table <- reactive({
    validate(need(!is.null(input$input_otu_table) && !is.null(input$mapping_file), "Please load an OTU table and mapping file."))
    output <- semi_join(post_fmt_selected(), recipient_unique(), by = "OTU")
    return(output)
  })
  output$FMT_rec_table <- renderDataTable({FMT_rec_table()})

  # FMT_rec, the number of OTUs in the post-transplant samples that came from the recipient
  FMT_rec <- reactive({
    if (is.null(FMT_rec_table())) {return()}
    return(nrow(FMT_rec_table()))
  })
  output$FMT_rec <- renderText({FMT_rec()})

  # R_Frac_FMT the proportion of recipient OTUs that remained in the post-transplant samples
  R_Frac_FMT <- reactive({
    if (is.null(recipient_unique())) {return()}
    num_recipient <- nrow(recipient_unique())
    fraction <- FMT_rec()/num_recipient
    return(fraction)
  })
  output$R_Frac_FMT <- renderText({R_Frac_FMT()})

  # FMT_FracR, the proportion of OTUs in post-transplant samples that came from the recipient
  FMT_FracR <- reactive({
    if (is.null(post_fmt_selected())) {return()}
    num_post_fmt <- nrow(post_fmt_selected())
    fraction <- FMT_rec()/num_post_fmt
    return(fraction)
  })
  output$FMT_FracR <- renderText({FMT_FracR()})


  ##########################################################################################
  # Include these metrics into a table for download                                        #
  ##########################################################################################

  # Populate table for Metrics
  metric_table <- reactive({
    validate(need(!is.null(input$input_otu_table) && !is.null(input$mapping_file), "Please load an OTU table and mapping file."))
    Item <- c("Input OTU Table", "# OTUs in original table",
                 "Mapping File", "Metadata Category to Select FMT Details",
                 "Donor", "# Donor Samples", "Recipient", "# Recipient Samples", "Post-FMT Recipient", "# Post-FMT Recipient Samples", "Minimum Relative Abundance Filter", "# OTUs After Minimum Relative Abundance Filter",
                 "Fleeting OTU Filter Threshold", "# OTUs After Fleeting OTU Filter (Final # OTUs)", "Comparison Metric",
                 "Exclude OTUs Found Only in Post-FMT Samples or Shared",
              "N_don: # specific & unique to donor", "N_rec: # specific & unique to recipient", "N_post_fmt: # specific to post-FMT",
              "N_post_fmt_unique: # Unique to post-fmt (should be low/zero)", "N_shared: # shared across donor and recipient throughout",
              "FMT_don: # from donor in post-FMT", "D_Frac_FMT: proportion of donor in post-FMT", "FMT_FracD: proportion of post-FMT from donor",
              "FMT_rec: # from recipient in post-FMT", "R_Frac_FMT: proportion of recipient in post-FMT", "FMT_FracR: propotion of post-FMT from recipient")
    Value <- c(input$input_otu_table$name, num_otus_raw(),
               input$mapping_file$name, input$comparison,
               donor(), N_donor_samples(), recipient(), N_recipient_samples(), post_fmt(), N_post_fmt_samples(), input$min_OTU_fraction, num_otus_preprocess(),
               input$min_fraction, N_otus_after_nonzero_filter(), input$comparison_test,
               input$remove_OTUs_test_specific,
               N_otus_unique_donor(), N_otus_unique_recipient(), N_otus_post_fmt_selected(),
               N_otus_unique_post_fmt(), N_otus_shared_throughout(),
               FMT_don(), D_Frac_FMT(), FMT_FracD(),
               FMT_rec(), R_Frac_FMT(), FMT_FracR())
    df <- data.frame(Item, Value)
    return(df)
  })

  output$metrics <- renderTable({metric_table()})
  output$xmetrics <- downloadHandler(filename = function() {
    'transplant_metric_summary.csv'
  }, content = function(con) {
    write.csv(metric_table(), con, row.names = FALSE)
  })

  ##########################################################################################
  #                                                                                        #
  # The following generates a visualization of the FMT metrics                             #
  #                                                                                        #
  ##########################################################################################

  metric_vis <- reactive({
    if (input$remove_OTUs_test_specific) {
      N_unique_vis <- 0
      N_shared_vis <- 0
    }
    else {
      N_unique_vis <- N_otus_unique_post_fmt()
      N_shared_vis <- N_otus_shared_throughout()
    }
    p <- visualize_metrics(N_otus_unique_donor(), N_otus_unique_recipient(), N_otus_post_fmt_selected(), FMT_don(), FMT_rec(), N_unique_vis, N_shared_vis, post_fmt())
    return(p)

  })
  output$metric_visualization <- renderPlot({metric_vis()})
  output$metric_vis_download = downloadHandler(filename = 'metric_vis.png', content = function(file) {
    device <- function(..., width, height) {
      grDevices::png(..., width = 2*width, height = height,
                     res = 500, units = "in")
    }
    ggsave(file, plot = metric_vis(), device = device)
  })


  ##########################################################################################
  #                                                                                        #
  # The following generates a PCOA plot based on a given distance metric                   #
  #                                                                                        #
  ##########################################################################################

  samples_x_unique_otus <- reactive({
    output <- full_table()[,full_relative_abundance_unique()$OTU]
    output$X.SampleID <- full_table()$X.SampleID
    output$Compare <- full_table()$Compare
    output2 <- table_reorder_first(output, "X.SampleID")
    output3 <-table_reorder_first(output2, "Compare")
    output4 <- droplevels(output3)
    return(output4)
  })

  distance_metric_table <- reactive({
    otus_only <- samples_x_unique_otus()[,c(-1,-2)]
    dm <- as.matrix(dist(otus_only, diag = TRUE, upper = TRUE, method = input$distance_method))
    row.names(dm) <- samples_x_unique_otus()$X.SampleID
    colnames(dm) <- samples_x_unique_otus()$X.SampleID
    return(dm)
  })

  pc_table <- reactive({
    pc <- cmdscale(distance_metric_table(), k = nrow(samples_x_unique_otus())-1, eig = TRUE, add = TRUE)
    return(pc)
  })

  pc_values_plus_metadata <- reactive({
    mid <- bind_cols(as.data.frame(pc_table()$points[,1]), as.data.frame(pc_table()$points[,2]), as.data.frame(pc_table()$points[,3]), as.data.frame(pc_table()$points[,4]), as.data.frame(pc_table()$points[,5]))
    colnames(mid) <- c("PC1", "PC2", "PC3", "PC4", "PC5")
    mid$X.SampleID <- as.character(samples_x_unique_otus()$X.SampleID)
    map <- mapping()[,c(-2,-3)]
    map$X.SampleID <- as.character(map$X.SampleID)
    pre_output <- left_join(mid, map, by = "X.SampleID")
    output <- droplevels.data.frame(pre_output)
    return(output)
  })

  output$pc1 <- renderText({
    eig <- pc_table()$eig
    return(percent(eig[1]/sum(eig)))
  })
  output$pc2 <- renderText({
    eig <- pc_table()$eig
    return(percent(eig[2]/sum(eig)))
  })
  output$pc3 <- renderText({
    eig <- pc_table()$eig
    return(percent(eig[3]/sum(eig)))
  })
  output$pc4 <- renderText({
    eig <- pc_table()$eig
    return(percent(eig[4]/sum(eig)))
  })
  output$pc5 <- renderText({
    eig <- pc_table()$eig
    return(percent(eig[5]/sum(eig)))
  })

  output$plot_x <- renderUI({
    selectInput('plot_x', 'X', names(pc_values_plus_metadata()))
  })
  output$plot_y <- renderUI({
    selectInput('plot_y', 'Y', names(pc_values_plus_metadata()), names(pc_values_plus_metadata())[[2]])
  })
  output$plot_color_by <- renderUI({
    selectInput('plot_color_by', 'Color By', c('None', names(pc_values_plus_metadata())), selected = input$comparison)
  })
  output$plot_facet_row <- renderUI({
    selectInput('facet_row', 'Facet Row', c(None='.', names(pc_values_plus_metadata())))
  })
  output$plot_facet_col <- renderUI({
  selectInput('facet_col', 'Facet Column', c(None='.', names(pc_values_plus_metadata())))
  })

  plot_inputs <- reactive({
    if(!is.null(input$plot_x) && !is.null(input$plot_y)){
      output <- paste(input$plot_x, input$plot_y, input$jitter, input$smooth, input$plot_color_by, input$facet_row, input$facet_col, input$point_size)
    }
    else {
      output <- NULL
    }
    return(output)
  })
  plot_item <- eventReactive(plot_inputs(),{

    xval <- input$plot_x
    yval <- input$plot_y
    p <- ggplot(data = pc_values_plus_metadata(), aes_string(x=xval, y=yval)) + geom_point(size = input$point_size) +
      theme_classic() +
      theme(axis.text=element_text(size=14, face="bold"), axis.title=element_text(size=16,face="bold"), axis.line=element_line(size=1), axis.text.x=element_text(angle = 45,hjust = 1), axis.ticks=element_line(size=1.5))


    if (!is.null(input$plot_color_by)) {
      if (input$plot_color_by != 'None') {
        p <- p + aes_string(color=input$plot_color_by)
      }
    }

    if (!is.null(input$facet_row) & !is.null(input$facet_col)) {
      facets <- paste(input$facet_row, '~', input$facet_col, sep = " ")
      if (facets != '. ~ .') {
          p <- p + facet_grid(facets)
      }
    }

    if (input$jitter)
      p <- p + geom_jitter()
    if (input$smooth)
      p <- p + geom_smooth()

    print(p)

  })
  output$plot1 <- renderPlot({plot_item()})
  output$plot_download = downloadHandler(filename = 'plot.png', content = function(file) {
      device <- function(..., width, height) {
        grDevices::png(..., width = 2*width, height = height,
                       res = 500, units = "in")
      }
      ggsave(file, plot = plot_item(), device = device)
    })

  output$pc_table <- renderDataTable({(pc_values_plus_metadata())}, options = list(lengthMenu = c(50,100,200), pageLength = 50, orderClasses = TRUE))
  output$pc_table_download <- downloadHandler(filename = function() {
    'test.csv'
  }, content = function(con) {
    write.csv(pc_values_plus_metadata(), con, row.names = FALSE)
  })



  ##########################################################################################
  #                                                                                        #
  # The following section is for when taxonomy is to be added (Being Updated)              #
  #                                                                                        #
  ##########################################################################################

  # Load taxonomy map
  tax_map <- reactive ({
    if (is.null(input$id_tax_file)) {return()}
    withProgress(message = "Loading Taxonomy", value = 0, {
      raw <- input$id_tax_file
      table <- read.delim(raw$datapath, header = FALSE)
      incProgress(0.5, detail = "Cleaning Up")
      colnames(table) <- c("OTU", "Taxon", "Quality_Score")
      table$OTU <- paste("OTU", table$OTU, sep = "_")
    })
    return(table)
  })

  output$tax_loaded <- renderText({
    if (!is.null(tax_map())) {
      output <- "Taxonomy map loaded successfully..."
    }
  })


  donor_nonzero_taxa <- reactive({
    withProgress(value = 0.5, message = "Collapsing Taxonomy", {
    out <- collapse_taxonomy(inner_join(donor_nonzero(), tax_map(), by = "OTU"))
    })
    return(out)
  })
  recipient_nonzero_taxa <- reactive({
    withProgress(value = 0.5, message = "Collapsing Taxonomy", {
    out <- collapse_taxonomy(inner_join(recipient_nonzero(), tax_map(), by = "OTU"))
    })
    return(out)
  })
  post_fmt_nonzero_taxa <- reactive({
    withProgress(value = 0.5, message = "Collapsing Taxonomy", {
    out <- collapse_taxonomy(inner_join(post_fmt_nonzero(), tax_map(), by = "OTU"))
    })
    return(out)
  })



  donor_unique_taxa <- reactive({
    out <- anti_join(donor_nonzero_taxa(), recipient_nonzero_taxa(), by = "Taxon")
    return(out)
  })
  recipient_unique_taxa <- reactive({
    out <- anti_join(recipient_nonzero_taxa(), donor_nonzero_taxa(), by = "Taxon")
    return(out)
  })
  post_fmt_unique_taxa <- reactive({
    out <- anti_join(anti_join(post_fmt_nonzero_taxa(), donor_nonzero_taxa(), by = "Taxon"), recipient_nonzero_taxa(), by = "Taxon")
    return(out)
  })
  post_fmt_nonunique_taxa <- reactive({
    a <- semi_join(post_fmt_nonzero_taxa(), donor_unique_taxa(), by = "Taxon")
    b <- semi_join(post_fmt_nonzero_taxa(), recipient_unique_taxa(), by = "Taxon")
    out <- bind_rows(a, b)
    return(out)
  })
  shared_pre_taxa <- reactive({
    output <- semi_join(donor_nonzero_taxa(), recipient_nonzero_taxa(), by = "Taxon")
    return(output)
  })
  shared_throughout_taxa <- reactive({
    output <- semi_join(shared_pre_taxa(), post_fmt_nonzero_taxa(), by = "Taxon")
    return(output)
  })

  post_fmt_selected_taxa <- reactive({
    if(input$remove_OTUs_test_specific) {
      return(post_fmt_nonunique_taxa())
    }
    else {
      return(post_fmt_nonzero_taxa())
    }
  })

  # Return numbers of OTUs unique to each condition
  N_otus_unique_donor_taxa <- reactive({
    if (is.null(donor_unique_taxa())) {return()}
    return(nrow(donor_unique_taxa()))})
  output$num_otus_unique_donor_taxa <- renderText({N_otus_unique_donor_taxa()})

  N_otus_unique_recipient_taxa <- reactive({
    if (is.null(recipient_unique_taxa())) {return()}
    return(nrow(recipient_unique_taxa()))})
  output$num_otus_unique_recipient_taxa <- renderText({N_otus_unique_recipient_taxa()})

  N_otus_unique_post_fmt_taxa <- reactive({
    if (is.null(post_fmt_unique_taxa())) {return()}
    return(nrow(post_fmt_unique_taxa()))})
  output$num_otus_unique_post_fmt_taxa <- renderText({N_otus_unique_post_fmt_taxa()})

  N_otus_nonunique_post_fmt_taxa <- reactive({
    if (is.null(post_fmt_nonqunique_taxa())) {return()}
    return(nrow(post_fmt_nonunique_taxa()))})
  output$num_otus_nonunique_post_fmt_taxa <- renderText({N_otus_nonunique_post_fmt_taxa()})

  N_otus_shared_pre_taxa <- reactive({
    if (is.null(shared_pre_taxa())) {return()}
    return(nrow(shared_pre_taxa()))})
  output$num_otus_shared_pre_taxa <- renderText({N_otus_shared_pre_taxa()})

  N_otus_shared_throughout_taxa <- reactive({
    if (is.null(shared_throughout_taxa())) {return()}
    return(nrow(shared_throughout_taxa()))})
  output$num_otus_shared_throughout_taxa <- renderText({N_otus_shared_throughout_taxa()})

  N_otus_post_fmt_nonzero_taxa <- reactive({
    if (is.null(post_fmt_nonzero_taxa())) {return()}
    return(nrow(post_fmt_nonzero_taxa()))
  })
  output$num_otus_post_fmt_nonzero_taxa <- renderText({N_otus_post_fmt_nonzero_taxa()})

  N_otus_post_fmt_selected_taxa <- reactive({
    if (is.null(post_fmt_selected_taxa())) {return()}
    return(nrow(post_fmt_selected_taxa()))
  })
  output$num_otus_post_fmt_selected_taxa <- renderText({N_otus_post_fmt_selected_taxa()})

  output$remove_OTUs_taxa <- renderText({
    if (input$remove_OTUs_test_specific){
      "Yes"
    }
    else {
      "No"
    }
  })

  #Create a combined table of OTUs unique for the different conditions, including nonzero fraction. Send to UI and allow download
  OTUs_unique_combined_taxa <- reactive({
    # Merge these tables together
    if (input$remove_OTUs_test_specific) {
      merged <- bind_rows(donor_unique_taxa(), recipient_unique_taxa())
    } else {
      merged <- bind_rows(bind_rows(donor_unique_taxa(), recipient_unique_taxa()), post_fmt_unique_taxa())
    }
    output <- merged %>% arrange(Specificity)
    return(output)
  })
  output$otus_unique_to_condition_taxa <- renderDataTable({OTUs_unique_combined_taxa()}, options = list(lengthMenu = c(50,100,200), pageLength = 50, orderClasses = TRUE))
  output$x3 <- downloadHandler(filename = function() {
    'otus_unique_to_each_condition_combined_taxa.csv'
  }, content = function(con) {write.csv(OTUs_unique_combined_taxa(), con, row.names = FALSE)})





  # Merge taxonomy information with final output table
  taxonomy_added_final_table <- reactive({
    if (!is.null(tax_map())) {
    mid <- inner_join(filtered_relative_abundance_unique_plus_differences(), tax_map(), by = "OTU")
    output <- table_reorder_first(mid, "Taxon") %>% arrange(Specificity)
    }
    else {output <- filtered_relative_abundance_unique_plus_differences()}
    return(output)
  })

  output$taxonomy_added_final <- renderDataTable(taxonomy_added_final_table(), options = list(lengthMenu = c(50,100,200), pageLength = 50, orderClasses = TRUE))
  output$x7 <- downloadHandler(filename = function() {
    'final_table_plus_taxonomy.csv'
  }, content = function(con) {
    write.csv(taxonomy_added_final_table(), con, row.names = FALSE)
  })

  # Create a table grouped by condition specificity and produce a table giving the number of unique
  # taxonomic assignments for each condition
  distinct_taxa_summary <- reactive({
    if (!is.null(tax_map())) {
    grouping <- group_by(taxonomy_added_final_table(), Specificity)
    return(summarise(grouping, Distinct_Taxa = n_distinct(Taxon)))
    }
    else {return(taxonomy_added_final_table())}
  })

  output$distinct_taxa_summary_table <- renderDataTable(distinct_taxa_summary(), options = list(lengthMenu = c(50,100,200), pageLength = 50, orderClasses = TRUE))
  output$x8 <- downloadHandler(filename = function() {
    'summary_of_distinct_taxa.csv'
  }, content = function(con) {
    write.csv(distinct_taxa_summary(), con, row.names = FALSE)
  })

  # Collapse the final table by shared taxonomy and average the mean/median OTU count & differences
  # for OTUs that share taxonomy
  collapsed_taxonomy <- reactive({
    if (!is.null(tax_map())) {
    return(collapse_taxonomy(taxonomy_added_final_table()))
    }
    else {return(taxonomy_added_final_table())}
  })

  output$collapsed_taxonomy_table <- renderDataTable(collapsed_taxonomy(), options = list(lengthMenu = c(50,100,200), pageLength = 50, orderClasses = TRUE))
  output$x9 <- downloadHandler(filename = function() {
    'table_collapsed_by_taxonomy.csv'
  }, content = function(con) {
    write.csv(collapsed_taxonomy(), con, row.names = FALSE)
  })



  # # Define Transplantation metrics + Taxonomy


  # FMT_don_table, the table of OTUs in the post-transplant samples that came from the donor.
  FMT_don_table_taxa <- reactive({
    validate(need(!is.null(input$input_otu_table) && !is.null(input$mapping_file), "Please load an OTU table and mapping file."))
    output <- semi_join(post_fmt_selected_taxa(), donor_unique_taxa(), by = "Taxon")
    return(output)
  })
  output$FMT_don_table_taxa <- renderDataTable({FMT_don_table_taxa()})

  # FMT_don, the number of OTUs in the post-transplant samples that came from the donor
  FMT_don_taxa <- reactive({
    if (is.null(FMT_don_table_taxa())) {return()}
    return(nrow(FMT_don_table_taxa()))
  })
  output$FMT_don_taxa <- renderText({FMT_don_taxa()})

  # D_Frac_FMT the proportion of donor OTUs that made it into the post-transplant samples
  D_Frac_FMT_taxa <- reactive({
    if (is.null(donor_unique_taxa())) {return()}
    num_donor <- nrow(donor_unique_taxa())
    fraction <- FMT_don_taxa()/num_donor
    return(fraction)
  })
  output$D_Frac_FMT_taxa <- renderText({D_Frac_FMT_taxa()})

  # FMT_FracD, the proportion of OTUs in post-transplant samples that came from the donor
  FMT_FracD_taxa <- reactive({
    if (is.null(post_fmt_selected_taxa())) {return()}
    num_post_fmt <- nrow(post_fmt_selected_taxa())
    fraction <- FMT_don_taxa()/num_post_fmt
    return(fraction)
  })
  output$FMT_FracD_taxa <- renderText({FMT_FracD_taxa()})

  # FMT_rec_table, the table of OTUs in the post-transplant samples that came from the recipient
  FMT_rec_table_taxa <- reactive({
    validate(need(!is.null(input$input_otu_table) && !is.null(input$mapping_file), "Please load an OTU table and mapping file."))
    output <- semi_join(post_fmt_selected_taxa(), recipient_unique_taxa(), by = "Taxon")
    return(output)
  })
  output$FMT_rec_table_taxa <- renderDataTable({FMT_rec_table_taxa()})

  # FMT_rec, the number of OTUs in the post-transplant samples that came from the recipient
  FMT_rec_taxa <- reactive({
    if (is.null(FMT_rec_table_taxa())) {return()}
    return(nrow(FMT_rec_table_taxa()))
  })
  output$FMT_rec_taxa <- renderText({FMT_rec_taxa()})

  # R_Frac_FMT the proportion of recipient OTUs that remained in the post-transplant samples
  R_Frac_FMT_taxa <- reactive({
    if (is.null(recipient_unique_taxa())) {return()}
    num_recipient <- nrow(recipient_unique_taxa())
    fraction <- FMT_rec_taxa()/num_recipient
    return(fraction)
  })
  output$R_Frac_FMT_taxa <- renderText({R_Frac_FMT_taxa()})

  # FMT_FracR, the proportion of OTUs in post-transplant samples that came from the recipient
  FMT_FracR_taxa <- reactive({
    if (is.null(post_fmt_selected_taxa())) {return()}
    num_post_fmt <- nrow(post_fmt_selected_taxa())
    fraction <- FMT_rec_taxa()/num_post_fmt
    return(fraction)
  })
  output$FMT_FracR_taxa <- renderText({FMT_FracR_taxa()})

  ##########################################################################################
  #                                                                                        #
  # The following generates a visualization of the FMT metrics (Taxa)                      #
  #                                                                                        #
  ##########################################################################################

  metric_vis_taxa <- reactive({
    if (input$remove_OTUs_test_specific) {
      N_unique_vis_taxa <- 0
      N_shared_vis_taxa <- 0
    }
    else {
      N_unique_vis_taxa <- N_otus_unique_post_fmt_taxa()
      N_shared_vis_taxa <- N_otus_shared_throughout_taxa()
    }
    p <- visualize_metrics(N_otus_unique_donor_taxa(), N_otus_unique_recipient_taxa(), N_otus_post_fmt_selected_taxa(), FMT_don_taxa(), FMT_rec_taxa(), N_unique_vis_taxa, N_shared_vis_taxa, paste(post_fmt(), "(Taxa)", sep = " "))
    return(p)

  })
  output$metric_visualization_taxa <- renderPlot({metric_vis_taxa()})
  output$metric_vis_download_taxa = downloadHandler(filename = 'metric_vis_taxa.png', content = function(file) {
    device <- function(..., width, height) {
      grDevices::png(..., width = 2*width, height = height,
                     res = 500, units = "in")
    }
    ggsave(file, plot = metric_vis_taxa(), device = device)
  })

  ##########################################################################################
  # Include these metrics into a table for download (Taxa)                                 #
  ##########################################################################################

  # Populate table for Metrics
  metric_table_taxa <- reactive({
    validate(need(!is.null(input$input_otu_table) && !is.null(input$mapping_file), "Please load an OTU table and mapping file."))
    Item <- c("Input OTU Table", "# OTUs in raw table", "Minimum Absolute Count Filter",
              "# OTUs After Minimum Absolute Count Filter", "Mapping File", "Metadata Category to Select FMT Details",
              "Donor", "# Donor Samples", "Recipient", "# Recipient Samples", "Post-FMT Recipient", "# Post-FMT Recipient Samples", "Minimum Relative Abundance Filter", "# OTUs After Minimum Relative Abundance Filter",
              "Fleeting OTU Filter Threshold", "# OTUs After Fleeting OTU Filter (Final # OTUs)", "Comparison Metric",
              "Exclude OTUs Found Only in Post-FMT Samples or Shared",
              "N_don: # specific & unique to donor", "N_rec: # specific & unique to recipient", "N_post_fmt: # specific to post-FMT",
              "N_post_fmt_unique: # Unique to post-fmt (should be low/zero)", "N_shared: # shared across donor and recipient throughout",
              "FMT_don: # from donor in post-FMT", "D_Frac_FMT: proportion of donor in post-FMT", "FMT_FracD: proportion of post-FMT from donor",
              "FMT_rec: # from recipient in post-FMT", "R_Frac_FMT: proportion of recipient in post-FMT", "FMT_FracR: propotion of post-FMT from recipient")
    Value <- c(input$input_otu_table$name, num_otus_raw(), input$min_count,
               num_otus_after_abs_filter(), input$mapping_file$name, input$comparison,
               donor(), N_donor_samples(), recipient(), N_recipient_samples(), post_fmt(), N_post_fmt_samples(), input$min_OTU_fraction, num_otus_preprocess(),
               input$min_fraction, N_otus_after_nonzero_filter(), input$comparison_test,
               input$remove_OTUs_test_specific,
               N_otus_unique_donor_taxa(), N_otus_unique_recipient_taxa(), N_otus_post_fmt_selected_taxa(),
               N_otus_unique_post_fmt_taxa(), N_otus_shared_throughout_taxa(),
               FMT_don_taxa(), D_Frac_FMT_taxa(), FMT_FracD_taxa(),
               FMT_rec_taxa(), R_Frac_FMT_taxa(), FMT_FracR_taxa())
    df <- data.frame(Item, Value)
    return(df)
  })

  output$metrics_taxa <- renderTable({metric_table_taxa()})
  output$xmetrics_taxa <- downloadHandler(filename = function() {
    'transplant_metric_summary_taxa.csv'
  }, content = function(con) {
    write.csv(metric_table_taxa(), con, row.names = FALSE)
  })





  # Handle downloading of zip file of all tables
  output$xall <- downloadHandler(filename = function() {
    paste(input$output_name, ".zip", sep = ".")
  }, content = function(fname) {
    initial_dir <- getwd()
    tmpdir <- tempdir()
    setwd(tempdir())
    fs <- c("settings.csv", "raw_table.csv", "fraction_nonzero.csv", "otus_unique_to_each_condition.csv",
            "filtered_differences_of_metric_table.csv", "otu_taxonomy_map.csv", "final_table_plus_taxonomy.csv",
            "summary_of_distinct_taxa.csv", "table_collapsed_by_taxonomy.csv", "transplant_metric_summary.csv",
            "transplant_metric_summary_plus_taxa.csv", "OTUs_for_qiime.txt")
    write.csv(settings(), file = "settings.csv", row.names = FALSE)
    write.csv(full_table(), file = "raw_table.csv", row.names = FALSE)
    write.csv(full_nonzero(), file = "fraction_nonzero.csv", row.names = FALSE)
    write.csv(OTUs_unique_each(), file = "otus_unique_to_each_condition.csv", row.names = FALSE)
    write.csv(otu_taxonomy(), file = "otu_taxonomy_map.csv", row.names = FALSE)
    write.csv(taxonomy_added_final_table(), file = "final_table_plus_taxonomy.csv", row.names = FALSE)
    write.csv(distinct_taxa_summary(), file = "summary_of_distinct_taxa.csv", row.names = FALSE)
    write.csv(collapsed_taxonomy(), file = "table_collapsed_by_taxonomy.csv", row.names = FALSE)
    write.csv(metric_table(), file = "transplant_metric_summary.csv", row.names = FALSE)
    write.csv(metric_table_taxa(), file = "transplant_metric_summary_plus_taxa.csv", row.names = FALSE)
    write(otus_for_qiime(), sep = "\n", file = "OTUs_for_qiime.txt", row.names = FALSE)
    zip(zipfile=fname, files=fs)
    setwd(initial_dir)
    if(file.exists(paste0(fname, ".zip"))) {file.rename(paste0(fname, ".zip"), fname)}
  }, contentType = "application/zip")









})
