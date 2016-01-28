#' Launch shiny app to analyze FMT efficacy
#' @return nothing - will launch Shiny app in viewer
#' @export
#'
launchFEAT <- function() {
  appDir <- system.file("shiny-FEAT", package = "FEAT")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `FEAT`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal", launch.browser = TRUE)
}


#' Within Shiny environment, split a data table into an experiment-specific table given the input parameters
#' @return output table with only data from selected conditions
#' @export
#'
split_into_experiment <- function(biom_table, mapping, compare, donor, recipient, post_fmt) {
  withProgress(message = 'Selecting experiment samples', value = 0, {
  table_out <- biom_table[biom_table[,compare] %in% c(donor, recipient, post_fmt), ]
  })
  return(table_out)
}

#' Within Shiny environment, normalize OTU counts for your OTU table, and impose a minimum filter
#' @return output OTU table with relative abundances above the provided minimum
#' @export
#'
normalize_and_filter <- function(biom_table, min_fraction) {
  withProgress(message = "Normalizing OTU counts and Filtering", value = 0.2, {
  metadata <- biom_table[, !grepl("OTU_", names(biom_table))]
  table_only <- biom_table[, grepl("OTU_", names(biom_table))]
  incProgress(0.4, detail = "Normalizing")
  normalized <- table_only/rowSums(table_only)
  max_otu_fraction <- apply(normalized, 2, max)
  incProgress(0.2, detail = "Filtering")
  otu_filtered <- normalized[, max_otu_fraction > min_fraction]
  otu_filtered$X.SampleID <- biom_table$X.SampleID
  output <- inner_join(otu_filtered, metadata, by = "X.SampleID")
  })
  return(output)
}


#' Function to get tables that are specific for each given condition
#' @return data table with only one condition represented
#' @export
#'
get_condition_only_table_function <- function(table,
                                              category, value) {
  output_table <- table %>% filter(get(category, envir = as.environment(table)) ==
                                     value)
  return(output_table)
}

#' function to reorder tables with a given value first
#' @return table with selected column moved to the front
#' @export
#'
table_reorder_first <- function(table, value) {
  cols <- c(value, names(table[-which(names(table) ==
                                        value)]))
  output_table <- data.frame(table[cols], row.names = NULL)
  return(output_table)
}

#' function that returns a table of OTUs and the fraction of samples in which they have non-zero values
#' @return table of OTUs and the fraction of samples in which they have non-zero values
#' @export
#'
nonzero_fraction_getter <- function(table) {
  nonzero_fraction <- as.data.frame(t(table %>% summarise_each(funs(fraction = (sum(. >
                                                                                      0)/length(.))), matches("OTU"))))
  nonzero_fraction$OTU <- row.names(nonzero_fraction)
  colnames(nonzero_fraction) <- c("fraction", "OTU")
  nonzero_fraction <- nonzero_fraction %>% select(OTU,
                                                  fraction)
}

#'  Define a function to produce tables populated only by OTUs that have non-zero counts in a given fraction of all of the samples for that condition
#' @return tables populated only by OTUs that have non-zero counts in a given fraction of all of the samples for that condition
#' @export
filter_function <- function(table, min_fraction, condition) {
  output_table <- table %>% filter(fraction >= min_fraction) %>%
    mutate(Specificity = condition)
  return(output_table)
}

#' Define a function that returns the mean or median OTU count in samples of a given condition
#' @return the mean or median OTU count in samples of a given condition
#' @export
#'
metric_function <- function(table, label, condition,
                            metric) {

  if (metric == "Median") {
      output_table <- as.data.frame(t(table %>% summarise_each(funs(median), matches("OTU"))))
      colnames(output_table) <- label
  }
  else if (metric == "Mean") {
      output_table <- as.data.frame(t(table %>% summarise_each(funs(mean), matches("OTU"))))
      colnames(output_table) <- label
  }
  output_table <- add_rownames(output_table, var = "OTU")
  return(output_table)
}


#' Define a function that filters tables by a minimum difference between means/medians
#' @return table filtered for a minimum difference between categories
#' @export
#'
filter_by_metric <- function(table, min_diff, difference_labels) {
  criteria <- list(interp(~one >= d | two >= d | three >= d,
                          .values = list(one = as.name(difference_labels[[1]]),
                                              two = as.name(difference_labels[[2]]), three = as.name(difference_labels[[3]]),
                                              d = min_diff)))
  output_table <- table %>% filter_(.dots = criteria)
  return(output_table)
}

#' Define a function that will collapse a table by shared taxonomy and average the mean/median OTU count & differences for OTUs that share taxonomy
#' @return table collapsed by taxonomy
#' @export
#'
collapse_taxonomy <- function(table) {
  mixed_column <- table %>% mutate(Sample_x_Taxon = paste(Specificity,
                                                             Taxon, sep = "_"))
  collapse <- mixed_column %>% group_by(Sample_x_Taxon) %>%
    summarise_each(funs(means_of_numbers_only)) %>%
    select(-c(OTU, Sample_x_Taxon, fraction, Quality_Score)) %>%
    arrange(Specificity)
  output <- table_reorder_first(table_reorder_first(collapse,
                                                    "Specificity"), "Taxon")
  return(output)
}

#' Define a function that will take a vector and return the mean of the vector only if it is a vector of numbers but returns the first instance of the vector when it is not a number
#' @return an average of only the numbers of a given vector
#' @export
#'
means_of_numbers_only <- function(x) {
  if (class(x) == "numeric") {
    return(mean(x))
  }
  else {
    return(x[1])
  }
}

#' Make a percentage
#' @return a percentage from a decimal
#' @export
#'
percent <- function(x, digits = 2, format = "f", ...) {
  paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
}

#' Function to produce FMT metric visualizations
#' @return a plot that depicts FMT efficacy
#' @export
#'
visualize_metrics <- function(N_donor, N_recipient, N_postFMT, FMT_donor, FMT_recipient, N_unique_postfmt, N_shared_postfmt, title_text) {

  N_don <<- N_donor
  N_rec <<- N_recipient
  N_FMT <<- N_postFMT
  FMT_don <<- FMT_donor
  FMT_rec <<- FMT_recipient
  N_unique <<- N_unique_postfmt
  N_shared <<- N_shared_postfmt
  D_frac_FMT <<- round(FMT_don/N_don, 3)
  R_frac_FMT <<- round(FMT_rec/N_rec, 3)
  FMT_FracD <<- round(FMT_don/N_FMT, 3)
  FMT_FracR <<- round(FMT_rec/N_FMT, 3)
  FMT_unique <<- round(N_unique/N_FMT, 3)
  FMT_shared <<- round(N_shared/N_FMT, 3)


  title_text <<- paste("Donor:", "A", "; Recipient:", "B", sep = " ")

  maximum <<- max(c(N_don, N_rec, N_FMT))

  buffer <<- max(c(1, maximum/5))

  # Create empty plot
  p <- ggplot()

  # Add rectangles

  ## Donor
  p <- p + annotate("text", x = (N_don/2)+buffer, y = 48, label = paste("Donor", N_don, sep = ": "), size = 8) +
    geom_rect(aes(xmin= buffer, xmax= N_don+buffer, ymin= 30, ymax = 45), fill = 'blue3', color = 'black', size = 2) # Donor

  ## D_Frac_FMT
  p <- p + geom_rect(aes(xmin= buffer, xmax= FMT_don+buffer, ymin= 10, ymax = 25), fill = 'darkgreen', color = 'black', size = 2) +
    geom_rect(aes(xmin= FMT_don+buffer, xmax= N_don+buffer, ymin= 10, ymax = 25), fill = 'dimgrey', color = 'black', size = 2) +
    annotate("text", x = (N_don/2)+buffer, y = 7, label = paste("D[FracFMT]: ", D_frac_FMT, "(", FMT_don,")", sep = ""), parse = T, size = 8) # Annotate

  ## Recipient
  p <- p + annotate("text", x = (N_don+N_rec/2)+2*buffer, y = 48, label = paste("Recipient", N_rec, sep = ": "), parse = F, size = 8) +
    geom_rect(aes(xmin= N_don + 2*buffer, xmax= N_don + 2*buffer + N_rec, ymin= 30, ymax = 45), fill = 'firebrick3', color = 'black', size = 2) # Recipient

  ## R_Frac_FMT
  p <- p + geom_rect(aes(xmin= N_don+2*buffer, xmax= N_don+N_rec+2*buffer-FMT_rec, ymin= 10, ymax = 25), fill = 'dimgrey', color = 'black', size = 2) +
    geom_rect(aes(xmin= N_don+N_rec+2*buffer-FMT_rec, xmax= N_don+N_rec+2*buffer, ymin= 10, ymax = 25), fill = 'darkgreen', color = 'black', size = 2) +
    annotate("text", x = (N_don+N_rec/2)+2*buffer, y = 7, label = paste("R[FracFMT]: ", R_frac_FMT, " (", FMT_rec, ")", sep = ""), parse = T, size = 8) # Annotate

  ## FMT
  p <- p + annotate("text", x = N_don+N_rec+(FMT_don+FMT_rec+N_shared+N_unique)/2+3*buffer, y = 48, label = paste("Post-FMT", N_FMT, sep = ": "), parse = F, size = 8) +
    geom_rect(aes(xmin= N_don+N_rec+3*buffer, xmax=N_don+N_rec+FMT_don+FMT_rec+N_shared+N_unique+3*buffer, ymin= 30, ymax = 45), fill = 'darkgreen', color = 'black', size = 2) # FMT

  ## FMT_FracD & FMT_FracR & FMT_Shared, & FMT_Unique
  p <- p + geom_rect(aes(xmin= N_don+N_rec+3*buffer, xmax= N_don+N_rec+3*buffer+FMT_don, ymin= 10, ymax = 25), fill = 'blue3', color = 'black', size = 2) +
    annotate("text", x = (N_don+N_rec+FMT_don/2)+3*buffer, y = 7, label = paste("FMT[FracD]: ", FMT_FracD, "(", FMT_don, ")", sep = ""), parse = T, size = 6) +
    geom_rect(aes(xmin= N_don+N_rec+3*buffer+FMT_don, xmax= N_don+N_rec+3*buffer+FMT_don+FMT_rec, ymin= 10, ymax = 25), fill = 'firebrick3', color = 'black', size = 2) +
    annotate('segment', x = (N_don+N_rec+FMT_don+FMT_rec/2)+3*buffer, y = 3, xend = (N_don+N_rec+FMT_don+FMT_rec/2)+3*buffer, yend = 9, colour = "black", size = 2, arrow = arrow(length = unit(0.5, "cm"))) +
    annotate("text", x = (N_don+N_rec+FMT_don+FMT_rec/2)+3*buffer, y = 1, label =  paste("FMT[FracR]: ", FMT_FracR, "(", FMT_rec, ")", sep = ""), parse =  T, size = 6) +
    geom_rect(aes(xmin= N_don+N_rec+3*buffer+FMT_don+FMT_rec, xmax= N_don+N_rec+3*buffer+FMT_don+FMT_rec+N_unique, ymin= 10, ymax = 25), fill = 'darkorchid3', color = 'black', size = 2) +
    annotate("text", x = (N_don+N_rec+FMT_don+FMT_rec+N_unique/2)+3*buffer, y = 7, label = paste("FMT[Unique]: ", FMT_unique, "(",N_unique,")", sep = ""), parse = T, size = 6) +
    geom_rect(aes(xmin= N_don+N_rec+3*buffer+FMT_don+FMT_rec+N_unique, xmax= N_don+N_rec+3*buffer+FMT_don+FMT_rec+N_unique+N_shared, ymin= 10, ymax = 25), fill = 'chartreuse3', color = 'black', size = 2) +
    annotate('segment', x = (N_don+N_rec+FMT_don+FMT_rec+N_unique+N_shared/2)+3*buffer, y = 3, xend = (N_don+N_rec+FMT_don+FMT_rec+N_unique+N_shared/2)+3*buffer, yend = 9, colour = "black", size = 2, arrow = arrow(length = unit(0.5, "cm"))) +
    annotate("text", x = (N_don+N_rec+FMT_don+FMT_rec+N_unique+N_shared/2)+3*buffer, y = 1, label = paste("FMT[Shared]: ", FMT_shared, "(", N_shared, ")",  sep = ""), parse = T, size = 6)

  p <- p + ggtitle(paste("FMT Metric Visualization:", title_text, sep = "\n")) + theme(axis.line=element_blank(),
                                                                                       axis.text.x=element_blank(),
                                                                                       axis.text.y=element_blank(),
                                                                                       axis.ticks=element_blank(),
                                                                                       axis.title.x=element_blank(),
                                                                                       axis.title.y=element_blank(),
                                                                                       legend.position="none",
                                                                                       panel.background=element_blank(),
                                                                                       panel.border=element_blank(),
                                                                                       panel.grid.major=element_blank(),
                                                                                       panel.grid.minor=element_blank(),
                                                                                       plot.background=element_blank(),
                                                                                       plot.title = element_text(face="bold", size = 36)) +
    coord_cartesian(ylim = c(0,55))



  return(p)
}


