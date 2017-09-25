require(shiny)
require(DT)
options(shiny.maxRequestSize = 200 * 1024 ^ 2)
require(shinyBS)
require(shinyjs)

shinyUI(
  navbarPage(
    title = "FMT Efficacy Analysis Toolkit v3.0",
    windowTitle = "FMT Efficacy Analysis Toolkit (FEAT)",
    position = 'static-top',

    tabPanel(id = 'test', title = "Home",
      h2("Home"),
      tags$style(
        type = "text/css",
        ".shiny-output-error { visibility: hidden;
        }",
               ".shiny-output-error:before { visibility: hidden; }"
      ),
      withMathJax(),
      useShinyjs(),
      h3("Description"),
      p(
        "This tool is designed to help quickly and interactively analyze microbiome sequencing data from FMT experiments. The three main inputs that you need are an OTU table, a mapping file, and a taxonomy map."
      ),
      h3("Inputs"),
      p(
        strong("OTU Table:"),
        "This tool will work best with manageably-sized OTU tables in the hdf5 biom format. The preferred method of getting a table in this format is to use the sequence data output from split_libraries.py (from QIIME), and dereplicate sequences to generate OTUs. Then, filter out singleton OTUs to drastically reduce the size of the OTU table. For an example on how to do this, and retain sequence information, see ",
        a("here", href = "https://gist.github.com/gregcaporaso/f3c042e5eb806349fa18", target =
            "_blank")
      ),
      p(
        strong("Metadata File:"),
        "A mapping file that contains one metadata category in which the donor, pre-transplant recipient, and post-transplant recipient can be identified. It is important that the SampleID field (from your standard QIIME mapping file) remains unchanged, as it is the field used to add metadata to the OTU table. Also, it is important that the file contain biomass data in a column titled 'biomass_ratio'. If you do not have biomass data, samples will be given a biomass of 1."
      ),
      p(
        strong('Taxonomy Map:'),
        "A .txt file that contains the taxonomy associated with each OTU ID. This can be from the greengenes reference, for example, or from the output of the assign_taxonomy.py function from QIIME."
      ),
      h3("Overview"),
      p(
        "This tool will split your OTU table into your experiment-specific samples. Next, a filtering step occurs after the OTU counts are normalized, that will remove OTUs that are less abundant than the given threshold. Subsequently, the tool will identify OTUs that are unique and specific to the donor and to the recipient samples, and perform analyses on the post-transplant samples based only on these information-rich OTUs. This includes calculating a distance matrix of the samples based on these OTUs, and allowing for the interactive visualization of this information plus any additional metadata."
      ),
      p(
        "Moreover, we have developed metrics to help quantify the ability of a microbiota to engraft in a recipient and for a microbiota to be pushed out by another (see figure)."
      ),
      img(src = 'FEAT-Metrics.png'),
      p(
        "Click through the sidebar to add inputs and change the parameters for analyzing the data."
      )
      ),
    tabPanel(title = "Inputs/Settings",
      h2('Select Inputs/Settings'),
      fluidRow(
        column(3, offset = 1,
               wellPanel(
                 fileInput(
                   "input_otu_table",
                   label = h5("Input OTU Table"),
                   accept = 'Document'
                 ),
                 bsPopover(
                   'input_otu_table',
                   "OTU Table",
                   "Absolute-filtered OTU table in hdf5 BIOM (.biom) format."
                 )
               )),
        column(3, wellPanel(
          fileInput(
            "metadata_file",
            label = h5("Metadata File"),
            accept = c('txt', 'text/plain')
          ),
          bsPopover(
            "metadata_file",
            "Metadata File",
            "This should be the same mapping file used to create your OTU table, and should contain at least one metadata category to select the details of the transplant for your experiment of interest."
          )
        )),
        column(3, wellPanel(
          fileInput(
            "id_tax_file",
            label = h5("Input Taxonomy Map"),
            accept = c('txt', 'text/plain')
          ),
          bsPopover(
            'id_tax_file',
            "OTU ID-to-Taxonomy Map",
            "A .txt file that contains the taxonomy associated with each OTU ID. This can be from the greengenes reference, for example, or from the output of the assign_taxonomy.py function in QIIME."
          )
        ))

      ),
      hr(),
      fluidRow(
        column(
          3,
          conditionalPanel(
            condition = "output.files_loaded",
            wellPanel(
              h3('Select FMT Details'),
              uiOutput("compare_options"),
              uiOutput("donor"),
              uiOutput("recipient"),
              uiOutput("post_fmt")
            )
          )
        ),
        column(
          5,
          conditionalPanel(condition = "output.files_loaded",
                           wellPanel(
                             h3('Other Settings'),
                             fluidRow(
                               column(
                                 6,
                                 numericInput(
                                   "min_OTU_fraction",
                                   label = h5("Minimum Relative Abundance Filter"),
                                   min = 0.0001,
                                   max = 0.1,
                                   value = 0.001
                                 ),
                                 bsPopover(
                                   'min_OTU_fraction',
                                   "Relative Abundance Filter",
                                   "The minimum fraction that an OTU must represent in a sample (in the sample it is most abundant), in order to be kept. (i.e. an OTU must represent at least this fraction of the microbiota of any one sample in order to be kept)"
                                 ),
                                 sliderInput(
                                   "min_fraction",
                                   label = h5("Minimum fraction of samples with a non-zero abundance"),
                                   min = 0,
                                   max = 1.0,
                                   value = 2 / 3
                                 ),
                                 bsPopover(
                                   "min_fraction",
                                   'Minimum Fraction',
                                   "The minimum fraction of samples within a condition that must have a non-zero abundance of an OTU in order for it to be retained as reliably present in that condition."
                                 )
                               ),
                               column(
                                 6,
                                 selectInput(
                                   'abundance_type',
                                   label = "Select Abundance Type",
                                   choices = c(
                                     "Absolute Abundance" = 'absolute',
                                     "Relative Abundance" = 'relative'
                                   ),
                                   selected = 'relative'
                                 ),
                                 bsPopover(
                                   'abundance_type',
                                   'Select Abundance Type',
                                   'Indicate whether to work with the relative or absolute abundances of taxa in your sample. If using absolute abundances, your metadata file MUST contain a column named "biomass_ratio", that contains the biomass of the sample (ug DNA per mg sample)',
                                   placement = 'top'
                                 ),
                                 uiOutput('check_taxonomy'),
                                 bsPopover(
                                   'toggle_taxonomy',
                                   'Add Taxonomy?',
                                   'Indicate whether or not to incorporate taxonomic information to the analysis. Essentially, do you want to work with OTUs or Taxonomic Groups)',
                                   placement = 'top'
                                 ),
                                 uiOutput('depth'),
                                 bsPopover(
                                   'depth',
                                   'Select Taxonomic Depth',
                                   'Select the depth with which to describe and collapse taxa.',
                                   placement = 'top'
                                 )
                               )
                             )
                           ))
        ),
        column(
          3,
          conditionalPanel(
            condition = "output.files_loaded",
            wellPanel(
              h3('Begin Analysis'),
              helpText('This may take several minutes.'),
              actionButton(
                "go",
                "Go!",
                icon = icon('ok', lib = 'glyphicon'),
                width = '100%',
                style = "color: #fff; background-color: #046E00; border-color: #0A3D08"
              ),
              bsPopover(
                'go',
                "Go!",
                "Click to select or update the Donor, Recipient, and Post-FMT metadata categories from the selected column of the mapping file, or the minimum relative abundance filter.",
                placement = 'top'
              )
            )
          ),
          conditionalPanel(
            condition = "output.data_analyzed",
            wellPanel(
              h3("Quick Summary"),
              span(strong(textOutput(
                'num_donor_samples'
              )), style = "color:blue"),
              span(strong(textOutput(
                'num_recipient_samples'
              )), style = "color:red"),
              span(strong(textOutput(
                'num_post_fmt_samples'
              )), style = "color:green"),
              hr(),
              span(strong(textOutput('raw_otu_count'))),
              span(strong(
                textOutput('num_otus_after_relative_filter')
              )),
              span(strong(
                textOutput('num_otus_after_fleeting_filter')
              ))
            )
          )
        )
      )
    )
    ,
    tabPanel(title = 'Transplant Metrics',
      fluidRow(
        column(3, h4("Donor"),
               textOutput("donor_id")),
        column(3, h4("Recipient"),
               textOutput("recipient_id")),
        column(3, h4("Post-FMT"),
               textOutput("post_fmt_id")),
        column(3, h4("Other"),
               textOutput('other_details'))
      ),
      hr(),
      fluidRow(column(
        5,
        wellPanel(
          style = 'overflow:hidden;background-color:white',
          h4("Numbers of OTUs/Taxa in each sample set"),
          br(),
          uiOutput("N_Donor"),
          bsPopover(
            'N_Donor',
            title = HTML(paste("N", tags$sub("Donor"), sep = "")),
            "The number of OTUs/Taxa that are reliably and exclusively in the donor.",
            placement = 'top'
          ),
          br(),
          uiOutput("N_Recipient"),
          bsPopover(
            'N_Recipient',
            title = HTML(paste("N", tags$sub("Recipient"), sep = "")),
            "The number of OTUs/Taxa that are reliably and exclusively in the recipient.",
            placement = 'top'
          ),
          br(),
          uiOutput("N_P_Total"),
          bsPopover(
            'N_P_Total',
            title = HTML(paste("N", tags$sub("Post-FMT"), sep = "")),
            "The number of OTUs/Taxa that are reliably in the post-transplant samples (minus number that are unique to the post-transplant samples and the number that are shared throughout, if this is selected for).",
            placement = 'top'
          ),
          br(),
          uiOutput("N_P_Unique"),
          bsPopover(
            'N_P_Unique',
            title = HTML(paste("N", tags$sub("Post-FMT (Unique)"), sep = "")),
            "Number of OTUs/Taxa that are reliably and exclusively in the post-transplant samples. This number should theoretically be low or zero, but can serve as a useful pseudo-negative control. These OTUs/Taxa were in NEITHER the Donor NOR the Recipient microbiotas.",
            placement = 'top'
          ),
          br(),
          uiOutput("N_P_Shared"),
          br(),
          uiOutput("N_shared_lost"),
          bsPopover(
            "N_shared_lost",
            title = HTML(paste("N", tags$sub("Lost | Shared"), sep = "")),
            content = "The number of OTUs/Taxa that were initially shared between the donor and recipient that were subsequently lost.",
            placement = 'top'
          ),
          bsPopover(
            'N_P_Shared',
            title = HTML(paste("N", tags$sub("Post-FMT (Shared)"), sep = "")),
            "Number of OTUs/Taxa that are reliably present in both donor AND recipient samples throughout.",
            placement = 'top'
          ),
          br(),
          uiOutput("N_P_Donor"),
          bsPopover(
            "N_P_Donor",
            title = HTML(paste("N", tags$sub("Donor | Post-FMT"), sep = "")),
            content = "The number of OTUs/Taxa in the post-transplant samples that came from the donor.",
            placement = 'top'
          ),
          br(),
          uiOutput("N_P_Recipient"),
          bsPopover(
            "N_P_Recipient",
            title = HTML(paste(
              "N", tags$sub("Recipient | Post-FMT"), sep = ""
            )),
            content = "The number of OTUs/Taxa in the post-transplant samples that came from the recipient.",
            placement = 'top'
          )
        )
      ),
      column(
        5,
        offset = 1,
        wellPanel(
          style = 'overflow:hidden;background-color:white',
          h4("FMT Metrics"),
          br(),
          uiOutput("D_Engraft"),
          bsPopover(
            'D_Engraft',
            title = HTML(paste("D", tags$sub("Engraft"), sep = "")),
            content = HTML(
              paste(
                "The proportion of the donor microbiota that engrafts into the recipient.",
                tags$img(src = 'D-engraft.png'),
                sep = ""
              )
            ),
            placement = 'top'
          ),
          br(),
          uiOutput("R_Persist"),
          bsPopover(
            'R_Persist',
            title = HTML(paste("R", tags$sub("Persist"), sep = "")),
            content = HTML(
              paste(
                "The proportion of the recipient microbiota that persists after the transplant.",
                tags$img(src = 'R-persist.png'),
                sep = ""
              )
            ),
            placement = 'top'
          ),
          br(),
          uiOutput("P_Donor"),
          bsPopover(
            "P_Donor",
            title = HTML(paste("P", tags$sub("Donor"), sep = "")),
            content = HTML(
              paste(
                "The proportion of OTUs/Taxa in the post-transplant samples that is derived from the donor.",
                tags$img(src = 'P-donor.png'),
                sep = ""
              )
            ),
            placement = 'top'
          ),
          br(),
          uiOutput("P_Recipient"),
          bsPopover(
            "P_Recipient",
            title = HTML(paste("P", tags$sub("Recipient"), sep = "")),
            content = HTML(
              paste(
                "The proportion of OTUs/Taxa in the post-transplant samples that is derived from the recipient.",
                tags$img(src = 'P-recipient.png'),
                sep = ""
              )
            ),
            placement = 'top'
          ),
          br(),
          uiOutput("P_Unique"),
          bsPopover(
            "P_Unique",
            title = HTML(paste("P", tags$sub("Unique"), sep = "")),
            content = HTML(
              paste(
                "The proportion of OTUs/Taxa in the post-transplant samples that are neither in the donor NOR in the recipient.",
                tags$img(src = 'P-unique.png'),
                sep = ""
              )
            ),
            placement = 'top'
          ),
          br(),
          uiOutput("P_Shared"),
          bsPopover(
            "P_Shared",
            title = HTML(paste("P", tags$sub("Shared"), sep = "")),
            content = HTML(
              paste(
                "The proportion of OTUs/Taxa in the post-transplant samples that are shared in both the donor AND the recipient.",
                tags$img(src = 'P-shared.png'),
                sep = ""
              )
            ),
            placement = 'top'
          ),
          hr(),
          uiOutput("Engraftment"),
          bsPopover(
            "Engraftment",
            title = "Engraftment",
            content = HTML(
              paste(
                "Engraftment metric to capture proportion of Post-FMT Recipient that is derived from the Donor.",
                tags$img(src = 'Engraftment.png'),
                sep = ""
              )
            ),
            placement = 'top'
          )
        )
      )),
      hr(),
      h5("Transplant Vizualization"),
      div(
        style = "position:relative; margin-bottom:50px",
        plotOutput("metric_visualization",
                   click = 'plot_click')
      ),
      fluidRow(column(10, offset = 1,
                      plotOutput('click_plot')
      )),
      hr(),
      dataTableOutput("click_info")

    ),
    tabPanel(title = 'Data Visualization',
      hr(),
      fluidRow(column(12,
                      wellPanel(
                        fluidRow(
                          column(
                            3,
                            selectInput(
                              'distance_method',
                              label = 'Select Distance Method',
                              choices = c(
                                'Euclidean' = "euclidean",
                                'Jaccard' = "jaccard",
                                'Bray-Curtis' = "bray",
                                'Chao' = "chao",
                                'Canberra' = "canberra",
                                'Manhattan' = "manhattan",
                                'Gower' = "gower",
                                'AltGower' = "altGower",
                                'Kulczynski' = "kulczynski",
                                'Morisita' = "morisita",
                                'Horn' = "horn",
                                'Binomial' = "binomial",
                                'Cao' = "cao",
                                'Mountford' = "mountford",
                                'Roup' = "roup"
                              ),
                              selected = "euclidean"
                            ),
                            uiOutput('plot_color_by')
                          ),
                          column(3,
                                 uiOutput('plot_x'),
                                 uiOutput('plot_y')),
                          column(3,
                                 uiOutput('plot_facet_row'),
                                 uiOutput('plot_facet_col')),
                          column(
                            3,
                            sliderInput(
                              "point_size",
                              "Point size:",
                              min = 2,
                              max = 10,
                              value = 6,
                              step = 1
                            ),
                            checkboxGroupInput(
                              inputId = 'plot_type',
                              label = 'Plot Type',
                              choices = c('Scatter', 'Boxplot'),
                              selected = 'Scatter',
                              inline = T
                            )
                          )
                        )
                      ))),
      hr(),
      fluidRow(column(6,
                      plotOutput("data_plot")),
               column(6,
                      plotOutput('spree_plot'))),

      hr(),
      fluidRow(
        column(2, offset = 1,
               h4("PC1 Percent Explained"),
               textOutput("pc1")),
        column(2,
               h4("PC2 Percent Explained"),
               textOutput("pc2")),
        column(2,
               h4("PC3 Percent Explained"),
               textOutput("pc3")),
        column(2,
               h4("PC4 Percent Explained"),
               textOutput("pc4")),
        column(2,
               h4("PC5 Percent Explained"),
               textOutput("pc5"))
      ),
      hr(),
       dataTableOutput('donor_unique_table')

      # verbatimTextOutput('test2')
    ),
    tabPanel(title = "Downloads",
             fluidRow(column(8,
                             wellPanel(
                               fluidRow(
                                 column(
                                   12,
                                   h3("Download Data"),
                                   helpText("The 'Download Data' button will download:", tags$ul(tags$li("A summary of your input parameters and files."), tags$li("A summary of the computed FMT metrics"), tags$li('Tables:'), tags$ul(tags$li("A table listing the taxa and their abundances for taxa unique to the Donor samples"), tags$li("A table listing the taxa and their abundances for taxa unique to the Recipient samples"), tags$li("A table listing the taxa and their abundances for taxa in the Post-FMT samples"), tags$li("A table listing the taxa and their abundances for taxa unique to the Donor samples that engrafted into the recipient"), tags$li("A table listing the taxa and their abundances for taxa unique to the Recipient samples that persisted following the transplant"), tags$li("A table listing the taxa and their abundances for taxa unique to the Post-FMT samples"),  tags$li("A table listing the taxa and their abundances for taxa that were shared between the Donor and Recipient samples throughout")), tags$li("A table listing the percent explained of the first five principal components"), tags$li("A data table with all of the metadata, computed diversity, PCs of distances between samples, and abundances for each sample")), "The 'Download Figures' button will download:", tags$ul(tags$li("A visual representaiton of the Transplant"), tags$li("A PCA plot of the first two principal components"), tags$li("The plot currently displayed in the Data Visualization tab"))),
                                   downloadButton('downloadData', 'Download Data'),
								   downloadButton('downloadFigures', 'Download Figures')
                                 )))))
             )
  )
)
