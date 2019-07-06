#if(!require(shiny)){install.packages("shiny")}
#if(!require(shiny)){install.packages("shinydashboard")}
#if(!require(shinythemes){install.packages("shinythemes")}
options(encoding = "UTF-8")
require(shiny)
require(shinydashboard)
library(shinythemes)
dashHeader <- dashboardHeader(title="Bioinformatic tool")

dashSidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem('Home',
             tabName='HomeTab',
             icon=icon("fas fa-church")
    ),
    menuItem("Quantification",
             tabName = "Quantification_Tab",
             icon=icon('bar-chart-o')
    ),
    menuItem("Descriptive statistics",
             tabName = "DescriptiveStatisticsTab",
             icon = icon("fas fa-file-alt")),
    menuItem('Normality Tests',
             tabName = 'NormalityTestTab',
             icon =  icon("fas fa-eye")
    ),
    menuItem('Heteroscedasticity test',
             tabName = "HeteroscedasticityTestTab",
             icon = icon("fas fa-venus-mars")
    ),
    menuItem("Mulitiple comparasion",
             tabName = "MulitipleComparasionTab",
             icon =  icon("list-alt")
    ),
    menuItem('Non-Parametric Tests',
             tabName = 'Non-ParametricTestsTab',
             icon =  icon("fas fa-align-right")
    ),
    menuItem("Correlation analysis",
             tabName = "CorrelationTab",
             icon = icon("fas fa-sitemap")
    ),
    menuItem("Cluster analysis",
             tabName = "ClusterTab",
             icon = icon("far fa-object-group")#"fas fa-bars"
    ),
    menuItem('Gene ID conversion',
             tabName = 'GeneIDconversionTab',
             icon = icon("fas fa-user-check")
    ),
    menuItem('GO classification',
             tabName = 'GOgroupTab',
             icon = icon("fas fa-sort-alpha-down")
    ),
    menuItem('GO enrichment',
             tabName = "EnrichgoTab",
             icon = icon("fas fa-project-diagram")
    ),
    menuItem('KEGG enrichment',
             tabName = "EnrichkeggTab",
             icon = icon("fab fa-connectdevelop")
    )#,
   # menuItem("DAVID functional analysis",
          #   tabName = "DAVIDTab",
          #   icon = icon("fas fa-dragon")
   # )
  )
)
#icon: https://fontawesome.com/icons?d=gallery&m=free
dashBody <- dashboardBody(
  fluidRow(
    tabItems(
      tabItem(tabName='HomeTab',
              h1("HUNGLIN's Bioinformatic tool"),
              p('If you have any further question, please contact hunglin59638@gmail.com'),
              p('The introduction of the application will be updated later.'),
              strong("If you use GO classification, GO enrichment and KEGG 
                     enrichment in published research, please cite:"),
              p(strong("Guangchuang Yu, Li-Gen Wang, Yanyan Han, Qing-Yu He. 
                       clusterProfiler: an R package for comparing biological themes 
                       among gene clusters. OMICS: A Journal of Integrative Biology. 
                       2012, 16(5):284-287.")),
              hr()
      ),
      tabItem(tabName = "Quantification_Tab",
              h1("Protein quantification"),
              textInput("stan_pos", 
                        label = h3("Where your standards start in the plate?"),
                        value = "A1"),
              numericInput(inputId= "stan_replication",
                           label = 'How many is the replication of the standards?',
                           3),
              textInput("stan_concentration", 
                        label = h3("What are the concentrations of standards?"),
                        value = "0,0.2,0.4,0.6,0.8,1.0,1.2"),
              hr(),
              textInput("sam_pos", 
                        label = h3("Where your samples start in the plate"),
                        value = "C1"),
              numericInput(inputId= "samplesize",
                           label = 'How many are your samples?',
                           6),
              numericInput(inputId= "sam_replication",
                           label = 'How many is the replication of the samples?',
                           3),
              numericInput(inputId = "sam_dilution_times",
                           label = "What is the dilution times",
                           10),
              numericInput("weight_need",
                           "The weight of protein you need.",
                           10),
              submitButton("Confirm"),
              hr(),
              fileInput(inputId= "quantificationfile",
                        label = "Please convert excel to csv",
                        multiple = F,
                        accept = c("text/csv",
                                   "text/comma-separated-values,text/plain",
                                   ".csv")
              ),
              mainPanel(
                verbatimTextOutput("quantitative_results"),
                hr("Experential concentration"),
                verbatimTextOutput("quantitative_prediction"),
                hr("Actual concentration"),
                verbatimTextOutput("dilution_results"),
                hr("The volume of the sample you need"),
                verbatimTextOutput("volume_need"),
                plotOutput("quantitative_plot")#,
                #verbatimTextOutput("testtable")
              )
      ),
      tabItem(tabName = "DescriptiveStatisticsTab",
              h1("Descriptive Statistics"),
              fileInput(inputId= "descriptivefile",
                        label = "Please choose CSV file.Column1 is group,column2 is value",
                        multiple = F,
                        accept = c("text/csv",
                                   "text/comma-separated-values,text/plain",
                                   ".csv")
              ),
              submitButton("Confirm"),
              hr(),
              mainPanel(
                verbatimTextOutput("descriptive_summary"),
                plotOutput("descriptive_barplot"),
                hr(),
                hr(),
                plotOutput("descriptive_boxplot")
              )
              
      ),
      tabItem(tabName='NormalityTestTab',
              h1('Normality Tests'),
              selectInput(inputId='Normalitytestmethod',
                          label='Choose stastics methods',
                          choices=c('Shapiro-Wilk test (3<n<50)', 'Kolmogorov-Smirnov test'),
                          selected='Shapiro-Wilk test (3<n<50)'),
              submitButton("Confirm"),
              hr(),
              fileInput(inputId= "normalityfile",
                        label = "Please choose CSV file. Column1 is group,Column2 is value",
                        multiple = F,
                        accept = c("text/csv",
                                   "text/comma-separated-values,text/plain",
                                   ".csv")
              ),
              hr(),
              mainPanel(
                # tableOutput("normalitytest_table"),
                verbatimTextOutput("normalitytest_summary"),
                plotOutput("graphical_test")
              )
      ),
      tabItem(tabName = "HeteroscedasticityTestTab",
              h1("Heteroscedasticity test"),
              fileInput(inputId= "heterofile",
                        label = "Please choose CSV file.Column1 is group,Column2 is value",
                        multiple = F,
                        accept = c("text/csv",
                                   "text/comma-separated-values,text/plain",
                                   ".csv")
              ),
              hr(),
              submitButton("Confirm"),
              hr(),
              mainPanel(
                verbatimTextOutput("hetero_result")
              )
      ),
      tabItem(tabName = "MulitipleComparasionTab",
              h1("Mulitiple Comparasion"),
              selectInput(inputId='mcomparemethod',
                          label='Choose stastics methods',
                          choices=c('Least Squares Means'),
                          selected='Least Squares Means'),
              fileInput(inputId= "mcomparefile",
                        label = "Please choose CSV file.Column1 is group,Column2 is value",
                        multiple = F,
                        accept = c("text/csv",
                                   "text/comma-separated-values,text/plain",
                                   ".csv")
              ),
              hr(),
              submitButton("Confirm"),
              hr(),
              mainPanel(
                verbatimTextOutput("mcompare_summary"),
                plotOutput("mcompare_plot"),
                hr(),hr(),
                verbatimTextOutput("mcompare_hoctest")
              )
      ),
      tabItem(tabName='Non-ParametricTestsTab',
              h1('Non-Parametric Tests'),
              selectInput(inputId='nonparmethod',
                          label='Choose stastics methods',
                          choices=c('Kruskal-Wallis test for equal variances', 
                                    "Welch's anova for unequal variances"),
                          selected='Kruskal-Wallis test for equal variances'),
              hr(),
              submitButton("Confirm"),
              hr(),
              strong("Please choose CSV file."),
              fileInput(inputId= "nonparfile",
                        label = "Column1 is group,column2 is value.",
                        multiple = F,
                        accept = c("text/csv",
                                   "text/comma-separated-values,text/plain",
                                   ".csv")
              ),
              hr(),
              mainPanel(
                verbatimTextOutput("NonParametricTests_summary"),
                verbatimTextOutput("NonParametricTests_hoctest")
              )
      ),
      tabItem(tabName = "CorrelationTab",
              h1("Correlation analysis"),
              selectInput(inputId = "method_cor",
                          label = "Choose method of correlaiton analysis",
                          choices = c("Pearson's product-moment correltion",
                                      "Spearman's rank order correlation"),
                          selected = "Pearson's product-moment correltion"),
              hr(),
              submitButton("Confirm"),
              hr(),
              strong("Please choose CSV file."),
              fileInput(inputId= "correlationfile",
                        label = "columns are variables",
                        multiple = F,
                        accept = c("text/csv",
                                   "text/comma-separated-values,text/plain",
                                   ".csv")
              ),
              mainPanel(
                verbatimTextOutput("cor_test"),
                tableOutput("cor_table"),
                plotOutput("cor_plot")
              ),
              hr(),
              downloadButton(outputId = "Correlation.csv", label = "Download result"),
              hr(),
              downloadButton(outputId = "Correlation_plot.pdf", label = "Download plot")
              
      ),
      tabItem(tabName = "ClusterTab",
              h1("Cluster analysis"),
              selectInput(inputId = "cluster_method",
                          label = "Choose a method of clustering",
                          choices = c("one","two"),
                          selected = "one"),
              hr(),
              strong("Please choose CSV file."),
              fileInput(inputId= "clusterfile",
                        label = "columns are variables",
                        multiple = F,
                        accept = c("text/csv",
                                   "text/comma-separated-values,text/plain",
                                   ".csv")
              )
              ),
      tabItem(tabName = 'GeneIDconversionTab',
              h1('Gene ID conversion'),
              selectInput(inputId = "species_choice",
                          label = "Choose species",
                          choices = c("Gallus gallus", "Sus scrofa"),
                          selected = "Gallus gallus"),
              selectInput(inputId = 'typefrom',
                          label = 'Current gene ID type',
                          choices = c('Accession','Gene symbol', 'Entrez ID', 'UniProt','Ensembl','GO','Gene name'),
                          selected = 'Accession'),
              selectInput(inputId = 'typeto',
                          label = 'Convert gene ID type to',
                          choices = c('Accession','Gene symbol', 'Entrez ID', 'UniProt','Ensembl','GO','Gene name'),
                          selected = 'Entrez ID'),
              hr(), 
              submitButton("Confirm"),
              hr(), 
              fileInput(inputId= "convertidfile",
                        label = "Please choose CSV file.First column is gene ID",
                        multiple = F,
                        accept = c("text/csv",
                                   "text/comma-separated-values,text/plain",
                                   ".csv")
              ),
              
              mainPanel(
                #textOutput("idtest"), 
                tableOutput("geneconvert_table")
              ),
              hr(),
              downloadButton(outputId = "Convert_gene_id.csv", label = "Download")
      ),
      tabItem(tabName = 'GOgroupTab',
              h1('GO classification'),
              selectInput(inputId = "species_choice",
                          label = "Choose species",
                          choices = c("Gallus gallus", "Sus scrofa"),
                          selected = "Gallus gallus"),
              selectInput(inputId = 'fromtype_ggo',
                          label= 'Current gene ID type',
                          choices= c('Accession','Gene symbol', 'Entrez ID'),
                          selected= 'Accession'
              ),
              checkboxGroupInput("regulated_gene_ggo", 
                                 label = h3("up or down-regulation"), 
                                 choices = list("up-regulated genes" = "up_regulated",
                                                "down-regulated genes" = "down_regulated",
                                                "total genes" = "total"),
                                 selected = list("up_regulated" = "up_regulated",
                                                 "down_regulated" = "down_regulated",
                                                 "total" = "total")
              ),
              checkboxGroupInput("GO_ontology",
                                 label = h3("GO Ontology"),
                                 choices = list("Biological Process" = "BP",
                                                "Cellular Component" = "CC",
                                                "Molecular Function" = "MF"),
                                 selected = list("Biological Process" = "BP",
                                                 "Cellular Component" = "CC",
                                                 "Molecular Function" = "MF")
                                 ),
              numericInput("level_ggo", 
                           label = h4("Specific GO level."),
                           value = "2"),
              submitButton("Confirm"),
              hr(),
              fileInput(inputId= "file_ggo",
                        label = "Please choose CSV file.First column is gene ID; sencond is ratio",
                        multiple = F,
                        accept = c(
                          "text/csv", "text/comma-separated-values,text/plain",".csv"
                        )
              ),
              hr(),
              downloadButton(outputId = "ggo_results.zip",
                             label = "Download Results"),
              hr(),
              downloadButton(outputId = "ggo_plots.zip",
                             label = "Download Plots"),
              hr(),
              mainPanel(
                #textOutput("idtest"), 
                textOutput("ggo_test")
              )
      ),
      tabItem(tabName = "EnrichgoTab",
              h1('GO over-representation test'),
              selectInput(inputId = "species_choice_ego",
                          label = h4("Choose species"),
                          choices = c("Gallus gallus", "Sus scrofa"),
                          selected = "Gallus gallus"),
              selectInput(inputId = 'fromtype_ego',
                          label= h4('Current gene ID type'),
                          choices= c('Accession','Gene symbol', 'Entrez ID'),
                          selected= 'Accession'
              ),
              selectInput("padjust_ego",
                          label = h4("Adjust p-values for multiple comparisons"),
                          choices = c("Bonferroni","Holm","Hochberg", "Hommel",
                                      "Benjamini & Hochberg","Benjamini & Yekutieli",
                                      "None"),
                          selected = "Benjamini & Hochberg"),
              numericInput("pvalueCutoff_ego", 
                           label = h4("Cutoff value of pvalue."),
                           value = "0.05"),
              numericInput("qvalueCutoff_ego",
                           label = h4("Cutoff value of qvalue."),
                           value = "0.2"),
              hr(),
              checkboxGroupInput("regulated_gene_ego", 
                                 label = h3("up or down-regulation"), 
                                 choices = list("up-regulated genes" = "up_regulated",
                                                "down-regulated genes" = "down_regulated",
                                                "total genes" = "total"),
                                 selected = list("up_regulated" = "up_regulated",
                                                 "down_regulated" = "down_regulated",
                                                 "total" = "total")
              ),
              hr(),
              checkboxGroupInput("ego_options",
                                 label = h3('Go enrichment results'),
                                 choices = list("GO dotplot" = "dotplot_ego",
                                                "GO DAG graph" = "dag_ego",
                                                "Gene-Concept Network for GO" = "Gene-Concept-Network_ego",
                                                "Enrichment Map for GO" = "emapplot_ego"
                                                #"KEGG enrichment" = "kegg_enrich",
                                                #"Gene-Concept Network for KEGG" = "cnetplot_kegg",
                                                #"Enrichment Map for KEGG" = "emapplot_kegg"
                                 ),
                                 selected = list("GO dotplot" = "dotplot_ego",
                                                 "GO DAG graph" = "dag_ego",
                                                 "Gene-Concept Network for GO" = "Gene-Concept-Network_ego",
                                                 "Enrichment Map for GO" = "emapplot_ego"
                                 )
              ),
              submitButton("Confirm"),
              hr(),
              fileInput(inputId= "file_ego",
                        label = "Please choose CSV file.First column is geneid; sencond is ratio",
                        multiple = F,
                        accept = c(
                          "text/csv", "text/comma-separated-values,text/plain",".csv"
                        )
              ),
              mainPanel(
                textOutput("test_ego"),
                #textOutput("test"),
                tableOutput("ego_table")
              ),
              hr(),
              hr(),
              downloadButton(outputId = "ego_results.zip",
                             label = "Download Results"),
              hr(),
              downloadButton(outputId = "ego_plots.zip",
                             label = "Download Plots")
      ),
      tabItem(tabName = "EnrichkeggTab",
              h1("KEGG over-representation test"),
              selectInput(inputId = "species_choice_ekegg",
                          label = h4("Choose species"),
                          choices = c("Gallus gallus", "Sus scrofa"),
                          selected = "Gallus gallus"),
              selectInput(inputId = 'fromtype_ekegg',
                          label= h4('Current gene ID type'),
                          choices= c('Accession','Gene symbol', 'Entrez ID'),
                          selected= 'Accession'
              ),
              selectInput("padjust_ekegg",
                          label = h4("Adjust p-values for multiple comparisons"),
                          choices = c("Bonferroni","Holm","Hochberg", "Hommel",
                                      "Benjamini & Hochberg","Benjamini & Yekutieli",
                                      "None"),
                          selected = "Benjamini & Hochberg"),
              numericInput("pvalueCutoff_ekegg", 
                           label = h4("Cutoff value of pvalue."),
                           value = "0.05"),
              numericInput("qvalueCutoff_ekegg",
                           label = h4("Cutoff value of qvalue."),
                           value = "0.2"),
              hr(),
              checkboxGroupInput("regulated_gene_ekegg", 
                                 label = h4("up or down-regulation"), 
                                 choices = list("up-regulated genes" = "up_regulated",
                                                "down-regulated genes" = "down_regulated",
                                                "total genes" = "total"),
                                 selected = list("up_regulated" = "up_regulated",
                                                 "down_regulated" = "down_regulated",
                                                 "total" = "total")
              ),
              hr(),
              checkboxGroupInput("ekegg_options",
                                 label = h4('KEGG enrichment results'),
                                 choices = list("KEGG dotplot" = "dotplot_ekegg",
                                                "Gene-Concept Network for KEGG" = "Gene-Concept-Network_ekegg",
                                                "Enrichment Map for KEGG" = "emapplot_ekegg"
                                 ),
                                 selected = list("KEGG dotplot" = "dotplot_ekegg",
                                                 "Gene-Concept Network for KEGG" = "Gene-Concept-Network_ekegg",
                                                 "Enrichment Map for KEGG" = "emapplot_ekegg"
                                 )
              ),
              submitButton("Confirm"),
              hr(),
              fileInput(inputId= "file_ekegg",
                        label = "Please choose CSV file.First column is geneid; sencond is ratio",
                        multiple = F,
                        accept = c(
                          "text/csv", "text/comma-separated-values,text/plain",".csv"
                        )
              ),
              mainPanel(
                textOutput("test_ekegg")#,
                #textOutput("test")
              ),
              hr(),
              hr(),
              downloadButton(outputId = "ekegg_results.zip",
                             label = "Download results"),
              hr(),
              downloadButton(outputId = "ekegg_plots.zip",
                             label = "Download plots"),
              hr(),
              strong("GO to Pathview Web to render pathway graphs"),
              actionButton("pathview",label = "Pathview",
                           onclick ="window.open('https://pathview.uncc.edu/', '_blank')")
              
              
      )
    )
  )
)


dashboardPage(
  header=dashHeader,
  sidebar=dashSidebar,
  body=dashBody,
  skin = "black",
  title="Bioimformatic tool"
)
