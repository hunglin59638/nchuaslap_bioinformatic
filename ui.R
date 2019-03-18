#if(!require(shiny)){install.packages("shiny")}
#if(!require(shiny)){install.packages("shinydashboard")}
options(encoding = "UTF-8")
require(shiny)
require(shinydashboard)

dashHeader <- dashboardHeader(title="Bioinformatic tool")

dashSidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem('Home',
             tabName='HomeTab',
             icon=icon('dashboard')
    ),
    menuItem("Quantification",
             tabName = "Quantification_Tab",
             icon=icon('bar-chart-o')
    ),
    menuItem("Descriptive statistics",
             tabName = "DescriptiveStatisticsTab",
             icon = icon("list-alt")),
    menuItem('Normality Tests',
             tabName = 'NormalityTestTab',
             icon =  icon("list-alt")
    ),
    menuItem('Heteroscedasticity test',
             tabName = "HeteroscedasticityTestTab",
             icon = icon("list-alt")
    ),
    menuItem("Mulitiple comparasion",
             tabName = "MulitipleComparasionTab",
             icon =  icon("list-alt")
    ),
    menuItem('Non-Parametric Tests',
             tabName = 'Non-ParametricTestsTab',
             icon =  icon("fas fa-calculator")
    ),
    menuItem("Correlation analysis",
             tabName = "CorrelationTab",
             icon = icon("fas fa-calculator")
    ),
    menuItem('Gene ID conversion',
             tabName = 'GeneIDconversionTab',
             icon = icon('cog')
    ),
    menuItem('GO classification',
             tabName = 'GOgroupTab',
             icon = icon("fas fa-chart-bar")
    ),
    menuItem('GO enrichment',
             tabName = "EnrichgoTab",
             icon = icon("fas fa-sitemap")
    ),
    menuItem('KEGG enrichment',
             tabName = "EnrichkeggTab",
             icon = icon("fas fa-allergies")),
    menuItem('Cross match',
             tabName = "CrossMatchTab",
             icon = icon("fas fa-allergies"))
  )
)

dashBody <- dashboardBody(
  fluidRow(
    tabItems(
      tabItem(tabName='HomeTab',
              h1("HUNGLIN's Bioinformatic"),
              p('If you have any further question, please contact hunglin59638@gmail.com.'),
              em('The introduction of the application will be updated later.'),
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
                        label = "Please choose CSV file.First column is group,Second is value",
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
                        label = "Please choose CSV file.First column is group,Second is value",
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
                        label = "Please choose CSV file.First column is group,Second is value",
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
                        label = "Please choose CSV file.First column is group,Second is value",
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
              fileInput(inputId= "nonparfile",
                        label = "Please choose CSV file.First column is group,Second is value",
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
              fileInput(inputId= "correlationfile",
                        label = "Please choose CSV file. columns are variables",
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
                tableOutput("ggo_table")
              )
      ),
      tabItem(tabName = "EnrichgoTab",
              h1('GO enrichment'),
              selectInput(inputId = "species_choice_ego",
                          label = "Choose species",
                          choices = c("Gallus gallus", "Sus scrofa"),
                          selected = "Gallus gallus"),
              selectInput(inputId = 'fromtype_ego',
                          label= 'Current gene ID type',
                          choices= c('Accession','Gene symbol', 'Entrez ID'),
                          selected= 'Accession'
              ),
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
                                 label = h3('Go enrichment results options'),
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
              h1("KEGG enrichment"),
              selectInput(inputId = "species_choice",
                          label = "Choose species",
                          choices = c("Gallus gallus", "Sus scrofa"),
                          selected = "Gallus gallus"),
              selectInput(inputId = 'fromtype_ekegg',
                          label= 'Current gene ID type',
                          choices= c('Accession','Gene symbol', 'Entrez ID'),
                          selected= 'Accession'
              ),
              checkboxGroupInput("regulated_gene_ekegg", 
                                 label = h3("up or down-regulation"), 
                                 choices = list("up-regulated genes" = "up_regulated",
                                                "down-regulated genes" = "down_regulated",
                                                "total genes" = "total"),
                                 selected = list("up_regulated" = "up_regulated",
                                                 "down_regulated" = "down_regulated",
                                                 "total" = "total")
              ),
              hr(),
              checkboxGroupInput("ekegg_options",
                                 label = h3('KEGG enrichment results options'),
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
              hr()
              
              
      ),
      tabItem("CrossMatchTab",
              h1("Cross match"),
              selectInput(inputId = "species_choice",
                          label = "Choose species",
                          choices = c("Gallus gallus", "Sus scrofa"),
                          selected = "Gallus gallus"),
              fileInput(inputId= "crossmatch1file",
                        label = paste("Please choose CSV file.",
                                      "Columns must contain 'Accession'",
                                      sep= "\n"),
                        multiple = F,
                        accept = c(
                          "text/csv", "text/comma-separated-values,text/plain",".csv"
                        )
              ),
              fileInput(inputId= "crossmatch2file",
                        label = "Please choose txt file.",
                        multiple = F,
                        accept = c(
                          "text/csv", "text/comma-separated-values,text/plain",".csv"
                        )
              ),
              downloadButton(label = "Download",
                             outputId = "crossmatch_result.csv")
              
      )
    )
  )
)


dashboardPage(
  header=dashHeader,
  sidebar=dashSidebar,
  body=dashBody,
  title="Bioimformatic tool"
)
