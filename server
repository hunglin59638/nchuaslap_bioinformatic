#if(!require(shiny)){install.packages("shiny")}
#if(!require(BiocManager)){install.packages("BiocManager")}
#if(!require(AnnotationHub)){BiocManager::install(c("AnnotationHub","BiocGenerics","parallel"))
#  }
#if(!require(clusterProfiler)){BiocManager::install(c("clusterProfiler","fgsea","enrichGO","DOSE"), version = "3.8")}
#if(!require(ggplot2)){install.packages("ggplot2")}
#if(!require(stringr)){install.packages("stringr")}
#if(!require(dplyr)){install.packages("dplyr")}
#if(!require(magrittr)){install.packages("magrittr")}
#if(!require(dplyr)){install.packages("dplyr")}
#if(!require(plyr)){install.packages("plyr")}
#if(!require(purrr)){install.packages("purrr")}
#if(!require(userfriendlyscience)){
#install.packages("userfriendlyscience")
#install.packages("httpuv")
#install.packages("rlang")
#install.packages("tibble")
#}
#if(!require(car)){install.packages("car")}

options(encoding = "UTF-8")
require(shiny)
require(AnnotationHub)
require(clusterProfiler)
require(AnnotationDbi)
require(ggplot2)
require(stringr)
require(magrittr)
require(dplyr)
require(plyr)
require(purrr)
library(userfriendlyscience)
library(FSA)
library(lattice)
library(car)
library(emmeans)

#donwload gallus gallus database
#ah <- AnnotationHub()
#databaseid<- query(ah, c("Gallus gallus","OrgDb")) 
#gallus<- ah[[databaseid$ah_id]]

shinyServer( function(input, output, session) {
  
  #protein quantification
  quan_data <- reactive({
    infile <- input$quantificationfile
    if (is.null(infile)) {
      return(NULL)
    }
    data<- read.csv(infile$datapath)
    rownames(data) <- data[ ,1]
    data <- data[ , -1]
    names(data) <- c(1:12)
    data
  })
  stan_concentration <- reactive({
    stan_string<- str_split(input$stan_concentration, pattern = ",")
    stan_string[[1]]
  })
  stan_row_col <- reactive({
    str_split(input$stan_pos, pattern = "")
  })
  stan_replication <- reactive(input$stan_replication)
  
  word2num_fun <- function(word) {
    words <- c("A","B","C","D","E","F","G","H")
    num<- seq(1,8)
    word2num <- t(data.frame(num)) %>% data.frame()
    names(word2num) <- words
    as.numeric(word2num[word])
  } #convert "A" to 1 ...
  
  #standard
  stan_lm <- reactive({
    data <- quan_data()
    stan_row_start <- stan_row_col()[[1]][1]
    stan_col_start <- stan_row_col()[[1]][2] %>% as.numeric()
    standardsize <- stan_concentration() %>% length()
    stan_replication <- stan_replication()
    stan_df <- data.frame(rep(NA, stan_replication))
    stan_times <- 0
    words <- c("A","B","C","D","E","F","G","H")
    
    while(stan_times < standardsize*stan_replication) {
      stan_df <- data.frame(stan_df, 
                            t(data[word2num_fun(stan_row_start),
                                   (stan_col_start:(stan_col_start+stan_replication-1))]))
      stan_col_start <- stan_col_start + stan_replication
      while ((stan_col_start -1) >= 12) {
        stan_col_start <- 1
        stan_row_start <- words[word2num_fun(stan_row_start)+1]
      }
      stan_times <- stan_times + stan_replication
    }
    
    stan_df<- stan_df[ ,-1]
    stan_mean <- sapply(stan_df, function(x){mean(x,na.rm=T, digit = 6)})
    stan_con_od <- data.frame(OD = stan_mean,
                              concentration = (stan_concentration() %>% as.numeric())
    )
    lm(data = stan_con_od, concentration ~ OD)
  })
  #sample
  sam_row_col <- reactive({
    str_split(input$sam_pos, pattern = "")
  })
  sam_replication <- reactive(input$sam_replication)
  samplesize <- reactive(input$samplesize)
  sam_mean <- reactive({
    data <- quan_data()
    sam_row_start <- sam_row_col()[[1]][1]
    sam_col_start <- sam_row_col()[[1]][2] %>% as.numeric()
    samplesize <- samplesize()
    sam_replication <- sam_replication()
    sam_df <- data.frame(rep(NA, sam_replication))
    sam_times <- 0
    words <- c("A","B","C","D","E","F","G","H")
    
    while(sam_times < samplesize*sam_replication) {
      sam_df <- data.frame(sam_df, 
                           t(data[word2num_fun(sam_row_start),
                                  (sam_col_start:(sam_col_start+sam_replication-1))]))
      sam_col_start <- sam_col_start + sam_replication
      while ((sam_col_start -1) >= 12) {
        sam_col_start <- 1
        sam_row_start <- words[word2num_fun(sam_row_start)+1]
      }
      sam_times <- sam_times + sam_replication
    }
    
    sam_df<- sam_df[ ,-1]
    sapply(sam_df, function(x){mean(x,na.rm=T, digits = 6)})
    
  })
  
  output$quantitative_results <- renderPrint({
    if (is.null(quan_data())) {
      cat("",sep = "")
    } else {
      stan_sum <- stan_lm() %>% summary()
      r_squared <- stan_sum$r.squared %>% as.character()
      a <- stan_sum$coefficients[2,1] %>% as.character()
      b <- stan_sum$coefficients[1,1] %>% as.character()
      if (b %>% as.numeric() > 0) {
        cat("R-squared: ",r_squared,"\n","Formula: ", "y = ", a,"x +",b,sep = "")
      } else {
        cat("R-squared: ",r_squared,"\n","Formula: ", "y = ", a,"x ",b,sep = "")
      }
    }
  })
  
  output$quantitative_prediction <- renderPrint({
    if (is.null(quan_data())) { cat("",sep = "")
    } else {
      fitted_value <- predict(object= stan_lm() ,newdata = data.frame(OD= sam_mean()))
      fitted_value
    }
  })
  
  output$dilution_results <- renderPrint({
    if (is.null(quan_data())) { cat("",sep = "")
    } else {
      fitted_value <- predict(object= stan_lm() ,newdata = data.frame(OD= sam_mean()))
      fitted_value * input$sam_dilution_times
    }
  })
  
  output$volume_need <- renderPrint({
    if (is.null(quan_data())) { cat("",sep = "")
    } else {
      fitted_value <- predict(object= stan_lm() ,newdata = data.frame(OD= sam_mean()))
      input$weight_need / (fitted_value * input$sam_dilution_times) 
    }
  })
  output$quantitative_plot <- renderPlot({
    if (is.null(quan_data())) { return(NULL) 
    } else {
      ggplot(stan_lm()$model, aes(OD, concentration)) +
        geom_point()+
        geom_smooth(method = "lm", se = F) +
        scale_y_continuous(breaks=seq(0, 1.2, 0.2))
    }
  })
  # output$testtable <- renderPrint({
  #   if (is.null(quan_data())) {return()
  #  } else {
  #     summary(stan_lm())
  #   }
  # })
  
  
  #Descriptive Statistics
  descriptivedata <- reactive({
    infile <- input$descriptivefile
    if (is.null(infile)) {return(NULL)
    } else {
      df <- read.csv(infile$datapath, header = T)
      names(df) <- c("group", "value")
      df
    }
  })
  
  output$descriptive_summary <- renderPrint(
    if (is.null(descriptivedata())) {
      cat("",sep = "")    } else {
        Summarize(value ~ group, data = descriptivedata())
      }
  )
  
  output$descriptive_barplot <- renderPlot(
    width = 450, height = 450,
    if (is.null(descriptivedata())) {return(NULL)
    } else {
      histogram(~ value | group,  data=descriptivedata(), 
                layout=c(1,length(levels(descriptivedata()$group))), 
                panel=function(x, ...) {                             
                  panel.histogram(x, ...) 
                  panel.mathdensity(dmath = dnorm, col = "black",
                                    args = list(mean=mean(x),sd=sd(x)))          
                })
    }
  )
  
  output$descriptive_boxplot <- renderPlot(
    width = 450, height = 450,
    if (is.null(descriptivedata())) {return(NULL)
    } else {
      boxplot(value ~ group,
              data = descriptivedata(),
              ylab="Value",
              xlab="Group")
    }
  )
  #Normality Test
  normalitydata <- reactive({
    infile <- input$normalityfile
    if (is.null(infile)) {
      # User has not uploaded a file yet
      return(NULL)
    }
    df <- read.csv(infile$datapath, header = T)
    names(df) <- c("group", "value")
    df
  })
  normalitytestmethod<- reactive(input$Normalitytestmethod)
  
  output$normalitytest_summary<- renderPrint({
    if (is.null(normalitydata())) { cat("",sep = "")
    } else {
      data_group <- base::split(normalitydata(), normalitydata()$group)
      result<- paste0("Normality test", "\n")
      for(i in seq(length(data_group))) {
        group.name<- data_group[[i]]$group[1]
        result <- paste0(result,"Group: ", group.name)
        if (normalitytestmethod()== 'Shapiro-Wilk test (3<n<50)') {
          shapiro.result <- shapiro.test(data_group[[i]]$value)
          words <- paste0(", Shapiro-Wilk normality test W = ",
                          shapiro.result$statistic,
                          " p-value = ", shapiro.result$p.value, "\n")
          result <- paste0(result, words)
          
        } else if (normalitytestmethod()=='Kolmogorov-Smirnov test') {
          ks.result <- ks.test(data_group[[i]]$value, pnorm,
                               mean = mean(data_group[[i]]$value),
                               sd = sd(data_group[[i]]$value))
          words <- paste0(", Kolmogorov-Smirnov normality test D = ",
                          ks.result$statistic,
                          " p-value = ", ks.result$p.value, "\n")
          result <- paste0(result, words)
        }
      }
      cat(result,sep = "\n", append = T)
    }
    
  })
  qqplot_line<- function(x) {
    par(pin=c(4,4))
    qqnorm(x)
    qqline(x)
  }
  output$graphical_test<- renderPlot({
    normaltest<- normalitydata()
    if (is.null(normaltest)) return(NULL)
    qqplot_line(normaltest[ ,2])
  })
  #  output$normalitytest_table<- renderTable({
  #    normalitydata()
  #  })
  
  #Heteroscedasticity test
  heterodata <- reactive({
    infile <- input$heterofile
    if (is.null(infile)) {
      return(NULL)
    } else {
      df <- read.csv(infile$datapath, header = T)
      names(df) <- c("group", "value")
      df    }
  })
  output$hetero_result <- renderPrint(
    if (is.null(heterodata())) { cat("",sep = "")
    } else {
      leveneTest(value~group, heterodata())
    }
    
  )
  
  #multiple comparasion
  mcompardata <- reactive({
    infile <- input$mcomparefile
    if (is.null(infile)) {
      # User has not uploaded a file yet
      return(NULL)
    }
    data <- read.csv(infile$datapath, header = T)
    names(data) <- c("group", "value")
    df
    mutate(data, group = factor(group, levels=unique(group)))
  })
  mcomparemethod <- reactive(input$mcomparemethod)
  mcompareresult <- reactive({
    data<- mcompardata()
    lm(value ~ group -1, data= data)
  })
  output$mcompare_summary <- renderPrint(
    if (is.null(mcompardata())) { cat("",sep = "")
    } else if ( mcomparemethod() == 'Least Squares Means') {
      data <- mcompareresult()
      summary(data)
    }
  )
  output$mcompare_plot <- renderPlot(
    width = 450, height = 450,
    if (is.null(mcompardata())) {return(NULL)
    } else {
      Info<- summary(mcompareresult())
      Coef <- as.data.frame(Info$coefficients[, 1:2])
      Coef <- within(Coef, {
        Lower <- Estimate - qt(p=0.90, df=Info$df[2]) * `Std. Error`
        Upper <- Estimate + qt(p=0.90, df=Info$df[2]) * `Std. Error`
        group <- rownames(Coef)
      })
      ggplot(Coef, aes(x=Estimate, y=group)) + geom_point() +
        geom_errorbarh(aes(xmin=Lower, xmax=Upper), height=.3) +
        ggtitle("value by group calculated from regression model")
    }
  )
  output$mcompare_hoctest <- renderPrint(
    if (is.null(mcompardata())) { cat("",sep = "")
    } else {
      lsm<- lsmeans(mcompareresult(), "group")
      contrast(lsm, method = "pairwise", interaction = T)
    }
  )
  
  #Non-Parametric Tests
  nonpardata <- reactive({
    infile <- input$nonparfile
    if (is.null(infile)) {
      # User has not uploaded a file yet
      return(NULL)
    }
    data<- read.csv(infile$datapath)
    names(data) <- c("group", "value")
    mutate(data, group = factor(group, levels=unique(group)))
  })
  nonparmethod<- reactive(input$nonparmethod)
  
  output$NonParametricTests_summary<- renderPrint({
    if (is.null(nonpardata())) { cat("",sep = "")
    } else if (nonparmethod()== 'Kruskal-Wallis test for equal variances') {
      kruskal.result<- kruskal.test(value ~ group, data = nonpardata())
      kruskal.result
    } else if (nonparmethod()== "Welch's anova for unequal variances") {
      Welch.result<- oneway.test(value ~ group, data = nonpardata(),var.equal=FALSE)
      Welch.result
    }
  })
  
  output$NonParametricTests_hoctest <- renderPrint({
    if (is.null(nonpardata())) { cat("",sep = "")
    } else if (nonparmethod()== 'Kruskal-Wallis test for equal variances') {
      #g <- factor(nonpardata()$group, levels=levels(nonpardata()$group))
      dunntest <- dunnTest(value ~ group,data=nonpardata(),method="bh")
      out <- capture.output(dunntest)
      cat("### Dunn test\n", out, sep = "\n")
      #Dunn (1964) Kruskal-Wallis multiple comparison
      #p-values adjusted with the Benjamini-Hochberg method.
    } else if (nonparmethod()== "Welch's anova for unequal variances") {
      userfriendlyscience::oneway(nonpardata()$group, 
                                  y = nonpardata()$value, posthoc = 'games-howell')
    }
  })
  
  #Correlation analysis
  cordata <- reactive({
    infile <- input$correlationfile
    if (is.null(infile)) {
      return(NULL)
    }
    read.csv(infile$datapath)
  })
  cor.method <- reactive({
    if (input$method_cor == "Pearson's product-moment correltion") {
      returnValue("pearson")
    } else if (input$method_cor == "Spearman's rank order correlation") {
      returnValue("spearman")
    }
  })
  cor.result <- reactive({
    data <- cordata()
    method <- cor.method()
    if (ncol(data) == 2) {
      var1<- data[,1]
      var2<- data[,2]
      cor.test(var1,var2, method = method)
    } else {
      cor(data, method = method)
    }
  })
  
  output$cor_test <- renderPrint({
    if (is.null(cordata())) {
      cat("","")
    } else if (ncol(cordata()) !=2) {
      cat("","")
    } else {
      cor.result()
    }
  })
  
  output$cor_table <- renderTable({
    if (is.null(cordata())) {
      return(NULL)
    } else if (ncol(cordata()) ==2) {
      return(NULL)
    } else {
      cor.result()
    }
  })
  
  output$Correlation.csv <- downloadHandler(
    filename = function() {
      paste("Correlation",".csv",sep = "")
    },
    content = function(fname) {
      cor.result<- as.data.frame(cor(cordata(), method = cor.method()))
      write.csv(cor.result, file = fname)
    },
    contentType = NA
  )
  output$cor_plot <- renderPlot({
    if (is.null(cordata())) { return(NULL) 
    } else {
      corrplot.mixed(cor(cordata(), method = cor.method()),
                     lower = "pie", upper = "circle")
    }
    
  })
  output$Correlation_plot.pdf <- downloadHandler(
    filename = function() {
      paste("Correlation_plot",".pdf",sep = "")
    },
    content = function(fname) {
      cor.result<- cor(cordata(), method = cor.method())
      pdf(file = fname)
      corrplot.mixed(cor.result, lower = "pie", upper = "circle")
      dev.off()
    },
    contentType = NA
  )
  
  # Setting gene database
  gene_db <- reactive({
    infile <- input$species_choice
    ah <- AnnotationHub()
    databaseid<- query(ah, c(infile,"OrgDb")) 
    ah[[databaseid$ah_id]]
    
  })
  
  #Gene ID coversion
  gene_label2value <- function(label) {
    if(label == 'Accession') {
      returnValue("ACCNUM")
    }
    else if (label =='Gene symbol') {
      returnValue("SYMBOL")
    }
    else if (label =='Entrez ID') {
      returnValue('ENTREZID')
    }
    else if (label == 'UniProt'){
      returnValue('UNIPROT')
    }
    else if (label == 'Ensembl') {
      returnValue("ENSEMBL")
    }
    else if (label == 'GO') {
      returnValue("GO")
    }
    else if (label == 'Gene name') {
      returnValue("GENENAME")
    }
  }
  current_gene_type<- reactive({
    infile<- input$typefrom
    gene_label2value(infile)
  })
  convert_gene_type<- reactive({
    infile<- input$typeto
    gene_label2value(infile)
  })
  convertiddata <- reactive({
    infile<- input$convertidfile
    if (is.null(infile)) {
      return(NULL)
    }
    data_convertid<- read.csv(infile$datapath)
    if (current_gene_type()=='ACCNUM'){
      str_extract(pattern = "^[A-Za-z0-9_]{6}[0-9]*",data_convertid[ ,1])
    } else if (current_gene_type()=='SYMBOL'){
      str_extract(pattern = "[A-Za-z0-9]+",data_convertid[ ,1])
    } else if (current_gene_type()=='ENTREZID') {
      str_extract(pattern = "[0-9]+",data_convertid[ ,1])
    } else {data_convertid[ ,1]}
  })
  
  convertid <- reactive({
    if (is.null(convertiddata())) {
      return(NULL)
    }
    else {
      id<- mapIds(x = gene_db() ,keys = convertiddata(),
                  keytype = current_gene_type(),
                  column = convert_gene_type()
      )
      df<- data.frame(names(id),id)
      names(df) <- c(input$typefrom, input$typeto)
      df
    }
  })
  #output$idtest <- renderPrint({convert_gene_type()})
  output$geneconvert_table<- renderTable({
    convertid()
  })
  output$Convert_gene_id.csv <- downloadHandler(
    filename = function() {
      paste("Converting gene id",".csv",sep = "")
    },
    content = function(fname) {
      write.csv(convertid(), file = fname,
                row.names = F)
    },
    contentType = NA
    
    
  )
  
  #GO classification
  current_gene_type_ggo <- reactive({
    infile<- input$fromtype_ggo
    if (infile =='Accession'){
      returnValue("ACCNUM")
    }
    else if (infile =='Gene symbol') {
      returnValue("SYMBOL")
    }
    else if (infile =='Entrez ID') {
      returnValue('ENTREZID')
    }
  })
  regulated_options_ggo <- reactive(input$regulated_gene_ggo)
  genelist <- reactive({
    infile_ggo <- input$file_ggo
    if (is.null(infile_ggo)) {
      return(NULL)
    }
    ggodata <- read.csv(infile_ggo$datapath)
    if(current_gene_type_ggo() == 'ACCNUM'){
      id <- str_extract(pattern = "^[A-Za-z0-9_]{6}[0-9]*",ggodata[ ,1])
      id <- mapIds(x = gene_db() ,keys = id,
                   keytype = "ACCNUM" ,column = "ENTREZID")
    } else if(current_gene_type_ggo() == 'SYMBOL'){
      id <- str_extract(pattern = "[A-Za-z0-9]+",ggodata[ ,1])
      id <- mapIds(x = gene_db() ,keys = id,
                   keytype = "SYMBOL" ,column = "ENTREZID")
    } else if(current_gene_type_ggo() == 'ENTREZID') {
      id <- str_extract(pattern = "[0-9]+",ggodata[ ,1]) 
      id <- mapIds(x = gene_db() ,keys = id,
                   keytype = 'ENTREZID' ,column = "ENTREZID") 
    }
    ggodata[ ,2] %<>% as.character() %>% 
      str_extract_all("^[0-9]*([0-9]|\\.)([0-9]*)") %>% as.numeric()
    genelist_total <- data.frame(geneid = id, 
                                 ratio = ggodata[ ,2],
                                 stringsAsFactors = F) 
    genelist_total%<>% filter(!is.na(genelist_total$geneid))
    genelist <- {}
    if ("up_regulated" %in% regulated_options_ggo() ) {
      genelist_up <- genelist_total %>% filter(ratio > 1)
      genelist_up$regulation <- rep('up',nrow(genelist_up))
      genelist[["up_regulated"]] <- genelist_up
      
    }
    if ("down_regulated" %in% regulated_options_ggo() ) {
      genelist_down <- genelist_total %>% filter(ratio < 1)
      genelist_down$regulation <- rep('down',nrow(genelist_down))
      genelist[["down_regulated"]] <- genelist_down
    }
    if ("total" %in% regulated_options_ggo() ) {
      genelist_total$regulation <- rep('total',nrow(genelist_total))
      genelist[["total"]] <- genelist_total
    }
    genelist
  })
  ggo_results <- reactive({
    ontology<- c("BP","CC","MF")
    if ("total" %in% regulated_options_ggo()){
      ggo_list_total<- {}
    }
    if ("up_regulated" %in% regulated_options_ggo() ) {
      ggo_list_up<- {}
    }
    if ("down_regulated" %in% regulated_options_ggo() ) {
      ggo_list_down<- {}
    }
    for(geneid in genelist()) {
      for(onts in ontology) {
        ggo <- groupGO(gene= as.character(geneid[,1]), 
                       OrgDb = gene_db(), keyType = "ENTREZID",
                       ont= onts, level=2, readable=TRUE) %>% as.data.frame()
        if (onts == "BP") {
          ggo$ontology<- rep("Biological process" ,nrow(ggo))
        } else if (onts == "CC") {
          ggo$ontology<- rep("Cellular componenent" ,nrow(ggo))
        } else {
          ggo$ontology<- rep("Molecular function" ,nrow(ggo))
        }
        ggo %<>% filter(Count > 0)
        if(TRUE %in% (geneid[,2] >1) && FALSE %in% (geneid[,2] >1)) {
          name_ggo<- paste0('total',"_",onts)
          ggo_list_total[[name_ggo]]<- ggo
        } else if(FALSE %in% (geneid[,2] >1)) {
          name_ggo<- paste0('down_regulated',"_",onts)
          ggo_list_down[[name_ggo]] <- ggo
        } else {
          name_ggo<- paste0('up_regulated',"_",onts)
          ggo_list_up[[name_ggo]] <- ggo
        }
      }
    }
    ggo_list <- list(groupGO_total_lv2 = ggo_list_total,
                     groupGO_up_lv2 = ggo_list_up,
                     groupGO_down_lv2 = ggo_list_down)
    rbind_list <- function(x) {ldply(x,rbind)}
    ggo_list %>% llply(rbind_list)
  })
  
  output$ggo_results.zip <- downloadHandler(
    filename =  function() {
      paste("ggo_results","zip",sep = ".")
    },
    contentType = "application/zip",
    content = function(fname) {
      paste_csv <- function(x) { paste0(x, '.csv')}
      fs <- c()
      tmpdir <- tempdir()
      setwd(tempdir())
      progress <- shiny::Progress$new()
      on.exit(progress$close())
      progress$set(message = "Begin to download results：", value = 0)
      
      for (filenames in names(ggo_results())) {
        fs <- c(fs,paste_csv(filenames))
        write.csv(eval(parse(text = paste0("ggo_results()[['",filenames,"']]"))),
                  file = paste0(filenames,".csv") ,row.names = F)
        progress$inc(1/length(names(ggo_results())), detail = "Please wait...")
      }
      zip::zip(zipfile= fname, files=fs)
      #https://tw.saowen.com/a/502ce2258e11f41406f355c727683774d48772c07d42e1adce948c9bc0d73c36
    }
  )
  output$ggo_plots.zip <- downloadHandler(
    filename =  function() {
      paste("ggo_plots","zip",sep = ".")
    },
    content = function(fname) {
      fs <- c()
      tmpdir <- tempdir()
      setwd(tempdir())
      progress <- shiny::Progress$new()
      on.exit(progress$close())
      progress$set(message = "Begin to download results：", value = 0)
      for (sort_df in ggo_results()) {
        sort_df %<>% dplyr::select(ontology,Description,Count, .id)
        sort_df <- sort_df[order(sort_df[['ontology']],-sort_df[['Count']]), ]
        sort_df$Description <- factor(sort_df$Description,
                                      levels= sort_df$Description)
        sort_df_BP <- filter(sort_df, ontology == "Biological process")
        sort_df_CC <- filter(sort_df, ontology == "Cellular componenent")
        sort_df_MF <- filter(sort_df, ontology == "Molecular function")
        if (nrow(sort_df_BP) >10) {
          sort_df_BP <- sort_df_BP[1:10, ]
        }
        if (nrow(sort_df_CC) >10) {
          sort_df_CC <- sort_df_CC[1:10, ]
        }
        if (nrow(sort_df_MF) >10) {
          sort_df_MF <- sort_df_MF[1:10]
        }
        fs <- c(fs,paste('groupgo_',
                         str_extract_all(sort_df$.id[1],pattern = "^[a-z]*"),
                         '.png',sep = ''))
        sort_df <- rbind(sort_df_BP,sort_df_CC,sort_df_MF)
        sort_df$Description <- factor(sort_df$Description,levels=sort_df$Description)
        
        ggplot(data = sort_df, mapping = aes(x=Description ,y= Count,fill=ontology))+
          geom_bar(stat="identity") +
          coord_flip()+
          scale_x_discrete(limits=rev(levels(sort_df$Description)))+
          labs( x='Term', fill ='Ontology') +
          theme_bw(base_size = 10)
        #ggplot(data = sort_df, mapping = aes(x=rev(Description) ,y= Count,fill=ontology))+
         # geom_bar(stat="identity",position = "dodge") +
        #  scale_y_continuous(labels = function (x) floor(x)) +
        #  coord_flip()+
        #  labs( x='Term', fill ='Ontology') +
        #  theme_bw(base_size = 10)
        ggsave(paste('groupgo_',
                     str_extract_all(sort_df$.id[1],pattern = "^[a-z]*"),
                     '.png',sep = ''))
        progress$inc(1/length(names(ggo_results())), detail = "Please wait...")
      }
      zip::zip(zipfile= fname, files=fs)
      
    }
    
  )
  output$ggo_table <- renderTable({
    if (is.null(genelist())) {
      return(NULL)
    } else {
      genelist()[["total"]]
    }
    
  })
  #GO enrichment
  gene_db_ego <- reactive({
      infile <- input$species_choice_ego
      ah <- AnnotationHub()
      databaseid<- query(ah, c(infile,"OrgDb")) 
      ah[[databaseid$ah_id]]
    })
  
  current_gene_type_ego <- reactive({
    infile<- input$fromtype_ego
    if (infile =='Accession'){
      returnValue("ACCNUM")
    }
    else if (infile =='Gene symbol') {
      returnValue("SYMBOL")
    }
    else if (infile =='Entrez ID') {
      returnValue('ENTREZID')
    }
  })
  
  regulated_options_ego <- reactive(input$regulated_gene_ego)
  
  ego_result_options <- reactive(input$ego_options)
  egodata <- reactive({
    infile_ego <- input$file_ego
    if (is.null(infile_ego)) {
      return(NULL)
    }
    read.csv(infile_ego$datapath)
  })
  
  genelist_ego <- reactive({ #make genelist
    if (is.null(egodata())) {
      return(NULL)
    } else {
      genedata <- egodata() 
      current <- current_gene_type_ego()
      regulated <- regulated_options_ego()
      if(current == 'ACCNUM'){
        id <- str_extract(pattern = "^[A-Za-z0-9_]{6}[0-9]*",genedata[ ,1])
        id <- mapIds(x = gene_db_ego() ,keys = id,
                     keytype = "ACCNUM" ,column = "ENTREZID")
      } else if(current == 'SYMBOL'){
        id <- str_extract(pattern = "[A-Za-z0-9]+",genedata[ ,1])
        id <- mapIds(x = gene_db_ego() ,keys = id,
                     keytype = "SYMBOL" ,column = "ENTREZID")
      } else if(current == 'ENTREZID') {
        id <- str_extract(pattern = "[0-9]+",genedata[ ,1])
        id <- mapIds(x = gene_db_ego() ,keys = id,
                     keytype = 'ENTREZID' ,column = "ENTREZID") 
      }
      genedata[ ,2] %<>% as.character() %>% 
        str_extract_all("^[0-9]*([0-9]|\\.)([0-9]*)") %>% as.numeric()
      genelist_total <- data.frame(geneid = id, 
                                   ratio = genedata[ ,2],
                                   stringsAsFactors = F) 
      genelist_total%<>% filter(!is.na(genelist_total$geneid))
      genelist <- {}
      if ("up_regulated" %in% regulated ) {
        genelist_up <- genelist_total %>% filter(ratio > 1)
        genelist_up$regulation <- rep('up',nrow(genelist_up))
        genelist[["up_regulated"]] <- genelist_up
        
      }
      if ("down_regulated" %in% regulated ) {
        genelist_down <- genelist_total %>% filter(ratio < 1)
        genelist_down$regulation <- rep('down',nrow(genelist_down))
        genelist[["down_regulated"]] <- genelist_down
      }
      if ("total" %in% regulated ) {
        genelist_total$regulation <- rep('total',nrow(genelist_total))
        genelist[["total"]] <- genelist_total
      }
      genelist
    }
    }) 
  
  ego_results <- reactive({
    ego_results<- list()
    ontology<- c("BP","CC","MF")
    
    for(geneid_ego in genelist_ego()) {
      for(onts in ontology)
      {
        ego<- enrichGO(gene         = geneid_ego[ ,1], 
                       OrgDb         = gene_db_ego(), 
                       keyType = "ENTREZID",
                       ont = onts, #'ALL','BP',"MF','CC'
                       pAdjustMethod = "BH", # BH is the other name of fdr
                       pvalueCutoff = 0.05,
                       qvalueCutoff  = 0.2,
                       minGSSize = 10, 
                       maxGSSize = 500, 
                       readable = T, #"ENTREZID"convert'SYMBOL'
                       pool = F)
        ego@result$regulation <- rep(geneid_ego$regulation[1], 
                                     nrow(ego@result))
        assign(paste0("ego_",geneid_ego$regulation[1],'_',onts), ego)
      }
      assign(value = list(BP = eval(parse(text = paste0("ego_",geneid_ego$regulation[1],'_',"BP"))),
                          CC = eval(parse(text = paste0("ego_",geneid_ego$regulation[1],'_',"CC"))),
                          MF = eval(parse(text = paste0("ego_",geneid_ego$regulation[1],'_',"MF")))
      ), x = paste0("ego_",geneid_ego$regulation[1]))
      ego_results[[paste0("ego_",geneid_ego$regulation[1])]] <- eval(parse(text = paste0("ego_",geneid_ego$regulation[1])))
    }
    ego_results
  })
  
  output$ego_results.zip  <- downloadHandler(
    filename = function() {
      paste0("ego_results",".zip")
    },
    content = function(fname){
      fs <- c()
      tmpdir <- tempdir()
      setwd(tempdir())
      progress <- shiny::Progress$new()
      on.exit(progress$close())
      progress$set(message = "Begin to download results：", value = 0)
      
      ego_list_to_df <- function(x) {
        ldply(x , data.frame)
      }
      ego_list_df<- llply(ego_results(), ego_list_to_df)
      for (result in ego_list_df) {
        write.csv(result, file = paste0("ego_",result$regulation[1],".csv"),
                  row.names =FALSE)
        fs <- c(fs, paste0("ego_",result$regulation,".csv"))
        progress$inc(1/length(ego_list_df), detail = "Please wait...")
      }
      
      zip::zip(zipfile = fname, files = fs)
    }
  )
  
  output$ego_plots.zip <- downloadHandler(
    filename = function() {
      paste0("ego_plots",".zip")
    },
    content = function(fname){
      fs <- c()
      tmpdir <- tempdir()
      setwd(tempdir())
      progress <- shiny::Progress$new()
      on.exit(progress$close())
      progress$set(message = "Begin to download plots：", value = 0)
      
      ego_plot<- function(object, plot.type) {
        if(plot.type == "dotplot_ego") {
          object %>% clusterProfiler::dotplot(x= "GeneRatio",
                                              font.size= 10,
                                              color ="p.adjust", 
                                              showCategory = 20#, #list how much items
                                              #title= 'up-regulated'
          ) +  
            ggplot2::scale_color_continuous(low='blue', high='red') +
            scale_size(range=c(2, 7))
        }
        if(plot.type == "dag_ego") {
          pdf(file = paste('enrichgo_',object@result$regulation[1],
                           '_',object@ontology,'_',
                           str_extract_all(plot.type,pattern = "^[a-z|A-Z]*[^_]*"),
                           '.pdf',sep = ''))
          plotGOgraph(object)
          dev.off()
        }
        if(plot.type == "Gene-Concept-Network_ego") {
          object %>% cnetplot(categorySize = "p.adjust", showCategory = 5)
        }
        if(plot.type == "emapplot_ego") {
          object %>% emapplot(showCategory = 20, color = "p.adjust", layout = "kk",vertex.label.cex=1.2)
        }
        if(plot.type != "dag_ego") {
          ggsave(paste('enrichgo_',object@result$regulation[1],
                       '_',object@ontology,'_',
                       str_extract_all(plot.type,pattern = "^[a-z|A-Z]*[^_]*"),
                       '.png',sep = ''))
        }
        
      }
      
      if("total" %in% regulated_options_ego()) {
        for (result in ego_results()[["ego_total"]]) {
          for ( plot_type in ego_result_options() ) {
            ego_plot(result, plot_type)
            if(nrow(data.frame(result)) != 0) {
              ego_plot(result, plot_type)
              if (plot_type != "dag_ego") {
                fs <- c(fs,paste('enrichgo_',result@result$regulation[1],
                                 '_',result@ontology,'_',
                                 str_extract_all(plot_type,pattern = "^[a-z|A-Z]*[^_]*"),
                                 '.png',sep = ''))
              } else {
                fs <- c(fs,paste('enrichgo_',result@result$regulation[1],
                                 '_',result@ontology,'_',
                                 str_extract_all(plot_type,pattern = "^[a-z|A-Z]*[^_]*"),
                                 '.pdf',sep = ''))
              }
            } else {next}
            progress$inc(1/(length(ego_results())*length(ego_result_options())*3), detail = "Please wait...")
          }
        }
      }
      if("up_regulated" %in% regulated_options_ego()) {
        for (result in ego_results()[["ego_up"]]) {
          for ( plot_type in ego_result_options() ) {
            ego_plot(result, plot_type)
            if(nrow(data.frame(result)) != 0) {
              ego_plot(result, plot_type)
              if (plot_type != "dag_ego") {
                fs <- c(fs,paste('enrichgo_',result@result$regulation[1],
                                 '_',result@ontology,'_',
                                 str_extract_all(plot_type,pattern = "^[a-z|A-Z]*[^_]*"),
                                 '.png',sep = ''))
              } else {
                fs <- c(fs,paste('enrichgo_',result@result$regulation[1],
                                 '_',result@ontology,'_',
                                 str_extract_all(plot_type,pattern = "^[a-z|A-Z]*[^_]*"),
                                 '.pdf',sep = ''))
              }
            } else {next}
            progress$inc(1/(length(ego_results())*length(ego_result_options())*3), detail = "Please wait...")
            
          }
        }
      }
      if("down_regulated" %in% regulated_options_ego()) {
        for (result in ego_results()[["ego_down"]]) {
          for ( plot_type in ego_result_options() ) {
            if(nrow(data.frame(result)) != 0) {
              ego_plot(result, plot_type)
              if (plot_type != "dag_ego") {
                fs <- c(fs,paste('enrichgo_',result@result$regulation[1],
                                 '_',result@ontology,'_',
                                 str_extract_all(plot_type,pattern = "^[a-z|A-Z]*[^_]*"),
                                 '.png',sep = ''))
              } else {
                fs <- c(fs,paste('enrichgo_',result@result$regulation[1],
                                 '_',result@ontology,'_',
                                 str_extract_all(plot_type,pattern = "^[a-z|A-Z]*[^_]*"),
                                 '.pdf',sep = ''))
              }
            } else {next}
            progress$inc(1/(length(ego_results())*length(ego_result_options())*3), detail = "Please wait...")
          }
        }
      }
      
      zip::zip(zipfile = fname, files = fs)
    }
    
  )
  
  output$test_ego<- renderPrint({
    if(length(ego_results()) == 0) {
      cat("Not yet" ,sep = "")
    } else {
      cat("Ready!" ,sep = "")
    }
  })
  output$ego_table <- renderTable({
    if (length(ego_results()) == 0) {
      return(NULL)
    } else {
      ego<- ego_results()[["ego_total"]]
      ego <- ego[["ego_total_BP"]]
      data.frame(ego)
    }
  })
  
#  output$test <- renderPrint({
 #   input$species_choice_ego
  #})
  
  
  #KEGG enrichment
  current_gene_type_ekegg <- reactive({
    infile<- input$fromtype_ekegg
    gene_label2value(infile)
  })
  regulated_options_ekegg <- reactive(input$regulated_gene_ego)
  result_options_ekegg <- reactive(input$ekegg_options)
  ekeggdata <- reactive({
    infile <- input$file_ekegg
    if (is.null(infile)) {
      return(NULL)
    }
    read.csv(infile$datapath)
  })
  
  genelist_ekegg <- reactive({
    if (is.null(ekeggdata())) {
      return(NULL)
    } else {
      genelist_fun(genedata = ekeggdata(), current = current_gene_type_ekegg(),
                   regulated = regulated_options_ekegg())
    }
  })
  kegg_species<- reactive({
    infile <- input$species_choice
    if(infile == "Gallus gallus") {
      returnValue("gga")
    }
    else if (infile== "Sus scrofa") {
      returnValue("ssc")
    }
  })
  ekegg_results <- reactive({
    ekegg_results<- list()
    # kegg_id<- bitr_kegg(genelist_down$ENTREZID,
    #                    fromType = "ncbi-geneid",
    #                     toType = "kegg",
    #                   organism='gga')
    for(geneid_ekegg in genelist_ekegg()) {
      ekegg<- enrichKEGG(geneid_ekegg$geneid , organism= kegg_species(),
                         keyType =  "ncbi-geneid",
                         pvalueCutoff=0.05, pAdjustMethod="BH",
                         qvalueCutoff=0.1,
                         use_internal_data= F) %>% 
        setReadable(., OrgDb = gene_db(),keytype = "ENTREZID")
      ekegg@result$regulation <- rep(geneid_ekegg$regulation[1], 
                                     nrow(ekegg@result))
      ekegg_results[[paste0("ekegg_",geneid_ekegg$regulation[1])]] <- ekegg
    }
    ekegg_results
  })
  
  output$ekegg_results.zip <- downloadHandler(
    filename = function() {
      paste("ekegg_results",".zip",sep = "")
    },
    content = function(fname) {
      fs<- c()
      tmpdir <- tempdir()
      setwd(tempdir())
      progress <- shiny::Progress$new()
      on.exit(progress$close())
      progress$set(message = "Begin to download results：", value = 0)
      
      for (result in ekegg_results()) {
        write.csv(data.frame(result), 
                  file = paste("enrichkegg_",result@result$regulation[1],
                               ".csv", sep = ""),
                  row.names = F)
        fs <- c(fs,paste("enrichkegg_",result@result$regulation[1],
                         ".csv", sep = ""))
        progress$inc(1/length(ekegg_results()), detail = "Please wait...")
      }
      zip::zip(zipfile = fname, files = fs)
    }
  )
  
  output$ekegg_plots.zip <- downloadHandler(
    filename = function() {
      paste("ekegg_plots",".zip",sep = "")
    },
    content = function(fname) {
      fs<- c()
      tmpdir <- tempdir()
      setwd(tempdir())
      progress <- shiny::Progress$new()
      on.exit(progress$close())
      progress$set(message = "Begin to download plots：", value = 0)
      
      ekegg_plot<- function(object, plot.type) {
        if(plot.type == "dotplot_ekegg") {
          object %>% clusterProfiler::dotplot(x= "GeneRatio",
                                              font.size= 10,
                                              color ="p.adjust", 
                                              showCategory = 20) 
        }
        if(plot.type == "Gene-Concept-Network_ekegg") {
          object %>% cnetplot(categorySize = "p.adjust", showCategory = 5)
        }
        if(plot.type == "emapplot_ekegg") {
          object %>% emapplot(showCategory = 20, color = "p.adjust", layout = "kk",vertex.label.cex=1.2)
        }
        ggsave(paste('enrichkegg_',result@result$regulation[1],'_',
                     str_extract_all(plot_type,pattern = "^[a-z|A-Z]*[^_]*"),
                     '.png',sep = ''))
      }
      
      for (result in ekegg_results()) {
        for (plot_type in result_options_ekegg()) {
          if(nrow(data.frame(result)) != 0) {
            ekegg_plot(result, plot_type)
            fs <- c(fs, paste('enrichkegg_',result@result$regulation[1],'_',
                              str_extract_all(plot_type,pattern = "^[a-z|A-Z]*[^_]*"),
                              '.png',sep = ''))
          } else {next}
          progress$inc(1/(length(ekegg_results())*length(result_options_ekegg())),
                       detail = "Please wait...")
        }
      }
      zip::zip(zipfile = fname, files = fs)
      
    }
  )
  output$test_ekegg<- renderPrint({
    result_options_ekegg()
    if(length(ekegg_results()) == 0) {
      cat("Not yet" ,sep = "")} else {
        cat("Ready!" ,sep = "")
      }
  })
  
  #Cross match
  crossmatchdata1 <- reactive({
    infile <- input$crossmatch1file
    if (is.null(infile)) {return(NULL)
    } else {
      read.csv(infile$datapath, header = T)
    }
  })
  crossmatchdata2 <- reactive({
    infile <- input$crossmatch2file
    if (is.null(infile)) {return(NULL)
    } else {
      read.csv(infile$datapath, header = T, sep = "\t")
    }
  })
  
  crossmatchresult <- reactive({
    if (is.null(crossmatchdata1())) {
      return(NULL)
    } else if (is.null(crossmatchdata1()$Accession)) {
      cat("Please rename the title of columns with 'Accession' ")
    }  else {
      data1 <- crossmatchdata1()
      data2 <- crossmatchdata2()
      data1$Accession %<>% sub("..$","", . ) #remove . (version)
      data1$EntrezID <- mapIds(x = gene_db() ,keys = data1$Accession ,
                               keytype = "ACCNUM" ,column = "ENTREZID")
      data1 %>% filter(EntrezID %in% data2$GeneID | EntrezID %in% NA)
      
    }
  })
  
  output$crossmatch_result.csv <- downloadHandler(
    filename = function() {
      paste0("matchacross_result", ".csv")
    },
    content = function(fname) {
      write.csv(crossmatchresult(), 
                file = fname,
                row.names = F)
    },
    contentType = NA
  )
  
})
