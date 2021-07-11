library(shiny)
library(tidyverse)
# library(viridis)
library(plotly)
library(data.table)
library(MASS)
library(shinythemes)
# library(shinyjs)
library(latex2exp)

ui <- navbarPage("Power Caculation in SW-CRT",

                 theme = shinytheme("cosmo"),

                 tabPanel("Instruction", includeMarkdown("introduction.rmd"), style = 'width:1000px;'),
                 

                 tabPanel("App",
                          style = 'width:3000px;',
                          sidebarLayout(
                            sidebarPanel(
                              radioButtons("model", "fitted model",
                                           choices = list("Model 1: Fixed treatment effect with no secular trend" = 1, 
                                                          "Model 2: Fixed treatment effect with linear secular trend" = 2,
                                                          "Model 3: Fixed treatment effect with non-linear secular trend" = 3,
                                                          "Model 4: Random treatment effect with no secular trend" = 4, 
                                                          "Model 5: Random treatment effect with linear secular trend" = 5,
                                                          "Model 6: Random treatment effect with non-linear secular trend" = 6), 
                                           selected = 1),
                              
                              
                              numericInput("seed"
                                           , "seed"
                                           , value = 1),
                              
                              
                              selectizeInput(
                                "cellSize"
                                , "cell size"
                                , choices = NULL
                                , multiple = TRUE
                                , options = list(create = TRUE)
                              ),
                              
                              selectizeInput(
                                "clusterSize"
                                , "cluster size"
                                , choices = NULL
                                , multiple = TRUE
                                , options = list(create = TRUE)
                              ),
                              
                              selectizeInput(
                                "groupSize"
                                , "group size"
                                , choices = NULL
                                , multiple = TRUE
                                , options = list(create = TRUE)
                              ),
                              
                              selectizeInput(
                                "nTransition"
                                , "transition length"
                                , choices = NULL
                                , multiple = TRUE
                                , options = list(create = TRUE)
                              ),

                              
                              selectizeInput(
                                "b0"
                                , "beta0"
                                , choices = NULL
                                , multiple = TRUE
                                , options = list(create = TRUE)
                              ),
                              
                              selectizeInput(
                                "b4"
                                , "beta4"
                                , choices = NULL
                                , multiple = TRUE
                                , options = list(create = TRUE)
                              ),
                              
                              selectizeInput(
                                "rhoG"
                                , "rho34"
                                , choices = NULL
                                , multiple = TRUE
                                , options = list(create = TRUE)
                              ),

                              selectizeInput(
                                "rhoO"
                                , "rho13"
                                , choices = NULL
                                , multiple = TRUE
                                , options = list(create = TRUE)
                              ),
                              
                              selectizeInput(
                                "rhoY"
                                , "rho24"
                                , choices = NULL
                                , multiple = TRUE
                                , options = list(create = TRUE)
                              ),
                              
                              radioButtons("secularTrend", "secular trend for true model",
                                           choices = list("no secular trend" = 0,
                                                          "linear secular trend" = 1,
                                                          "non-linear secular trend" = 2),
                                           selected = 0),
                              
  
                              selectizeInput(
                                "t"
                                , "time effect size"
                                , choices = NULL
                                , multiple = TRUE
                                , options = list(create = TRUE)
                              ),
                              

                              selectizeInput(
                                "SIGMA"
                                , "s.d. of random cluster effect"
                                , choices = NULL
                                , multiple = TRUE
                                , options = list(create = TRUE)
                              ),
                              
                              radioButtons("randomTrt", "treatment effect for true model",
                                           choices = list("fixed treatment effect only" = 0,
                                                          "include random treatment effect" = 1),
                                           selected = 0),
                              
                              selectizeInput(
                                "SIGMA2"
                                , "s.d. of random treatment effect"
                                , choices = NULL
                                , multiple = TRUE
                                , options = list(create = TRUE)
                              ),
 
                              numericInput(
                                "nIteration"
                                , "iteration"
                                , value = 10
                              ),
                              
                              selectInput(
                                "x1",
                                "val1",
                                choices = list("Beta1" = "trt_indtrt1",
                                               "Beta2" = "trt_indtrt2",
                                               "Beta3" = "trt_indtrt3",
                                               "Beta4" = "trt_indtrt4")),
                              
                              selectInput(
                                "x2",
                                "val2",
                                choices = list("Beta1" = "trt_indtrt1",
                                               "Beta2" = "trt_indtrt2",
                                               "Beta3" = "trt_indtrt3",
                                               "Beta4" = "trt_indtrt4")),
                              
                              selectizeInput(
                                "lincom"
                                , "linear combination"
                                , choices = NULL
                                , multiple = TRUE
                                , options = list(create = TRUE)
                              ),
                              
                              actionButton("go", "Submit"),
                              
                              selectInput("dataset", "Choose a dataset:",
                                          choices = c("Model Fit", "Power")),
                              
                              downloadButton("downloadData", "Download")
                              
                              ),
                            
                            mainPanel(
                              tabsetPanel(
                                tabPanel("Power", tableOutput("power")),
                                tabPanel("Model Fit", tableOutput("modelfit")),
                                tabPanel("Plot", plotlyOutput("plot"))
                              )
                            )
                          )
                        )
                 )




server= function(input, output) {
  
  source("./fittedModelDetermine.function.R")
  source("./gData.function.R")
  source("./test.function.R")
  source("./ModelFit.function.R")
  
  sim <- eventReactive(input$go, {
    
    nIteration = as.numeric(input$nIteration)
    
    cellSize = as.numeric(input$cellSize)
    clusterSize = as.numeric(input$clusterSize)
    groupSize = as.numeric(input$groupSize)
    nTransition = as.numeric(input$nTransition)
    b0 = as.numeric(input$b0)
    b4 = as.numeric(input$b4)
    rhoO = as.numeric(input$rhoO)
    rhoY = as.numeric(input$rhoY)
    rhoG = as.numeric(input$rhoG)
    
    SIGMA = as.numeric(input$SIGMA)
    SIGMA2 = as.numeric(input$SIGMA2)
    
    
    t = as.numeric(input$t)
    x1 = input$x1
    x2 = input$x2
    lincom = as.numeric(input$lincom)
    
    glm.model = fittedModelDetermine(input$model)
    
    secularTrend = as.numeric(input$secularTrend)
    randomTrt = as.numeric(input$randomTrt)
    
    
    
    result = data.frame(iteration=numeric(),
                        N=numeric(),
                        m=numeric(),
                        nCluster=numeric(),
                        nGroup=numeric(),
                        nTransit=numeric(),
                        beta0=numeric(), 
                        beta1=numeric(),
                        beta2=numeric(), 
                        beta3=numeric(), 
                        beta4=numeric(), 
                        rho13=numeric(), 
                        rho24=numeric(),
                        rho34=numeric(),
                        sigma=numeric(),
                        sigma2=numeric(),
                        coef_trt1=numeric(),
                        lci_trt1=numeric(),
                        uci_trt1=numeric(),
                        pvalue_trt1=numeric(),
                        coef_trt2=numeric(),
                        lci_trt2=numeric(),
                        uci_trt2=numeric(),
                        pvalue_trt2=numeric(),
                        coef_trt3=numeric(),
                        lci_trt3=numeric(),
                        uci_trt3=numeric(),
                        pvalue_trt3=numeric(),
                        coef_trt4=numeric(),
                        lci_trt4=numeric(),
                        uci_trt4=numeric(),
                        pvalue_trt4=numeric(),
                        coef_time=numeric(),
                        lci_time=numeric(),
                        uci_time=numeric(),
                        pvalue_time=numeric(),
                        test1_lci_x1=numeric(),
                        test1_uci_x1=numeric(),
                        test1_lci_x2=numeric(),
                        test1_uci_x2=numeric(),
                        test1_pvalue=numeric(),
                        test2_lci=numeric(),
                        test2_uci=numeric(),
                        test2_pvalue=numeric())
    
    
    for (val1 in 1:length(cellSize)){
      m = cellSize[val1]
      
      for (val2 in 1:length(clusterSize)){
        nCluster = clusterSize[val2]
        
        for (val3 in 1:length(groupSize)){
          nGroup = groupSize[val3]
          
          for (val4 in 1:length(nTransition)){
            nTransit = nTransition[val4]
            nPeriod = nCluster + 1 + nTransit
            
            for (val5 in 1:length(b0)){
              beta0 = b0[val5]
              
              for (val6 in 1:length(b4)){
                beta4 = b4[val6]
                
                for (val7 in 1:length(rhoO)){
                  rho13 = rhoY[val7]
                  
                  for (val8 in 1:length(rhoY)){
                    rho24 = rhoY[val8]
                    
                    for (val9 in 1:length(rhoG)){
                      rho34 = rhoG[val9]
                      
                      for (val10 in 1:length(SIGMA)){
                        sigma = SIGMA[val10]
                        
                        for (val11 in 1:length(SIGMA2)){
                          sigma2 = SIGMA2[val11]
                          if(randomTrt == 0 & is.na(sigma2) == F){
                            sigma2 = NA
                            warning('When randomTrt option is FALSE (default), the input of sigma2 is ignored.')
                          }
                          
                          set.seed(input$seed)
                          i=1
                          while (i <= nIteration) {
                            
                            
                            
                            
                            result_vec = ModelFit(glmFun=glm.model,
                                                  nCluster, nGroup, nTransit, m, 
                                                  beta0, beta4, rho34, rho24, rho13, t, 
                                                  sigma, sigma2,
                                                  secularTrend, randomTrt,
                                                  x1, x2, lincom)
                            
                            
                            if (!is.null(result_vec)){
                              
                              result[nrow(result)+1, ] <- append(i, result_vec)
                              i = i + 1
                              
                            }
                            
                            
                          } # i iteration
                        } # SIGMA2
                      } # SIGMA
                    } # rhoG
                  } # rhoY
                } # rhoO
              } # b4
            } # b0
          } # nTransit
        } # nGroup
      } # nCluster
    } # m
    return(result)
    })


  power_tb <- eventReactive(input$go, {

    result = sim()
    power_data = result %>%
      mutate(sig_trt1 = ifelse(pvalue_trt1 < 0.05, 1, 0),
             sig_trt2 = ifelse(pvalue_trt2 < 0.05, 1, 0),
             sig_trt3 = ifelse(pvalue_trt3 < 0.05, 1, 0),
             sig_trt4 = ifelse(pvalue_trt4 < 0.05, 1, 0),
             sig_test1 = ifelse(test1_pvalue < 0.05, 1, 0),
             sig_test2 = ifelse(test2_pvalue < 0.05, 1, 0)
      ) %>%
      group_by(
        m,
        nCluster,
        nGroup,
        nTransit,
        beta0,
        beta4,
        rho13,
        rho24,
        rho34,
        sigma,
        sigma2,
      ) %>%
      summarise(power_trt1 = mean(sig_trt1),
                power_trt2 = mean(sig_trt2),
                power_trt3 = mean(sig_trt3),
                power_trt4 = mean(sig_trt4),
                test1_power = mean(sig_test1),
                test2_power = mean(sig_test2)
      )
    return(power_data)
  })
  
  output$modelfit <- renderTable(
    sim(),
    options = list(scrollX = TRUE,
                   scrollY = 800))
  
  output$power <- renderTable(
    power_tb(),
    options = list(scrollX = TRUE,
                   scrollY = 800))



  datasetInput <- reactive({
    switch(input$dataset,
           "Model Fit" = sim(),
           "Power" = power_tb())
  })

  output$table <- renderTable({
    datasetInput()
  })

  output$downloadData <- downloadHandler(filename = function() {
    return(paste(input$dataset, '.csv', sep=''))

  }, content = function(file) {
    write.csv(datasetInput(), file)
  })
  
  
  output$plot <- renderPlotly({
    plot_ly(power_tb(), x = ~N, y = ~power_trt_gy, type = 'scatter', mode = 'lines')
  })

  
}




shinyApp(ui = ui, server = server)
