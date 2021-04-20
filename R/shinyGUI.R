#' A shiny app for the package GFDsurv
#'
#' This function provides a shiny app for calculating
#' CASANOVA, medSANOVA and copSANOVA test statistics and respective p-values.
#'
#'
#'
#'
#' @aliases GFDsurvGUI
#'
#'
#' @import shiny
#' @import utils
#' @importFrom shinyjs useShinyjs hide hidden show
#' @importFrom shinyWidgets sliderTextInput
#' @importFrom tippy tippy_this
#' @importFrom shinythemes shinytheme
#'
#' @export

GFDsurvGUI <- function() {
  requireNamespace("shiny", quietly = TRUE)

  if (!("package:shiny" %in% search())) {
    attachNamespace("shiny")
  }
  requireNamespace("shinyWidgets", quietly = TRUE)
  if (!("package:shinyWidgets" %in% search())) {
    attachNamespace("shinyWidgets")
  }
  requireNamespace("tippy", quietly = TRUE)

  if (!("package:tippy" %in% search())) {
    attachNamespace("tippy")
  }


      ui <- fluidPage(theme = shinythemes::shinytheme("cerulean"),
                      shinyjs::useShinyjs(),
                      titlePanel("Tests for GFDsurv"),
                      sidebarLayout(
                        sidebarPanel(
                          splitLayout(
                            fileInput("infile", "Choose CSV File",
                                      accept = c(
                                        "text/csv",
                                        "text/comma-separated-values,text/plain",
                                        ".csv")),
                            checkboxInput("header", "Header", TRUE),
                            selectInput("sep","Seperator in csv", c(",",
                                                                    ";",
                                                                    ".",
                                                                    "|"))
                          ),

                          tags$head(tags$style(HTML("
                              .shiny-split-layout > div {
                                overflow: visible;
                              }
                              "))), #for selectinput in splitlayout with full dropdown view
                          tags$style(HTML("
                                 input[type=number] {
                                                              -moz-appearance:textfield;
                                                    }
                                  input[type=number]::{
                                                  -moz-appearance:textfield;
                                                    }
                        input[type=number]::-webkit-outer-spin-button,
                        input[type=number]::-webkit-inner-spin-button {
                        -webkit-appearance: none;
                        margin: 0;
                        }
                        ")),


                        h3(id="titleLoadData","Load dataset first!", style = "color:red"),

                        shinyjs::hidden(
                          selectInput("Method", "Select Testing Method:",
                                      c("CASANOVA: Cumulative Aalen survival analyis-of-variance" = "casanova",
                                        "MedSANOVA: Median survival analyis-of-variance"= "medSANOVA",
                                        "CopSANOVA: Concordance probability survival analyis-of-variance"="copSANOVA"))
                        ),

                          splitLayout(cellWidths = c("40%","40%"),
                            uiOutput(outputId = 'dynamicInput'),
                            uiOutput(outputId = 'dynamicInput2')

                          ),

                        splitLayout(cellWidths = c("50%","5%","45%"),
                          shinyjs::hidden(
                            textInput("formula", "Formula ", "timeFactor ~ FactorA * FactorB")
                          ),
                          shinyjs::hidden(
                            actionButton("infoButton", "", icon = icon("info-circle"))
                          ),

                            tippy::tippy_this("infoButton", "Interaction effects need to be specified! Example:<br><br>
                                              - 1 Factor: <br>
                                              timefactor ~ factorA <br><br>
                                              - 2 Factors without interactions:
                                              <br> timefactor ~ factorA + factorB
                                              <br>  <br>- 2 Factors with interactions:
                                              <br> timefactor ~ factorA + factorB + factorA:factorB
                                              <br> or
                                              <br> timefactor ~ factorA * factorB "
                                       , placement = "right")
                        ),

                        splitLayout(cellWidths = c("60%"),
                          shinyjs::hidden(
                            checkboxInput("nested", "Are the levels of nested factors the same for each level main factor?", FALSE)
                          )
                        ),


                        shinyjs::hidden(
                          h5(id="titleWeights",strong("Which weight functions w should be combined?"), style = "color:grey")
                        ),


                          splitLayout(cellWidths = c("35%","60%","5%"),
                                      shinyjs::hidden(
                            checkboxGroupInput("Weights", "Pre-specified weights",selected = c("crossing","proportional"),
                                      choiceNames = list("Crossing", "Proportional"),
                                      choiceValues = list("crossing", "proportional"))
                                      ),

                            shinyjs::hidden(
                            selectInput("weights1","Specify the exponents (r,g) of weights  w(x) = x^r(1-x)^g ",
                                        paste0("(",expand.grid(0:10,0:10)[,1],",",expand.grid(0:10,0:10)[,2],")"),
                                        multiple=TRUE,selectize = TRUE)
                            )

                          ),


                          splitLayout(cellWidths = c("30%","70%"),
                                      shinyjs::hidden(
                            radioButtons("variante", "Variance estimation based on",
                                       c("one.sided"= "onesided",
                                         "two.sided" = "twosided"), inline = TRUE)
                                      ),
                                  shinyjs::hidden(
                            numericInput("var_level", "confidence interval with level",
                                                   min = 0, max = 1,
                                                   value = 0.9, width = "20%")
                                  )
                          ),



                          splitLayout(cellWidths = c("15%"),
                            shinyjs::hidden(
                            numericInput("nperm", "Number of permutations", value = 1999)
                            )
                          ),

                          splitLayout(cellWidths = c("20%","10%","20%","10%","30%"),
                                      shinyjs::hidden(
                          shinyWidgets::sliderTextInput(
                              inputId = "sliderBoot",
                              label = "Bootstrap type:",
                              choices = c("wild","weird"),
                              selected = "wild")
                                      ),
                            shinyjs::hidden(
                              selectInput("Platz1", "Distribution:",
                                          c("Poisson"="pois","Normal"="norm"))
                            ),
                            shinyjs::hidden(
                            selectInput("methodBoot", "Distribution:",
                                        c("Poisson"="pois","Normal"="norm"))
                            ),
                            shinyjs::hidden(
                              selectInput("Platz2", "Distribution:",
                                          c("Poisson"="pois","Normal"="norm"))
                            ),
                            shinyjs::hidden(
                            checkboxInput("correction", "correction for liberal test", TRUE)
                            )

                          ),

                          splitLayout(cellWidths = c("15%","85%"),
                                      shinyjs::hidden(
                          numericInput("nboot", "Number of bootstrap iterations", value = 99)
                                      )
                          ),

                        shinyjs::hidden(
                          checkboxInput("plots", "Plot the surival curves", TRUE)
                        ),

                        splitLayout(cellWidths = c("15%","85%"),
                            shinyjs::hidden(
                            numericInput("tau", "Endpoint tau of the relevant time window [0,tau]",NULL)
                            )),

                        shinyjs::hidden(
                          actionButton("tau_suggest", "Specify tau automatically", class = "btn-primary")
                        ),

                        shinyjs::hidden(
                          actionButton("process", "Calculate", class = "btn-primary")
                        )
                        , width = 6
                        ),




                        mainPanel(

                          verbatimTextOutput("result"),
                          plotOutput("result_plot"),
                          plotOutput("result_plot2"),
                          width = 6

                        )
                      )
      )


       server <- function(input, output,session) {

         datasetInput <- reactive({

           req(input$infile)

           if (is.null(input$infile))
             return(NULL)
           read.csv(input$infile$datapath, header = input$header, sep = as.character(input$sep))
         })


         observeEvent(input$infile, {

           if(is.null(input$infile)){
             shinyjs::hide(id = "Method")
             shinyjs::hide(id = "weights1")
             shinyjs::hide(id = "Weights")
             shinyjs::hide(id = "formula")
             shinyjs::hide(id = "infoButton")
             shinyjs::hide(id = "variante")
             shinyjs::hide(id = "var_level")
             shinyjs::hide(id = "nperm")
             shinyjs::hide(id = "sliderBoot")
             shinyjs::hide(id = "methodBoot")
             shinyjs::hide(id = "correction")
             shinyjs::hide(id = "nboot")
             shinyjs::hide(id = "tau")
             shinyjs::hide(id = "tau_suggest")
             shinyjs::hide(id = "process")
             shinyjs::hide(id = "titleWeights")
             shinyjs::hide(id = "Platz1")
             shinyjs::hide(id = "Platz2")
             shinyjs::hide(id = "nested")
             shinyjs::hide(id = "plots")



           }else{
             shinyjs::show(id = "Method")
             shinyjs::show(id = "formula")
             shinyjs::show(id = "infoButton")
             shinyjs::show(id = "nested")
             shinyjs::show(id = "process")
             shinyjs::show(id = "plots")
             shinyjs::hide(id = "titleLoadData")


             observeEvent(input$Method, {

               if (input$Method != "casanova") {

                 shinyjs::hide(id = "titleWeights")
                 shinyjs::hide("Weights")
                 shinyjs::hide("weights1")

               }


               if (input$Method == "casanova") {

                 shinyjs::show(id = "titleWeights")
                 shinyjs::show("Weights")
                 shinyjs::show("weights1")

               }


               if (input$Method != "medSANOVA") {

                 shinyjs::hide("variante")
                 shinyjs::hide("var_level")

               }

               if (input$Method == "medSANOVA") {

                 shinyjs::show("variante")
                 shinyjs::show("var_level")

               }

               if (input$Method != "copSANOVA") {
                 # data  <- as.data.frame(datasetInput())
                 shinyjs::hide("methodBoot")
                 shinyjs::hide("nboot")
                 shinyjs::hide("sliderBoot")
                 shinyjs::hide("correction")
                 shinyjs::hide("tau")
                 shinyjs::hide("tau_suggest")
                 shinyjs::show("nperm")

               }
               if (input$Method == "copSANOVA") {
                 # data  <- as.data.frame(datasetInput())
                 shinyjs::show("methodBoot")
                 shinyjs::show("nboot")
                 shinyjs::show("sliderBoot")
                 shinyjs::show("correction")
                 shinyjs::show("tau")
                 shinyjs::show("tau_suggest")
                 shinyjs::hide("nperm")

               }




             })## observeevent




             observeEvent(input$sliderBoot, {

               if (input$sliderBoot != "weird" && input$Method == "copSANOVA") {
                 # data  <- as.data.frame(datasetInput())
                 shinyjs::show("methodBoot")

               }
               if (input$sliderBoot == "weird") {
                 # data  <- as.data.frame(datasetInput())
                 shinyjs::hide("methodBoot")

               }


             })

           }
         })












         values <- reactiveValues()

         output$dynamicInput <- renderUI({

           # This input exists if the `static`
           # one is equal to `A` only
           if (input$Method == "casanova" || input$Method == "medSANOVA"|| input$Method == "copSANOVA") {
             selectInput(inputId = 'dynamic',
                         label = "Name of event variable",
                         choices = colnames(datasetInput()))
           } else {
             return(NULL)
           }

         })

         ## this bit fixes the issue
         observe({
           if (input$Method == "casanova" || input$Method == "medSANOVA"|| input$Method == "copSANOVA") {
             values$dyn <- input$dynamic
           } else {
             values$dyn <- NULL
           }
         })

         values2 <- reactiveValues()

         output$dynamicInput2 <- renderUI({
           if (input$Method == "casanova" || input$Method == "medSANOVA"|| input$Method == "copSANOVA") {

             # This input exists if the `static`
             # one is equal to `A` only
             selectInput(inputId = 'dynamic2',
                         label = "Label of censored variable",
                         choices = unique(datasetInput()[,values$dyn]))
           } else {
             return(NULL)}

         })
         ## this bit fixes the issue

         observe({
           if (input$Method == "casanova" || input$Method == "medSANOVA"|| input$Method == "copSANOVA") {
             values2$dyn <- input$dynamic2
           } else {
             values2$dyn <- NULL
           }
         })









         observeEvent(input$tau_suggest, {
           if (input$formula == "~ + ") {

             output$result <- renderPrint({
               "'formula' missing or invalid"
             })
           }
           data <- as.data.frame(datasetInput())
           updateNumericInput(session, "tau", value = copsanova_tau(formula= isolate(input$formula),
                                                                    event = input$dynamic,
                                                                    data = isolate(data)))

         }
         )


         observeEvent(input$process, {

          if (input$formula == "timeFactor ~ FactorA*FactorB") {

            output$result <- renderPrint({
              "'formula' missing or invalid"
            })

          } else {

            if(length(unique(as.data.frame(datasetInput())[,input$dynamic])) != 2){
              output$result <- renderPrint({
                "ERROR: more or less then two censcoring types"
              })

            } else {

            if (input$Method == "casanova" ){
            data <- as.data.frame(datasetInput())
            event <- data[,input$dynamic]
            data[,input$dynamic] <- ifelse(event == input$dynamic2,1,0)

            rg <- list()
            givenWeights <- isolate(input$Weights)
            if("proportional" %in% givenWeights){
              rg[[length(rg)+1]] <- c(0,0)
            }
            if("crossing" %in% givenWeights){
              crossing <- TRUE
            } else {
              crossing <- FALSE
            }

            inputWeights <- isolate(input$weights1)

            if(is.null(inputWeights)){}else{
              expand.grid(0:10,0:10)
              kombiAll <- paste0("(",expand.grid(0:10,0:10)[,1],",",expand.grid(0:10,0:10)[,2],")")
              kombi <- expand.grid(0:10,0:10)[which(kombiAll %in% inputWeights),]

              for (i in 1:(dim(kombi)[1])){
                rg[[length(rg)+1]] <- as.numeric(kombi[i,])
              }
            }
            if(length(rg)==0){
              rg = list(c(0,0))
            }

            showModal(modalDialog("Calculating!"))
            output_cas <- casanova(formula= isolate(input$formula),
                     event = input$dynamic,
                     data = isolate(data),
                     nperm = isolate(input$nperm),
                     cross = crossing,
                     nested.levels.unique = isolate(input$nested),
                     rg = isolate(rg))
            removeModal()

            output$result <- renderPrint({
              output_cas
            })


            if(input$plots){
                if(length(all.vars(as.formula(input$formula)[[3]]))!=2){
                  output$result_plot <- renderPlot({
                    plot(output_cas,  direction = "vertical")
                  })
                  output$result_plot2 <- renderPlot({

                  })
                } else{
                  output$result_plot <- renderPlot({
                    plot(output_cas, by.group = TRUE, nr.group = 1)
                  })
                  output$result_plot2 <- renderPlot({
                    plot(output_cas,  by.group = TRUE, nr.group = 2)
                  })
                }
              }
            }

            if (input$Method == "medSANOVA" ){

              data <- as.data.frame(datasetInput())
              event <- data[,input$dynamic]
              data[,input$dynamic] <- ifelse(event == input$dynamic2,1,0)

              showModal(modalDialog("Calculating!"))

              output_med <-  medsanova(formula= isolate(input$formula),
                         event = input$dynamic,
                         data = isolate(data),
                         var_method = isolate(input$variante),
                         nperm = isolate(input$nperm),
                         nested.levels.unique = isolate(input$nested)
                         )
              removeModal()

              output$result <- renderPrint({
                output_med
              })

              if(input$plots){
                if(length(all.vars(as.formula(input$formula)[[3]]))!=2){
                  output$result_plot <- renderPlot({
                    plot(output_med,  direction = "vertical")
                  })
                  output$result_plot2 <- renderPlot({

                  })
                } else{
                  output$result_plot <- renderPlot({
                    plot(output_med, by.group = TRUE, nr.group = 1)
                  })
                  output$result_plot2 <- renderPlot({
                    plot(output_med,  by.group = TRUE, nr.group = 2)
                  })
                }
              }



            }
            if (input$Method == "copSANOVA" ){
              data <- as.data.frame(datasetInput())
              event <- data[,input$dynamic]
              data[,input$dynamic] <- ifelse(event == input$dynamic2,1,0)

              bootstrapMethod <- isolate(input$sliderBoot)
              bootstrapDis <- isolate(input$methodBoot)
              correction <- isolate(input$correction)
             if(bootstrapMethod == "wild"){
               if(bootstrapDis =="pois"){
                 if(correction == TRUE){
                   weights <- "corrLibPois"
                 } else {
                   weights <- "pois"
                 }
               } else {
                 if(correction == TRUE){
                  weights <- "corrLibNorm"
                }else{
                  weights <- "norm"
               }
               }
             } else {
               if(correction == TRUE){
                 weights <- "corrLibWeird"
               }else{
                 weights <- "weird"
               }

             }
              showModal(modalDialog("Calculating!"))

              output_cop <- copsanova(formula= isolate(input$formula),
                          event = input$dynamic,
                          data = isolate(data),
                          BSiter = isolate(input$nboot),
                          weights = weights,
                          tau = isolate(input$tau),
                          nested.levels.unique = isolate(input$nested)
                )
              removeModal()

              output$result <- renderPrint({
                output_cop
              })

              if(input$plots){
                if(length(all.vars(as.formula(input$formula)[[3]]))!=2){
                  output$result_plot <- renderPlot({
                    plot(output_cop,  direction = "vertical")
                  })
                  output$result_plot2 <- renderPlot({

                  })
                } else{
                  output$result_plot <- renderPlot({
                    plot(output_cop, by.group = TRUE, nr.group = 1)
                  })
                  output$result_plot2 <- renderPlot({
                    plot(output_cop,  by.group = TRUE, nr.group = 2)
                  })
                }
              }


            }

          }

          }
         }
         ) #end of observeEvent(input$process

      }


    shinyApp(ui = ui, server = server)

}

