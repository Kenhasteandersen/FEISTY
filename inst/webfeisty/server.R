#===============================================================================
# Server for the Shiny web application. 
# You can run the application by clicking  the 'Run App' button above, OR:.
# use webFEISTY() to trigger it from anywhere (if the FEISTYR package is loaded)
#===============================================================================

server <- function(input, output) {
  
  observeEvent(input$Setup,if(input$Setup=="setupBasic"){ 
    shinyjs::hide("region")
    shinyjs::hide("bottom")
    shinyjs::hide("photic")
    shinyjs::hide("nSizeGroups")
    shinyjs::hide("etaMature")
    shinyjs::hide("dfpho")
    shinyjs::show("szprod")
    shinyjs::show("lzprod")
    shinyjs::show("bprodin")
    shinyjs::show("Tp")
    shinyjs::hide("Tm")
    shinyjs::show("Tb")
    shinyjs::show("sh_de")
    shinyjs::hide("F")
    shinyjs::hide("etaF")
  }else if (input$Setup=="setupBasic2") {
    shinyjs::hide("region")
    shinyjs::hide("bottom")
    shinyjs::hide("photic")
    shinyjs::hide("dfpho")
    shinyjs::show("sh_de")
    shinyjs::show("etaMature")
    shinyjs::show("nSizeGroups")
    shinyjs::show("szprod")
    shinyjs::show("lzprod")
    shinyjs::show("bprodin")
    shinyjs::show("Tp")
    shinyjs::hide("Tm")
    shinyjs::show("Tb")
    shinyjs::show("F")
    shinyjs::show("etaF")
  }else if (input$Setup=="setupVertical") {
    shinyjs::show("region")
    shinyjs::show("bottom")
    shinyjs::show("photic")
    shinyjs::hide("nSizeGroups")
    shinyjs::show("szprod")
    shinyjs::show("lzprod")
    shinyjs::hide("etaMature")
    shinyjs::show("dfpho")
    shinyjs::hide("bprodin")
    shinyjs::hide("Tp")
    shinyjs::hide("Tm")
    shinyjs::hide("Tb")
    shinyjs::hide("sh_de")
    shinyjs::hide("F")
    shinyjs::hide("etaF")
  }else if (input$Setup=="setupVertical2") {
    shinyjs::hide("region")
    shinyjs::show("bottom")
    shinyjs::show("photic")
    shinyjs::show("nSizeGroups")
    shinyjs::show("szprod")
    shinyjs::show("lzprod")
    shinyjs::show("etaMature")
    shinyjs::show("dfpho")
    shinyjs::hide("bprodin")
    shinyjs::show("Tp")
    shinyjs::show("Tm")
    shinyjs::show("Tb")
    shinyjs::hide("sh_de")
    shinyjs::show("F")
    shinyjs::show("etaF")
  }
  )
  
  sim <- eventReactive(
    c(
      input$szprod,input$lzprod,input$bprodin,input$USEdll,input$Setup,input$nSizeGroups,
      input$Tp,input$Tm,input$Tb,input$region, input$bottom,input$photic,input$etaMature,input$dfpho,
      input$sh_de,input$F,input$etaF
    ),
    {
      # setup simulation
      if (input$Setup == "setupBasic") {
        p = setupBasic(szprod = input$szprod, lzprod = input$lzprod, bprodin=input$bprodin,depth=input$sh_de,Tp=input$Tp,Tb=input$Tb)
        #setupini = c(input$szprod,input$lzprod,input$bprod,input$sh_de,input$Tp,input$Tb)
      }else if (input$Setup == "setupBasic2") {
        p = setupBasic2(szprod = input$szprod, lzprod = input$lzprod, bprodin=input$bprodin,depth=input$sh_de,Tp=input$Tp,Tb=input$Tb,
                        nStages =input$nSizeGroups, # Number of size groups
                        etaMature=input$etaMature,
                        F=input$F,
                        etaF=input$etaF)
        #setupini = c(input$szprod,input$lzprod,input$bprod,input$nSizeGroups,depth=input$sh_de,input$Tp,input$Tb,input$etaMature,input$F,input$etaF)
      }else if (input$Setup == "setupVertical") {
        p = setupVertical(szprod = input$szprod, lzprod = input$lzprod, dfpho=input$dfpho,
                          #nStages  = input$nSizeGroups, # Number of size groups
                          region   = as.integer(input$region),
                          depth    = input$bottom,
                          photic   = input$photic)
        #setupini = c(input$szprod,input$lzprod,input$bent,input$nSizeGroups,input$region, input$bottom, input$photic)
        
      }else if (input$Setup == "setupVertical2") {
        p = setupVertical2(szprod = input$szprod, lzprod = input$lzprod, dfpho=input$dfpho,
                           nStages  = input$nSizeGroups, # Number of size groups
                           Tp = input$Tp,
                           Tm = input$Tm,
                           Tb = input$Tb,
                           depth    = input$bottom,
                           photic   = input$photic,
                           shelfdepth=250,
                           visual=1.5,
                           etaMature=input$etaMature,
                           F=input$F,
                           etaF=input$etaF)
        #setupini = c(input$szprod,input$lzprod,input$bent,input$nSizeGroups,input$region,input$bottom,input$photic,input$etaMature,input$F,input$etaF)
        
      }
      
      # Simulate
      return( simulateFEISTY(bCust  = TRUE,
                             p      = p, 
                             tEnd   = 500,
                             tStep  = 1,
                             times  = seq(from=0, to=500, by=1), #to=tEnd but must give a number directly
                             yini   = p$u0,  
                             USEdll = input$USEdll,
                             Rmodel = derivativesFEISTYR) )
    })
  
  # Make plots
  output$plotSimulationShiny <- renderPlot( plotSimulationShiny(sim()) )
  output$plotNetwork <- renderPlot( plotNetwork(sim(), ref_b=5) )
  output$plotDiet <- renderPlot( plotDiet(sim()) )
}
