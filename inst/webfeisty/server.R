#===============================================================================
# Server for the Shiny web application. 
# You can run the application by clicking  the 'Run App' button above, OR:.
# use webFEISTY() to trigger it from anywhere (if the FEISTYR package is loaded)
#===============================================================================

server <- function(input, output) {
  
  observeEvent(input$Setup,if(input$Setup=="setupBasic"){ 
    hide("region")
    hide("bottom")
    hide("photic")
    hide("nSizeGroups")
    hide("etaMature")
    hide("dfpho")
    show("szprod")
    show("lzprod")
    show("bprodin")
    show("Tp")
    hide("Tm")
    show("Tb")
    show("sh_de")
    hide("F")
    hide("etaF")
  }else if (input$Setup=="setupBasic2") {
    hide("region")
    hide("bottom")
    hide("photic")
    hide("dfpho")
    show("sh_de")
    show("etaMature")
    show("nSizeGroups")
    show("szprod")
    show("lzprod")
    show("bprodin")
    show("Tp")
    hide("Tm")
    show("Tb")
    show("F")
    show("etaF")
  }else if (input$Setup=="setupVertical") {
    show("region")
    show("bottom")
    show("photic")
    hide("nSizeGroups")
    show("szprod")
    show("lzprod")
    hide("etaMature")
    show("dfpho")
    hide("bprodin")
    hide("Tp")
    hide("Tm")
    hide("Tb")
    hide("sh_de")
    hide("F")
    hide("etaF")
  }else if (input$Setup=="setupVertical2") {
    show("region")
    show("bottom")
    show("photic")
    show("nSizeGroups")
    show("szprod")
    show("lzprod")
    show("etaMature")
    show("dfpho")
    hide("bprodin")
    show("Tp")
    show("Tm")
    show("Tb")
    hide("sh_de")
    show("F")
    show("etaF")
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
                        mesopelagic=250,
                        visual=1.5,
                        etaMature=input$etaMature,
                        F=input$F,
                        etaF=input$etaF)
      #setupini = c(input$szprod,input$lzprod,input$bent,input$nSizeGroups,input$region,input$bottom,input$photic,input$etaMature,input$F,input$etaF)
      
    }

    # Simulate
    return( simulateFEISTY(bCust    = FALSE,
                           p, 
                           tEnd   = 100,
                           tStep  = 0.1,
                           times  = seq(from=0, to=100, by=0.1), #to=tEnd but must give a number directly
                           yini   = p$u0,  
                           USEdll = input$USEdll,
                           Rmodel = derivativesFEISTYR) )
  })
  
  # Make plots
  output$plotSimulation <- renderPlot( plotSimulation(sim()) )
  output$plotNetwork <- renderPlot( plotNetwork(sim()) )
  output$plotDiet <- renderPlot( plotDiet(sim()) )
}
