#===============================================================================
# Server for the Shiny web application. 
# You can run the application by clicking  the 'Run App' button above, OR:.
# use webFeisty() to trigger it from anywhere (if the FeistyR package is loaded)
#===============================================================================

server <- function(input, output) {
  
  observeEvent(input$Setup,if(input$Setup==1){ 
    hide("region")
    hide("bottom")
    hide("photic")
    hide("nSizeGroups")
    hide("etaMature")
    hide("bent")
    show("szprod")
    show("lzprod")
    show("bprod")
    show("temps")
    show("tempb")
    show("sh_de")
  }else if (input$Setup==2) {
    hide("region")
    hide("bottom")
    hide("photic")
    hide("bent")
    show("sh_de")
    show("etaMature")
    show("nSizeGroups")
    show("szprod")
    show("lzprod")
    show("bprod")
    show("temps")
    show("tempb")
  }else if (input$Setup==3) {
    show("region")
    show("bottom")
    show("photic")
    show("nSizeGroups")
    show("szprod")
    show("lzprod")
    show("etaMature")
    show("bent")
    hide("bprod")
    hide("Tp")
    hide("Tb")
    hide("sh_de")
  }
  )
  
  sim <- eventReactive(c(
    input$szprod,input$lzprod,input$bprod,input$USEdll,input$Setup,input$nSizeGroups,
    input$Tp,input$Tb,input$region, input$bottom,input$photic,input$etaMature,input$bent,
    input$sh_de
  ),
  {
    # setup simulation
    if (input$Setup == 1) {
      p = setupBasic(szprod = input$szprod, lzprod = input$lzprod, bprod=input$bprod,depth=input$sh_de,Tp=input$Tp,Tb=input$Tb)
      setupini = c(input$szprod,input$lzprod,input$bprod,input$sh_de,input$Tp,input$Tb)
    }else if (input$Setup == 2) {
      p = setupBasic2(szprod = input$szprod, lzprod = input$lzprod, bprod=input$bprod,depth=input$sh_de,Tp=input$Tp,Tb=input$Tb,
                      nStages =input$nSizeGroups, # Number of size groups
                      etaMature=input$etaMature)
      setupini = c(input$szprod,input$lzprod,input$bprod,input$nSizeGroups,depth=input$sh_de,input$Tp,input$Tb,input$etaMature)
    }else if (input$Setup == 3) {
      p = setupVertical(szprod = input$szprod, lzprod = input$lzprod, bent=input$bent,
                      nStages  = input$nSizeGroups, # Number of size groups
                      region   = as.integer(input$region),
                      depth    = input$bottom,
                      photic   = input$photic,
                      mesopelagic=250,
                      visual=1.5,
                      etaMature=input$etaMature)
      setupini = c(input$szprod,input$lzprod,input$bent,input$nSizeGroups,input$region,input$bottom,input$photic,input$etaMature)
      
    }

    # Simulate
    return( simulateFeisty(bUseRDerivative    = FALSE,
                           setup  = input$Setup,
                           setupini = setupini,
                           p, 
                           tEnd   = 100,
                           times  = seq(from=0, to=100, by=1), #to=tEnd but must give a number directly
                           yini   = p$u0,  
                           USEdll = input$USEdll,
                           Rmodel = derivativesFeistyR,
                           simpleOutput = TRUE) )
  })
  
  # Make plots
  output$plotSimulation <- renderPlot( plotSimulation(sim()) ) 
}
