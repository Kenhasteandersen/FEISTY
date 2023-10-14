#===============================================================================
# Server for the Shiny web application. 
# You can run the application by clicking  the 'Run App' button above, OR:.
# use webFeisty() to trigger it from anywhere (if the FeistyR package is loaded)
#===============================================================================

server <- function(input, output) {

  sim <- eventReactive(c(
    input$pprod,input$bprod,input$USEdll,input$Setup
  ),
  {
    # setup simulation
    if (input$Setup == 1) {
      p = setupBasic(pprod = input$pprod, bprod=input$bprod)
    }else if (input$Setup == 2) {
      p = setupBasic2(pprod = input$pprod, bprod=input$bprod, nStages=9)   
    }else if (input$Setup == 3) {
      p = setupVertical(pprod = input$pprod)
    }

    # Simulate
    return( simulateFeisty(cus    = FALSE,
                           setup  = 1,
                           setupini = c(input$pprod,input$pprod,input$bprod,10,8),# setupbasic(smzprod,lgzprod,bprod,Ts,Tb)
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
