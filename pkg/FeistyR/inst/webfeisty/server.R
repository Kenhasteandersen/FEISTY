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
    return( simulateFeisty(p, tEnd = 100, USEdll=input$USEdll, simpleOutput = TRUE) )
  })
  
  set_param <- eventReactive(c(
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
    return( p)
  })
  
  sim2 <- eventReactive(c(
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
    return( simulateFeisty(p, tEnd = 100, USEdll=input$USEdll, simpleOutput = FALSE) )
  })
  
  # Make plots
  #output$plotSimulation <- renderPlot( plotSimulation(sim()) ) 
 output$plotNetwork <- renderPlot( plot_network(set_param(), sim2()) )
}
