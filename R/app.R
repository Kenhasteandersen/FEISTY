#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#

library(shiny)
source("Plots.R")

#
# Define UI
#
ui <- fluidPage(

    # Application title
    titlePanel("FEISTY"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput("pprod",
                        "Primary prod. (1/yr):",
                        min = 0,
                        max = 500,
                        step = 5,
                        value = 100)
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("plotSimulation")
        )
    )
)
#
# Define server logic
#
server <- function(input, output) {

  sim <- eventReactive({
    input$pprod
  },
  {
    # setup simulation
    p = setupBasic(depth = 500, pprod = input$pprod)
    
    # Simulate
    return( simulate(p) )
  })
  
  # Make plots
  output$plotSimulation <- renderPlot( plotSimulation(sim()) ) 
}
#
# Run the application 
#
shinyApp(ui = ui, server = server)
