#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#

# need this package to run : latex2exp
# install.packages("shiny")
# install.packages("latex2exp")

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
                        min = 1,
                        max = 500,
                        step = 5,
                        value = 100)
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("plotSimulation", height="600px")
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
