#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# button added for selection 
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
          
          
            radioButtons("USEdll", 
                         h5("ODE solved by"),
                         choices = list("Fortran dll" = TRUE, "R" = FALSE),
                         selected = TRUE),
          
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

  sim <- eventReactive(c(
    input$pprod,input$USEdll
  ),
  {
    # setup simulation
    p = setupBasic(pprod = input$pprod,bprod=5)
    
    # Simulate
    return( simulate(p, tEnd = 100, USEdll=input$USEdll) )
  })
  
  # Make plots
  output$plotSimulation <- renderPlot( plotSimulation(sim()) ) 
}
#
# Run the application 
#
shinyApp(ui = ui, server = server)
