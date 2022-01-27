#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("FEISTY"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput("pprod",
                        "Primary prod. (1/yr):",
                        min = 0,
                        max = 5,
                        step = 0.1,
                        value = 3)
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("plotSSB")
        )
    )
)

# Define server logic required to draw a histogram
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
  output$plotSSB <- renderPlot( plotSSBtime(sim()) ) 
}

# Run the application 
shinyApp(ui = ui, server = server)
