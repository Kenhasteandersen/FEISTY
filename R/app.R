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
                         ("Run by"),
                         choices = list("Fortran dll" = TRUE, "R" = FALSE),
                         selected = TRUE,inline = TRUE),
            
            radioButtons("Setup", 
                         ("Setup"),
                         choices = list("SetupBasic" = 1, "SetupBasic2" = 2, "SetupVertical (no bprod input)" = 3),
                         selected = 1 ,inline = TRUE),
            radioButtons("region", 
                         ("Region (temp. profile, only for vertical ver.)"),
                         choices = list("Default (10 Celsius)" = 4, "Tropical" = 1, "Temperate" = 2, "Boreal"=3),
                         selected = 4 ,inline = TRUE),
            
            sliderInput("pprod",
                        "Primary prod. (1/yr):",
                        min = 1,
                        max = 500,
                        step = 5,
                        value = 100),
            
            sliderInput("bprod",
                        "Small benthos prod. (1/yr):",
                        min = 1,
                        max = 50,
                        step = 1,
                        value = 5),
            sliderInput("nSizeGroups",
                        "Fish stage number (not for setupbasic):",
                        min = 3,
                        max = 45,
                        step = 3,
                        value = 3),
            sliderInput("temp",
                        "Temperature (not for vertical ver. ):",
                        min = 0,
                        max = 28,
                        step = 0.1,
                        value = 10)
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
    input$pprod,input$bprod,input$USEdll,input$Setup,input$nSizeGroups,input$temp,input$region
  ),
  {
    # setup simulation
    if (input$Setup == 1) {
      p = setupBasic(pprod = input$pprod,
                     bprod = input$bprod,
                     temp = input$temp)
    } else if (input$Setup == 2) {
      p = setupBasic2(
        pprod = input$pprod,
        bprod = input$bprod,
        nSizeGroups = input$nSizeGroups,
        temp = input$temp
      )
    } else if (input$Setup == 3) {
      p = setupVertical(pprod = input$pprod,
                        nSizeGroups = input$nSizeGroups,
                        region =as.integer(input$region))
    }

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
