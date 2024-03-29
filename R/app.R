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
library(shinyjs)
library(latex2exp)
source("Plots.R")

#
# Define UI
#
ui <- fluidPage(
    useShinyjs(),
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
                         ("Region (Temp. profile)"),
                         choices = list("Default (the whole water column is 10 \u00B0C)" = 4,
                                        "Tropical" = 1, "Temperate" = 2, "Boreal"=3),
                         selected = 4 ,inline = TRUE),
            sliderInput("bottom",
                        "Seafloor depth (m):",
                        min = 100,
                        max = 5500,
                        step = 50,
                        value = 1500),
            sliderInput("photic",
                        "Euphotic zone depth (m):",
                        min = 50,
                        max = 700,
                        step = 10,
                        value = 150),
            sliderInput("szprod",
                        "Small zoop. carrying capacity (g/m2):",
                        min = 1,
                        max = 500,
                        step = 1,
                        value = 100),
            sliderInput("lzprod",
                        "Large zoop. carrying capacity (g/m2):",
                        min = 1,
                        max = 500,
                        step = 1,
                        value = 100),
            sliderInput("bprod",
                        "Small benthos carrying capacity (g/m2):",
                        min = 1,
                        max = 50,
                        step = 1,
                        value = 5),
           sliderInput("bent",
                       "Detrital flux out of photic zone (g/m2):",
                       min = 100,
                       max = 1000,
                       step = 5,
                       value = 150),
            sliderInput("nSizeGroups",
                        "Fish stage number:",
                        min = 3,
                        max = 45,
                        step = 3,
                        value = 6),
           sliderInput("etaMature",
                       "Mature size relative to asymptotic size:",
                       min = 0.002,
                       max = 0.3,
                       step = 0.002,
                       value = 0.25),
            sliderInput("temps",
                        "Surface temperature (top 100m, \u00B0C):",
                        min = 0,
                        max = 28,
                        step = 0.1,
                        value = 10),
            sliderInput("tempb",
                        "Bottom temperature (\u00B0C):",
                        min = 0,
                        max = 28,
                        step = 0.1,
                        value = 8)
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
    }else if (input$Setup==2) {
      hide("region")
      hide("bottom")
      hide("photic")
      hide("bent")
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
      hide("temps")
      hide("tempb") 
      }
    )

  sim <- eventReactive(c(
    input$szprod,input$lzprod,input$bprod,input$USEdll,input$Setup,input$nSizeGroups,
    input$temps,input$tempb,input$region, input$bottom,input$photic,input$etaMature,input$bent
  ),
  {
    # setup simulation
    if (input$Setup == 1) {
      p = setupBasic(szprod = input$szprod,
                     lzprod = input$lzprod,
                     bprod = input$bprod,
                     temps = input$temps,
                     tempb = input$tempb)
    } else if (input$Setup == 2) {
      p = setupBasic2(
        szprod = input$szprod,
        lzprod = input$lzprod,
        bprod = input$bprod,
        nSizeGroups = input$nSizeGroups,
        temps = input$temps,
        tempb = input$tempb,
        etaMature= input$etaMature
      )
    } else if (input$Setup == 3) {
      p = setupVertical(szprod = input$szprod,
                        lzprod = input$lzprod,
                        nSizeGroups = input$nSizeGroups,
                        region =as.integer(input$region),
                        bottom=input$bottom,
                        photic=input$photic,
                        etaMature= input$etaMature,
                        bent=input$bent)
    }

    # Simulate
    return( simulate(p, tEnd = 200, USEdll=input$USEdll) )
  })
  
  # Make plots
  output$plotSimulation <- renderPlot( plotSimulation(sim()) ) 
}
#
# Run the application 
#
shinyApp(ui = ui, server = server)
