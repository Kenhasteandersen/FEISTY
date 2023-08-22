#===============================================================================
# User interface for the Shiny web application. 
# You can run the application by clicking  the 'Run App' button above, OR:.
# use webFeisty() to trigger it from anywhere (if the FeistyR package is loaded)
#===============================================================================

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
            
            radioButtons("Setup", 
                         h5("Setup"),
                         choices = list("SetupBasic" = 1, "SetupBasic2" = 2, "SetupVertical (no bprod input)" = 3),
                         selected = 1),
          
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
                        value = 5)
        ),

        # Show a plot of the generated distribution
        mainPanel(
          
           tabsetPanel(
             tabPanel("Rates", plotOutput(outputId = "plotSimulation", height="600px")),
             
             tabPanel("Network", plotOutput(outputId = "plotNetwork", height="600px")),
           
           tabPanel("Diet", plotOutput(outputId = "plotDiet", height="600px")))
    )
    ))
