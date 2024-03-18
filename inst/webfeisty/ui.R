#===============================================================================
# User interface for the Shiny web application. 
# You can run the application by clicking  the 'Run App' button above, OR:.
# use webFEISTY() to trigger it from anywhere (if the FEISTYR package is loaded)
#===============================================================================

ui <- fluidPage(
  useShinyjs(),
    # Application title
    titlePanel("FEISTY"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
          
          h3('Setup:'),
          radioButtons("Setup", 
                       ("Setup"),
                       choices = list("SetupBasic" = "setupBasic", 
                                      "SetupBasic2" = "setupBasic2", 
                                      "SetupVertical" = "setupVertical",
                                      "SetupVertical2" = "setupVertical2"),
                       selected = "setupVertical2" ,inline = TRUE),
          radioButtons("sh_de", 
                       ("Region (shallow or deep)"),
                       choices = list("Continental shelf (< 200m)" = 100, "Deep water (>= 200m)" = 300),
                       selected = 100 ,inline = TRUE),
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
          
          hr(),
          h3('Production:'),
          
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
          sliderInput("bprodin",
                      "Benthos carrying capacity (g/m2):",
                      min = 1,
                      max = 50,
                      step = 1,
                      value = 5),
          sliderInput("dfpho",
                      "Detrital flux out of photic zone (g/m2):",
                      min = 100,
                      max = 1000,
                      step = 5,
                      value = 150),
          
          hr(),
          h3('Temperature:'),
          
          radioButtons("region", 
                       ("Region (Temp. profile)"),
                       choices = list("Default (the whole water column is 10 \u00B0C)" = 4,
                                      "Tropical" = 1, "Temperate" = 2, "Boreal"=3),
                       selected = 4 ,inline = TRUE),
          sliderInput("Tp",
                      "Surface temperature (top 100m, \u00B0C):",
                      min = 0,
                      max = 28,
                      step = 0.1,
                      value = 10),
          sliderInput("Tm",
                      "Mid-water temperature (500m - up to 1500m, \u00B0C):",
                      min = 0,
                      max = 28,
                      step = 0.1,
                      value = 10),          
          sliderInput("Tb",
                      "Bottom temperature (\u00B0C):",
                      min = 0,
                      max = 28,
                      step = 0.1,
                      value = 10),
          
          hr(),
          h3('Fishing parameters:'),
          sliderInput("F",
                      "Fishing mortality (1/yr):",
                      min = 0,
                      max = 3,
                      step = 0.1,
                      value = 0),
          sliderInput("etaF",
                      "fish size with 50% selectivity relative to asymptotic size:",
                      min = 0.002,
                      max = 0.3,
                      step = 0.01,
                      value = 0.05),
          
          hr(),
          h3('Numerical solver:'),
          
          sliderInput("nSizeGroups",
                      "Number of fish stages:",
                      min = 3,
                      max = 45,
                      step = 3,
                      value = 9),
          radioButtons("USEdll", 
                       ("Run by"),
                       choices = list("Fortran dll" = TRUE, "R" = FALSE),
                       selected = TRUE,inline = TRUE),
          sliderInput("etaMature",
                      "Mature size relative to asymptotic size:",
                      min = 0.002,
                      max = 0.3,
                      step = 0.002,
                      value = 0.25),
          
        ),

        # Show a plot of the generated distribution
        mainPanel(
          tabsetPanel(
            tabPanel("Rates", plotOutput(outputId = "plotSimulation", height="600px")),
            
            tabPanel("Network", plotOutput(outputId = "plotNetwork", height="600px")),
            
            tabPanel("Diet", plotOutput(outputId = "plotDiet", height="600px")))
        )
    )
)
