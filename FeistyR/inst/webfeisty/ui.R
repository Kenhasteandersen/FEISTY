#===============================================================================
# User interface for the Shiny web application. 
# You can run the application by clicking  the 'Run App' button above, OR:.
# use webFeisty() to trigger it from anywhere (if the FeistyR package is loaded)
#===============================================================================

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
          sliderInput("Tp",
                      "Surface temperature (top 100m, \u00B0C):",
                      min = 0,
                      max = 28,
                      step = 0.1,
                      value = 10),
          sliderInput("Tb",
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
