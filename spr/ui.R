
library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
    
    # Application title
    titlePanel("CSEO Female BRF biomass (lbs)"),
    
    
    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            sliderInput("wt",
                        "Avg fish wt (lbs):",
                        min = 1,
                        max = 3.5,
                        value = 1.75,
                        step = 0.05),
            sliderInput("wt_se",
                        "Avg wt se:",
                        min = 0,
                        max = 0.5,
                        value = 0.02255,
                        step = 0.02),
            sliderInput("M",
                        "Natural mortality (M):",
                        min = 0.12,
                        max = 0.25,
                        value = 0.19,
                        step = 0.01),
            sliderInput("Fc",
                        "Fishing mortality (F):",
                        min = 0.0,
                        max = 0.3,
                        value = 0.22,
                        step = 0.01),
            sliderInput("density",
                        "Density (#/km2):",
                        min = 1000,
                        max = 30000,
                        value = 3000,
                        step = 250),
            sliderInput("density_se",
                        "Density se:",
                        min = 0,
                        max = 1000,
                        value = 400,
                        step = 50),
            sliderInput("ratio",
                        "Proportion female:",
                        min = 0.30,
                        max = 0.6,
                        value = 0.5,
                        step = 0.01)
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
            tabsetPanel(
                type = "tabs", 
                tabPanel("Plot",plotOutput("distPlot")),
                tabPanel("SPR", plotOutput("sprPlot"),
                         plotOutput("sprPlot2"))
            ),
            hr(),
            print("Vertical line is the GIS estimated area of quality black rockfish habitat (1781 km2)")
        )
    )))
