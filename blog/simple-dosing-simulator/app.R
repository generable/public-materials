#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#
source(here::here('blog', 'simple-dosing-simulator', 'app-functions.R'))
library(shiny)

addResourcePath('www', here::here('blog', 'simple-dosing-simulator', 'www'))
ui <- fluidPage(
  titlePanel("Dose–Risk Simulation"),
  tabsetPanel(
    # --- Tab 1: Simulation ---
    tabPanel("Simulation",
             sidebarLayout(
               sidebarPanel(
                 numericInput("stroke_rate",
                              "Stroke baseline risk: ",
                              min = 0, max = 1, step = 0.005, value = 0.08),
                 numericInput("bleed_rate",
                              "Bleed baseline risk: ",
                              min = 0, max = 1, step = 0.005, value = 0.03),
                 numericInput("bleed_weight",
                              "Bleed harmfulness: ",
                              min = 0, max = 10, step = 0.05, value = 1),
                 numericInput("stroke_weight",
                              "Stroke harmfulness: ",
                              min = 0, max = 10, step = 0.05, value = 1)
               ),
               mainPanel(
                 plotOutput("riskPlot"),
                 plotOutput("utilityPlot")
               )
             )
    ),
    
    # --- Tab 2: Overview (markdown help text) ---
    tabPanel("Simulation overview",
             includeHTML("www/about-simulation.html")
    )
  )
)

server <- function(input, output) {
  results <- reactive({sim_risks(weights = c(Bleed = input$bleed_weight, Stroke = input$stroke_weight), p0_bleed = input$bleed_rate, p0_stroke = input$stroke_rate)})
  output$riskPlot <- renderPlot({cowplot::plot_grid(results()$risks, results()$composite, nrow = 1)})
  output$utilityPlot <- renderPlot({results()$utility})
}

shinyApp(ui = ui, server = server)