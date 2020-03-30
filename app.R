# load libraries
library(shiny); library(ggplot2); library(viridis); library(mgcv); 
library(boot); library(magrittr); library(tidymv)

#setwd("~/Dropbox/hzam_shiny/")

# define UI
ui <- fluidPage(

    # Application title
    titlePanel("HZAM: simulating hybrid zones with assortative mating"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            actionButton("go","Run Simulation",width="100%"),
            br(),br(),
            numericInput("max_generations","Number of Generations",10,min=1,max=500),
            numericInput("K_half","Per-species carrying capacity",800,min=1,max=2000),
            selectInput('mating_trait_dominance', 'Mating Trait Expression:', c("Dominant / Partially Dominant"=1,
                                                                                 "Additive"=2)),
            sliderInput("dom_coefficient","Partial Dominance Coefficient",value=0.75,min=0,max=1),
            div(helpText("Indicates proportion of homozygous dominant phenotype expressed by heterozygotes"),style="font-size:75%"),
            sliderInput("hybrid_fitness","Hybrid Fitness",value=0.75,min=0,max=1),
            sliderInput("pref_ratio","Assortative Mating Preference Ratio",value=0.75,min=0.00001,max=0.99999),
            div(helpText("Indicates ratio of heterospecific : homospecific matings"),style="font-size:75%"),
            sliderInput("male_trait_loci","Number of Trait Loci",value=1,min=1,max=5),
            sliderInput("neutral_loci","Number of Neutral Loci",value=3,min=1,max=10),
            sliderInput("meandispersal", "Average Dispersal Radius",value=0.01,min=0.0001,max=0.25),
            sliderInput("growth_rate", "Population Growth Rate",value=1.05,min=1,max=1.5),
            checkboxInput("fit_neutral","Fit Neutral Cline",value = F),
            checkboxInput("fit_trait","Fit Trait Cline",value = F)
        ),
            mainPanel(#style="position:fixed;margin-left:32vw;",
                plotOutput("plot"),
                div(style="height:5vh;"),
                tableOutput("table"),
                div(style="height:5vh;"),
                #img(src="melidectes.jpeg"),
                div(style="height:5vh;"),
                h4("About HZAM"),
                div(p("HZAM, originally written by", a(href="https://www.zoology.ubc.ca/~irwin/irwinlab/",'Darren Irwin', .noWS = "outside"), 
                "and adapted for Shiny by", a(href="https://elinck.org/",'Ethan Linck,', .noWS = "outside"), "is a program to
                simulate the influence of assortative mating in a one-dimensional hybrid zone in concert with the 
                      effects of underdominance, variation in dispersal, and other parameters. Code is available ",
                      a(href="http://github.com/elinck/hzam_shiny",'here.', .noWS = "outside"), 
                      "If you found this program useful, please cite the following publication:",
                a(href="https://doi.org/10.1086/708529",
                         'Irwin, D.E. In press. Assortative mating in hybrid zones is remarkably ineffective in promoting speciation. 
                         American Naturalist.', .noWS = c("after-begin", "before-end"))))
                )
    )
)


# define server
server <- function(input, output) {
    
    source("./dev.R")
    
    sim.data <- eventReactive(input$go, ignoreNULL = F,{

        tmp <- run_simulation(max_generations = input$max_generations, hybrid_fitness = input$hybrid_fitness, pref_ratio = input$pref_ratio,
                              male_trait_loci = input$male_trait_loci,
                              neutral_loci = input$neutral_loci, K_half = input$K_half,
                              meandispersal = input$meandispersal, growth_rate = input$growth_rate, dom_coefficient=input$dom_coefficient,
                              mating_trait_dominance = input$mating_trait_dominance)
        
        tmp
    })
    
    plot.data <- eventReactive(sim.data(),{
        format_plot_data(sim_output = sim.data(), num_trait_loci = input$male_trait_loci, neutral_loci = input$neutral_loci)
    })
    
    output$plot <- renderPlot({
        plot_freqs(plot.data(), fit_neutral=input$fit_neutral, fit_trait=input$fit_trait)
    })
    
    cline.data <- eventReactive(sim.data(),{
        calc_clines(format_output = plot.data())
    })
    
    output$table <- renderTable({
        cline.data()
    })
}

# run app
shinyApp(ui = ui, server = server)
