# load libraries
library(shiny); library(ggplot2); library(viridis); library(mgcv); 
library(boot); library(magrittr); library(tidymv)

# define UI
ui <- fluidPage(
    
    # # title panel
    titlePanel("HZAM: simulating hybrid zones with assortative mating"),
    
    # control panel 
    sidebarLayout(
        sidebarPanel(
            actionButton("new","Reset Simulation",width="100%"),
            div(helpText("Press to initiate or reset a simulation"),style="font-size:75%"),
            actionButton("run","Run Simulation",width="100%"),
            div(helpText("Press to run simulation"),style="font-size:75%"),
            numericInput("max_generations","Number of Generations",20,min=1,max=500),
            numericInput("update_interval","Plot Update Interval",10,min=1,max=500),
            selectInput('mating_trait_dominance', 'Mating Trait Expression:', c("Dominant / Partially Dominant"=1,
                                                                                "Additive"=2)),
            sliderInput("dom_coefficient","Partial Dominance Coefficient",value=0.75,min=0,max=1),
            div(helpText("Proportion of homozygous dominant phenotype expressed by heterozygotes"),style="font-size:75%"),
            sliderInput("K_half","Per-species carrying capacity",800,min=400,max=1600,step = 200),
            sliderInput("hybrid_fitness","Hybrid Fitness",value=0.98,min=0,max=1),
            sliderInput("pref_ratio","Assortative Mating Preference Ratio",value=0.75,min=0.00001,max=0.99999),
            div(helpText("Indicates ratio of heterospecific : homospecific matings"),style="font-size:75%"),
            sliderInput("male_trait_loci","Number of Trait Loci",value=3,min=1,max=5),
            sliderInput("neutral_loci","Number of Neutral Loci",value=3,min=1,max=5),
            sliderInput("meandispersal", "Average Dispersal Radius",value=0.01,min=0.001,max=0.1,step=0.001),
            sliderInput("growth_rate", "Population Growth Rate",value=1.05,min=1,max=1.5),
            checkboxInput("fit_neutral","Fit Neutral Cline",value = F),
            checkboxInput("fit_trait","Fit Trait Cline",value = F)
        ),
        
        # plot and table panel (with explainer)
        mainPanel(#style="position:fixed;margin-left:32vw;",
            plotOutput("plot"),
            div(style="height:5vh;"),
            textOutput("gencounter"),
            div(style="height:5vh;"),
            tableOutput("table"),
            div(style="height:5vh;"),
            h4("About HZAM"),
            div(p("HZAM, originally written by", a(href="https://www.zoology.ubc.ca/~irwin/irwinlab/",'Darren Irwin', .noWS = "outside"), 
                  "and adapted for Shiny by", a(href="https://elinck.org/",'Ethan Linck,', .noWS = "outside"), "is a program to
                simulate the influence of assortative mating in a one-dimensional hybrid zone in concert with the 
                      effects of underdominance, variation in dispersal, and other parameters. Code for the program 
                  (and a version that can be run locally) is available ",
                  a(href="http://github.com/elinck/hzam_shiny",'here.', .noWS = "outside"),
                "If you find bugs or have feature suggestions, please write me at ",
                a(href="mailto:ethanblinck@gmail.com","ethanblinck@gmail.com.", 
                  .noWS = c("after-begin", "before-end")))), 
            div(style="height:1.5vh;"),
            h4("Citing HZAM"),
            div(p("If you found this program useful, please cite the following publication:",
                    a(href="https://doi.org/10.1086/708529",
                      'Irwin, D.E. In press. Assortative mating in hybrid zones is remarkably ineffective in promoting speciation. 
                         American Naturalist.', .noWS = c("after-begin", "before-end"))
                    ))
        )
    )
)


# define server
server <- function(input, output, session) {
    
    source("./dev.R")
    
    # set up reactiveVal with population information
    matrix <- reactiveVal()
    vals <- reactiveValues(counter = 0, done=0, initial=0)
    
    # run when reset button is pushed to generate blank simulation 
    observeEvent(input$new, ignoreNULL = F,{
        
        # resets simulation completion tracker
        vals$done <- 0
        
        # resets iteration tracker
        vals$counter <- 0
        
        # indicates program shouldn't fit cline
        vals$initial <- 0
        
        # calls function to make population matrix
        pop_matrix <- make_matrix(K_half = input$K_half, neutral_loci = input$neutral_loci, 
                                  male_trait_loci = input$male_trait_loci)
        
        # weird reactive value assignment thing
        matrix(pop_matrix)
    })
    
    # run simulation
    observeEvent(input$run, ignoreNULL = F,{
        
        observe({
            
            isolate({
                
                # variable to indicate program should fit cline
                vals$initial <- 1
                
                # validation if program has run
                if(vals$done==1){
                    matrix(NULL)
                    validate(need(is.null(matrix())==FALSE,""))
                }
                
                # function to run simulation for plot update interval
                new_pop <- run_simulation(pop_matrix=matrix(),hybrid_fitness = input$hybrid_fitness, pref_ratio = input$pref_ratio,
                                          meandispersal = input$meandispersal, growth_rate = input$growth_rate, dom_coefficient=input$dom_coefficient, 
                                          mating_trait_dominance = input$mating_trait_dominance, K_half = input$K_half, input$update_interval)
                
                # assign last state to population matrix
                matrix(new_pop)
                
                # update iteration counter 
                vals$counter <- vals$counter + 1 
                
            })
            
            # this updates the plot at the proper interval
            if (isolate(vals$counter) < (input$max_generations / input$update_interval)){
                invalidateLater(0, session)
            }
            
            # this indicates the simulation has finished running
            if (isolate(vals$counter) >= (input$max_generations / input$update_interval)){
                vals$done <- 1
            }
            
        })
        
    })
    
    
    # format plot data every time pop matrix is updated
    plot.data <- reactive({
        
        format_plot_data(sim_output = matrix(), num_trait_loci = input$male_trait_loci, neutral_loci = input$neutral_loci)
        
    })
    
    # plot data every time pop matrix is updated
    output$plot <- renderPlot({
        
        validate(need(is.null(matrix())==FALSE,"Please reset the simulation"))
        
        plot_freqs(plot.data(), fit_neutral=input$fit_neutral, fit_trait=input$fit_trait)
        
    })
    
    # fit cline every time pop matrix is updated
    cline.data <- reactive({
        
        validate(need(is.null(matrix())==FALSE,""))
        
        if(vals$initial==1){
            calc_clines(format_output = plot.data())
        }
    })
    
    # fit cline table time pop matrix is updated
    output$table <- renderTable({
        
        if(vals$initial==1){
            cline.data()
        }
    })
    
    # generation counter
    output$gencounter <- renderText({
        
        validate(need(is.null(matrix())==FALSE,""))
        
        paste0("Generation ",(vals$counter*input$update_interval))
    })
}

# run app
shinyApp(ui = ui, server = server)

