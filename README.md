# HZAM

HZAM is a program to simulate hybrid zone dynamics in the presence of assortative mating, underdominance, and other variables. Originally written by [Darren Irwin](https://www.zoology.ubc.ca/~irwin/irwinlab/), I adapted it for Shiny with minimal changes—namely a simplified set of parameters, and support for partial dominance. To run the app on the web, use the following link: https://elinck.shinyapps.io/hzam_shiny/

Note that the server may time out at 100 active hours per month. To rehost locally, implement the following steps: 


1. Install [R](https://www.r-project.org/) and / or [RStudio](https://www.rstudio.com/products/rstudio/download/).

2. Paste the following code block into an R console: 

        pkgs <- c("shiny","mgcv","ggplot2","magrittr","viridis","boot","tidymv")
        dl_pkgs <- subset(pkgs,!pkgs %in% rownames(installed.packages()))
        if(length(dl_pkgs)!=0){
          for(i in dl_pkgs) install.packages(i)
        }
        library(shiny)
        runGitHub(username="elinck",repo="hzam_shiny")

If you find this tool useful, please cite the following paper: 

[Irwin, D.E. *In press.* Assortative mating in hybrid zones is remarkably ineffective in promoting speciation. American Naturalist.](https://doi.org/10.1086/708529)

(Thanks to [C.J. Battey](http://cjbattey.com/) for the procrastination inspiration and local hosting directions.)