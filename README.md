---
title: "Repo Overview"
author: "bjsugg --- @ncsu.edu"
date: "April 2021"
---

## Purpose and Primary Files  

This repo contains all relevant files used for a class project in *ST502 - Fundamentals of Statistical Inference II* at NC State University. A **complete end-to-end, in-depth overview** of the purpose of this project, the inference being explored, various methods of calculating confidence intervals, and GIF creation process can be found in the `project1gifGuide.RMD` file, with a HTML version also available for download under `project1gifGuide.html`.  

## Other Files  

All other files used to create this project are included as well, generally falling into the following categories:  

### Draft R Script  

These files were used for drafting the R code behind this project:  

* `project1codeFinal - GIF FUN.R` - The heart of this project, with all code used in the Monte Carlo simulation of data that is behind the final GIF output.  

* `gifBuildScript` - Exploratory code that was drafted for the initial GIF creation with the - `gganimate` package. Final version of this code is contained within the primary `.RMD` file.  

### The Data  

* `binom95CI.RData` - The output of the Monte Carlo simulation for each sample size n, from 1 to 200, which took 10-11 hours to generate. More details on the tables within this file can be found in the primary `.RMD` file.  

### GIF Output  

Four GIF files are included, which are generated with the code chunks within the primary `.RMD` file.

## Final Comments

Anyone can follow the below steps to recreate the GIF files contained within this repo:  

1. Download the `project1gifGuide.RMD` and `binom95CI.RData` files, and preferably store them together in the same directory.  

2. Open the `.RMD` file in R and install any packages in the `setup` code chunk that are not installed already.  

3. Navigate to the *Data Load* section of the `.RMD` file, and ensure the path of the current working directory is set to the location of the downloaded `.RData` file.  

4. Knit the `.RMD` file and subsequently find the four generated GIFs in the current working directory.  

Enjoy!  
