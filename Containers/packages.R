ipak <- function(pkg){
		    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
	    	            install.packages(new.pkg, dependencies = TRUE, repos = "http://cran.us.r-project.org")
            sapply(pkg, require, character.only = TRUE)
}

# usage
packages <- c("Rcpp", "doParallel", "foreach", "ggplot2", "data.table", "RcppArmadillo",
	      "glue", "doRNG", "iterators", "HDInterval")
ipak(packages)

if (!require("BiocManager", quietly = TRUE))
	    install.packages("BiocManager")

BiocManager::install("SummarizedExperiment")
