ipak <- function(pkg){
		    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
	    	            install.packages(new.pkg, dependencies = TRUE, repos = "http://cran.us.r-project.org")
            sapply(pkg, require, character.only = TRUE)
}

# usage
packages <- c("ggplot2", "readr")
ipak(packages)

