ipak <- function(pkg){
		    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
	    	            install.packages(new.pkg, dependencies = TRUE, repos = "http://cran.us.r-project.org")
            sapply(pkg, require, character.only = TRUE)
}

# usage
packages <- c("stringi", "ggplot2", "readr", "precrec", "pROC", "ggpubr", "latex2exp")
ipak(packages)

