list.of.packages <- c('knitr', 'optparse', 'kableExtra', 'VennDiagram')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
