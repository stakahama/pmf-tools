ftirpath <- "~/Documents/Work/UCSD/projects/CalMex/programs/ftir"
pjinfo <- read.delim(file.path(ftirpath,"userinputs/projectinfo.txt"))
filters <- subset(pjinfo,Type=="S" &
                  !FilterID %in% sprintf("TJ%03d",c(14,123,126,128,129,166)),
                  FilterID)
writeLines(unlist(filters),"samples.txt")
