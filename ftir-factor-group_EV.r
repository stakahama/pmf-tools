
## inputs:
input <- commandArgs()
pattern <- "--file=(.+)"
srcpath <- gsub('~+~'," ",dirname(sub(pattern,"\\1",input[grepl(pattern,input)])),fixed=TRUE)

argv <- tail(input,-grep("--args",input,fixed=TRUE))
filename <- argv[1]

## contents:
source(file.path(srcpath,"functions/io.R"))
args <- readArgs(filename)
for(p in names(args))
  assign(p,args[[p]])

## libraries
library(plyr)
library(dplyr)
library(reshape2)
library(ggplot2)
theme_set(theme_bw())

## files
file.clust <- list.files(file.path(FOLDER, "Allplots"),
                        "factors-[0-9]+-clusters\\.txt", full=TRUE)
files.ev <- list.files(FOLDER,"ExplainedVariation.txt", full=TRUE, rec=TRUE)


## read

readMelt <- function(x) {
  df <- add_rownames(read.delim(x, check.names=FALSE), "sample")
  melt(df, id.vars="sample", 
       variable.name="comp",
       value.name="ev")
}

pattern.comp <- "(soln_[0-9]{3})-([0-9]{2})"

clust <- with(read.delim(file.clust), setNames(grp, comp))
simgrid <- add_rownames(read.delim(file.path(FOLDER, "simgrid.txt"), row.names=1), "soln")

table.ev <- ldply(files.ev, readMelt) %>% filter(comp!="Resid")
table.ev$grp <- clust[as.character(table.ev$comp)]
table.ev$soln <- sub(pattern.comp,"\\1",table.ev$comp)

## if similar factors exist within the same solution set,
##   aggregate them
table.ev <- table.ev %>% group_by(soln, grp, sample) %>%
  summarize(ev = sum(ev)) %>% ungroup()

table.ev <- inner_join(table.ev, select(simgrid, soln, nFactors), by="soln")

## table.ev <- table.ev %>% group_by(nFactors, grp) %>%
##   mutate(soln.index=factor(unclass(factor(soln)))) %>% ungroup()

## plot and export

outfile <- c(
  pdf=sub("\\.txt$", "_EV.pdf", basename(file.clust)),
  txt=sub("\\.txt$", "_EV.txt", basename(file.clust))
  )


library(lattice)
library(latticeExtra)


## for strip text labeling
table.ev <- table.ev %>% 
  mutate(nFactors=factor(nFactors),
         grp=factor(grp))

mycolors <- with(list(n=nlevels(table.ev$grp)),
                 setNames(colorRampPalette(c("steelblue","darkred"))(n),seq(n)))
out <- useOuterStrips(
  densityplot(~ev| grp*nFactors,
              groups = soln,
              data=table.ev,
              plot.points=FALSE, 
              as.table=TRUE,
              par.settings=list(superpose.line=list(col=mycolors)),              
              xlab="Explained Variation (%)")
  )

pdf(file.path(FOLDER, "Allplots", outfile["pdf"]),
    width=8,height=5)
print(out)
dev.off()

write.table(table.ev %>% select(soln, grp, sample, ev) %>%
              convertFactor("grp") %>% dplyr::rename(EV=ev),
            file.path(FOLDER, "Allplots", outfile["txt"]),
            sep="\t", row.names=FALSE)
