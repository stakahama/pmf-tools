## for coloring dendrogram branches and leaves.
## example usage:

## hc <- hclust(dist(normorg[use,]))
## d <- as.dendrogram(hc)
## ngroups <- 10
## grps <- cutree(hc,k = ngroups)
## brew <- rainbow(ngroups)
## p <- dendrapply(d,colLab,grps,brew)
## plot(p)

colLab <- function(n,g,mycols) {
  allchildren <- sapply(split(names(g),g),function(a,b)
                        all(b %in% a), labels(n) )
  if(any(allchildren))
    attr(n, "edgePar") <-
      c(attr(n,"edgePar"), list(col = mycols[which(allchildren)]))
  if(is.leaf(n))
    attr(n, "nodePar") <-
      c(attr(n,"nodePar"),list(pch=NA,
                       lab.col = mycols[g[attr(n,"label")]],
                       lab.font= 1, lab.cex=0.6))
  return(n)
}

attr(colLab,"example") <- 
  "hc <- hclust(dist(X))
d <- as.dendrogram(hc)
ngroups <- 10
grps <- cutree(hc,k = ngroups)
brew <- rainbow(ngroups)
p <- dendrapply(d,colLab,grps,brew)
plot(p)"
