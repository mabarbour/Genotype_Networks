adonis.table <- function(adonis.object, data.frame){
  x <- adonis.object
  mean.sample.size <- mean(table(data.frame$Genotype))
  genetic.variance <- (x$aov.tab$MeanSqs[1] - x$aov.tab$MeanSqs[2])/mean.sample.size
  heritability <- genetic.variance/(genetic.variance + x$aov.tab$MeanSqs[2])
  
  table.out <- data.frame(df = x$aov.tab$Df[1:2],
                          SS = x$aov.tab$SumsOfSqs[1:2],
                          MS = x$aov.tab$MeanSqs[1:2],
                          Heritability = c(heritability,""),
                          P = c(x$aov.tab$'Pr(>F)'[1],""))
  table.out
}

