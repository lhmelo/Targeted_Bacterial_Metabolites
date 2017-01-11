require(dplyr)

ANOVAfunction <- function(data) {

model <- function(data)
{
  aov(Area ~ experiment, data = data)
  
}


anova.output <- dlply (data, .(Compound.Name), model)

juicy <- function(x) {
  c(summary(x)[[1]][["F value"]][[1]],
    summary(x)[[1]][["Pr(>F)"]][[1]])
}


anova.summary <- ldply(anova.output, juicy)
anova.summary <- data.frame(anova.summary)

colnames(anova.summary)[2] <- 'Fvalue'
colnames(anova.summary)[3] <- 'Pvalue'

return(invisible(anova.summary))

}
