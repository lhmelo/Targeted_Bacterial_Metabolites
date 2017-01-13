
BacteriaTargeted_test <- function(data, type, extraction) {
  
Unique.names <-unique(data$Compound.Name)

Next.row <- 1

Check.data = data.frame(matrix(vector(), 0, 3,
                               dimnames=list(c(), c("Compound.Name", "AnyNA", "ALLNA"))),
                        stringsAsFactors=F)

Names <- colnames(data)

temp.good = data.frame(matrix(vector(), 0, length(Names),
                                        dimnames=list(c(), c(Names))),
                                 stringsAsFactors=F)

#i = loop for each compound
for (i in 1:length(Unique.names))
{Compound_set <- subset(data, Compound.Name == Unique.names[i])

Check.data[i,1] <- as.character(Unique.names[i])

if (anyNA(Compound_set$Area)==TRUE)
{
  Check.data[i,2] <-TRUE
  {
    if (all(is.na(Compound_set$Area))==TRUE)
      Check.data[i,3] <-TRUE
  }
}
}

#Make a list of compounds we detected
Check.data_trim <- Check.data[is.na(Check.data$ALLNA),]

for (i in 1:length(Check.data_trim$Compound.Name)) {
  
  temp <- subset(data, Compound.Name == Check.data_trim[i,1])
  
  temp.good<- rbind(temp.good, temp)
  


temp.good <- replace_na(temp.good, list(Area=0))
}
write_csv(temp.good, paste(type, extraction, "csv", sep = "."))
print("done")

test.array <- list(df=temp.good, names=Unique.names)
return(invisible (test.array))

}


