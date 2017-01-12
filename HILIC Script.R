require(dplyr)
require(reshape2)

setwd("~/Desktop/Google Drive/3. Ingalls Lab/Summer2016Bacteria")
source("BacteriaTargeted.R")

#Neg HILIC-----
data.neg <- read.csv("QCd_QE_Bacteria_HILIC_neg.csv", skip=1)

data.neg <- data.neg[!grepl("Heavy", data.neg$Compound.Name),]



data.neg$column <- "HILIC"

data.neg$experiment <- NA
data.neg$experiment <- ifelse(grepl("DSS3", data.neg$Replicate.Name, ignore.case = T), "DSS3", 
                  ifelse(grepl("SA11", data.neg$Replicate.Name, ignore.case = T), "SA11", 
                  ifelse(grepl("Glu", data.neg$Replicate.Name, ignore.case = T), "VGlu",
                  ifelse(grepl("gly", data.neg$Replicate.Name, ignore.case = T), "VGly",
                  ifelse(grepl("SWT-1|SWT-2|SWT-3", data.neg$Replicate.Name, ignore.case = T), "VSWT","other")))))

data.neg$extraction <-NA
data.neg$extraction <-ifelse(grepl("AQ", data.neg$Replicate.Name, ignore.case = T), "AQ", "DCM")

Vibrio.data.neg <- filter(data.neg, experiment=="VGlu" | experiment=="VGly" | experiment=="VSWT")

#Pos HILIC-----
#
data.pos <- read.csv("QCd_QE_Bacteria_HILIC_pos.csv", skip=1)

data.pos <- data.pos[!grepl("Heavy", data.pos$Compound.Name),]




data.pos$column <- "HILIC"

data.pos$experiment <- NA
data.pos$experiment <- ifelse(grepl("DSS3", data.pos$Replicate.Name, ignore.case = T), "DSS3", 
                       ifelse(grepl("SA11", data.pos$Replicate.Name, ignore.case = T), "SA11",                         ifelse(grepl("Glu", data.pos$Replicate.Name, ignore.case = T), "VGlu",
                       ifelse(grepl("gly", data.pos$Replicate.Name, ignore.case = T), "VGly",
                      ifelse(grepl("SWT-1|SWT-2|SWT-3", data.pos$Replicate.Name, ignore.case = T), "VSWT","other")))))

data.pos$extraction <-NA
data.pos$extraction <-ifelse(grepl("AQ", data.pos$Replicate.Name, ignore.case = T), "AQ", "DCM")

Vibrio.data.pos <- filter(data.pos, experiment=="VGlu" | experiment=="VGly" | experiment=="VSWT")



#Run the function to get rid of NAs and empty analysis----

VibrioNeg <- BacteriaTargeted(Vibrio.data.neg, "VibrioHILIC", "neg")
VibrioPos <- BacteriaTargeted(Vibrio.data.pos, "VibrioHILIC", "pos")

   
   
VibrioNeg <- read_csv("VibrioHILIC.neg.csv")
VibrioPos <- read_csv("VibrioHILIC.pos.csv")
Vibrio.HILIC <- rbind(VibrioNeg, VibrioPos)
write_csv(Vibrio.HILIC, "Vibrio.HILIC.csv")


#PCA----
source("PCAfunction.R")

PCAfunction(VibrioNeg, "Vibrio HILIC Neg", "AQ")
PCAfunction(VibrioPos, "Vibrio HILIC Pos", "AQ")

#ANOVA----

source("ANOVAfunction.R")
VibrioHILICpos_ANOVAsumm <- ANOVAfunction(VibrioPos)
VibrioHILICneg_ANOVAsumm <- ANOVAfunction(VibrioNeg)

#Break up HILIC Pos---

Unique.Names <- unique.names(VibrioPos)
length(Unique.Names)
Unique.Names_1 <-Unique.Names[1:12]
Unique.Names_2 <-Unique.Names[13:24]
Unique.Names_3 <-Unique.Names[25:36]
Unique.Names_4 <-Unique.Names[37:48]
Unique.Names_5 <-Unique.Names[49:60]
Unique.Names_6 <-Unique.Names[61:63]


VibrioPos1 <- subset(VibrioPos, Compound.Name %in% Unique.Names_1)
VibrioPos2 <- subset(VibrioPos, Compound.Name %in% Unique.Names_2)
VibrioPos3 <- subset(VibrioPos, Compound.Name %in% Unique.Names_3)
VibrioPos4 <- subset(VibrioPos, Compound.Name %in% Unique.Names_4)
VibrioPos5 <- subset(VibrioPos, Compound.Name %in% Unique.Names_5)
VibrioPos6 <- subset(VibrioPos, Compound.Name %in% Unique.Names_6)


#Faceted bargraphs
```{r}

#set file path for plots
plotsavepath <- "~/Desktop/Google Drive/3. Ingalls Lab/Summer2016Bacteria/CyanoPlots"

source("FacetGraphs.R")

FacetGraphs(VibrioNeg, "Vibrio", "HILIC", "Neg")

FacetGraphs(VibrioPos1, "Vibrio", "HILIC", "Pos")
FacetGraphs(VibrioPos2, "Vibrio", "HILIC", "Pos")
FacetGraphs(VibrioPos3, "Vibrio", "HILIC", "Pos")
FacetGraphs(VibrioPos4, "Vibrio", "HILIC", "Pos")
FacetGraphs(VibrioPos5, "Vibrio", "HILIC", "Pos")
FacetGraphs(VibrioPos6, "Vibrio", "HILIC", "Pos")


```


