require(dplyr)
setwd("~/Desktop/Google Drive/3. Ingalls Lab/Summer2016Bacteria")
source("BacteriaTargeted.R")

data <- read.csv("QC_outputBacteria_cyano_160727_akb.csv")

data <- filter(data,!Compound.Name== "13C Acetyl CoA",
               !Compound.Name== "13C-Isethionic Acid",
               !Compound.Name== "13C-Sulfoacetic Acid",
               !Compound.Name== "13C-Sulfolactic Acid",
               !Compound.Name== "d3-Cysteic Acid",
               !Compound.Name== "d4 Succinic Acid",
               !Compound.Name== "d4 Taurine",
               !Compound.Name== "D-Methionine",
               !Compound.Name== "Heavy Alanine",
               !Compound.Name== "d4-Tryptamine",
               !Compound.Name== "d-IAA",
               !Compound.Name== "D-Phenylalanine",
               !Compound.Name== "D-Tryptophan",
               !Compound.Name== "Vitamin B2_IS",
               !Compound.Name== "D8-Arachidonic Acid",
               !Compound.Name== "Heavy Histidine",
               !Compound.Name== "Heavy Isoleucine")

data$column <- "cyano"

data$experiment <- NA
data$experiment <- ifelse(grepl("DSS3", data$Replicate.Name, ignore.case = T), "DSS3", 
                  ifelse(grepl("SA11", data$Replicate.Name, ignore.case = T), "SA11", 
                  ifelse(grepl("Glu", data$Replicate.Name, ignore.case = T), "VGlu",
                  ifelse(grepl("gly", data$Replicate.Name, ignore.case = T), "VGly",
                  ifelse(grepl("SWT-1|SWT-2|SWT-3", data$Replicate.Name, ignore.case = T), "VSWT","other")))))

data$extraction <-NA
data$extraction <-ifelse(grepl("AQ", data$Replicate.Name, ignore.case = T), "AQ", "DCM")

Vibrio.data <- filter(data, experiment=="VGlu" | experiment=="VGly" | experiment=="VSWT")
                      
Vibrio.data.AQ <- filter(Vibrio.data, extraction=="AQ")
Vibrio.data.DCM <- filter(Vibrio.data, extraction=="DCM")

VibrioAQ <- BacteriaTargeted(Vibrio.data.AQ, "Vibrio", "AQ")
VibrioDCM <- BacteriaTargeted(Vibrio.data.DCM, "Vibrio", "DCM")

#VibrioAQ <- read_csv("Vibrio.AQ.csv")
#VibrioDCM <- read_csv("Vibrio.DCM.csv")
Vibrio.Cyano <- rbind(VibrioAQ, VibrioDCM)
write_csv(Vibrio.Cyano, "Vibrio.Cyano.csv")

drops <- c("X")
Vibrio.Cyano <- Vibrio.Cyano[ , !(names(Vibrio.Cyano) %in% drops)]

#PCA----
source("PCAfunction.R")

PCAfunction(VibrioAQ, "Vibrio Cyano", "AQ")
PCAfunction(VibrioDCM, "Vibrio Cyano", "DCM")

#ANOVA----
#
#

source("ANOVAfunction.R")
VibrioCyanoAQ_ANOVAsumm <- ANOVAfunction(VibrioAQ)
VibrioCyanoDCM_ANOVAsumm <- ANOVAfunction(VibrioDCM)

model <- function(data)
{
  aov(Area ~ experiment, data = data)
  
}


anova.output <- dlply (VibrioAQ, .(Compound.Name), model)

juicy <- function(x) {
  c(summary(x)[[1]][["F value"]][[1]],
    summary(x)[[1]][["Pr(>F)"]][[1]])
}


anova.summary <- ldply(anova.output, juicy)
anova.summary <- data.frame(anova.summary)

colnames(anova.summary)[2] <- 'Fvalue'
colnames(anova.summary)[3] <- 'Pvalue'

