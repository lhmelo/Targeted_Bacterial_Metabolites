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

BacteriaTargeted(Vibrio.data.neg, "VibrioHILIC", "neg")
source("BacteriaTargeted_test.R")
VibrioHILIC.array <- BacteriaTargeted_test(Vibrio.data.pos, "VibrioHILIC", "pos")

VibrioPos <- VibrioHILIC.array$df   
   
VibrioNeg <- read_csv("VibrioHILIC.neg.csv")
VibrioPos <- read_csv("VibrioHILIC.pos.csv")
Vibrio.HILIC <- rbind(VibrioNeg, VibrioPos)
write_csv(Vibrio.HILIC, "Vibrio.HILIC.csv")


#PCA----
source("PCAfunction.R")

PCAfunction(VibrioNeg, "Vibrio HILIC Neg", "AQ")
PCAfunction(VibrioPos, "Vibrio HILIC Pos", "AQ")





