#EddyAllData is normalized to internal standards, not to TIC
#This is updated to reflect that Area column in Eddy_All_Data.csv was not normalized to internal standard.  Replacing Area with pooPlus.Model.Norm column.
#This is plotting bar graphs of Enrichment factors, normalized to Station 7, for each compound
# _1451 Added estimate of standard deviation on the ratio of means, the enrichment factor
# _1941 Return bar graphs for Area by Station

setwd("~/Desktop/Google Drive/3. Ingalls Lab/Summer2016Bacteria/SRH")
#dir()

require(stringr)
require(plyr)
require(cowplot)
require(dplyr)

rm(list=ls())

#data1to2----

data <- read.csv("BacterialCultures_NoStandards_181116.csv", stringsAsFactors = FALSE)
data <- filter(data,Precursor.Ion.Name== "S-ribosylhomocysteine")
data <- data[grepl("1to2", data$Replicate.Name),]
data <- data[grepl("SWT|Glu", data$Replicate.Name),]
data$experiment <- NA
data$experiment <- ifelse(grepl("Glu", data$Replicate.Name, ignore.case = T), "VGlu",
                      ifelse(grepl("gly", data$Replicate.Name, ignore.case = T), "VGly",
                      ifelse(grepl("SWT-1|SWT-2|SWT-3", data$Replicate.Name, ignore.case = T), "VSWT","other")))



data$experiment <- as.factor(str_sub(data$experiment))
data$Area <- as.numeric(data$Area)


AreabyStation <- ddply(data, "experiment", summarize, area.mean=mean(Area), area.sd=sd(Area), Length=NROW(Area), tfrac=qt(p=.90, df=Length-1), Lower=area.mean-tfrac*area.sd/sqrt(Length), Upper=area.mean + tfrac*area.sd/sqrt(Length))
AreabyStation$experiment <- factor(AreabyStation$experiment, levels=c("VGlu", "VGly", "VSWT"))

p <- ggplot(AreabyStation, aes(y=area.mean, x=experiment)) + 
  geom_bar(stat="identity", fill="white", colour="black", width = 0.5, position = "dodge") +
  geom_errorbar(aes(ymin=area.mean-area.sd, ymax=area.mean+area.sd), width=0.2) +
  labs(title="SRH", y="Area", x="Station")+
  scale_y_continuous(expand = c(0,0))
print(p)
save_plot("SRH_VibrioCultures",".pdf",sep=""), p,  base_aspect_ratio = 2)

#Data1to100----

data2 <- read.csv("BacterialCultures_NoStandards_181116.csv", stringsAsFactors = FALSE)
data2 <- filter(data2,Precursor.Ion.Name== "S-ribosylhomocysteine")
data2 <- data2[grepl("1to100", data2$Replicate.Name),]
data2 <- data2[grepl("SWT|Glu", data2$Replicate.Name),]
data2$experiment <- NA
data2$experiment <- ifelse(grepl("Glu", data2$Replicate.Name, ignore.case = T), "VGlu",
                          ifelse(grepl("gly", data2$Replicate.Name, ignore.case = T), "VGly",
                                 ifelse(grepl("SWT-1|SWT-2|SWT-3", data2$Replicate.Name, ignore.case = T), "VSWT","other")))



data2$experiment <- as.factor(str_sub(data2$experiment))
data2$Area <- as.numeric(data2$Area)

AreabyStation2 <- ddply(data2, "experiment", summarize, area.mean=mean(Area), area.sd=sd(Area), Length=NROW(Area), tfrac=qt(p=.90, df=Length-1), Lower=area.mean-tfrac*area.sd/sqrt(Length), Upper=area.mean + tfrac*area.sd/sqrt(Length))
AreabyStation2$experiment <- factor(AreabyStation2$experiment, levels=c("VGlu", "VGly", "VSWT"))

p <- ggplot(AreabyStation2, aes(y=area.mean, x=experiment)) + 
  geom_bar(stat="identity", fill="white", colour="black", width = 0.5, position = "dodge") +
  geom_errorbar(aes(ymin=area.mean-area.sd, ymax=area.mean+area.sd), width=0.2) +
  labs(title="SRH", y="Area", x="Station")+
  scale_y_continuous(expand = c(0,0))
print(p)
save_plot("SRH_VibrioCultures_1to100",".pdf",sep=""), p,  base_aspect_ratio = 2)

#adapt everything below here to plot bar graph of enrichment factors by compound

bp2 <- ggplot(EF.summary, aes(x=Station.list, y=EF)) + 
  geom_hline(aes(yintercept=1), linetype="dashed", colour="grey54") + 
  geom_bar(stat="identity", fill="white", colour="black", width = 0.5, position = "dodge") +
  xlab("Station") +
  ylab("Enrichment Factor; normalized to Station 7") +
  labs(title = (paste(as.character(LM_Summary[i,1])))) + 
  scale_y_continuous(expand = c(0,0)) +
  geom_errorbar(aes(ymin=EF-StdDev, ymax=EF+StdDev), width=0.2)
bp2 <- bp2 + theme(axis.text.x = element_text(vjust=1)) 
#print (bp2)
save_plot(paste(as.character(LM_Summary[i, 1 ]),"_b",".pdf",sep=""), bp2,  base_aspect_ratio = 2)



}

print(i)
#All compound by compound loops are done


LM_Summary_trim <- LM_Summary[!is.na(LM_Summary$pANOVA),]

#Make a list of samples to with some, but not all NA data, to check
Check_Data = data.frame(matrix(vector(), 0, 1,
                               dimnames=list(c(), c("Compounds"))),
                        stringsAsFactors=F)
next.row2 <- 1
for (m in 1:length(LM_Summary$allNA))
{
  if(is.na(LM_Summary[m,2])==TRUE)
  {
    next
  }
  else
  {
    Check_Data[next.row2, 1] <- LM_Summary[m,1]
    next.row2 <- next.row2+1
  }
}


#sort data.frame according to Enrichment factors
#test2 <- test[order(test$EF_10_7),]


#save(LM_Summary, file="LM_Summary_HILIC_TQS.RData")
#save(LM_Summary_trim, file="LM_Summary_trim_HILIC_TQS.RData")
#save(Summary_stats, file="Summary_stats_HILIC_TQS.RData")
#save(Check_Data, file="Check_Data_Cyano_HILIC_TQS.RData")
#write.csv(LM_Summary, file="LM_Summary_HILIC_TQS.csv")
#write.csv(LM_Summary_trim, file="LM_Summary_trim_HILIC_TQS.csv")
#write.csv(Summary_stats, file="Summary_stats_HILIC_TQS.csv")
#write.csv(Check_Data, file="Check_Data_HILIC_TQS.csv")