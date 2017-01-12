require(cowplot)
require(ggplot2)
require(reshape2)
require(ggbiplot)

PCAfunction <- function(data, type, extraction) {

keeps <- c("Replicate.Name", "Compound.Name", "Area", "column", "experiment", "extraction")
data <- data[keeps]

casted.data <- dcast(data, Replicate.Name + column + experiment + extraction ~ Compound.Name, value.var = "Area")

casted.data$experiment <- as.factor(str_sub(casted.data$experiment))
experiment <- casted.data$experiment

data.pca <- prcomp(casted.data[, 5:length(casted.data)], center = TRUE, scale. = TRUE) 

#print(data.pca)
#ids=paste0(data$experiment, data$injReplicate)
g <- ggbiplot(data.pca, obs.scale = 1, var.scale = 1, 
              groups = experiment, ellipse = TRUE, 
              circle = FALSE, labels=experiment, var.axes = TRUE)
g <-g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.position = 'none') + labs(title = paste("Vibrio Catabolote Repression:", type, extraction, sep = " "))
g <- g + scale_x_continuous(limits = c(-8, 10)) + scale_y_continuous(limits = c(-8, 8))
print(g)

save_plot(paste(type, extraction, "PCA", ".pdf",sep=" "), g,  base_aspect_ratio = 2, path = "~/Desktop/Google Drive/3. Ingalls Lab/Summer2016Bacteria/PCAplots")

}