#Barplot
library(ggplot2)
library(reshape)
x1 <- read.table("bargraph_percentage.txt", header=TRUE)
x1$Sample <- reorder(x1$Sample, rowSums(x1[-1]))
mdf1=melt(x1,id.vars="Sample")
pdf("barplot.pdf", width = 30, height = 30)
ggplot(mdf1, aes(Sample, value, fill=variable)) + geom_bar(stat="identity") + ylab("Percentage of Filtered Positions") + xlab("Samples") + theme(text = element_text(size=9)) + scale_fill_manual(name="Reason for filtered out positions", values=c("#08306b", "black", "orange", "darkgrey", "#fdd0a2", "#7f2704")) + ggtitle("title_here") + ylim(0, 100) + theme(text = element_text(size=10), panel.background = element_rect(fill = 'white', colour = 'white'), plot.title = element_text(size=20, face="bold", margin = margin(10, 0, 10, 0)), axis.ticks.y = element_blank(), axis.ticks.x = element_blank(),  axis.text.x = element_text(colour = "black", face= "bold.italic", angle = 90)) + theme(legend.position = c(0.6, 0.7), legend.direction = "horizontal")
dev.off()

#FQ heatmap
library(GMD)
x1 <- read.table("temp_Only_filtered_positions_for_closely_matrix_FQ.txt", header=TRUE)
y1 <- as.matrix(x1)
yt1 <- t(y1)
range<-colorRampPalette(c("white", "black", "yellow", "#a63603", "blue"), space="rgb")
pdf("temp_Only_filtered_positions_for_closely_matrix_FQ.pdf", width = 30, height = 30)
heatmap.3(yt1, dendrogram="none", Rowv=TRUE, Colv=TRUE, color.FUN=range, main="title_here (Red = Positions filtered out due to Low FQ)", key=FALSE)
legend("topleft",legend=c("Unmapped","Reference", "Variant","LowFQ","Other"),fill=c("white", "black", "yellow", "#a63603", "blue"), border=TRUE, bty="n", y.intersp = 0.8, cex=0.9)
dev.off()

#Dp Heatmap
library(GMD)
x1 <- read.table("temp_Only_filtered_positions_for_closely_matrix_DP.txt", header=TRUE)
y1 <- as.matrix(x1)
yt1 <- t(y1)
range<-colorRampPalette(c("white", "black", "yellow", "#a63603", "blue"), space="rgb")
pdf("temp_Only_filtered_positions_for_closely_matrix_DP.pdf", width = 30, height = 30)
heatmap.3(yt1, dendrogram="none", Rowv=TRUE, Colv=TRUE, color.FUN=range, main="title_here (Red = Positions filtered out due to Low DP)", key=FALSE)
legend("topleft",legend=c("Unmapped","Reference", "Variant","LowDP","Other"),fill=c("white", "black", "yellow", "#a63603", "blue"), border=TRUE, bty="n", y.intersp = 0.8, cex=0.9)
dev.off()


#Depth Barplot
library(GMD)
library(reshape)
library(ggplot2)
x1 <- read.table("DP_bargraph_counts.txt", header=TRUE)
mdf1=melt(x1,id.vars="Sample")
pdf("barplot_DP.pdf", width = 30, height = 30)
ggplot(mdf1, aes(Sample, value, fill=variable)) + geom_bar(stat="identity") + ylab("No of Filtered Positions(DP)") + xlab("Samples") + theme(text = element_text(size=9)) + scale_fill_manual(name="Reason for filtered out positions", values=c("black", "orange", "darkgrey", "#fdd0a2", "#7f2704")) + ggtitle("title_here") + theme(text = element_text(size=15), panel.background = element_rect(fill = 'white', colour = 'white'), plot.title = element_text(size=20, face="bold", margin = margin(10, 0, 10, 0)), axis.ticks.y = element_blank(), axis.ticks.x = element_blank(),  axis.text.x = element_text(colour = "black", face= "bold.italic", angle = 90)) + theme(legend.position = c(0.7, 1), legend.direction = "horizontal")
dev.off()

#Depth Heatmap
library(GMD)
x1 <- read.table("filtered_DP_values.txt", header=TRUE)
y1 <- as.matrix(x1)
y1[is.na(y1)] <- 0
y1[y1 >= 10] <- NA
yt1 <- t(y1)
matrix.min  <- min(yt1[ yt1!=0 ], na.rm=TRUE)
matrix.min
matrix.max <- max( yt1, na.rm = TRUE )
matrix.max
pairs.breaks <- c(0, seq( 1, matrix.max, 0.5))
range<-colorRampPalette(c("#084594", "#4292c6", "#9ecae1", "#deebf7", "#ffffe5","#fee391","#fe9929","#cc4c02", "#cc4c02"))
pdf("DP_position_analysis.pdf", width = 15, height = 15)
heatmap.3(yt1, dendrogram="none", Rowv=FALSE, Colv=FALSE, main="title_here Dp position analysis", sub = "Black = Variant position with depth >= 10(Passed), Dark Blue = No variant/Reference Allele", key=TRUE, na.color="black", breaks = pairs.breaks, color.FUN=range)
dev.off()
