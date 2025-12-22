library(ggplot2)
library(dplyr)
library(cowplot)
library(svglite)


storage_fastq<-read.csv("~/Desktop/Subramaniam_Lab/SEA_AD_Project/meta_data/storage_data_fastq.csv")
storage_h5<-read.csv("~/Desktop/Subramaniam_Lab/SEA_AD_Project/meta_data/storage_data_h5.csv")



#### For FASTQ ###############################################################################################

# Reorder Tissue levels
storage_fastq$Tissue <- factor(storage_fastq$Tissue, levels = c("MTG", "DLPFC", "MEC", "PVC", "STG"))

# Pre-compute totals per Tissue for the top labels
totals <- storage_fastq %>%
  group_by(Tissue) %>%
  summarise(total = sum(Storage), .groups = "drop")

# Build the stacked bar plot 
fastq <- ggplot(storage_fastq, aes(x = Tissue, y = Storage, fill = Type)) +
  geom_col() +
  geom_text(
    data = totals,
    aes(x = Tissue, y = total, label = total),
    fontface = "bold",      # bold text
    vjust = -0.5,           # lift the label above the bar
    inherit.aes = FALSE) +
  labs(title='Storage Needs for Each Brain Region (FASTQ)',x = "Tissue", y = "Storage (TB)") +
  scale_fill_manual(values = c(snRNAseq = "#704D9E",snATACseq = "#A9A9A9" ))+
  theme_minimal()+
  theme(
    axis.text.y  = element_text(size = 16),  
    axis.title.y = element_text(size = 17),
    plot.title = element_text(size = 20, face="bold"))


fastq

#### For H5 ##################################################################################################

# Reorder Tissue levels
storage_h5$Tissue <- factor(storage_h5$Tissue, levels = c("ROSMAP-DLPFC","MTG", "DLPFC", "MEC", "PVC", "STG"))

# Pre-compute totals per Tissue for the top labels
totals <- storage_h5%>%
  group_by(Tissue) %>%
  summarise(total = sum(Storage), .groups = "drop")

# Build the stacked bar plot 
h5 <- ggplot(storage_h5, aes(x = Tissue, y = Storage, fill = Type)) +
  geom_col() +
  geom_text(
    data = totals,
    aes(x = Tissue, y = total, label = total),
    fontface = "bold",      # bold text
    vjust = -0.5,           # lift the label above the bar
    inherit.aes = FALSE) +
  labs(title='Storage Needs for Each Brain Region (H5, Peak, Fragment)',x = "Tissue", y = "Storage (TB)") +
  scale_fill_manual(values = c(snRNAseq = "#5D4A72",snATACseq = "#A9A9A9"))+
  theme_minimal()+
  theme(
    axis.text.y  = element_text(size = 16),  
    axis.title.y = element_text(size = 17),
    plot.title = element_text(size = 20, face="bold"))


h5



#### Set same Y-axis and get rid of double legends ###########################################################
y_min <- 0
y_max <- max(c(ggplot_build(h5)$data[[1]]$y,
               ggplot_build(fastq)$data[[1]]$y))

y_max <- ceiling(y_max)

h5 <- h5 + ylim(y_min, y_max)
fastq <- fastq + ylim(y_min, y_max)

 theme(legend.position = "none")

h5 <- h5 + theme(legend.position = "none") # removes legend

combined<-plot_grid(fastq, h5, ncol = 2, align = 'v',rel_widths = c(1, 1), scale = 1  )
combined


ggsave("~/Desktop/Subramaniam_Lab/SEA_AD_Project/meta_data/storage_plot.svg", 
       plot = combined, width = 17, height = 5, units = "in")






