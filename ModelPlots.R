#Plot's for Will's models 
#Canan Karakoc
#30 April 2024

#Plotting stuff
library(tidyverse)

mytheme <- theme_bw()+
  theme(axis.ticks.length = unit(.25, "cm"))+
  theme(legend.text = element_text(size=16))+
  theme(axis.text = element_text(size = 18, color = "black"), axis.title = element_text(size = 18))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(fill=NA, colour = "black", 
                                    size=1))+
  theme(strip.text.x = element_text(size = 18))+
  theme(legend.title=element_blank())+
  theme(panel.border = element_rect(fill=NA, colour = "black", 
                                    linewidth=1)) +
  theme(axis.text.x.top = element_blank(), axis.title.x.top = element_blank(),
        axis.text.y.right = element_blank(), axis.title.y.right = element_blank())+
  theme(axis.title.x = element_text(margin=margin(10,0,0)),
        axis.title.y = element_text(margin=margin(0,10,0,0)),
        axis.text.x = element_text(margin=margin(10,0,0,0)),
        axis.text.y = element_text(margin=margin(0,10,0,0)))+
  theme(legend.background=element_blank())

# Color blind palette
cbpalette <- c("#0072B2", "#D55E00","#009E73", "#CC79A7", "#56B4E9", "#999999", "#F0E442", "#000000")

# Consumer-resource model 
consumerResource <- read.table("~/GitHub/SporeCosts/data/consumer_resource_spore_model.tsv", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
consumerResource$values <- strsplit(consumerResource$V2, ",")

max_values <- max(sapply(consumerResource$values, length))
value_cols <- paste0("value", 1:max_values)
values_df  <- as.data.frame(do.call(rbind, lapply(consumerResource$values, `length<-`, max_values)))
colnames(values_df) <- value_cols

new_data <- cbind(consumerResource[,1, drop=FALSE], values_df)

new_data_matrix <- as.matrix(new_data[, -1])
rownames(new_data_matrix) <- new_data[, 1]
new_data_transposed <- t(new_data_matrix)

new_data_final <- as.data.frame(new_data_transposed, stringsAsFactors = FALSE)
for (i in 1:ncol(new_data_final)) {
  new_data_final[, i] <- as.numeric(new_data_final[, i])
}

legend_names <- c("resource concentration", "cell density", "spore density")  # Replace with your desired legend names

plotModel <- new_data_final%>%
  pivot_longer(-hours, names_to = "population", values_to = "abundance")

plotModel$population <- factor(plotModel$population, levels = c("resource_concentration", "cell_concentration", "spore_concentration"))
  
annotation_positions <- data.frame(
  population = c("resource", "vegetative", "spore"),
  hours = c(1, 30, 34),  # Example positions for each level
  abundance = c(10^3.5, 10^9.5, 10^8.5)  # Example positions for each level
)

plot1 <- ggplot(plotModel, aes(x = hours, y = abundance, color = population))+
  geom_line(size = 1.6)+
  mytheme+
  labs(x = "Time (h)", y = "Concentration")+
  theme(legend.position = c(0.25, 0.88))+
  scale_color_manual(labels = c("resource concentration", "cell density", "spore density"), 
                       values = c("#CC79A7", "#0072B2", "#D55E00"))+
  scale_y_continuous(limits = c(10^0, 10^10), trans = "log10", breaks = c(10^1,10^3,10^5,10^7,10^9), 
                     labels = function(x) parse(text = paste0("10^", format(log10(x), digits = 2))), sec.axis = dup_axis())+
    scale_x_continuous(sec.axis = dup_axis())+
  theme(legend.position = "none")+
  geom_text(data = annotation_positions,
            aes(label = population), color = c("#CC79A7", "#0072B2", "#D55E00"),
            x = annotation_positions$hours, y = log10(annotation_positions$abundance),
            hjust = -0.1, vjust = 0.5, size = 6)   

ggsave("~/GitHub/SporeCosts/figures/modelResults.png", plot = plot1, width = 5.6, height = 4.4, dpi = 300)
ggsave("~/GitHub/SporeCosts/figures/modelResults.pdf", plot = plot1, width = 5.6, height = 4.4, dpi = 300)
#############################################################################

# Spore efficiency 
efficiency <- read.table("~/GitHub/SporeCosts/data/efficiency_data.tsv", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
efficiency$values <- strsplit(efficiency$V2, ",")

max_values_eff <- max(sapply(efficiency$values, length))
value_cols_eff <- paste0("value", 1:max_values_eff)
values_df_eff  <- as.data.frame(do.call(rbind, lapply(efficiency$values, `length<-`, max_values_eff)))
colnames(values_df_eff) <- value_cols_eff

new_data_eff <- cbind(efficiency[,1, drop=FALSE], values_df_eff)

new_data_matrix_eff <- as.matrix(new_data_eff[, -1])
rownames(new_data_matrix_eff) <- new_data_eff[, 1]
new_data_transposed_eff <- t(new_data_matrix_eff)

new_data_final_eff <- as.data.frame(new_data_transposed_eff, stringsAsFactors = FALSE)
for (i in 1:ncol(new_data_final_eff)) {
  new_data_final_eff[, i] <- as.numeric(new_data_final_eff[, i])
}

efficiencyP <-  ggplot(new_data_final_eff, aes(x = spore_vs_cell_atp_ratio, y = sporulation_efficiency))+
  geom_vline(xintercept = 0.2668621700879765, linetype = "dashed")+
    geom_line(size = 1.3, color = "#0072B2")+
  mytheme+
    labs(x = "Ratio of energetic costs, spore/cell", y = "Sporulation efficiency")+ 
  scale_y_continuous(trans = "log10", breaks = c(10^-3, 10^-2, 10^-1), 
                     labels = function(x) parse(text = paste0("10^", format(log10(x), digits = 2))), sec.axis = dup_axis())+
    scale_x_continuous(trans = "log10", labels = function(x) parse(text = paste0("10^", format(log10(x), digits = 2))),sec.axis = dup_axis())+
    theme(legend.position = "none")+
    geom_segment(aes(x = 10^-4, y = 10^-2.8, xend = 10^-3.5, yend = 10^-2.8), linetype = "dashed") +
    annotate("text", x = 10^-2.2, y = 10^-2.89, label = "Emprical estimate", vjust = -0.5, size = 5)

  ggsave("~/GitHub/SporeCosts/figures/efficiency.png", plot = efficiencyP, width = 5.6, height = 4.4, dpi = 300)
  ggsave("~/GitHub/SporeCosts/figures/efficiency.pdf", plot = efficiencyP, width = 5.6, height = 4.4, dpi = 300)
##########################################################################################################
  
evoratio <- read.table("~/GitHub/SporeCosts/data/evo_ratio_conditional.tsv", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  evoratio$values <- strsplit(evoratio$V2, ",")
  
max_values_evo <- max(sapply(evoratio$values, length))
value_cols_evo <- paste0("value", 1:max_values_evo)
values_df_evo  <- as.data.frame(do.call(rbind, lapply(evoratio$values, `length<-`, max_values_evo)))
colnames(values_df_evo) <- value_cols_evo
  
new_data_evo <- cbind(evoratio[,1, drop=FALSE], values_df_evo)
  
new_data_matrix_evo <- as.matrix(new_data_evo[, -1])
rownames(new_data_matrix_evo) <- new_data_evo[, 1]
new_data_transposed_evo <- t(new_data_matrix_evo)
  
new_data_final_evo <- as.data.frame(new_data_transposed_evo, stringsAsFactors = FALSE)
for (i in 1:ncol(new_data_final_evo)) {
    new_data_final_evo[, i] <- as.numeric(new_data_final_evo[, i])
  }
  

evorates <- new_data_final_evo %>%
pivot_longer(-deletion_size, names_to = "rates", values_to = "values") %>%
ggplot( aes(x = log10(deletion_size), y = log10(values), color = rates))+
  geom_point(size = 2)+
    geom_line(size = 1.2, aes(linetype = rates))+
    scale_linetype_manual(values = c("dashed", "solid", "solid"))+
    xlab(expression("Deletion size (bp),"~Delta))+
  ylab(expression(atop("Deletion vs. subtition fixation rate under", paste("relaxed selection for spore formation, ",
                                       italic(d[del])(Delta)/italic(d[sub])))))+
  mytheme+
  theme(axis.title.y = element_text(size = 11))+
  theme(legend.position = c(0.3, 0.85))+
  theme(legend.text = element_text(size = 14))+
  scale_color_manual(values = c("#CC79A7", "#0072B2", "#D55E00"))
  
  
ggsave("~/GitHub/SporeCosts/figures/evorates.png", plot = evorates, width = 5.6, height = 4.4, dpi = 300)
ggsave("~/GitHub/SporeCosts/figures/evorates.pdf", plot = evorates, width = 5.6, height = 4.4, dpi = 300)

############################################################################################################

evoratioN <- read.table("~/GitHub/SporeCosts/data/evo_ratio.tsv", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
evoratioN$values <- strsplit(evoratioN$V2, ",")

max_values_evoN <- max(sapply(evoratioN$values, length))
value_cols_evoN <- paste0("value", 1:max_values_evoN)
values_df_evoN  <- as.data.frame(do.call(rbind, lapply(evoratioN$values, `length<-`, max_values_evoN)))
colnames(values_df_evoN) <- value_cols_evoN

new_data_evoN <- cbind(evoratioN[,1, drop=FALSE], values_df_evoN)

new_data_matrix_evoN <- as.matrix(new_data_evoN[, -1])
rownames(new_data_matrix_evoN) <- new_data_evoN[, 1]
new_data_transposed_evoN <- t(new_data_matrix_evoN)

new_data_final_evoN <- as.data.frame(new_data_transposed_evoN, stringsAsFactors = FALSE)
for (i in 1:ncol(new_data_final_evoN)) {
  new_data_final_evoN[, i] <- as.numeric(new_data_final_evoN[, i])
}


evoratesN <- new_data_final_evoN %>%
  pivot_longer(-deletion_size, names_to = "rates", values_to = "values") %>%
  ggplot( aes(x = log10(deletion_size), y = log10(values), color = rates))+
  geom_hline (yintercept = -0.6, linetype = "dashed")+
  geom_line(size = 1.2)+
  xlab(expression("Deletion size (bp),"~Delta))+
  ylab(expression(atop("Deletion vs. subtition fixation rate under", paste("relaxed selection for spore formation, ",                                                                        italic(d[del])(Delta)/italic(d[sub])))))+
  mytheme+
  theme(axis.title.y = element_text(size = 11))+
  theme(legend.position = c(0.3, 0.85))+
  theme(legend.text = element_text(size = 14))+
  scale_color_manual(values = c("#0072B2", "#D55E00"))+
  scale_x_continuous(sec.axis = dup_axis())+
  scale_y_continuous(sec.axis = dup_axis())

ggsave("~/GitHub/SporeCosts/figures/evoratesN.png", plot = evoratesN, width = 5.6, height = 4.4, dpi = 300)
ggsave("~/GitHub/SporeCosts/figures/evoratesN.pdf", plot = evoratesN, width = 5.6, height = 4.4, dpi = 300)

############################################################################################################

masize    <- read.table("~/GitHub/SporeCosts/data/ma_size.csv", header = FALSE, sep = ",", stringsAsFactors = FALSE)
masummary <- read.table("~/GitHub/SporeCosts/data/ma_summary.csv", header = FALSE, sep = ",", stringsAsFactors = FALSE)

### Combine with empirical - Other code ###

library(ggpubr)
library(grid)
figure <- ggarrange(efficiencyPlot, plot1, efficiencyP,
                    labels = c("A", "B", "C"),
                    ncol = 1, nrow = 3)
ggsave("~/GitHub/SporeCosts/figures/figure4.pdf", plot = figure, width = 6, height = 12, dpi = 300)
