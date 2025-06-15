#setwd() to relevant folder getwd() to check (CAN DO THIS VIA SESSION TAB)
getwd()
rm(list=ls())
setwd("/Users/zhangshengyuan/Library/CloudStorage/OneDrive-UniversityCollegeLondon/05 PIV/share_code")
###################FORMAT THE DATA for the drug combinations#####################
extract.2combtherapy <- function(mergedcsv, drugno_A, drugno_B, ConcUnit){
  # mergedcsv <- "merged.csv"
  # drugno_A <- 1
  # drugno_B <- 2
  platedf <- read.csv(mergedcsv)
  #trying to sort error for duplicate rows(didn't work)
  #platedf <- read.csv("merged.csv", row.names = NULL)
  # work out which of the 4 drugs the user has NOT chosen
  exclude <- c(1:4)
  exclude <- exclude[exclude != drugno_A]
  exclude <- exclude[exclude != drugno_B]
  # remove wells with other drugs
  platedf <- platedf[platedf[paste0("drug", exclude[1], "_conc")] == 0, ]
  platedf <- platedf[platedf[paste0("drug", exclude[2], "_conc")] == 0, ]
  # select nonzero drug wells
  platedf$conc1 <- platedf[paste0("drug", drugno_A, "_conc")]
  platedf$conc2 <- platedf[paste0("drug", drugno_B, "_conc")]
  analysisdf <- platedf[platedf$conc1 != 0, ]
  analysisdf <- analysisdf[analysisdf$conc2 != 0, ]
  # include the data with concentration of 0
  analysisdf <- rbind(analysisdf, platedf[platedf$conc1 == unique(analysisdf$conc1)[1,] & platedf$conc2 == 0,])
  analysisdf <- rbind(analysisdf, platedf[platedf$conc1 == unique(analysisdf$conc1)[2,] & platedf$conc2 == 0,])
  analysisdf <- rbind(analysisdf, platedf[platedf$conc1 == unique(analysisdf$conc1)[3,] & platedf$conc2 == 0,])
  analysisdf <- rbind(analysisdf, platedf[platedf$conc1 == unique(analysisdf$conc1)[4,] & platedf$conc2 == 0,])
  analysisdf <- rbind(analysisdf, platedf[platedf$conc1 == 0 & platedf$conc2 == unique(analysisdf$conc2)[1,],])
  analysisdf <- rbind(analysisdf, platedf[platedf$conc1 == 0 & platedf$conc2 == unique(analysisdf$conc2)[2,],])
  analysisdf <- rbind(analysisdf, platedf[platedf$conc1 == 0 & platedf$conc2 == unique(analysisdf$conc2)[3,],])
  analysisdf <- rbind(analysisdf, platedf[platedf$conc1 == 0 & platedf$conc2 == unique(analysisdf$conc2)[4,],])
  # include the response of 0 when no drug is added
  df0 <- platedf[platedf$conc1 == 0 & platedf$conc2 ==0,]
  df0$viral_percent_inh <- platedf[1,]$virus_control_od
  analysisdf <- rbind(analysisdf, df0) 
  
  # make a clean analysis dataframe
  df <- data.frame(matrix(nrow = nrow(analysisdf), ncol = 7))
  colnames(df) <- c("block_id", "drug1", "drug2", "conc1", "conc2", "response", "conc_unit") 
  
  
  #PairIndex
  df$block_id <- 1
  
  # drug name and concentration
  df$drug1 <- analysisdf[[paste0("drug", drugno_A, "_name")]]
  df$drug2 <- analysisdf[[paste0("drug", drugno_B, "_name")]]
  df$conc1 <- analysisdf[[paste0("drug", drugno_A, "_conc")]]
  df$conc2 <- analysisdf[[paste0("drug", drugno_B, "_conc")]]
  
  # drug response
  df$response<- analysisdf$viral_percent_inh
  
  # unit of concentration
  df$conc_unit <- ConcUnit
  
  # trying to sort row problem by Ensuring unique row names in the resulting data frame didn't work
  #row.names(df) <- make.names(seq_len(nrow(df)), unique = TRUE)
  
  df
}

#####DRUG COMBINATION AND ASSAYS TO COMBINE AND TO RUN SYNERGY TEST BETWEEN#####

#1=rem,2=favi,3=moln,4=riba

dfcombo <- extract.2combtherapy(mergedcsv="merged.csv", drugno_A = 1, drugno_B = 2, ConcUnit = "mg/L")
dfcombo2 <- extract.2combtherapy(mergedcsv="merged.csv", drugno_A = 2, drugno_B = 3, ConcUnit = "mg/L")
dfcombo3 <- extract.2combtherapy(mergedcsv="merged.csv", drugno_A = 1, drugno_B = 3, ConcUnit = "mg/L")
dfcombo4 <- extract.2combtherapy(mergedcsv="merged.csv", drugno_A = 1, drugno_B = 4, ConcUnit = "mg/L")
dfcombo5 <- extract.2combtherapy(mergedcsv="merged.csv", drugno_A = 2, drugno_B = 4, ConcUnit = "mg/L")
dfcombo6 <- extract.2combtherapy(mergedcsv="merged.csv", drugno_A = 3, drugno_B = 4, ConcUnit = "mg/L")

#### RUN THE DOSE/SYNERGY PLOTS ####

library(synergyfinder)
#res <- ReshapeData(data = df, data_type = "inhibition", impute = T,impute_method = NULL,noise = F)
res <- ReshapeData(data = dfcombo, data_type = "inhibition", impute = F,impute_method = NULL,noise = F)
res2 <- ReshapeData(data = dfcombo2, data_type = "inhibition", impute = F,impute_method = NULL,noise = F)
res3 <- ReshapeData(data = dfcombo3, data_type = "inhibition", impute = F,impute_method = NULL,noise = F)
res4 <- ReshapeData(data = dfcombo4, data_type = "inhibition", impute = F,impute_method = NULL,noise = F)
res5 <- ReshapeData(data = dfcombo5, data_type = "inhibition", impute = F,impute_method = NULL,noise = F)
res6 <- ReshapeData(data = dfcombo6, data_type = "inhibition", impute = F,impute_method = NULL,noise = F)
# check the missing value in raw data, usually the result without drug treatment
which(is.na(res$response))
# Calculate the synergy score
res <- CalculateSynergy(data = res, method = c("ZIP", "HSA", "Bliss", "Loewe"), 
                        Emin = NA, Emax = NA, correct_baseline = "part")
res2 <- CalculateSynergy(data = res2, method = c("ZIP", "HSA", "Bliss", "Loewe"), 
                         Emin = NA, Emax = NA, correct_baseline = "part")
res3 <- CalculateSynergy(data = res3, method = c("ZIP", "HSA", "Bliss", "Loewe"), 
                         Emin = NA, Emax = NA, correct_baseline = "part")
res4 <- CalculateSynergy(data = res4, method = c("ZIP", "HSA", "Bliss", "Loewe"), 
                         Emin = NA, Emax = NA, correct_baseline = "part")
res5 <- CalculateSynergy(data = res5, method = c("ZIP", "HSA", "Bliss", "Loewe"), 
                         Emin = NA, Emax = NA, correct_baseline = "part")
res6 <- CalculateSynergy(data = res6, method = c("ZIP", "HSA", "Bliss", "Loewe"), 
                         Emin = NA, Emax = NA, correct_baseline = "part")

# save the results
saveRDS(res, "rem-favi.rds")
saveRDS(res2, "favi-mol.rds")
saveRDS(res3, "rem-mol.rds")
saveRDS(res4, "rem-riba.rds")
saveRDS(res5, "favi-riba.rds")
saveRDS(res6, "mol-riba.rds")


######### PLOT RESPONSE
# read in results
res <- readRDS("rem-favi.rds")
res2 <- readRDS("favi-mol.rds")
res3 <- readRDS("rem-mol.rds")
res4 <- readRDS("rem-riba.rds")
res5 <- readRDS("favi-riba.rds")
res6 <- readRDS("mol-riba.rds")

# plot the response-heatmap
library(ggplot2)
library(RColorBrewer)

# plot with ggplot
plot.response <- function(res){
  # prepare the dataset before plotting
  df <- data.frame(matrix(nrow = length(res$response_statistics$response_mean), ncol = 7))
  colnames(df) <- c("drug1","drug2","unit1","unit2","conc1","conc2","response")
  df$drug1 <- res$drug_pairs$drug1
  df$drug2 <- res$drug_pairs$drug2
  df$unit1 <- res$drug_pairs$conc_unit1
  df$unit2 <- res$drug_pairs$conc_unit2
  df$conc1 <- res$response_statistics$conc1
  df$conc2 <- res$response_statistics$conc2
  df$response <- round(res$response_statistics$response_mean)
  df$x <- 0
  df$y <- 0
  # give the concentrations a new label for plot
  for(i in 1:nrow(df)){
    for(j in 1:length(unique(df$conc1))){
      if(df[i,]$conc1==unique(df$conc1)[j]){df[i,]$x <- j}
    }
  }
  for(i in 1:nrow(df)){
    for(j in 1:length(unique(df$conc2))){
      if(df[i,]$conc2==unique(df$conc2)[j]){df[i,]$y <- j}
    }
  }
  
  library(ggplot2)
  ggplot(df,aes(x=x, y=y, fill= response)) + # plot content
    # theme_classic()+
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          text = element_text(colour = "black",size = 15),
          axis.title.x = element_text(size = 15),    # X-axis title text size
          axis.title.y = element_text(size = 15),    # Y-axis title text size
          axis.text.x = element_text(size = 15),     # X-axis tick text size
          axis.text.y = element_text(size = 15),     # Y-axis tick text size
          legend.title = element_text(size = 14),
          legend.position = "bottom")+ # custom the plot
    labs(x = paste0(df$drug1," (",df$unit1,")"), y = paste0(df$drug2," (",df$unit2,")"),
         fill = "Inhibition") + # change labels of the plots
    geom_tile() +
    # geom_text(aes(label = response), color = "black", size = 4) +
    scale_fill_gradient(low = "#f8d86a", high = "#757CBB",
                        limits = c(0, 100), oob = scales::squish) + # cunstom legend
    scale_x_continuous(breaks = c(1:length(unique(df$conc1))), labels = c(sort(round(unique(df$conc1),1))))+
    scale_y_continuous(breaks = c(1:length(unique(df$conc2))), labels = c(sort(round(unique(df$conc2),1))))
}

library(ggpubr)
p1 <- plot.response(res)
p1
p2 <- plot.response(res2)
p2
p3 <- plot.response(res3)
p3
p4 <- plot.response(res4)
p4
p5 <- plot.response(res5)
p5
p6 <- plot.response(res6)
p6

# plot the plots one in a seperate page
pdf("PIV_response_plots_text.pdf",width = 5.4, height = 4.7)
p1
p2
p3
p4
p5
p6
dev.off()

# combine the plots together
pdf("response_plots.pdf", width = 3, height = 18)
ggarrange(p1,p2,p3,p4,p5,p6,# objects to combine together
          ncol = 1, # plot in one column
          common.legend = TRUE, legend = "bottom" # shared legend at bottom
)%>%
  annotate_figure(
    top = text_grob("Dose response", size = 15, face = "bold")
  ) # add figure plot
dev.off()

pl_HSA <- function(res){
  p <- Plot2DrugSurface(data = res, plot_block = 1, drugs = c(1, 2), # select data and drug
                        plot_value = "HSA_synergy", # Specifies the synergy model 
                        dynamic = F, # A static plot will be generated.
                        high_value_color = "#757CBB", # choose the colour for high value
                        low_value_color = "#f8d86a", # choose the colour for low value
                        color_range = c(-100,100), # Sets the color scale range
                        z_range= c(-50,50), # Sets the vertical axis (z-axis) range for the 3D surface plot.
                        summary_statistic = NULL, # select summary statistic in title, e.g.c("mean",  "median")
                        plot_title = NULL, # plot title
                        text_size_scale = 1) # text size
}

# pl_HSA(res)
#to create and save as pdf: 
pdf("PIV_score_plots.pdf", width = 3.6, height = 3.6)
pl_HSA(res)
pl_HSA(res2)
pl_HSA(res3)
pl_HSA(res4)
pl_HSA(res5)
pl_HSA(res6)
dev.off()

############cytotoxic 3D plot ############
dfcyto <- extract.2combtherapy(mergedcsv="cytotoxic_merged.csv", drugno_A = 1, drugno_B = 2, ConcUnit = "mg/L")
dfcyto2 <- extract.2combtherapy(mergedcsv="cytotoxic_merged.csv", drugno_A = 2, drugno_B = 3, ConcUnit = "mg/L")
dfcyto3 <- extract.2combtherapy(mergedcsv="cytotoxic_merged.csv", drugno_A = 1, drugno_B = 3, ConcUnit = "mg/L")
dfcyto4 <- extract.2combtherapy(mergedcsv="cytotoxic_merged.csv", drugno_A = 1, drugno_B = 4, ConcUnit = "mg/L")
dfcyto5 <- extract.2combtherapy(mergedcsv="cytotoxic_merged.csv", drugno_A = 2, drugno_B = 4, ConcUnit = "mg/L")
dfcyto6 <- extract.2combtherapy(mergedcsv="cytotoxic_merged.csv", drugno_A = 3, drugno_B = 4, ConcUnit = "mg/L")

res_cyto <- ReshapeData(data = dfcyto, data_type = "inhibition", impute = F,impute_method = NULL,noise = F)
res_cyto2 <- ReshapeData(data = dfcyto2, data_type = "inhibition", impute = F,impute_method = NULL,noise = F)
res_cyto3 <- ReshapeData(data = dfcyto3, data_type = "inhibition", impute = F,impute_method = NULL,noise = F)
res_cyto4 <- ReshapeData(data = dfcyto4, data_type = "inhibition", impute = F,impute_method = NULL,noise = F)
res_cyto5 <- ReshapeData(data = dfcyto5, data_type = "inhibition", impute = F,impute_method = NULL,noise = F)
res_cyto6 <- ReshapeData(data = dfcyto6, data_type = "inhibition", impute = F,impute_method = NULL,noise = F)

# check the missing value in raw data, usually the res_cytoult without drug treatment
library(lattice)
pl_cyto <- function(res, main = NULL){
  Plot2DrugSurface(data = res, plot_block = 1, drugs = c(1, 2),
                          plot_value = "response", dynamic = F,
                          high_value_color = "#ED7170",
                          low_value_color = "#757CBB",
                          color_range = c(-100,100),
                          z_range= c(-100,100),
                          summary_statistic = NULL,
                          plot_title = NULL,
                          text_size_scale = 1)
}

#to create and save as pdf: 
pdf("PIV_cyto_plots.pdf", width = 4, height = 4)
pl_cyto(res = res_cyto)
pl_cyto(res_cyto2)
pl_cyto(res_cyto3)
pl_cyto(res_cyto4)
pl_cyto(res_cyto5)
# add a legend below the last plot
pl_cyto(res_cyto6)
dev.off()

