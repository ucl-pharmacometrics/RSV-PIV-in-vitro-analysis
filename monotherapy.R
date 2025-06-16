rm(list=ls())
getwd()
####### STEP 1 Prepare the formatted data #########
# The following code uses the output from STEP 1 data formatting (see e.g. merged.csv)
# The formatting of cytotoxicity keeps the same as viral results
# the unit of concentration in the file is mg/L
# If just only replot is needed, skip to STEP4
####### function to analyse monotherapy data #########
#   For monotherapies we will only extract the data for each drug where it was given without
#   any of the other drugs
extract.monotherapy <- function(mergedcsv, drugname){
  # mergedcsv = "merged.csv"
  # drugname = "Remdesivir" 
  platedf <- read.csv(mergedcsv)
  # match drugname with drugno
  for (i in 1:4){
    if(platedf[paste0("drug", i, "_name")][1,] == drugname){drugno = i}
  }
  # work out which of the 4 drugs the user has NOT chosen
  exclude <- which(c(1:4) != drugno)
  # select nonzero drug wells
  platedf <- platedf[platedf[paste0("drug", drugno, "_conc")] != 0, ]
  # remove nonzero combination wells
  platedf <- platedf[platedf[paste0("drug", exclude[1], "_conc")] == 0, ]
  platedf <- platedf[platedf[paste0("drug", exclude[2], "_conc")] == 0, ]
  platedf <- platedf[platedf[paste0("drug", exclude[3], "_conc")] == 0, ]
  # make a clean analysis dataframe
  analysisdf <- data.frame(matrix(nrow = nrow(platedf), ncol = 7))
  colnames(analysisdf) <- c("drugno", "drug", "conc_uM", "conc","viral_percent_inh", 
                            "mean_percent_inh", "sem_percent_inh" ) 
  # add the concentration and inhibition results 
  analysisdf$conc <-  platedf[[paste0("drug", drugno, "_conc")]]
  analysisdf$viral_percent_inh <-  platedf$viral_percent_inh
  #define standard error of mean function
  std.error <- function(x) sd(x)/sqrt(length(x))
  # calculate the mean and sem for each concentration
  analysisdf$mean_percent_inh <- -999
  analysisdf$sem_percent_inh <- -999
  for(i in unique(analysisdf$conc)){
    analysisdf$mean_percent_inh[analysisdf$conc == i] <- mean(analysisdf$viral_percent_inh[analysisdf$conc == i])
    analysisdf$sem_percent_inh[analysisdf$conc == i] <- std.error(analysisdf$viral_percent_inh[analysisdf$conc == i]) 
  }
  # include other information
  analysisdf$drug <- platedf[[paste0("drug", drugno, "_name")]]
  analysisdf$drugno <- drugno
  analysisdf
}

####### STEP 2 format the data per drug######
# for each of 1-4 drugs makes a dataset of just the monotherapy data
# use the extract.monotherapy function to format the data:
# mergedcsv - the name of the formatted file

#### Monotherapy results
df1 <- extract.monotherapy(mergedcsv = "merged.csv", drugname = "Remdesivir")
df2 <- extract.monotherapy(mergedcsv = "merged.csv", drugname = "Favipiravir")
df3 <- extract.monotherapy(mergedcsv = "merged.csv", drugname = "EIDD-1931")
df4 <- extract.monotherapy(mergedcsv = "merged.csv", drugname = "Ribavirin")
# combine the 4 extracted monotherapy data
dfmono <- rbind(df1, df2, df3, df4)
# save the new monotherapy dataset
write.csv(dfmono, "monotherapy_data.csv", row.names = F)

#### cytotoxicity results
df1.c <- extract.monotherapy(mergedcsv = "cytotoxic_merged.csv", drugname = "Remdesivir")
#df1.c <- df1[df1$conc != 32.00,]
df2.c <- extract.monotherapy(mergedcsv = "cytotoxic_merged.csv", drugname = "Favipiravir")
df3.c <- extract.monotherapy(mergedcsv = "cytotoxic_merged.csv", drugname = "EIDD-1931")
df4.c <- extract.monotherapy(mergedcsv = "cytotoxic_merged.csv", drugname = "Ribavirin")

dfcyto <- rbind(df1.c, df2.c, df3.c, df4.c)
# save the new monotherapy dataset
write.csv(dfcyto, "cytotoxicity_monotherapy_data.csv", row.names = F)

######## plot raw monotherapy data ########
plot.raw <- function(dfmono, drugname){
  library(ggplot2)
  # select the data with specified drugname
  df.v <- dfmono[dfmono$drug == drugname, ]
  pl.drug <- ggplot() +
    theme_bw() +
    scale_x_log10() +
    labs(x = "Drug concentration (mg/L)", y = "% viral inhibition",
         title = drugname) +
    geom_point(aes(x = conc, y = viral_percent_inh), data = df.v) +
    geom_line(aes(x = conc, y = mean_percent_inh), data = df.v)+
    geom_errorbar(aes(x = conc, y = mean_percent_inh, ymin = mean_percent_inh - sem_percent_inh,
                      ymax = mean_percent_inh + sem_percent_inh), data = df.v, width=0.1) + 
    ylim(-50, 150) #ylim sets y scale limits so can make the different drugs on same scale, but be careful of missing off extreme data
  pl.drug
}
# plots the data
# dfmono - monotherapy dataset defined above
pl.1 <- plot.raw(dfmono, "Remdesivir")
pl.1
pl.2 <- plot.raw(dfmono, "Favipiravir")
pl.2
pl.3 <- plot.raw(dfmono, "EIDD-1931")
pl.3
pl.4 <- plot.raw(dfmono, "Ribavirin")
pl.4

# combine the plots together
library(ggpubr)
ggarrange(pl.1,pl.2,pl.3,pl.4, ncol = 2, nrow = 2)

####### STEP 3 estimate EC50 ########
# This function estimates the EC50 and optionally other parameters
# Model parameters:
# - EC50: estimated
# - gamma: the Hill slope, can be fixed to 1/estimated
# - Emax: can be fixed to 100/estimated
# - E0: can be fixed to 0/estimated

# Different combinations of fixed/estimated parameters (m1-8) were tested in the function below
est.hill <- function(dfmono, drugname, gam){
  # specify the 8 models:
  m1 <- y ~  100 * x / (EC50 + x)
  m2 <- y ~  100 * x ^ gam / (EC50 ^ gam + x ^ gam)
  m3 <- y ~  Emax * x / (EC50  + x)
  m4 <- y ~  E0 + 100 * x / (EC50  + x)
  m5 <- y ~  E0 + Emax * x / (EC50  + x)
  m6 <- y ~  Emax * x ^ gam / (EC50 ^ gam  + x ^gam)
  m7 <- y ~  E0 + 100 * x ^ gam / (EC50 ^ gam  + x ^ gam)
  m8 <- y ~  E0 + Emax * x ^ gam / (EC50 ^ gam  + x ^ gam)
  
  df <- dfmono[dfmono$drug == drugname, ]
  # 1. take EC50 initial est as half of the conc range
  EC50_init1 <- median(df$conc)
  x <- df$conc
  y <- df$viral_percent_inh
  r1 <- try(nls(formula = m1, start = c(EC50 = EC50_init1)),
            silent = TRUE)
  # update the initial estimate for model 2
  EC50_init2 <- summary(r1)$coefficients[1]
  r2 <- try(nls(formula = m2, start = c(EC50 = EC50_init2, gam = gam), control = nls.control(maxiter = 500)),
            silent = TRUE)
  r3 <- try(nls(formula = m3, start = c(EC50 = EC50_init2, Emax = 100), control = nls.control(maxiter = 500)),
            silent = TRUE)
  r4 <- try(nls(formula = m4, start = c(EC50 = EC50_init2, E0 = 0), control = nls.control(maxiter = 500)),
            silent = TRUE)
  r5 <- try(nls(formula = m5, start = c(EC50 = EC50_init2, E0 = 0, Emax = 100), control = nls.control(maxiter = 500)),
            silent = TRUE)
  r6 <- try(nls(formula = m6, start = c(EC50 = EC50_init2, Emax = 100, gam = gam), control = nls.control(maxiter = 500)),
            silent = TRUE)
  r7 <- try(nls(formula = m7, start = c(EC50 = EC50_init2, E0 = 0, gam = gam), control = nls.control(maxiter = 500)),
            silent = TRUE)
  r8 <- try(nls(formula = m8, start = c(EC50 = EC50_init2, Emax = 100, E0 = 0, gam = gam)),
            silent = TRUE)
  library(nlstools)
  # make a dataframe for the output
  out <- data.frame(matrix(nrow = 8, ncol = 10))
  colnames(out) <- c("drug","model", "AIC", "EC50","EC50_se","EC50_low", "EC50_high",
                     "gam","Emax","E0")
  out$drug <- df[1,]$drug
  out$model <- c("m1", "m2", "m3", "m4", "m5", "m6", "m7", "m8")
  out$AIC[1] <- ifelse(is.numeric(try(AIC(r1), silent = TRUE)) == TRUE, AIC(r1), "fail")
  out$AIC[2] <- ifelse(is.numeric(try(AIC(r2), silent = TRUE)) == TRUE, AIC(r2), "fail")
  out$AIC[3] <- ifelse(is.numeric(try(AIC(r3), silent = TRUE)) == TRUE, AIC(r3), "fail")
  out$AIC[4] <- ifelse(is.numeric(try(AIC(r4), silent = TRUE)) == TRUE, AIC(r4), "fail")
  out$AIC[5] <- ifelse(is.numeric(try(AIC(r5), silent = TRUE)) == TRUE, AIC(r5), "fail")
  out$AIC[6] <- ifelse(is.numeric(try(AIC(r6), silent = TRUE)) == TRUE, AIC(r6), "fail")
  out$AIC[7] <- ifelse(is.numeric(try(AIC(r7), silent = TRUE)) == TRUE, AIC(r7), "fail")
  out$AIC[8] <- ifelse(is.numeric(try(AIC(r8), silent = TRUE)) == TRUE, AIC(r8), "fail")
  # EC50, with se, 95%CI(low,high)
  out$EC50[1] <- ifelse(out$AIC[1] == "fail", "fail", summary(r1)$coefficients[1])
  out$EC50_se[1] <- ifelse(out$AIC[1] == "fail", "fail", summary(r1)$coefficients[1,2])
  out$EC50_low[1] <- ifelse(is.numeric(try(confint2(r1), silent = TRUE)) == TRUE, confint2(r1)[1], "fail")
  out$EC50_high[1] <- ifelse(is.numeric(try(confint2(r1), silent = TRUE)) == TRUE, confint2(r1)[2], "fail")
  out$EC50[2] <- ifelse(out$AIC[2] == "fail", "fail", summary(r2)$coefficients[1])
  out$EC50_se[2] <- ifelse(out$AIC[2] == "fail", "fail", summary(r2)$coefficients[1,2])
  out$EC50_low[2] <- ifelse(is.numeric(try(confint2(r2), silent = TRUE)) == TRUE, confint2(r2)[1], "fail")
  out$EC50_high[2] <- ifelse(is.numeric(try(confint2(r2), silent = TRUE)) == TRUE, confint2(r2)[3], "fail")
  out$EC50[3] <- ifelse(out$AIC[3] == "fail", "fail", summary(r3)$coefficients[1])
  out$EC50_se[3] <- ifelse(out$AIC[3] == "fail", "fail", summary(r3)$coefficients[1,2])
  out$EC50_low[3] <- ifelse(is.numeric(try(confint2(r3), silent = TRUE)) == TRUE, confint2(r3)[1], "fail")
  out$EC50_high[3] <- ifelse(is.numeric(try(confint2(r3), silent = TRUE)) == TRUE, confint2(r3)[3], "fail")
  out$EC50[4] <- ifelse(out$AIC[4] == "fail", "fail", summary(r4)$coefficients[1])
  out$EC50_se[4] <- ifelse(out$AIC[4] == "fail", "fail", summary(r4)$coefficients[1,2])
  out$EC50_low[4] <- ifelse(is.numeric(try(confint2(r4), silent = TRUE)) == TRUE, confint2(r4)[1], "fail")
  out$EC50_high[4] <- ifelse(is.numeric(try(confint2(r4), silent = TRUE)) == TRUE, confint2(r4)[3], "fail")
  out$EC50[5] <- ifelse(out$AIC[5] == "fail", "fail", summary(r5)$coefficients[1])
  out$EC50_se[5] <- ifelse(out$AIC[5] == "fail", "fail", summary(r5)$coefficients[1,2])
  out$EC50_low[5] <- ifelse(is.numeric(try(confint2(r5), silent = TRUE)) == TRUE, confint2(r5)[1], "fail")
  out$EC50_high[5] <- ifelse(is.numeric(try(confint2(r5), silent = TRUE)) == TRUE, confint2(r5)[4], "fail")
  out$EC50[6] <- ifelse(out$AIC[6] == "fail", "fail", summary(r6)$coefficients[1])
  out$EC50_se[6] <- ifelse(out$AIC[6] == "fail", "fail", summary(r6)$coefficients[1,2])
  out$EC50_low[6] <- ifelse(is.numeric(try(confint2(r6), silent = TRUE)) == TRUE, confint2(r6)[1], "fail")
  out$EC50_high[6] <- ifelse(is.numeric(try(confint2(r6), silent = TRUE)) == TRUE, confint2(r6)[4], "fail")
  out$EC50[7] <- ifelse(out$AIC[7] == "fail", "fail", summary(r7)$coefficients[1])
  out$EC50_se[7] <- ifelse(out$AIC[7] == "fail", "fail", summary(r7)$coefficients[1,2])
  out$EC50_low[7] <- ifelse(is.numeric(try(confint2(r7), silent = TRUE)) == TRUE, confint2(r7)[1], "fail")
  out$EC50_high[7] <- ifelse(is.numeric(try(confint2(r7), silent = TRUE)) == TRUE, confint2(r7)[4], "fail")
  out$EC50[8] <- ifelse(out$AIC[8] == "fail", "fail", summary(r8)$coefficients[1])
  out$EC50_se[8] <- ifelse(out$AIC[8] == "fail", "fail", summary(r8)$coefficients[1,2])
  out$EC50_low[8] <- ifelse(is.numeric(try(confint2(r8), silent = TRUE)) == TRUE, confint2(r8)[1], "fail")
  out$EC50_high[8] <- ifelse(is.numeric(try(confint2(r8), silent = TRUE)) == TRUE, confint2(r8)[5], "fail")
  # gam
  out$gam[1] <- ifelse(out$AIC[1] == "fail", "fail", 1)
  out$gam[2] <- ifelse(out$AIC[2] == "fail", "fail", summary(r2)$coefficients[2])
  out$gam[3] <- ifelse(out$AIC[3] == "fail", "fail", 1)
  out$gam[4] <- ifelse(out$AIC[4] == "fail", "fail", 1)
  out$gam[5] <- ifelse(out$AIC[5] == "fail", "fail", 1)
  out$gam[6] <- ifelse(out$AIC[6] == "fail", "fail", summary(r6)$coefficients[3])
  out$gam[7] <- ifelse(out$AIC[7] == "fail", "fail", summary(r7)$coefficients[3])
  out$gam[8] <- ifelse(out$AIC[8] == "fail", "fail", summary(r8)$coefficients[4])
  
  # Emax
  out$Emax[1] <- ifelse(out$AIC[1] == "fail", "fail", 100)
  out$Emax[2] <- ifelse(out$AIC[2] == "fail", "fail", 100)
  out$Emax[3] <- ifelse(out$AIC[3] == "fail", "fail", summary(r3)$coefficients[2])
  out$Emax[4] <- ifelse(out$AIC[4] == "fail", "fail", 100)
  out$Emax[5] <- ifelse(out$AIC[5] == "fail", "fail", summary(r5)$coefficients[3])
  out$Emax[6] <- ifelse(out$AIC[6] == "fail", "fail", summary(r6)$coefficients[2])
  out$Emax[7] <- ifelse(out$AIC[7] == "fail", "fail", 100)
  out$Emax[8] <- ifelse(out$AIC[8] == "fail", "fail", summary(r8)$coefficients[2])
  
  # E0
  out$E0[1] <- ifelse(out$AIC[1] == "fail", "fail", 0)
  out$E0[2] <- ifelse(out$AIC[2] == "fail", "fail", 0)
  out$E0[3] <- ifelse(out$AIC[3] == "fail", "fail", 0)
  out$E0[4] <- ifelse(out$AIC[4] == "fail", "fail", summary(r4)$coefficients[2])
  out$E0[5] <- ifelse(out$AIC[5] == "fail", "fail", summary(r5)$coefficients[2])
  out$E0[6] <- ifelse(out$AIC[6] == "fail", "fail", 0)
  out$E0[7] <- ifelse(out$AIC[7] == "fail", "fail", summary(r7)$coefficients[2])
  out$E0[8] <- ifelse(out$AIC[8] == "fail", "fail", summary(r8)$coefficients[3])
  out
}

##### tables for esthill ####
#  Now we can apply the est.hill function for each drug and look at the output.
# dfmono - the monotherapy data we just define in STEP 2
# gam - initial estimate for the Hill coefficient, bigger number leads to steeper slope
#       usually start with 1

table1 = est.hill(dfmono = dfmono, drugname = "Remdesivir", gam = 5)
table1

# Output: xx_high, xx_low - 95% confidence interval of the parameter,
#         EC50_se - standard error of the EC50,
#         AIC - The best fitting model has the lowest AIC ,
#               BUT we need to check if the parameters make sense.
#  
#  Only fit the model to data that makes sense. e.g. Emax around 100
#  Extract the best fit model estimations: lowest AIC/ specific model
est1 <- table1[table1$AIC == min(table1$AIC),]
# est1 <- table1[table1$model == "m2",]

table2 = est.hill(dfmono = dfmono, drugname = "Favipiravir", gam = 1)
table2
est2 <- table2[table2$AIC == min(table2$AIC), ]

table3 = est.hill(dfmono = dfmono, drugname = "EIDD-1931", gam = 1)
table3
est3 <- table3[table3$AIC == min(table3$AIC), ]

table4 = est.hill(dfmono = dfmono, drugname = "Ribavirin", gam = 1)
table4
est4 <- table4[table4$AIC == min(table4$AIC), ]

# combine the estimated results together
param <- rbind(table1, table2, table3, table4)
est <- rbind(est1, est2, est3, est4)

# save the estimation for later use, no need to run this part each time
write.csv(param, "EC50_table.csv")
write.csv(est, "EC50_est.csv")

###### STEP 4 plot the results #####
pl.mono <- function(drugname, table, ymin = 0, ymax = 100, l = 0.03){
  # select the data by drugname
  dfpred <- dfmono[dfmono$drug == drugname, ]
  dftoxi <- dfcyto[dfcyto$drug == drugname, ]
  dftoxi$drug <- "Cytotoxicity"
  dfplot <- rbind(dfpred, dftoxi)
  # use the drug name from the dataset as title for the plot
  drugname <- dfpred$drug[1]
  #  Generate a range of concentrations starting lower and ending higher than those observed
  cmin <- min(dfpred$conc) 
  cmax <- max(dfpred$conc) 
  cpreds <- seq(cmin/2, cmax*2, 0.1)
  crange <- log2(cmax/cmin)
  # select the model with lowest AIC
  estparam <- table[table$drug == drugname,]
  EC50 <- as.numeric(estparam$EC50)
  E0 <- as.numeric(estparam$E0)
  Emax <- as.numeric(estparam$Emax)
  gam <- as.numeric(estparam$gam)
  # make a vector of predictions
  hill_pred <- E0 + Emax * cpreds ^ gam / (EC50 ^ gam + cpreds ^ gam)
  # prepare the annotation
  ec50 <- round(as.numeric(estparam$EC50),2)
  
  library(ggplot2)
  library(ggprism)
  pl.fit <- ggplot() +
    theme_prism() + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
          text = element_text(size = 15, face = "bold"),
          axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "none") + 
    labs(y = "% viral inhibition",title = drugname) +
    geom_point(aes(x = conc, y = mean_percent_inh, colour = drug), data = dfplot, size=1.5)+
    geom_line(aes(x = conc, y = mean_percent_inh), data = dftoxi, colour = "#ED7170", linewidth=1)+
    geom_line(aes(x = cpreds, y = hill_pred), colour = "#757CBB", linewidth = 1)+
    geom_errorbar(aes(x = conc, y = mean_percent_inh, ymin = mean_percent_inh - sem_percent_inh, 
                      ymax = mean_percent_inh + sem_percent_inh, colour = drug), data = dfplot, width=.1) +
    scale_colour_manual(name = NULL, values = c("#ED7170", "#757CBB"))+
    scale_x_continuous(trans = 'log2',
                       breaks = sort(unique(dfpred$conc)),
                       labels = ~paste(round(sort(unique(dfpred$conc)),2)),
                       name = "Drug concentration(mg/L)")+
    scale_y_continuous(breaks = seq(0,100,50), limits = c(ymin, ymax))+
    geom_hline(aes(yintercept = 50), colour="grey", linetype="dashed")+
    annotate("text", x = l*crange, y = 55, size = 5,
             label = paste0("EC50 = ",round(as.numeric(estparam$EC50),1),"mg/L"))

}

# read in the results from STEP 3
est <- read.csv("EC50_est.csv")
dfmono <- read.csv("monotherapy_data.csv")
dfcyto <- read.csv("cytotoxicity_monotherapy_data.csv")

# table - table of previously selected best fit model parameters
# ymin, ymax - y-axis scale, default 0-100
# l - EC50 label position, scaled by concentration range
pl.1 <- pl.mono(drugname = "Remdesivir", table = est, ymin = -10 , ymax = 110, l = 0.035)
pl.1
pl.2 <- pl.mono(drugname = "Favipiravir", table = est, ymin = -10 , ymax = 110, l = 0.3)
pl.2
pl.3 <- pl.mono(drugname = "EIDD-1931", table = est, ymin = -10 , ymax = 110, l = 0.018)
pl.3
pl.4 <- pl.mono(drugname = "Ribavirin", table = est, ymin = -10 , ymax = 110, l = 0.07)
pl.4

pdf("monotherapy_plots.pdf",width = 4, height = 4)
pl.1
pl.2
pl.3
pl.4
dev.off()

###### STEP 5 Output the EC50 results #####
# read in parameter estimation table
tab1 <- read.csv("EC50_est.csv")

df1 <- data.frame(matrix(nrow = 4, ncol = 4))
colnames(df1) <- c("virus","drug","EC50","95%CI")
df1$virus <- "PIV" 
df1$drug <- tab1$drug
df1$EC50 <- round(tab1$EC50,1)
df1$`95%CI` <- paste0("(",round(tab1$EC50_low,1),",",round(tab1$EC50_high,1),")")

# df <- rbind(df1, df2)
write.csv(df1, "EC50_ci_est.csv", row.names = F)