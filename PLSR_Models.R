### Updated on 04.04.2025
### Written by Nargiz Rüter

### PLSR for multi year IS data ###
### The code contains loops to build PLSR model using randomly generated sample data for following trait: "Biom_wet"

#### Load libraries needed to run the script #####
list.of.packages <- c("pls","dplyr","reshape2","here","plotrix","ggplot2","gridExtra","spectratrait","reshape","agricolae",
                      "baseline","spectrolab","caret","tidyverse","readr","tidyr","ggrepel","directlabels","ggpubr",
                      "animation","pdftools","gghighlight","stringr","gtools","gsheet","GGally","rgdal","magick",
                      "pavo","mdatools","vip","pdp","magrittr","enpls","crayon", "highlight","seecolor")
invisible(lapply(list.of.packages, library, character.only = TRUE))
# install.packages("remotes")
# remotes::install_github("TESTgroup-BNL/spectratrait")

# Avoid scientific notation globally
options(scipen = 999)

# setwd("C:/Users/.../Desktop/...")

# [COLOR BLIND FRIENDLY PALETTE]
# palette_name <- paste("COLOR_BLIND_FRIENDLY_PALETTE.pdf", sep = "")
# pdf(palette_name,  width=15, height=10)
palette <- c( "#5566AA", "#117733", "#44AA66", "#55AA22", "#668822", "#99BB55", "#558877", "#88BBAA",
              "#AADDCC", "#44AA88", "#DDCC66", "#FFDD44", "#FFEE88", "#BB0011", "#000000", "#E69F00",
              "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

seecolor::print_color(c( "#5566AA", "#117733", "#44AA66", "#55AA22", "#668822", "#99BB55", "#558877", "#88BBAA",
                         "#AADDCC", "#44AA88", "#DDCC66", "#FFDD44", "#FFEE88", "#BB0011", "#000000", "#E69F00",
                         "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"), type = "r")

columns_of_int <- c("Biom_wet_g")

columns_of_interest <- columns_of_int[1]
    
g=1
REF_DF <- read.csv("./r_model_input/Randomly_generated_sample_data_structure.csv")

df <- as.data.frame(REF_DF)
df$nRow <- 1:nrow(df)
df <- df %>% relocate(nRow, .before = 1)
df$Unique_ID <- paste(df$nRow, df$SAMPLE_ID, sep = "_")
df <- df %>% relocate(Unique_ID, .before = 1)
df$Year <- REF_DF$Year[g]
df <- df %>% relocate(Year, .before = 1)
df$Unique_ID <- paste(df$Unique_ID, df$Year, sep = "_")

# Assigning wavelength names into wvl variable 
names(df)[which(names(df) == "X410"):length(df)]
wvl <- names(df)[which(names(df) == "X410"):length(df)][grep("X", names(df)[which(names(df) == "X410"):length(df)])] # wavelengths

# Optimal number of iterations per year is 1000
inIter <- 10 #1000    #inIter <- 10

inBands <- wvl

# Naming the output folder according to the years
inYear <- REF_DF$Year[g]

inProcs <- REF_DF$Processing[g]

for (col_name in columns_of_interest) { # col_name = columns_of_interest[1]
  
  if (col_name %in% colnames(df)) {
    
    # Variable(s) to be predicted, example: Wet Biomass
    #inVar <- names(df)[which(names(df) == "Biom_wet")]
    #inVar <- "Biom_wet"
    inVar <- col_name
    
    inYear <- REF_DF$Year[g]
    
    inProcs <- REF_DF$Processing[g]
    
    # Keep only complete rows of inVar and spec data before fitting
    df<-df[complete.cases(df[,(colnames(df) %in% inVar)]),] 
    df<-df[complete.cases(df[,(colnames(df) %in% inBands)]),]  
    
    # it creates a matrix to see which model had what specific ncomps, R2, rmsep etc.
    nComps_iter_df <- data.frame(Year = character(inIter), #inVarTrait = character(inIter), #Model = numeric(inIter), 
                                 nComps = numeric(inIter), R2 = numeric(inIter), RMSEP = numeric(inIter), PRESS = numeric(inIter))
    
    # It creates an empty dataframe to save all coefficients from all inIter models
    coefs_iter_df <- data.frame(matrix(ncol = length(wvl)+1 , nrow = inIter))
    names(coefs_iter_df)[1] <- "Intercept"
    names(coefs_iter_df)[2:(length(wvl)+1)] <- as.numeric(gsub('X', '', wvl))
    
    # It creates an empty dataframe to save all vips from all inIter models
    vips_iter_df <- data.frame(matrix(ncol = length(wvl), nrow = inIter))
    names(vips_iter_df) <- as.numeric(gsub('X', '', wvl))
    
    # It creates an empty dataframe to save all X loadings from all inIter models' selected component
    loads_iter_df <- data.frame(matrix(ncol = length(wvl), nrow = inIter))
    names(loads_iter_df) <- as.numeric(gsub('X', '', wvl))
    
    # Create Output Directory
    dir.create(paste("./r_model_output/", "PPR1_", inYear, "_", inProcs, sep = ""), recursive = T) #Single
    pati_main <- paste("./r_model_output/", "PPR1_", inYear, "_", inProcs, sep = "")
    pati <- paste(pati_main, "/", "PPR1_", inYear, "_", inProcs, sep = "")
    
    ##### Start of Calculations ####
    set.seed(6955866) 
    for (iter in seq(inIter)){ #iter=1
      
      inYear <- REF_DF$Year[g]
      
      inProcs <- REF_DF$Processing[g]
      
      print(paste0("##################"))
      print(paste0("  Model --->  ", iter ,"  <---"))
      
      # Splitting data into Validation (test) and Calibration (train) data # inIter times partitioning (cal.data.1; val.data.1) and inIter times model
      data.source <- spectratrait::create_data_split(dataset=df, approach="dplyr", split_seed=(iter+pi*(88033)), prop=0.8)
      cal.data <- data.source$cal_data # CAL calibration data
      assign(paste0("cal.data", iter), cal.data)
      #cal.data[iter] <- cal.data
      val.data <- data.source$val_data # VAL external validation data
      assign(paste0("val.data", iter), val.data)
      #cal.data[iter] <- cal.data
      
      # Structuring formula to be applied on all bands
      form <- paste(inBands,collapse="+")
      form <- as.formula(paste(inVar,"~",form))
      
      ##### Part 1: #####
      #Find number of components 
      set.seed(iter+pi*(69558))
      nComps <- if (0.8 * as.numeric(dim(df)[1]) <= 25) {
        round((as.numeric(dim(df)[1])* 0.8) * 0.8, 0)
      } else {
        25
      } #25 
      nCompsi <- if (0.8 * as.numeric(dim(df)[1]) <= 25) {
        round((as.numeric(dim(df)[1])* 0.8) * 0.8, 0)
      } else {
        25
      } #25
      i.no <- 1 # number of iterations
      
      outMatR2 <- matrix(data=NA,nrow=i.no,ncol=nComps)
      outMatRMSEP <- matrix(data=NA,nrow=i.no,ncol=nComps)
      outMatPRESS <- matrix(data=NA,nrow=i.no,ncol=nComps)
      
      model.set <- plsr(form,data=cal.data,ncomp=nComps,
                        validation="LOO", 
                        method="oscorespls")
      resR2 <- pls::R2(model.set, intercept = F)[[1]]
      
      # Calculate the adjusted R²
      nw <- nrow(cal.data)  # Number of observations
      pw <- nComps  # Number of components used in the model
      
      adjusted_R2 <- 1 - ((1 - resR2) * (nw - 1) / (nw - pw - 1))
      
      # Print the adjusted R²
      #adjusted_R2
      
      # Check if all R2 values are negative in this iteration
      if (all(resR2 < 0)) {
        print("All R2 values are negative in this iteration")
        #break  # Move to the next iteration of the outer loop # this will stop further iterations
        next # this will go back to the next iteration
      }else{
        
        outMatR2[1,seq(model.set$validation$ncomp)] <- resR2
        resRMSEP <- as.numeric(RMSEP(model.set,estimate="CV",intercept=F)$val)
        outMatRMSEP[1,seq(model.set$validation$ncomp)] <-resRMSEP
        # Lower values of PRESS indicate better predictive power
        resPRESS <- as.vector(model.set$validation$PRESS)
        outMatPRESS[1,seq(model.set$validation$ncomp)] <-resPRESS

        summary(model.set)
        # explained variance per each component
        explvar(model.set)
        
        # Set up the plot area
        ### Here we decide the number of components by plotting average of R2 and PRESS of inIter model iterations.
        # Plot all inIter models' R2 & PRESS
        outMatPRESS_df_t <- as.data.frame(t(outMatPRESS))
        outMatPRESS_df_t$comps <- 1:nrow(outMatPRESS_df_t)
        outMatR2_df <- as.data.frame(t(outMatR2))
        outMatR2_df$comps <- 1:nrow(outMatR2_df)
        outMatPRESS_df_t <- merge(outMatPRESS_df_t,outMatR2_df, by="comps" )
        first_neg_index <- which(outMatPRESS_df_t$V1.y < 0)[1]
        
        # If there is at least one negative value in V1.y, filter out rows after the first negative value
        if (!is.na(first_neg_index)) {
          outMatPRESS_df_t <- outMatPRESS_df_t[1:(first_neg_index - 1), ]
        } else {
          # If there are no negative values in V1.y, do nothing
        }
        
        # Find the row where V1.y is positive and V1.x is the minimum
        min_V1x_row <- outMatPRESS_df_t[outMatPRESS_df_t$V1.y > 0 & outMatPRESS_df_t$V1.x == min(outMatPRESS_df_t$V1.x), ]
        
        # Check if there are rows that meet the criteria
        if (nrow(min_V1x_row) > 0) {
          nComps <- min_V1x_row$comps[1]
        } else {
          # Handle the case where no rows meet the criteria
          cat("No rows found where V1.y is positive and V1.x is minimum.\n")
          
          inIter <- inIter - 1
          # In case of no match, proceed with the next iteration
          next
        }
        
        max_R2_index <- which.max(outMatR2_df$V1)
        max_comps <- outMatR2_df$comps[max_R2_index]
        
        max_R2 <- round(as.numeric(max(outMatR2_df$V1)),2)
        min_R2 <- round(as.numeric(min(outMatR2_df$V1)),2)
        min_PRESS <- round(as.numeric(min(outMatPRESS_df_t$V1.x)),2)
        max_PRESS <- round(as.numeric(max(outMatPRESS_df_t$V1.x)),2)
        
        # save each nComps in R session
        nComps <- nComps
        #assign(paste0("nComps", iter), nComps)
        
        # Save Some Statistics at each iteration into a dataframe
        nComps_iter_df[iter, "nComps"] <- nComps
        nComps_iter_df[iter, "R2"] <- max_R2
        nComps_iter_df[iter, "RMSEP"] <- as.numeric(min(outMatRMSEP))
        nComps_iter_df[iter, "PRESS"] <- min_PRESS
        nComps_iter_df[iter, "Model"] <- iter
        nComps_iter_df[iter, "Year"] <- inYear
        nComps_iter_df[iter, "Processing"] <- inProcs

        nComps_iter_df <- na.omit(nComps_iter_df)
        
        ##### Part 2: #####
        set.seed(6955866)
        
        model.pro <- plsr(form,data=cal.data,ncomp=nComps,
                          validation= "LOO", 
                          method="oscorespls")
        assign(paste0("model.pro", iter), model.pro)
        
        summary(model.pro)
        # explained variance per each component
        explvar(model.pro)
        
        ### Model Fit
        pred_cal_df <- data.frame(model.pro$fitted.values[,,nComps]) #11 #9
        pred_cal_df$measured <- model.pro$model[[inVar]]
        names(pred_cal_df)[1] <- "predicted"
        pred_cal_df <- pred_cal_df
        pred_cal_df$cal_data_Unique_ID <- cal.data$Unique_ID
        pred_cal_df$cal_data_measured <- cal.data[[inVar]]
        #pred_cal_df$nRow <- 1:nrow(pred_cal_df)
        assign(paste0("pred_cal_df", iter), pred_cal_df)
        
        ### Model applied on external Val Data
        pred_val <- predict(model.pro, newdata = val.data, ncomp = nComps)
        # Combine predicted and actual values in a data frame
        pred_val_df <- data.frame(pred_val)
        colnames(pred_val_df)[1] <- "predicted"
        pred_val_df$measured <- val.data[[inVar]]
        pred_val_df <- pred_val_df
        pred_val_df$val_data_Unique_ID <- val.data$Unique_ID
        pred_val_df$val_data_measured <- val.data[[inVar]]
        #pred_val_df$nRow <- 1:nrow(pred_val_df)
        assign(paste0("pred_val_df", iter), pred_val_df)
        
        # save each Coefficients in a dataframe
        coefs <- as.data.frame(t(as.data.frame(coef(model.set,ncomp=nComps,intercept=TRUE))))
        # save each VIPs in a dataframe
        vips <- as.data.frame(t(as.data.frame(spectratrait::VIP(model.set)[nComps,])))
        
        # save the most contributed bands to the nth component at each iteration
        x_loadings <- as.data.frame(model.set$loadings[,nComps])
        names(x_loadings) <- "x_load"
        x_loadings <- as.data.frame(t(x_loadings))
        
      }
        
      # add variables into the dataframe
      coefs_iter_df[iter, ] <- coefs 
      vips_iter_df[iter, ] <- vips 
      loads_iter_df[iter, ] <- x_loadings 
        
    } 
    
    # save inIter models' statistics before averaging
    write.csv(coefs_iter_df, paste0(pati_main,"/", inYear, "_", inProcs, "_", "REF_DF_MODELs-n-iters-COEFS.csv"))
    write.csv(vips_iter_df, paste0(pati_main,"/", inYear, "_", inProcs, "_", "REF_DF_MODELs-n-iters-VIPS.csv"))
    write.csv(nComps_iter_df, paste0(pati_main,"/", inYear, "_", inProcs, "_", "REF_DF_MODELs-n-iters-STATISTICS.csv"))
    write.csv(loads_iter_df, paste0(pati_main,"/", inYear, "_", inProcs, "_", "REF_DF_MODEL-n-iters-X_Loadings.csv"))
    
    
    ##### Part 3: #####
    #### Ensemble PLSR Models
    
    # Find the most commonly selected nComps among inIter models
    # Create a frequency table of the values in the nComps column
    freq_table <- table(nComps_iter_df$nComps)
    # Find the nComps with the highest frequency
    common_nComp <- as.numeric(names(freq_table)[which.max(freq_table)])
    
    # Step 1: Find data frame names starting with "pred_cal_df"
    # Step 2: Count the number of data frames found after excluding "pred_cal_df" because "pred_cal_df"="pred_cal_df1"
    # Step 3: Check the count and take action accordingly
    if(length(ls(pattern = "^pred_cal_df\\d*$")[-1])< 2){
      cat("Less than 2 data frames starting with 'pred_cal_df'. Breaking the loop or taking other action.\n")
      # You can break a loop or take another action here
    } else {
      cat("2 or more data frames starting with 'pred_cal_df' found.\n")
      # Continue with your loop or other operations
      
      
      #### THE HOLY MODEL ####
      ### Final Model Fit -> Averages of all predicted and measured of all models'
      
      ### Take average of predicted and measured from cal and cal datasets
      set.seed(6955866)
      
      # This ensembles iterations that were not empty
      ensemble_cal_list <- lapply(ls(pattern = "^pred_cal_df\\d*$")[-1],get)
      ensemble_val_list <- lapply(ls(pattern = "^pred_val_df\\d*$")[-1],get)
      
      ensemble_cal_pre <- do.call(rbind, ensemble_cal_list)
      ensemble_val_pre <- do.call(rbind, ensemble_val_list)
      
      ensemble_cal_preds <- ensemble_cal_pre %>%
        group_by(cal_data_Unique_ID) %>%
        summarize(
          predicted = mean(predicted),
          measured = mean(measured),
          #cal_data_ID = mean(cal_data_ID),
          cal_data_measured = mean(cal_data_measured))
      
      ensemble_val_preds <- ensemble_val_pre %>%
        group_by(val_data_Unique_ID) %>%
        summarize(
          predicted = mean(predicted),
          measured = mean(measured),
          #val_data_ID = mean(val_data_ID),
          val_data_measured = mean(val_data_measured))    
      
      # Check if field mesaurements are the same (it must be the same)
      # View(cbind(sort(unique(ensemble_val_preds$measured)),sort(unique(df[[inVar]]))))
      # View(cbind(sort(unique(ensemble_cal_preds$measured)),sort(unique(df[[inVar]]))))
      
      # Model Fit -> Cal Data
      mod_fit <- lm(predicted~measured, ensemble_cal_preds)
      # Print the summary of the linear regression model
      summary(mod_fit)
      summary(mod_fit)$adj.r.squared
      R2_plot_cal <- round(summary(mod_fit)$adj.r.squared, 2)
      
      # Calculate RMSEP for calibration data
      residuals_cal <- ensemble_cal_preds$predicted - ensemble_cal_preds$measured
      mse_cal <- mean(residuals_cal^2)
      rmsep_cal <- round(sqrt(mse_cal),1)
      #print(paste("RMSEP for Calibration Data:", round(rmsep_cal, 2)))
      CI_pred_cal <- predict(mod_fit, newdata = ensemble_cal_preds, interval = "confidence") #Confidence Interval
      
      # Apply the model on val.data
      # Predict validation data -> take average of predicted and measured
      mod_fit_val <- lm(predicted~measured, data = ensemble_val_preds)
      # Print the summary of the linear regression model
      summary(mod_fit_val)
      summary(mod_fit_val)$adj.r.squared
      R2_plot_val <- round(summary(mod_fit_val)$adj.r.squared, 2)
      
      # Calculate RMSEP for validation data
      residuals_val <- ensemble_val_preds$predicted - ensemble_val_preds$measured
      mse_val <- mean(residuals_val^2)
      rmsep_val <- round(sqrt(mse_val),1)
      #print(paste("RMSEP for Validation Data:", round(rmsep_val, 2)))
      CI_pred_val <- predict(mod_fit_val, newdata = ensemble_val_preds, interval = "confidence") #Confidence Interval
      
      ########################################################
      
      
      ##### Part 4: #####
      # Calculated out of inIter models' average
      # Avg_Coefficents
      coefs_iter_avg_df <- as.data.frame(colMeans(na.omit(coefs_iter_df)))
      names(coefs_iter_avg_df)[1] <- "Avg_Coefficents"
      coefs_iter_avg_df$Wavelength <- rownames(coefs_iter_avg_df) #there is intercept here, so it can't be as.numeric
      
      # Avg_VIPs
      vips_iter_avg_df <- as.data.frame(colMeans(na.omit(vips_iter_df))) 
      names(vips_iter_avg_df)[1] <- "Avg_VIPs"
      vips_iter_avg_df$Wavelength <- as.numeric(rownames(vips_iter_avg_df))
      
      # Avg_X_Loadings
      loads_iter_avg_df <- as.data.frame(colMeans(na.omit(loads_iter_df)))  
      names(loads_iter_avg_df)[1] <- "Avg_X_Loadings"
      loads_iter_avg_df$Wavelength <- as.numeric(rownames(loads_iter_avg_df))
      
      # Avg_Models_Statistics
      #na.omit(nComps_iter_df)
      nComps_iter_avg_df <- colMeans(na.omit(nComps_iter_df[ ,3:5])) #5:ncol(nComps_iter_df)
      nComps_iter_avg_df <- as.data.frame(t(nComps_iter_avg_df))
      nComps_iter_avg_df$Year <- inYear
      nComps_iter_avg_df$Processing <- inProcs
      nComps_iter_avg_df$Model <- "Average"
      nComps_iter_avg_df$nComps <- common_nComp
      nComps_iter_avg_df$PRESS <- round(nComps_iter_avg_df$PRESS,1)
      nComps_iter_avg_df$RMSEP <- round(nComps_iter_avg_df$RMSEP,1)
      nComps_iter_avg_df$R2 <- round(nComps_iter_avg_df$R2,2)
      nComps_iter_avg_df <- nComps_iter_avg_df[ ,c("Year", "Processing", "Model", "nComps", "PRESS", "RMSEP", "R2")]
      nComps_iter_avg_df$Cal_R2 <- round(summary(mod_fit)$adj.r.squared, 2)
      nComps_iter_avg_df$Val_R2 <- round(summary(mod_fit_val)$adj.r.squared, 2)
      nComps_iter_avg_df$Cal_RMSEP <- rmsep_cal
      nComps_iter_avg_df$Val_RMSEP <- rmsep_val
      nComps_iter_avg_df
      
      
      ##### Part 5: #####
      # The Holy model results as graphs etc.
      
      # Define the file path where variables will be saved
      save_file_path <- paste0(pati_main, "/", inYear, "_", inProcs, "_", Sys.Date(), "_plot_data.RData")
      # Save all relevant variables to an RData file
      save(ensemble_cal_preds, CI_pred_cal, mod_fit, R2_plot_cal, rmsep_cal,
           ensemble_val_preds, CI_pred_val, mod_fit_val, R2_plot_val, rmsep_val, pati_main, 
           inYear, inProcs, file = save_file_path)
      
      # Plot predicted vs actual values
      png(file = paste0(pati_main,"/", inYear, "_", inProcs, "_", "cal_prediction.png"), width = 1250, height = 850, res = 300)
      par(mar = c(3, 3, 0.4, 0.25)+0.1)#par(mar = c(bottom, left, top, right) + additional)
      #main=paste0(inYear," ",inProcs," ", inIter," Models' Cal. Data Prediction"),
      plot(predicted ~ measured, ensemble_cal_preds, pch=15, cex=1, xlab = "", ylab = "", xaxt = "n", yaxt = "n",
           xlim=c(0,2500), ylim=c(0,2000)) 
      # Add axes and labels
      axis(1, mgp=c(0, 0.45, 0), cex.axis=0.875, cex.lab = 0.375, font.lab = 2)
      axis(2, mgp=c(0, 0.45, 0), cex.axis=0.875, cex.lab = 0.375, font.lab = 2)
      mtext(expression("Measured [g" ~ m^{-2} ~ "]"), side=1, line=2.1, cex=1.25, font=2)
      mtext(expression("Predicted [g" ~ m^{-2} ~ "]"), side=2, line=1.5, cex=1.25, font=2)
      abline(mod_fit)
      # add expressions to top left corner
      mtext(paste(" R² =", sprintf(as.character(R2_plot_cal))), side = 3, line = -1.6, at = par("usr")[1], adj = 0, cex = 1.5)
      mtext(paste(" RMSEP =", sprintf(as.character(rmsep_cal))), side = 3, line = -3, at = par("usr")[1], adj = 0, cex = 1.5)
      #mtext(paste(" nComps =", sprintf(as.character(common_nComp))), side = 3, line = -4.5, at = par("usr")[1], adj = 0, cex = 1.5)
      graphics.off()
      
      # Plot predicted vs actual values
      png(file = paste0(pati_main,"/", inYear, "_", inProcs, "_", "val_prediction.png"), width = 1250, height = 850, res = 300)
      par(mar = c(3, 3, 0.4, 0.25)+0.1)#par(mar = c(bottom, left, top, right) + additional)
      #main=paste0(inYear," ",inProcs," ",inIter," Models' Val. Data Prediction"), 
      plot(predicted ~ measured, data = ensemble_val_preds, pch=15, cex=1, xlab = "", ylab = "", xaxt = "n", yaxt = "n",
           xlim=c(0,2500), ylim=c(0,2000))  
      # Add axes and labels
      axis(1, mgp=c(0, 0.45, 0), cex.axis=0.875, cex.lab = 0.375, font.lab = 2)
      axis(2, mgp=c(0, 0.45, 0), cex.axis=0.875, cex.lab = 0.375, font.lab = 2)
      mtext(expression("Measured [g" ~ m^{-2} ~ "]"), side=1, line=2.1, cex=1.25, font=2)
      mtext(expression("Predicted [g" ~ m^{-2} ~ "]"), side=2, line=1.5, cex=1.25, font=2)
      abline(mod_fit_val)
      # add expressions to top left corner
      mtext(paste(" R² =", sprintf(as.character(R2_plot_val))), side = 3, line = -1.6, at = par("usr")[1], adj = 0, cex = 1.5)
      mtext(paste(" RMSEP =", sprintf(as.character(rmsep_val))), side = 3, line = -3, at = par("usr")[1], adj = 0, cex = 1.5)
      #mtext(paste(" nComps =", sprintf(as.character(common_nComp))), side = 3, line = -4.5, at = par("usr")[1], adj = 0, cex = 1.5)
      graphics.off()
      
      
      ### Plot
      # get the x value of the highest point
      
      plot_x_loadings <- as.data.frame(t(loads_iter_avg_df))
      plot_x_loadings <- plot_x_loadings[1,]
      plot_x_loadings <- melt(plot_x_loadings)
      
      top_15_positive <- plot_x_loadings %>%
        top_n(n = 30, wt = value)
      
      top_15_negative <- plot_x_loadings %>%
        arrange(desc(value)) %>%
        tail(n = 30)
      
      # Save positively and negatively correlated bands of the final model
      # which is the average of 30 models' x loadings. 
      Top_bands_pos <- as.data.frame(top_15_positive)
      names(Top_bands_pos) <- c("Wavelength", "top_30_positive")
      Top_bands_pos$Year <- inYear
      Top_bands_pos$Processing <- inProcs
      Top_bands_pos$Model <- "Average"
      
      Top_bands_neg <- as.data.frame(top_15_negative)
      names(Top_bands_neg) <- c("Wavelength", "top_30_negative")
      Top_bands_neg$Year <- inYear
      Top_bands_neg$Processing <- inProcs
      Top_bands_neg$Model <- "Average"
      
      # Define the file path where variables will be saved
      save_file_path <- paste0(pati_main, "/", inYear, "_", Sys.Date(),"_" , inProcs, "_plot_wvl_contrib_data.RData")
      # Save all relevant variables to an RData file
      save(top_15_positive, top_15_negative, plot_x_loadings, inProcs, pati_main, inYear, file = save_file_path)
      
      wvl_contrib <- ggplot() +
        geom_vline(data = top_15_positive, aes(xintercept = as.numeric(variable)), color = "#AADDCC", linetype = "solid", alpha = 0.85, size = 0.9) +
        geom_vline(data = top_15_negative, aes(xintercept = as.numeric(variable)), color = "#E69F00", linetype = "solid", alpha = 0.7, size = 0.9) +
        geom_line(data = plot_x_loadings, aes(x = factor(variable), y = `value`, group = inProcs, color = inProcs), size = 2.5, alpha = 0.9) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "#BB0011", size=1.5, alpha = 0.8) + # add a dashed red line at 0 on the y-axis
        labs( #subtitle = paste0(inYear," ", inProcs, " ", inProcs," PLSR models' top 30 bands"),
          x = "Wavelength [nm]", y = "Loadings") +
        scale_y_continuous(breaks = seq(round(min(top_15_negative$value),2), round(max(top_15_positive$value),2), 0.02), 
                           limits = c(min(top_15_negative$value),max(top_15_positive$value))) +
        scale_x_discrete(breaks = levels(factor(plot_x_loadings$variable))[seq(1, length(levels(factor(plot_x_loadings$variable))), 25)]
                         ,labels = c(410,600,680,775,920,1175,1415,1650,1860,2060,2235,2400)) +
        theme(axis.text.y = element_text(hjust = 1, size=18, face = "bold"),
              axis.text.x = element_text(angle = 45,hjust = 1, size=18, face = "bold"),
              #plot.subtitle=element_text(size=14,face="bold", color="black", hjust = 0.5),
              axis.title.x = element_text(size = 24, face = "bold"),
              axis.title.y = element_text(size = 24, face = "bold"),
              legend.position = "top", # move the legend to the top
              #legend.title = element_text(size = 18),
              #legend.text = element_text(size = 18),
              legend.margin = margin(t = 0, unit = "cm"), # adjust the top margin of the legend
              plot.margin = unit(c(0.5, 0.5, 0, 0.5), "cm")) +
        theme(panel.background = element_blank(),
              #panel.grid.major = element_line(color = "gray"),
              #panel.grid.minor = element_blank(),
              panel.border = element_rect(color = "black", fill = NA, size = 1),
              plot.margin = unit(c(0.15, 0.25, 0.15, 0.15), "cm"))+ # top, right, bottom, left
        scale_color_manual(name = "", #VIPs:
                           values = c( inProcs = "black"))# +
      ggsave(filename = file.path(paste0(pati_main,"/", inYear, "_", inProcs, "_", "common-comps-", "REF_DF_MODEL_X_Loadings.pdf")), plot=wvl_contrib,
             device="pdf", width = 24, height = 12, units = "cm", dpi = 300)
      ggsave(filename = file.path(paste0(pati_main,"/", inYear, "_", inProcs, "_", "common-comps-", "REF_DF_MODEL_X_Loadings.png")), plot=wvl_contrib,
             device="png", width = 24, height = 12, units = "cm", dpi = 300)
      wvl_contrib
      graphics.off()
      ###
      
      ##### Save Coefficients, VIPs, X_Loadings #####
      # save averaged model's statistics
      write.csv(coefs_iter_avg_df, paste0(pati_main,"/", inYear, "_", inProcs, "_", "MAIN_MODELs-avg-COEFS.csv"))
      write.csv(vips_iter_avg_df, paste0(pati_main,"/", inYear, "_", inProcs, "_", "MAIN_MODELs-avg-VIPS.csv"))
      write.csv(nComps_iter_avg_df, paste0(pati_main,"/", inYear, "_", inProcs, "_", "MAIN_MODELs-avg-STATISTICS.csv"))
      write.csv(loads_iter_avg_df, paste0(pati_main,"/", inYear, "_", inProcs, "_", "MAIN_MODEL-avg-X_Loadings.csv"))
      write.csv(Top_bands_pos, paste0(pati_main,"/", inYear, "_", inProcs, "_", "MAIN_MODEL-avg-top_30_POS_bands.csv"))
      write.csv(Top_bands_neg, paste0(pati_main,"/", inYear, "_", inProcs, "_", "MAIN_MODEL-avg-top_30_NEG_bands.csv"))
      # save to recreate final
      write.csv(ensemble_cal_preds, paste0(pati_main,"/", inYear, "_", inProcs, "_", "MAIN_MODELs-avg-CAL-PREDS.csv"))
      write.csv(ensemble_val_preds, paste0(pati_main,"/", inYear, "_", inProcs, "_", "MAIN_MODELs-avg-VAL-PREDS.csv"))

    }
    
  } else {
    print(paste("Column not found:", col_name))
    # Handle case where column is not found
    # For example, continue to next iteration of columns_of_interest
    next  # This will continue to the next iteration of 'col_name'
    
  }
  
}


### END ###

