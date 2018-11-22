#######################################################################################################################
# Analyse qPCR data from lightcycler 480
#
# - read in raw fluo values and sample information from txt file
# - txt file needs to be exported from lightcycler software
# - plot raw data
# - calculate background corrected/substracted data and plot them
# - calculate efficiencies from standard curve and plot together with amplification cycle for all standards
# - calcuate Ct values with cycle threshold method (Cq calculations are also available)
# --> completed using indicated packages
# 
# - normalize to housekeeping genes and compute mean +/- SD for technical triplicates
# - calculate mean +/- SD for "biological" replicates
# - normalize to control 
# - plot all target genes normalized to one housekeeping gene on one graph and save
# - save different results (Efficiencies, Ct-values, means, normalization) to excel
# --> completed in an grubby way --> should be optimized
#
# 24.10.2017
#######################################################################################################################


# ==================
# Load the packages:
# ==================

# for reading-in of data
if ( ! require(ReadqPCR) )     { source("https://bioconductor.org/biocLite.R") ; biocLite("ReadqPCR") ; library(ReadqPCR) }

# for preprocessing data
if ( ! require(chipPCR) )     { install.packages("chipPCR");     library(chipPCR) }

# for normalization of data 
if ( ! require(NormqPCR) )     { source("https://bioconductor.org/biocLite.R") ; biocLite("NormqPCR") ; library(NormqPCR) }

# Load package for table formatting
if ( ! require(xtable) )     { install.packages("xtable");     library(xtable) }

# for saving data to excel
if ( ! require(XLConnect) )      { install.packages("XLConnect");      library(XLConnect) }

# for calculating row means 
if ( ! require(dplyr) )      { install.packages("dplyr");      library(dplyr) }

# for normalizing to control
if ( ! require(data.table) )      { install.packages("data.table");      library(data.table) }

# for plotting
if ( ! require(ggplot2) )      { install.packages("ggplot2");      library(ggplot2) }
if ( ! require(gridExtra) )      { install.packages("gridExtra");      library(gridExtra) }
if ( ! require(grid) )      { install.packages("grid");      library(grid) }





RTgutGeneExpression_wrapper <- function(InName_data, InName_info, title_name,replicates, n_std, n_samples, n_controls, OutName){
  
  
  #===================================================================================================================
  # DATA IMPORT
  # package "ReadqPCR"
  #========================
  
  # 1. Step: Load RT-qPCR data into R using the "read.LC480" function 
  # With function read.LC480 raw fuoresence data of Roche LightCycler 480 may be read in. 
  # The result is saved in an object of class "CyclesSet". At the moment the function works only txt-file
  # uncommand get information about class "CyclesSet" , which is an extension of the class "eSet"
  # ?CyclesSet 
  # getClass("CyclesSet")
  
  conc_LPS <- read.LC480(file=InName_data)
  
  # read in an additional sample information file (needs to be exported from Lightcycler software from "Sample editor")
  conc_LPS_info <- read.LC480SampleInfo(file=InName_info)
  
  # these information can be merged into one data file
  conc_LPS_data <- merge(conc_LPS, conc_LPS_info)
  
  # fluorescence data can be accessed as follows - save in a data.frame
  fluo <- as.data.frame(exprs(conc_LPS_data))
  names <- conc_LPS_info@data$"Sample name" # extract the sample names from the information file 
  colnames(fluo) <- names # add the sample names to the data 
  
  # remove columns that do not contain any measurements
  # get index for all columns that contain the word "sample"
  # it is really important, that the column names of columns that include data have a different name
  count_sample <- grep(pattern="^Sample", colnames(fluo))
  # remove the columns from the data frame
  # "tail" returns the last part of a vector or data frame
  fluo[count_sample[1]:tail(count_sample,n=1)] <- list(NULL)
  
  # save the concentration of the standard curve
  conc_std <- as.numeric(conc_LPS_data@phenoData@data[9][1:(replicates*n_std),1])
  
  # save the sample information in a vector
  sample <- names
  # remove all unlabeled starting with "Sample"
  s <- grep("^Sample", sample) # show all indices with start with "Sample"
  sample <- sample[1:(s[1]-1)] # save everything which has a sample name different than "Sample"
  
  # save the target gene names in a vector
  target <- c(conc_LPS_data@phenoData@data$"Target name")
  # remove the repitition of names and additionally remove empty "strings" 
  # it is really important that each measured set has an associated target name
  target_names <- unique(target[target !=""]) 
  
  # calculate the number of wells measured per target gene
  n_wells <- (n_std+n_samples+n_controls)*replicates
  
  # save the numbers of cycles separatly in a vector
  cycles <- c(conc_LPS_data@featureData@data$"Cycle number")
  
  # combine the number of cycles and the fluo values in a dataframe using "cbind"
  fluo_data <- cbind(cycles, fluo)
  
  
  # plot raw data and automatically save them 
  png(filename = paste(title_name, "raw_data.png", sep="_"), width = 900, height = 700, res = 100)
  par(oma=c(0,0,2,0))
  matplot(fluo_data[,-1], type = "l", pch = 19, col = 1, lty = 1,
          xlab = "Cycle", ylab = "Raw fluorescence", main = "Raw fluorescence data")
  # abline(h = 5, lty = 2, col="red")
  title(paste(title_name), outer=TRUE)
  dev.off()
  
  
  #===================================================================================================================
  # BACKGROUND SUBSTRACTION
  # package "chipPCR"
  # using the CPP command
  #========================
  
  # Despite the growing interest in amplification-based techniques, there are only few open source
  # tools for pre-processing real-time amplification data.
  # chipPCR is a versatile R package for pre-processing (e.g., imputation of
  # missing values, smoothing) and quality analysis of raw data for amplification curves coming from conventional
  # quantitative polymerase chain reactions (qPCR), and quantitative isothermal amplification (qIA) reactions.
  
  # The CPP function is a wrapper for pre-processing algorithms, which includes normalization, background subtraction,
  # outlier removal using the (fixNA function) in the background range and tests for positive amplification
  # reactions.
  
  # it uses "lm.coefs" calculates the estimated background of the amplification curve
  # 4 methods are possible to chose (one and three robust methods - better in case of outliers):
  # "lmrob" - MM-type estimator for linear regressions - default
  # "rq" - quantile regression fit
  # "least" - least squares linear regression
  # "rfit" - rank-based estimates of regression coefficients
  
  # includes different regression types and the function bg.max to detect and correct background noise
  background_sub_CPP <- function(x){
    CPP(seq(length(x)),x, 
        method.reg="lmbrob",  # MM-type estimator for linear regressions - default
        method.norm="none", # No normalization is used, other possible inputs "minm", "max", "lugn", "zscore"
        bg.range=NULL)$y  #bg.range is numeric vector of length 2 indicating the background, 
    #if NULL, the background is calculated by "bg.max function"
  }
  
  
  fluo_data_bgr <- apply(X=fluo_data[,-1],
                         MARGIN=2,
                         FUN=background_sub_CPP)
  
  # plot background substracted data
  png(filename = paste(title_name, "background_sub.png", sep="_"), width = 900, height = 700, res = 100)
  par(oma=c(0,0,2,0))
  matplot(cycles, fluo_data_bgr, type = "l", pch = 19, col = 1, lty = 1, 
          xlab = "Cycle", ylab = "Raw fluorescence", main = "background substracted - CPP-method")
  title(paste(title_name), outer=TRUE)
  dev.off()
  
  
  
  #===================================================================================================================
  # Threshold Cycle Method
  # package "chipPCR"
  #==============
  # th.cyc calculates the number of cycles at which the fluorescence exceeds a defined threshold, 
  # called the threshold cycle (Ct).The th.cyc calculates the intersection of the user defined 
  # Ct value threshold (r) and a linear regression or quadratic polynomial in the range of the user defined Ct value. 
  # inputs as vectors
  # it returns a matrix with the calculated ct values and the background threshold applied
  
  
  ct_method <- function(x, cycles){
    th.cyc(x = cycles, # vector containing the time or cycle values
           y = x,  # vector containing the fluo values
           r=1) # fluorescence values which defines the threshold
  }
  
  
  ct_values <- apply(X = fluo_data_bgr, # background corrected data
                     MARGIN=2,  # by columns
                     FUN=ct_method,  # calculates the threshold
                     cycles)
  ct_values <- as.data.frame(ct_values) # save as data.frame
  
  
  #===================================================================================================================
  # AMPLIFICATION EFFICIENCY
  # package "chipPCR"
  #==============
  # effcalc is used to calculate the amplification efficiency of a dilution curve --> for standards!
  # needs processed data (=calculated Cq/Ct values)
  
  
  # Plot the Amplification curves of the Standard dilutions with replicates
  # count the column numbers in which the word "standard" appears
  # will be used for plotting only the standard without rearranging the data.frame with fluo values
  count_standard <- grep(pattern="^Standard", colnames(fluo_data_bgr))
  
  # calculate the number of wells measured per standard per target
  standard_wells <- length(count_standard)/length(target_names)
  
  # calculate the index which can be used to get the numbers from "count_standard"
  # these index numbers can be used to extract the standard fluo data from the data frame fluo_data_bgr
  index <- seq(1,standard_wells*length(target_names)+1,by=standard_wells)
  
  ampli_eff <- list() # empty list to store the amplification efficiencies for each target
  
  for(i in 1:length(target_names)){
    rainbowcols = rep(rainbow(5), each=3)
    png(filename = paste(title_name, target_names[i], "standard_efficiencies.png", sep="_"), width = 1200, height = 700, res = 100)
    par(mfrow=c(1,2), oma=c(0,0,2,0))
    # plot the amplification curves for the standards
    matplot(cycles, fluo_data_bgr[,count_standard[index[i]]:(count_standard[index[i]]+(standard_wells-1))], 
            type="l", col=rainbowcols, 
            xlab = "Cycle", ylab = "Raw fluorescence", main = paste("Standard dilution", target_names[i], sep=" - "))
    legend(0, max(fluo_data_bgr[,count_standard[index[i]]:(count_standard[index[i]]+(standard_wells-1))]),
           colnames(fluo_data_bgr[,count_standard[index[i]]:(count_standard[index[i]]+(standard_wells-1))]), 
           col = rainbowcols, pch = 19, bty = "n")
    
    # calculate and plot the efficiencies 
    ct <- ct_values[1,count_standard[index[i]]:(count_standard[index[i]]+(standard_wells-1))]
    res <- cbind(conc_std, t(ct)) # combine the concentrations of the standard with the ct values
    res <- res[res[,2]<=50,]  # remove all rows where the calculated ct values exceed 50 --> something wrong with pipetting
    eff <- effcalc(res[,1],res[,2]) # calculate the amplification efficiencies
    
        # save the results in a list
    ampli_eff[[i]] <- eff@amplification.efficiency # to access the amplification efficiencies
    plot(effcalc(res[,1],res[,2]), main=paste("Efficiency plot", target_names[i], sep=" - ")) # plot the efficiencies against the concentration
    
    title(paste(title_name, target_names[i], sep=" - "), outer=TRUE)
    dev.off()
    
  }
  
  #names(ampli_eff) <- target_names[1:length(target_names)]
  eff_values <- do.call("rbind", ampli_eff) #bind the different lists into one data frame
  #rownames(eff_values) <- c()
  
  eff <- data.frame(target_names, eff_values, eff_values/100*2)
  colnames(eff)<- c("target", "efficiencies [%]", "efficiencies")
  
  
  #===================================================================================================================
  # calculate with the efficiencies and Ct values
  
  ct_values_all <- data.frame(t(ct_values)) # convert dataframe (change rows with columns)
  ct_values_only <- ct_values_all[,1] # save only ct_values without sample information in a dataframe
  
  # calculate the Efficiency^(-CT)
  E_ct <- (rep(eff[,3], each=n_wells))^(-ct_values_only)
  
  # normalize based on "18S" housekeeping gene
  E_ct_18s <- E_ct[(n_wells*4+1):(n_wells*5)]
  E_ct_norm_18s <- E_ct/rep(E_ct_18s,6)
  
  # normalize based on EF1 housekeeping gene
  E_ct_EF1 <- E_ct[(n_wells*5+1):(n_wells*6)]
  E_ct_norm_EF1 <- E_ct/rep(E_ct_EF1,6)
  
  # combine data in a data frame
  ct_data <- cbind.data.frame(rep(target_names, each=n_wells), # target (=gene) names
                              rep(eff[,3], each=n_wells),      # efficiencies for each gene
                              sample,                          # name of samples
                              ct_values_only,                  # calculated ct values
                              E_ct,                            # calculated E_ct values
                              E_ct_norm_18s,                   # normalized on 18s
                              E_ct_norm_EF1)                   # normalized on EF1
  
  colnames(ct_data) <- c("target", "efficiency",  "sample", "Ct value", "E_ct", "norm_18s", "norm_EF1")
  
  
  # save data without standards and controls in separate data frames
  datNew1 <- ct_data[-which(ct_data$sample == "Standard"), ]
  datNew2 <- datNew1[-which(datNew1$sample == "H20"), ]
  ct_data_samples <- datNew2[-which(datNew2$sample == "RT3"), ]
  
  
  # calculate mean +/- SD out of the three technical replicates for each target
  # "summarize" will create a new data frame with mean +/- SD values
  # use "mutate" instead of "summarize" in case the calculated means and SD shall be added to the existing data frame
  
  ct_data_mean_18s_tech <- ct_data_samples %>% group_by(target, sample) %>% summarize(Mean = mean(norm_18s), std = sd(norm_18s))
  ct_data_mean_EF1_tech <- ct_data_samples %>% group_by(target, sample) %>% summarize(Mean = mean(norm_EF1), std = sd(norm_EF1))
  
  
  # calculate mean +/- SD out of the sample treatments for each target
  # for 18 s
  
  dat_18s <- data.frame(do.call(rbind, strsplit(as.vector(ct_data_mean_18s_tech$sample), split =" ")),stringsAsFactors = FALSE)
  dat_18s$X2 = as.integer(dat_18s$X2)
  ct_data_mean_18s_tech_1 <- cbind.data.frame(ct_data_mean_18s_tech, dat_18s$X2)
  
  ct_data_mean_18s <- ct_data_mean_18s_tech_1 %>% group_by(target, dat_18s$X2) %>% 
    summarize(mean = mean(Mean), std = sd(Mean))
  colnames(ct_data_mean_18s) <- c("target", "conc", "Mean", "SD")
  
  # for EF1
  ct_data_mean_EF1_tech_1 <- cbind.data.frame(ct_data_mean_EF1_tech, dat_18s$X2)
  
  ct_data_mean_EF1 <- ct_data_mean_EF1_tech_1 %>% group_by(target, dat_18s$X2) %>% 
    summarize(mean = mean(Mean), std = sd(Mean))
  colnames(ct_data_mean_EF1) <- c("target", "conc", "Mean", "SD")
  
  
  
  # normalize each target gene to its own control and plot the results
  # define a function that divides the mean and the SD by the first mean (=control mean)
  normByMean = function(x) {
    x$SD = x$SD / x$Mean[1]
    x$Mean = x$Mean / x$Mean[1]
    return(x)
  }
  
  
  # for 18s housekeeping gene
  # split the data frame into a list of data frames according to the target name
  t1 <- split(ct_data_mean_18s, f = ct_data_mean_18s$target)
  
  # normalize all values with the control mean
  t2 <- lapply(t1,normByMean)
  # combine the list back into one data frame
  norm_18s <- do.call("rbind",t2)
  
  
  
  # for EF1 housekeeping gene
  # split the data frame into a list of data frames according to the target name
  m1 <- split(ct_data_mean_EF1, f = ct_data_mean_EF1$target)
  # normalize all values with the control mean
  m2 <- lapply(m1,normByMean)
  # combine the list back into one data frame
  norm_EF1 <- do.call("rbind",m2)
  
  
  
  #===================================================================================================================
  # plot the data 
  #normalized to EF1
  
  png(filename = paste(title_name, "Normalized_EF1.png", sep="_"), width = 1000, height = 700, res = 100)
  plots <- list() # safe all plots in a list
  for(i in 1:length(target_names)){
    
    xx <- as.data.frame(m2[i])
    colnames(xx) <- c("target", "conc", "Mean", "SD")
    
    pl = ggplot(data=xx,aes(x=target,y=Mean,fill=factor(conc))) +
      geom_bar(stat="identity", position=position_dodge()) + labs(x="conc. [ug/ml",y="fold change") +
      coord_cartesian(ylim = c(0,max(xx$Mean+xx$SD)), expand=T)
    
    # for ggplot in a for loop to work: needs to be specifically printed with "print"
    plots[[i]] <- (pl + theme(legend.title=element_blank()) + geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD),
                                                                            size=.4,    # Thinner lines
                                                                            width=.4,
                                                                            position=position_dodge(.9)) +
                     ggtitle(levels(xx$target)[i]))
    
  }
  # use grid.arrange from gridExtra package to plot all plots in one
  grid.arrange(plots[[1]], plots[[2]], plots[[3]],plots[[4]], plots[[5]], plots[[6]], nrow=2, 
               top = textGrob("Normalized to EF1",gp=gpar(fontsize=20,font=2)))
  dev.off()
  
  
  
  #===================================================================================================================
  # plot the data 
  # normalized to 18s
  
  png(filename = paste(title_name, "Normalized_18s.png", sep="_"), width = 1000, height = 700, res = 100)
  plots <- list() # safe all plots in a list
  for(i in 1:length(target_names)){
    
    xx <- as.data.frame(t2[i])
    colnames(xx) <- c("target", "conc", "Mean", "SD")
    
    pl = ggplot(data=xx,aes(x=target,y=Mean,fill=factor(conc))) +
      geom_bar(stat="identity", position=position_dodge()) + labs(x="conc. [ug/ml",y="fold change") +
      coord_cartesian(ylim = c(0,max(xx$Mean+xx$SD)), expand=T)
    
    # for ggplot in a for loop to work: needs to be specifically printed with "print"
    plots[[i]] <- (pl + theme(legend.title=element_blank()) + geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD),
                                                                            size=.4,    # Thinner lines
                                                                            width=.4,
                                                                            position=position_dodge(.9)) +
                     ggtitle(levels(xx$target)[i]))
    
  }
  # use grid.arrange from gridExtra package to plot all plots in one
  grid.arrange(plots[[1]], plots[[2]], plots[[3]],plots[[4]], plots[[5]], plots[[6]], nrow=2, 
               top = textGrob("Normalized to 18s",gp=gpar(fontsize=20,font=2)))
  dev.off()
  
  
  #===================================================================================================================
  # EXCEL
  # export ct_data and mean +/- SD into excel
  sheets <- list(eff, ct_data, ct_data_mean_18s, ct_data_mean_EF1, norm_18s, norm_EF1 )
  names(sheets) <- c("efficiencies", "Ct_values", "Ct_values_mean_18s", "Ct_values_mean_EF1", "Norm_18s", "Norm_EF1")
  wb <- loadWorkbook(file.path(paste("Results", OutName,sep="_")), create = TRUE)
  createSheet(wb, name=names(sheets))
  writeWorksheet(wb, sheets, names(sheets), header=TRUE)
  saveWorkbook(wb)
  
}

