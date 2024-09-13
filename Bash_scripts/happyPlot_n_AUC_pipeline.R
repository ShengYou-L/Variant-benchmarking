library(happyR)
library(bayestestR)
library(ggplot2)
library(magrittr)
library(tibble)
theme_set(theme_minimal())

args <- commandArgs(trailingOnly = TRUE)
#typeof(args)
print(args)
argslem <- length(args)
#print(argslem)
outputPDF <- outputAUC <- "NONE"
Xmin <- Ymin <- 0
Xmax <- Ymax <- 1
fileArray <- c()
argument <- FALSE
for (i in 1:argslem) {
  if (args[i] == "--Xmin") {
    Xmin <- as.double(args[i+1])
    print(paste("Minimum value for x-axis: ", Xmin))
    argument <- TRUE
  } else if (args[i] == "--Xmax") {
    Xmax <- as.double(args[i+1])
    print(paste("Maximum value for x-axis: ", Xmax))
    argument <- TRUE
  } else if (args[i] == "--Ymin") {
    Ymin <- as.double(args[i+1])
    print(paste("Minimum value for y-axis: ", Ymin))
    argument <- TRUE
  } else if (args[i] == "--Ymax") {
    Ymax <- as.double(args[i+1])
    print(paste("Maximum value for y-axis: ", Ymax))
    argument <- TRUE
  } else if (args[i] == "--graph-output") {
    outputPDF <- args[i+1]
    print(paste("Output graph: ", outputPDF))
    argument <- TRUE
  } else if (args[i] == "--auc-output") {
    outputAUC <- args[i+1]
    print(paste("Output AUC: ", outputAUC))
    argument <- TRUE
  } else if (args[i] == "--sample") {
    Sample <- args[i+1]
    print(paste("Sample: ", Sample))
    argument <- TRUE
  } else if (args[i] == "--seq-type") {
    ngs <- args[i+1]
    print(paste("Sequencing type: ", ngs))
    argument <- TRUE
  } else {
    if (argument == TRUE) {
      argument <- FALSE
    } else if (argument == FALSE) {
      print(paste("File: ", args[i]))
      fileArray <- c(fileArray, args[i])
    } else {
      print("Unrecognised input")
      # stop()
    }
  }
}
print("Check inputs: ")
print(paste("Input files: ", fileArray))

print(paste("X-axis: ", Xmin, Xmax))
print(paste("Y-axis: ", Ymin, Ymax))

hr_variant_list <- c("SNP", "SNP_PASS", "SNP_SEL", 
                     "INDEL", "INDEL_PASS", "INDEL_SEL")
for (each in fileArray) {
  happy_input <- each
  happy_prefix <- sub(".summary.csv", "", happy_input)
  
  happy_data <- happyR::read_happy(happy_prefix)
  
  # some variables for looping and plots
  num.plots <- 3
  my.plots <- vector(num.plots, mode="list")
  counter <- 1
  
  # show data frame
  # knitr::kable(head(happy_data$pr_curve$all[,1:9]))
  if (outputPDF != "NONE") {
    print(sprintf("[Process]: Start plotting graphs for %s", happy_input))
    # Filters
    hapfilter <- c("ALL", "PASS", "SEL")
    for (f in hapfilter) {
      # Choose filter (ALL, PASS, SEL) *normally ALL or SEL as PASS almost identical with SEL
      pr <- pr_data(happy_data, filter = f)
      # Plot ROC/PR curve
      my.plots[[counter]] <- ggplot(pr, aes(x = METRIC.Recall, y = METRIC.Precision, col = Type)) +
        geom_line() + theme_minimal() +
        geom_point(data = happy_data$summary) +
        scale_x_continuous(limits = c(Xmin, Xmax)) +
        scale_y_continuous(limits = c(Ymin, Ymax)) +
        ggtitle(paste("ROC/PR curve for", Sample, ngs, "vs GIAB variants, Filter = ", f))
      
      counter = counter + 1
      graphics.off()
      print(sprintf("[Status]: Graph for filter %s has been plotted successfully", f))
    }
    # Save all plots in my.plots into a PDF
    pdf(outputPDF, onefile = TRUE)
    for (my.plot in my.plots) {
      print(my.plot)
    }
    graphics.off()
  }
  
  if (outputAUC != "NONE") {
    write(each, file = outputAUC, append = TRUE)
    Info <- sprintf("Sample: %s    Sequencing type: %s", Sample, ngs)
    write(Info, file = outputAUC, append = TRUE)
    print(sprintf("[Process]: Start calculating AUC values for %s", happy_input))
    df_list <- list()
    for (x in hr_variant_list) {
      df <- data.frame(
        METRIC.Recall = happy_data$pr_curve[[x]]$METRIC.Recall,
        METRIC.Precision = happy_data$pr_curve[[x]]$METRIC.Precision
      )
      df_list[[x]] <- df
      
      # Save in respective vectors
      recall <- df[['METRIC.Recall']]
      precision <- df[['METRIC.Precision']]
      
      # Replace N/A fields with 0
      recall <- replace(recall, is.na(recall), 0)
      precision <- replace(precision, is.na(precision), 0)
      
      # Call area under curve function
      auc_value <- area_under_curve(recall, precision)
      output_auc <- sprintf("The AUC value for %s is %s", x, auc_value)
  
      write(x, file = outputAUC, append = TRUE)
      write(output_auc, file = outputAUC, append = TRUE)
      write(strrep("*****", 5), file = outputAUC, append = TRUE)
      print("[Status]: Done calculating AUC values")
    }
  }
}