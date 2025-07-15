args <- commandArgs(trailingOnly=TRUE)

if(!length(args)){
  print("Usage:")
  print("Rscript HXMS_Spectra_Normalizer.R 'sim_input_data_path.txt' 'exp_input_data_path.txt' 'height'/'width', 'output_data_path.txt'")
  print("Example:")
  print("Rscript HXMS_Spectra_Normalizer.R example_normalizer_sim_input.txt example_normalizer_exp_input.txt height example_normalizer_output.txt")
  q()
}

library(tidyverse)
library(magrittr)
library(grDevices) 

input_sim_file_path <- args[1]
input_exp_file_path <- args[2]
norm_type <- args[3]
output_file_path <- args[4]
plot_file_path <- paste0('graphsim_', norm_type, '.png')


##Normalization functions
trapezoidal_integrator <- function(y, x){
  # definite integral of y with respect to x
  area <- sum(diff(x)*(head(y, -1) + tail(y, -1))/2)
  return(area)
}

normalize_by_area <- function(time_data){
  total_area <- trapezoidal_integrator(time_data$Intensity, time_data$mz)
  time_data <- time_data%>%
    mutate(Intensity = Intensity/total_area)
  return(time_data)
}

normalize_by_height <- function(time_data){
  max_height <- max(time_data$Intensity)
  time_data <- time_data%>%
    mutate(Intensity = Intensity/max_height)
  return(time_data)
}
##Goodness of Fit functions
tp_errors <- function(timedata){
  simdata <- filter(timedata, Source == "Simulation")
  expdata <- filter(timedata, Source == "Experiment")
  time <- as.numeric(as.character(timedata$Time[1]))
  avg <- 2/diff(range(expdata$mz))
  s <- trapezoidal_integrator((expdata$Intensity-simdata$Intensity)^2, expdata$mz)
  a <- trapezoidal_integrator((expdata$Intensity-avg)^2, expdata$mz)
  return(cbind.data.frame(s, a, time))
}

gof <- function(fulldata){
  errors <- fulldata%>%
    # split the data based on time and source
    split(list(.$Source, .$Time))%>% 
    # normalize data to area
    map_dfr(normalize_by_area)%>%
    # split the data based on time
    split(.$Time)%>%
    # find squared and average errors
    map_dfr(tp_errors)
  totalerror <- trapezoidal_integrator(errors$s, 
                                       as.numeric(as.character(errors$time)))
  #print(paste("The totalerror was", totalerror)) 
  totalavgerror <- trapezoidal_integrator(errors$a, 
                                          as.numeric(as.character(errors$time)))
  #print(paste("The totalavgerror was", totalavgerror))
  gof_value <- signif((1-(totalerror/totalavgerror)), 4)
  return(gof_value)
}

##Import data
simdata <- read.table(input_sim_file_path, header = TRUE)%>% 
  mutate(Source = 'Simulation') 
expdata <- read.table(input_exp_file_path, header = TRUE)%>% 
  mutate(Source = 'Experiment') 
coldata <- rbind(simdata, expdata)%>% 
  mutate(Time = as.factor(signif(Time, digits = 5)), 
         Source = as.factor(Source)) 

##Normalize
if(norm_type == "height"){
  coldata <- coldata%>%
    # split the data based on time and source
    split(list(.$Source, .$Time))%>% 
    # normalize data to area
    map_dfr(normalize_by_height)
}else if(norm_type == "area"){
  coldata <- coldata%>%
    # split the data based on time and source
    split(list(.$Source, .$Time))%>% 
    # normalize data to area
    map_dfr(normalize_by_area)
}
coldata <- coldata %>%
  arrange(Source, Time, mz)

##Create plot
simplot <- ggplot(coldata, aes(x = mz, y = Intensity, color = Source)) + 
  geom_line(size = 3) + 
  facet_wrap(~ Time) + 
  theme_bw() + 
  theme(axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.text.x = element_text(angle = 90), 
        panel.grid = element_blank(), 
        legend.position = 'top', 
        text = element_text(size = 60) 
  ) + 
  xlab('m/z') 

##Export data
write.table(coldata, 
            file = output_file_path,
            quote = FALSE,
            row.names = FALSE)

##Export plot
png(plot_file_path, width = 3840, height = 2160) 
simplot 
dev.off() 

gof_val <- gof(coldata)
summary_string <- paste0("The file '", output_file_path, 
                         "' contains the simulation and experimental data normalized to ", norm_type, 
                         ". The time courses are plotted in the file '",
                         plot_file_path, 
                         "'. The GoF was ", gof_val)
print(summary_string)