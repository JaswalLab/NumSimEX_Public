# Code to graph of C:\Users\jaspe\Documents\GitHub\Jaswal-NumericalSimulations\JaswalNS\Collation\SimOutput.txt 

library(tidyverse) 
library(grDevices) 
simdata <- read.table('C:\Users\jaspe\Documents\GitHub\Jaswal-NumericalSimulations\JaswalNS\Collation\SimOutput.txt', header = TRUE)%>% 
  mutate(Time = as.factor(Time))
simplot <- ggplot(simdata, aes(x = mz, y = Intensity, color = Time)) + 
  geom_line(size = 3) + 
  facet_wrap(~ Time) + 
  theme_bw() + 
  theme(axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.text.x = element_text(angle = 90), 
        panel.grid = element_blank(), 
        legend.position = 'none', 
        text = element_text(size = 60) 
  ) + 
  xlab('m/z') 
png('hxms_simulation.png', width = 3840, height = 2160) 
simplot 
dev.off() 
