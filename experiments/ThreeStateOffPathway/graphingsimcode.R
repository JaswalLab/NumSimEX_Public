# Code to graph best-sim-results-2024-12-11.txt 

library(tidyverse) 
library(grDevices) 
coldata <- read.table('best-sim-results-2024-12-11.txt', header = TRUE)%>% 
  mutate(Time = as.factor(signif(Time, digits = 5)), 
  Source = as.factor(Source)) 
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
png('graphsim.png', width = 3840, height = 2160) 
simplot 
dev.off() 
