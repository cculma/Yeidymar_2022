# install.packages("emmeans")
rm(list = ls())
library(emmeans)
library(ggplot2)

# Read raw data -----------------------------------------------------------
a1 <- read.csv("~/Documents/git/Yeidymar_2022/Y_Final.csv", header = T, ) 
head(a1)
colnames(a1)
# [1] "Population"            "GenotypeID"            "Resitance"            
# [4] "Population_GenotypeID" "LeafletID"             "DPI"                  
# [7] "Leaflet_Weight"        "CFU"                   "LOG_CFU"              
# [10] "Genome_copies"   

lev1 <- colnames(a1)[1:6]
lev1
# [1] "Population"            "GenotypeID"            "Resitance"            
# [4] "Population_GenotypeID" "LeafletID"             "DPI"  

# To convert the columns to factor
a1[,lev1] <- lapply(a1[,lev1], factor)
str(a1)

mod3 <- lm(LOG_CFU ~ DPI * Resitance, data = a1) # this model can be modified
aov3 <- anova(mod3)
summary(mod3)
summary(mod3)$r.squared 

# The coefficient of determination of the simple linear regression model for the data set is 0.5959747.

# up to here is the same code, now this is new to fin the interactions

mod3.1 <- emmeans(mod3, pairwise ~ Resitance | DPI)

# plot simple
p1 <- emmip(mod3, Resitance ~ DPI, CIs = T)
p1

# extract mean values
df1 <- as.data.frame(mod3.1$emmeans)
head(df1)

# save as tables ----------------------------------------------------------

setwd("~/Documents/git/Yeidymar_2022/")
write.csv(df1, "Resistance_DPI.csv", quote = F, row.names = F)

# plot your results -------------------------------------------------------

# plot 1

plot1 <- ggplot(data = df1, aes(x=DPI, y=emmean, group=Resitance))  + theme_classic(base_family = "Arial", base_size = 12) + geom_point(aes(shape = Resitance), size = 1, alpha = 0.6) + geom_line(aes(linetype = Resitance), alpha = 0.6) + theme(axis.text.x = element_text(angle = 0, hjust = 0.95, vjust = 0.2)) + geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.1, alpha = 0.6) + labs(y = "CFU_log", x = "DPI")

plot1
setwd("~/Documents/git/Yeidymar_2022/")
ggsave(filename = "CFU_DPI.png", plot = plot1, dpi = 300, width = 6, height = 4)

#end
############################################################