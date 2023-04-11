# Data analysis for Yeidymar

library(ggplot2)
library(tidyverse)
library(emmeans)
library(agricolae)
library(hrbrthemes)

a1 <- read.csv("~/Documents/git/Yeidymar_2022/Alf_03.csv", header = T)
head(a1)
tail(a1)
colnames(a1)
str(a1)

lev1 <- colnames(a1)[1:4]

a1[,lev1] <- lapply(a1[,lev1], factor)
str(a1)

# 2 * 65 * 3 * 4 * 3

head(a1)
a1 <- a1 %>% unite(col = "pop_gen", 1:2, sep = "_", remove = F)
str(a1)
a1$pop_gen <- as.factor(a1$pop_gen)
summary(a1$pop_gen)
cc <- count(a1, pop_gen)

colnames(a1)
head(a1)

a2 <- aggregate(LOG_CFU ~ pop_gen + LeafletID + DPI, data = a1, FUN = "mean")
a3 <- aggregate(Genome_copies ~  pop_gen + LeafletID + DPI, data = a1, FUN = "mean")

nrow(a3)
nrow(a2)
colnames(a2)

a4 <- full_join(a2, a3, by = c("pop_gen","LeafletID","DPI"))
str(a4)
# 65 * 3 * 4 = 780

summary(a4$pop_gen)
cc <- count(a4, pop_gen)

# fix repeated data e.g. 99, 90, 86


head(a4)
a5 <- a4 %>% dplyr::filter(DPI %in% c("9"))

#~~~~~~~~~~~~~~~~
# LOG_CFU Number of bacteria
# Genome_copies gene number by qPCR

cor(a4$LOG_CFU, a4$Genome_copies, use = "complete.obs")

ggplot(a4, aes(x=LOG_CFU, y=Genome_copies)) + geom_point() + theme_ipsum(base_family = "Arial", base_size = 12)

head(a4)
str(a4)

ggplot(a4, aes(x = DPI, y = Genome_copies, fill = DPI)) + geom_boxplot() 
ggplot(a4, aes(x = DPI, y = LOG_CFU, fill = DPI)) + geom_boxplot() 

ggplot(a4, aes(x = plant_id_name, y = log_CFU, fill = DPI)) + geom_boxplot() # too much

colnames(a4)

a5 <- a4 %>% dplyr::filter(DPI %in% "0")

a5 <- a4 %>% dplyr::filter(pop_gen %in% c("5425_150", "5425_103"))
head(a5)

ggplot(a5, aes(x = DPI, y = LOG_CFU, group = pop_gen)) + geom_line(aes(color = pop_gen)) + geom_point(aes(color = pop_gen))


# AOV

?lm
?aov

mod1 <- lm(LOG_CFU ~ DPI, data = a4)

mod1 <- lm(LOG_CFU ~ pop_gen, data = a5)
anova(mod1)

mod2 <- lm(log_CFU ~ DPI + plant_id_name, data = a4)
anova(mod2)

colnames(a4)
colnames(a1)
nlevels(a1$GenotypeID)
mod3 <- lm(LOG_CFU ~ DPI + GenotypeID + Population + DPI:GenotypeID + DPI:Population + GenotypeID:Population + DPI:GenotypeID:Population, data = a1)
anova(mod3)
# 100 = 99.5% diferente entre grupos (DPI), 
# 0.05 que las dif sean por azar.
# H0 DPI 0 =3 =6 = 9 
# H1 DPI 0 =3 = 6 dif 9 
# (0 vs 3) 0.05 error, (0 vs 6) 0.05 error, 




mod4 <- aov(mod3)
mod4 <- aov(mod1)
?TukeyHSD
# lsmeans(mod3, "plant_id_name")
tukey.test <- TukeyHSD(mod4)
?TukeyHSD

p_sig <- as.data.frame(tukey.test[["pop_gen"]])

colnames(p_sig)[4] <- "p_adj"
dim(p_sig)
head(p_sig)

p_sig1 <- p_sig %>% dplyr::filter(p_adj < 0.05)
dim(p_sig1)

p_0 <- p_sig1
p_9 <- p_sig1

#~~~~~~~~~~~~~~
1.287146 - 1.214723

0.1571975 - 0.1571961


plant.lm <- lm(weight ~ group, data = PlantGrowth)
plant.av <- aov(plant.lm)
summary(plant.av)

tukey.test <- TukeyHSD(plant.av)
tukey.test

# least squared means 
# arithmentic mean

library(lsmeans)
### Factorial experiment
head(warpbreaks)
str(warpbreaks)
dat1 <- warpbreaks
warp.lm <- lm(breaks ~ wool * tension, data = warpbreaks)
summary(warp.lm)
aov(warp.lm)
anova(warp.lm)

( warp.lsm <- lsmeans (warp.lm,  ~ wool | tension,
                       options = list(estName = "pred.breaks")) )

class(warp.lsm)
pairs(warp.lsm) # remembers 'by' structure
contrast(warp.lsm, method = "poly", by = "wool")

( warp.lsm1 <- lsmeans (warp.lm, "tension") )
contrast(warp.lsm1, "trt.vs.ctrlk")
