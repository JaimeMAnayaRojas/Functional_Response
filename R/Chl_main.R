# clean memory
rm(list=ls())
# load libraries

library(ggplot2)
library(dplyr)
library(readxl)

# read file
getwd()
data <- read_xlsx("data/Meso_data_S2.xlsx", sheet = "Measurements", skip = 4)

head(data)
names(data)
names(data)[13]
names(data)[13] <- "Cyano.ugL"
names(data)[14]
names(data)[14] <- "Total.ugL"


# get mean and sd of Total.ugL base on treatment
data$Treatment = factor(data$Treatment)

levels(data$Treatment)
levels(data$Treatment)[4] <- "P+"


levels(data$Treatment)
levels(data$Treatment)[c(1,3)] <- "P-"

levels(data$Treatment)
levels(data$Treatment)[c(2,3)] <- c("LowV", "HighV")


data %>% 
    group_by(Treatment) %>% summarise(mean = mean(Total.ugL), sd = sd(Total.ugL)) -> sumdata

sumdata

## fit model

library(brms)

# fit model

fit <- brm(log(Total.ugL) ~ Treatment + (1|Block), data = data, family = gaussian(), prior = c(set_prior("normal(0, 10)", class = "Intercept"), set_prior("normal(0, 10)", class = "b")), chains = 4, cores = 4, iter = 2000, warmup = 1000, control = list(adapt_delta = 0.99))
summary(fit)

postPhyto <- posterior_samples(fit)

# get posterior predictions for the treatment

pControl <- exp(postPhyto$b_Intercept)
pLV <- exp(postPhyto$b_Intercept + postPhyto$b_TreatmentLowV)
pHV <- exp(postPhyto$b_Intercept + postPhyto$b_TreatmentHighV)

## 

LOS <- function(x){

  p = 100*  length(which(x >0 )) / length(x)
  return(p)


}

# Level of support
LOS(pLV - pControl) # effect of low virulence


LOS(pHV- pControl) # effect of high virulence


mean(pHV/pControl) # the density of phytoplankton is 2.0 times higher in the high virulence treatment than in the control
hdi(pHV/pControl)


# put the posterior in a dataframe

predPhyto <- data.frame(value = c(pControl, pLV, pHV),
Treatment = rep(c("P-", "LowV", "HighV"), each = 4000))

# get the mean and HPDi of the posterior predictions

library(HDInterval)

ppp <- cbind(pControl, pLV, pHV)
data.frame(m = apply(ppp, 2, median),
l = t(apply(ppp, 2, hdi, prob = 0.95))[,1],
u = t(apply(ppp, 2, hdi, prob = 0.95))[,2], 
Treatment = c("P-", "LowV", "HighV")) -> sum_stat

sum_stat$Treatment <- factor(sum_stat$Treatment, levels = c("P-", "LowV", "HighV"))
data$Treatment <- factor(data$Treatment, levels = c("P-", "LowV", "HighV"))



##

p1 <- ggplot(sum_stat, aes(y = m, x = Treatment)) +
  geom_jitter(data = data, aes(y = Total.ugL, x = Treatment, 
  group = factor(Block), 
  colour = factor(Block)), width = 0.1) + 
  geom_point(size = 3, color = "black") +
  geom_errorbar(aes(ymin = l, ymax = u), 
  size = 1.25, color = "black", width = 0) +
  # geom_line(data = data, aes(y = Total.ugL, x = Treatment, group = factor(Block), 
  # colour = factor(Block)), alpha = 0.5) + 
  theme_bw() + theme_classic() + # remove legend for the geom_jitter
  theme(legend.position = "none") +
  ylab("Total (ugL)") + # add title and axis labels
  xlab("Treatment") +
  ggtitle("A) Phytoplankton biomass (ug/L)") 
  
# save plot to pdf file

#ggsave("plots/Phytoplankton.pdf", width = 10, height = 10, units = "cm")

### Do the same for cyanobacteria 

fit2 <- brm(log1p(Cyano.ugL) ~ Treatment + (1|Block), data = data, family = gaussian(), prior = c(set_prior("normal(0, 10)", class = "Intercept"), set_prior("normal(0, 10)", class = "b")), chains = 4, cores = 4, iter = 2000, warmup = 1000, control = list(adapt_delta = 0.99))
summary(fit2)

postCyano <- posterior_samples(fit2)

# get posterior predictions for the treatment

pControl <- exp(postCyano$b_Intercept) - 1 
pLV <- exp(postCyano$b_Intercept + postCyano$b_TreatmentLowV) -1 
pHV <- exp(postCyano$b_Intercept + postCyano$b_TreatmentHighV) -1 

# Level of support
LOS(pLV - pControl) # effect of low virulence

LOS(pHV- pControl) # effect of high virulence

mean(pHV/pControl) # the density of phytoplankton is 1.8 times higher in the high virulence treatment than in the control
hdi(pHV/pControl)



# put the posterior in a dataframe

predCyano <- data.frame(value = c(pControl, pLV, pHV),
Treatment = rep(c("P-", "LowV", "HighV"), each = 4000))

# get the mean and HPDi of the posterior predictions

library(HDInterval)

ppp <- cbind(pControl, pLV, pHV)
data.frame(m = apply(ppp, 2, median),
l = t(apply(ppp, 2, hdi, prob = 0.95))[,1],
u = t(apply(ppp, 2, hdi, prob = 0.95))[,2],
Treatment = c("P-", "LowV", "HighV")) -> sum_stat

sum_stat$Treatment <- factor(sum_stat$Treatment, levels = c("P-", "LowV", "HighV"))
data$Treatment <- factor(data$Treatment, levels = c("P-", "LowV", "HighV"))

##

p2 <- ggplot(sum_stat, aes(y = m, x = Treatment)) +
  geom_jitter(data = data, aes(y = Cyano.ugL, x = Treatment, 
  group = factor(Block), 
  colour = factor(Block)), width = 0.1) + 
  geom_point(size = 3, color = "black") +
  geom_errorbar(aes(ymin = l, ymax = u), 
  size = 1.25, color = "black", width = 0) +
  # geom_line(data = data, aes(y = Total.ugL, x = Treatment, group = factor(Block), 
  # colour = factor(Block)), alpha = 0.5) + 
  theme_bw() + theme_classic() + # remove legend for the geom_jitter
  theme(legend.position = "none") +
  ylab("Cyano (ugL)") + # add title and axis labels
  xlab("Treatment") +
  ggtitle("B) Cyanobacteria biomass (ug/L)")

# save plot to pdf file

# put plot p1 and p2 together

library(cowplot)

p12 <- plot_grid(p1, p2, ncol = 2, align = "v", rel_heights = c(1, 1))

# save cowplot to pdf file

ggsave("plots/Phyto_Cyano.pdf", width = 17, height = 10, units = "cm")









#