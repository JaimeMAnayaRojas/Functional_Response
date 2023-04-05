####  This is the script to analyses Nanda's project module data
library(lme4)
library(car)
library(lmerTest)
library("plyr")
library(data.table)

rm(list=ls(all=TRUE))



###########################################################################################################
# First get the data
getwd() # check your working directory

Data <- (read.csv("data/Data.csv"))



Data$N = 1
Data$Family = factor(Data$Family)
levels(Data$Family)[2:7]
levels(Data$Family)

Data$PGenotype = factor(Data$Family)

levels(Data$PGenotype) <- c(NA, "IBB", "IBB", "IBB", "LKZ", "LKZ", "LKZ")
levels(Data$PGenotype)


vars = names(Data)[c(1:5,7,10, 13, 11,12)]
vars
df = ddply(Data, vars ,summarise, Size= mean(Size), Parasites= mean(Parasites), N = sum(N))
df

df$N

# So, something we want to look at is if parasite exposure and infection alter life history traits such somathic growth and survival
# of copepods
# You have look at a set of copepods for four weeks, so we can compare the growth in those time intervals.
# we need to organize the data in a way that we can use, so I will split it into weeks for now.

w1 <- subset(df, Week == 1)
w2 <- subset(df, Week == 2)
w3 <- subset(df, Week == 3)
w4 <- subset(df, Week == 4)

# Something we need to do is to calculate how the size of an individual change from week 1 to 2, from week 2 to 3, and week 3 to 4,
# so, we need to add the size of the same individual from week 2 in front of week 1, and so oo.

# so, I can create a new data for each one

w2a <- w2
w3a <- w3
w4a <- w4

# now we change the "Size" variable to "Size2"
id = which(names(w2a) == "Size") 
names(w2a)[id] <- "Size2"
names(w3a)[id] <- "Size2"
names(w4a)[id] <- "Size2"

# Merge the week sub data sets ---------------------------------------------
names(Data)

w1$name = paste(w1$Plate, w1$Well, sep = "-")
w2a$name = paste(w2a$Plate, w2a$Well, sep = "-")
w1c = merge.data.table(w1, w2a[,c("name", "Size2")], by = "name" )
head(w1c)

w2$name = paste(w2$Plate, w2$Well, sep = "-")
w3a$name = paste(w3a$Plate, w3a$Well, sep = "-")
w2c = merge.data.table(w2, w3a[,c("name", "Size2")], by = "name" )
head(w2c)


w3$name = paste(w3$Plate, w3$Well, sep = "-")
w4a$name = paste(w4a$Plate, w4a$Well, sep = "-")
w3c = merge.data.table(w3, w4a[,c("name", "Size2")], by = "name" )
head(w3c)


## Now put the weeks c into a single dataframe

new_df = data.frame(rbind(w1c, w2c, w3c))
head(new_df)
range(new_df$N)
# calculate growth

new_df$Treatment = factor(new_df$Treatment)
levels(new_df$Treatment) # what does it mean "C/E in Treatment?
levels(new_df$Treatment)[2] = "Exposed"


# calculate growth
new_df$Gz = log(new_df$Size2/new_df$Size)

new_df$Infection = ifelse(new_df$Parasites > 0, "I+", "I-")

new_df$Infection[is.na(new_df$Infection)] = "I-"


new_df$Infection =  paste(new_df$Treatment, new_df$Infection, "-")
new_df$Infection = factor(new_df$Infection)
levels(new_df$Infection)
levels(new_df$Infection) = c("C", "I-", "I+")

new_df$W = factor(new_df$Week)

library(ggplot2)

str(new_df)


# Estimate how many individuals are in the Infection, Treatment, and Week groups
# This is to make sure that we have enough individuals to estimate the parameters
# of the model
# We can use the function "ddply" from the plyr package to do this

ddply(new_df, .(Treatment, W), summarise, N = sum(N))


ddply(new_df, .(Infection, Treatment, W), summarise, N = sum(N))

#

ggplot(new_df, aes(x = Size, y = Gz, colour = Infection)) + geom_point() +
  geom_smooth(method = "lm") + 
  scale_colour_manual(values = c("black", "orange", "red")) + xlab("Initial Size (mm)")+
  ylab("Individual growth ln(Final size/Initial size)")   + facet_grid(. ~ W)


new_df$z = log(new_df$Size) - log(500)
new_df$W = factor(new_df$Week)


ggplot(new_df, aes(x = z, y = Gz, colour = Infection)) + geom_point() +
  geom_smooth(method = "lm") + 
  scale_colour_manual(values = c("black", "orange", "red")) + xlab("Initial Size (mm)")+
  ylab("Individual growth ln(Final size/Initial size)")   + facet_grid(. ~ W)


# Model for all data ------------------------------------------------------
# Model for only individuals ----------------------------------------------

Ind_df = subset(new_df, N = 1)


png("plots/IndividualGrowth.png")
ggplot(Ind_df, aes(x = Size, y = Gz, colour = Infection)) + geom_point() +
  geom_smooth(method = "lm") + 
   scale_colour_manual(values = c("black", "orange", "red")) + xlab("Initial Size (mm)")+
   ylab("Individual growth ln(Final size/Initial size)")   + facet_grid(W ~ .)
graphics.off()

# mod2 = brm(Gz ~ 1 + z * Infection * W + (1|name/Block), Ind_df)
# summary(mod2)


mod3 = lmer(Gz ~ z * Infection * W + (1|name/Block), Ind_df)
summary(mod3)

Anova(mod3)

# For Nanda
# there is a significant 3 way interaction effect
# z:Infection:W  15.4976  4   0.003773 ** 

