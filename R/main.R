# test the effects of Schistocephalus infections on the functional response of copepods

# clean up the workspace
rm(list = ls())

# load the libraries
library(ggplot2)
library(dplyr)
library(cmdstanr)

library(rethinking)

# load the data from excel
library(readxl)

data <- read_xlsx("data/Batch_1_2.xlsx")

head(data)


dim(data)
unique(data$Treatment)
# Create unique individual IDs, and change the name of columns

data %>% 
  mutate(ID = paste0(Plate, Individual)) %>% 
  filter(Treatment != "infected LV", Treatment != "double infection") %>% 
  rename(R = prey, Infection = Treatment, Eaten = "prey consumed") %>%
  mutate(I = ifelse(Infection == "infected", 1,0), nI = ifelse(Infection == "exposed", 1,0)) -> data



dim(data)

summary(lm(Eaten ~ Infection, data = data))

head(model.matrix(~ Infection, data = data))

# plot Eaten ~ R, colored by Infection

ggplot(data, aes(x = R, y = Eaten, color = Infection)) +
  geom_point() +
  geom_smooth(method = "smooth", se = FALSE) +
  theme_bw() +  ylim(0, 20) + 
  labs(x = "Number of prey", y = "Number of prey consumed") -> p1

sort(unique(data$R))


# save the plot as a pdf file
ggsave("plots/Batch_1.pdf", p1, width = 5, height = 5)


# Fit a stan model using cmdstanr
# First, create a stan file
# Path: stan_models/functional_response.stan



factor(data$ID) -> data$ID

levels(data$ID) <- c(1:length(levels(data$ID)))

dlist = list(
    N = length(data$Eaten),
    Nid = length(unique(data$ID)),
    ID = as.integer(data$ID),
    R = data$R,
    Eaten = data$Eaten,
    I = data$I,
    i = data$nI 
    # pD = data$pD,
    # Egg = as.integer(data$Egg). 1, 0
    # cSize
    # ISize = data$I * data$cSize 
  )

# Now, fit the model
# Now, compile the stan model

model <- cmdstan_model("stan/FR_model.stan")


fit <- model$sample(
    data = dlist,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 2000,
    iter_sampling = 4000,
    refresh = 100,
    save_warmup = FALSE,
    show_messages = TRUE, 
    adapt_delta = 0.93,
    
    )

# get summary

precis(fit, depth = 1, digits = 3)

precis(fit, depth = 2)

# get posteriors from fit

post <- fit$draws(format = "df")

# estimate the level of support for the hypothesis.
LOS <- function(x){

    p = length(which(x> 0))/ length(x)

    return(p*100)


}

# is the effect of infection on atack rate positive?
LOS(post$a_I) # % that is positive
100 - LOS(post$a_I) # % that is negative

# is the effect of infection on handling time positive?
LOS(post$h_I) # % that is positive
100 - LOS(post$h_I) # % that is negative


# plot the functional response curves based on the posterior draws

# write predicted values function

p.link <- function(post, R, I, nI){
    
    A = exp(post$a + post$a_I*I + post$a_i*nI)

    H = exp(post$h + post$h_I*I + post$h_i*nI)

    B = exp(post$b)

    p <- (A * (R^B)) / (1 + A * H * (R^B))

    return(p)

}


# used sapply to apply the function to the posterior draws

R <- c(0, 2, 4, 8, 16, 32)
C <- sapply(1:length(R), function(i) p.link(post, R = R[i], I = 0, nI = 0))

I <- sapply(1:length(R), function(i) p.link(post, R = R[i], I = 1, nI = 0))

nI <- sapply(1:length(R), function(i) p.link(post, R = R[i], I = 0, nI = 1))


LOS(C - I)


# put the median and 95% HPD intervals in a data frame for each treatment (C, I, nI)


medians <- c(apply(C, 2, median), apply(I, 2, median), apply(nI, 2, median))
CI = rbind(t(apply(C, 2, HPDI, prob = 0.95)), t(apply(I, 2, HPDI, prob = 0.95)), t(apply(nI, 2, HPDI, prob = 0.95)))



df <- data.frame(
    R = rep(R, 3),
    median = medians,
    lower = CI[,1],
    upper = CI[,2],
    treatment = rep(c("C", "I", "nI"), each = length(R))

)


df

# plot the functional response curves using ggplot2, adding lower and upper 95% HPD intervals as ribbons


p2 <- ggplot(df, aes(x = R, y = median, color = treatment)) +
    geom_line() +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = treatment), alpha = 0.2) +
    theme_bw() +
    labs(x = "Number of prey", y = "Number of prey consumed") +
    ylim(0, 11) +
   # geom_smooth(method = "smooth", se = FALSE) +
    theme(legend.position = "none") 
  


p2

# save p2 as a pdf file

ggsave("plots/Batch_1.pdf", p2, width = 5, height = 5)
