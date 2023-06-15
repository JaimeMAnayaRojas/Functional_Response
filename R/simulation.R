## Here, I want to simulate the functional response i.e. comsumption rate of copepods
## I will use that data to test the stan model and check how many individuals do I need for the experiment

library(ggplot2)
library(dplyr)

# Now, let's create the simulation

# First, I will create a function that will simulate the functional response of copepods

FR_Holling <- function(a, h, R, b=2){
  
  FR <- (a * (R^b)) / (1 + a * h * (R^b))
  
  return(FR)
  
}

# test the function
R =  (seq(from = 0.1, to = 10, length = 100))
e <- FR_Holling(a = 10, h = 0.5, b= 1, R = R)

plot(R, e, type = "l", xlim = c(0, max(R)), ylim = c(0, 5))



# Load the MASS library for the mvrnorm function
library(MASS)

# Define the number of individuals per group
N <- 1000

# Define the means of a and h
mean_aC <- 10
mean_aI <- mean_aC* (1-0.04)
mean_ai <- mean_aC * (1-0.025)

mean_hC <- log(0.3)
mean_hI <- mean_hC * (1-0.025)
mean_hi <- mean_hC * (1-0.0.03)


# Define the covariance matrix for a and h
# The off-diagonal entries represent the covariance between a and h
# For example, a covariance of 0.05 will induce a correlation of 0.5 between a and h (assuming sd(a) = sd(h) = 0.1)
sigma <- matrix(c(0.1, 0.05, 0.05, 0.1), ncol=2)


# Generate the parameters for each group using a multivariate normal distribution
parameters_C <- mvrnorm(N, mu = c(mean_aC, mean_hC), Sigma = sigma) # control group
parameters_i <- mvrnorm(N, mu = c(mean_ai, mean_hi), Sigma = sigma) # I-
parameters_I <- mvrnorm(N, mu = c(mean_aI, mean_hI), Sigma = sigma) # I+

# Extract the a and h values from the generated parameters

dens(parameters_C[, 1], xlim = c(7, 12), ylim = c(0, 5), xlab = "a", ylab = "Density")
dens(parameters_i[, 1], add = TRUE)
dens(parameters_I[, 1], add = TRUE, col = "red")



a <- c(parameters_C[, 1], parameters_i[, 1], parameters_I[, 1])
h <- c(parameters_C[, 2], parameters_i[, 2], parameters_I[, 2])
h <- exp(h)

range(h)
# Put parameters in a data frame
df <- data.frame(ID = 1:length(a), a = a, h = h)

# Simulate values for R
R <- 0:25

# Apply the function to the data frame
df2 <- data.frame()
for (j in 1:length(R)){
  Eaten = (FR_Holling(a = df$a, h = df$h, b = 1, R = R[j]))
  ID = df$ID
  Exp = rep(c("C", "I-", "I+"), each = N)
  R2 = rep(R[j], length(ID))
  df2 <- rbind(df2, data.frame(ID, Exp, R2, Eaten))
}

# add random error to Eaten

df2$Eaten <- df2$Eaten + rnorm(nrow(df2), mean = 0, sd = 0.15)


df3 <- df2 %>% group_by(Exp, R2) %>% summarise(mean = mean(Eaten), sd = sd(Eaten))

df3
# plot the functional response for each group, add legend for each group, and sd for each group

ggplot(df3, aes(x = R2, y = mean, group = Exp, color = Exp)) + 
  geom_line() + geom_point() + 
  theme_bw() + 
  theme(legend.position = "none") + # now fix axis
  scale_x_continuous(breaks = seq(0, 25, 5), limits = c(0, 25)) +
  scale_y_continuous(breaks = seq(0, 12, 2), limits = c(0, 4)) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.1)




dim(df3)


sim_N <- function(df2, n){
      # which individuals are in each group?
    df2a <- df2 %>% group_by(ID, Exp) %>% summarise(Eaten = mean(Eaten))
    # filter the data for Exp == "C" and choose n individuals from df2a

    df2aC <- df2a %>% filter(Exp == "C")
    ID_C <- df2aC$ID[sample(nrow(df2aC), n)]

    # filter the data for Exp == "I-" and choose n individuals from df2a

    df2ai <- df2a %>% filter(Exp == "I-")
    ID_i <- df2ai$ID[sample(nrow(df2ai), n)]

    # filter the data for Exp == "I+" and choose n individuals from df2a

    df2aI <- df2a %>% filter(Exp == "I+")
    ID_I <- df2aI$ID[sample(nrow(df2aI), n)]

    # combine the IDs
    IDs <- c(ID_C, ID_i, ID_I)
    # filter the data for the selected IDs, and keep only R = c(0, 2, 4 , 8, 16, 32)

    df4 <- df2 %>% filter(ID %in% IDs, R2 %in% c(0, 2, 4 , 8, 16, 32))
    
    df4$ID <- as.factor(df4$ID)
    levels(df4$ID) <- c(1:length(levels(df4$ID)))
    return(df4)

}


## Now, let's test the stan model
# Load the rstan library
library(cmdstanr)
#df4 <- read.csv("data/15sim_data.csv")


###

# run stan in a loop for different values of n

n <- c(10, 15, 20, 25, 50, 150)

 # Compile the model
library(rethinking)

mod <- cmdstan_model("R/FR_modB2.stan", quiet = FALSE)
df_res <- data.frame()
i = 1

for( i in 1:length(n)){
    df <- sim_N(df2, n[i])
    df$ID <- as.factor(df$ID)

    dlist = list(
      N = length(df$Eaten),
      Nid = (length(unique(df$ID))),
      Eaten = df$Eaten,
      ID = as.integer(df$ID),
      nI = ifelse(df$Exp == "I-", 1, 0 ),
      I = ifelse(df$Exp == "I+", 1, 0 ),
      R = df$R2  
    )

   

    fit <- mod$sample(
      data = dlist, 
      seed = 123, 
      chains = 4, 
      parallel_chains = 4,
      refresh = 500 # print update every 500 iters
    )

    fit$summary()

    # draws x variables data frame
    draws_df <- fit$draws(format = "df")
    # Make a data frame with the estimated median 95% highest posterior density intervals, and the true values


    df_temp <- data.frame(
      Parameter = c("aC", "aI", "ai", "hC", "hI", "hi"),
      Median = exp(c(median(draws_df$a), median(draws_df$a + draws_df$a_I), median(draws_df$a + draws_df$a_i), 
          median(draws_df$h), median(draws_df$h + draws_df$h_I), median(draws_df$h + draws_df$h_i))),
      
      low = exp(c(HPDI(draws_df$a, prob = 0.95)[1], HPDI(draws_df$a + draws_df$a_I, prob = 0.95)[1], HPDI(draws_df$a + draws_df$a_i, prob = 0.95)[1],
      HPDI(draws_df$h, prob = 0.95)[1], HPDI(draws_df$h + draws_df$h_I, prob = 0.95)[1], HPDI(draws_df$h + draws_df$h_i, prob = 0.95)[1])),
      
      up = exp(c(HPDI(draws_df$a, prob = 0.95)[2], HPDI(draws_df$a + draws_df$a_I, prob = 0.95)[2], HPDI(draws_df$a + draws_df$a_i, prob = 0.95)[2],
      HPDI(draws_df$h, prob = 0.95)[2], HPDI(draws_df$h + draws_df$h_I, prob = 0.95)[2], HPDI(draws_df$h + draws_df$h_i, prob = 0.95)[2])),
      
      
      true = c(mean_aC, mean_aI, mean_ai, exp(mean_hC), exp(mean_hI), exp(mean_hi)),
      N = n[i]
      
    )

    df_res <- rbind(df_res, df_temp)

    
}

df_res


# plot Median and 95% HPDI for each parameter vs N, with true values, and parameters as facets

pt <- ggplot(df_res, aes(x = N, y = Median, ymin = low, ymax = up)) + 
  geom_point() + geom_errorbar() + facet_wrap(Parameter ~ ., scale = "free") + 
  geom_hline(aes(yintercept = true), color = "red") 

# save plot as pdf

ggsave("FR_modB2.pdf", pt, width = 20, height = 15, units = "cm")


