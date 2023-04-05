####  This is the script to analyses Nanda's project module data

# https://xcelab.net/rm/statistical-rethinking/

# Install the rethinking package
# https://github.com/rmcelreath/rethinking



library(rethinking)
rm(list=ls(all=TRUE))

data <- read.csv("data/FRdata_2.csv")
data$Eaten <- data$R - data$Left

# remove dead
str(data)
data <- data[-which(is.na(data$Eaten)),]

data$Expose = factor(data$Infection.status)

levels(data$Expose)
levels(data$Expose) <- c("E", "C","E","E") 


data$Infection = factor(data$Infection.status)
levels(data$Infection)
levels(data$Infection) <- c("I-", "C","I-","I+") 



# run the analyses in Batch 2

B4 <- subset(data, Batch == 4)



B4$Expose
B4$Exp = factor(B4$Expose)
levels(B4$Expose)
levels(B4$Exp) <- c(1,0) # recode the exposure to  exopse = 1 and non-expose or c = 0


B4$ID = factor(B4$Individual)

levels(B4$ID) <- 1:length(levels(B4$ID))


df <- B4[,c("Eaten", "R", "Exp", "ID")]

df$Exp = as.numeric(as.character(df$Exp))

# get formula for rethinking
#glimmer(Eaten ~ R + Exp + (1|ID), df)


dlist = list(
  N = length(df$Eaten),
  Nid = (length(levels(df$ID))),
  Eaten = df$Eaten,
  ID = as.integer(df$ID),
  Exp = df$Exp,
  R = df$R
  
)


mod2 <- stan("R/FR_mod.stan", data= dlist, iter = 1000,# warmup = 1000, 
             chains =4, cores = 4,
             control = list(adapt_delta = 0.9, max_treedepth = 13))


Females_Sum = precis(mod2, depth = 2, prob=.95)

write.csv(Females_Sum, "results/Females_Sum.csv")
# tracerplot(mod2)

post = extract.samples(mod2)

names(post)

length(post$a)

p_link <- function(post, Exp = 0, R){
  
  FR <- with(post, 
             
             ((10^(a + a_Exp*Exp)) * R) / (1 + (10^(a + a_Exp*Exp)) * (10^(h + h_Exp*Exp)) * R)
             
             )
  return(FR)
  
}

median(p_link(post, Exp= 1, 5))

mean(p_link(post, Exp= 1, 25))

xR = seq(from=0.1, to= 25, length= 100)


Exp0 <- sapply(1:length(xR), function(i) p_link(post, Exp = 0, xR[i]))
Exp1 <- sapply(1:length(xR), function(i) p_link(post, Exp = 1, xR[i]))


LOS <- function(x=NULL){
  
  out =   (length(which(x > 0)) / length(x))*100
  
  return(round(out, 3))
  
}


LOS(Exp1-Exp0)# 63.4% that the curve of exposure is higher than the control

apply(Exp1,2, mean)



sdf <- data.frame(R= c(xR,xR), 
                  Expose = c(rep("C",length(xR)),rep("E",length(xR))),
                  Eaten = c(apply(Exp0, 2, median),apply(Exp1, 2, median)),
                  lc = c(t(apply(Exp0, 2, HPDI, prob = .68))[,1], t(apply(Exp1, 2, HPDI, prob = .68))[,1]),
                  uc = c(t(apply(Exp0, 2, HPDI, prob = .68))[,2], t(apply(Exp1, 2, HPDI, prob = .68))[,2]))


C = ggplot(B4, aes(x = R, y = Eaten, group= Expose, fill = Expose)) + 
  geom_point(size = 1.5, shape = 21, position = position_jitterdodge(dodge.width=1.2), alpha = 0.2, show.legend = F) + 
  geom_line(data=sdf, aes(x=R, y=Eaten, colour = Expose), size = 1) + 
  geom_ribbon(data=sdf, aes(ymin = lc, ymax = uc, fill= Expose, group= Expose,), 
              alpha = 0.2, size = 0.01, show.legend = F) +
  scale_fill_manual(values=c("deepskyblue",
                             "black"))+
  scale_color_manual(values=c("deepskyblue",
                              "black")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ylab("Eaten prey \n (mean, 68% CI)") +
  xlab("Prey density") + ylim(c(0,13))


## Results females
LOS(10^(post$a + post$a_Exp) - 10^(post$a)) # 71.7% that expose copepods had a higher attack rate
LOS(10^(post$h + post$h_Exp) - 10^(post$h)) # 86.85% that exposed copepods have a higher handling time... 



#  then handling time of the expose copepods is a bit longer.

### Model for the males

B2 <- subset(data, Batch == 2)
B2$Infection


B2$ID = factor(B2$Individual)

levels(B2$ID) <- 1:length(levels(B2$ID))


df <- B2[,c("Eaten", "R", "Infection", "ID")]


# get formula for rethinking


dlist = list(
  N = length(df$Eaten),
  Nid = (length(levels(df$ID))),
  Eaten = df$Eaten,
  ID = as.integer(df$ID),
  E_nIf = ifelse(df$Infection == "I-", 1, 0 ),
  E_If = ifelse(B2$Infection == "I+", 1, 0 ),
  R = df$R
  
)


modB2 <- stan("R/FR_modB2.stan", data= dlist, iter = 1000,# warmup = 1000, 
             chains =4, cores = 4,
             control = list(adapt_delta = 0.9, max_treedepth = 13))


malesSum = precis(modB2, depth = 2, prob = .95)
write.csv(malesSum, "results/Males_vals.csv")
# tracerplot(mod2)

post = extract.samples(modB2)

names(post)

length(post$a)

p_linkINF <- function(post, nIf= 0, If = 0,  R){
  
  FR <- with(post, 
             
             ((10^(a + a_E_If*If + a_E_nIf*nIf)) * R) / (1 + (10^(a + a_E_If*If + a_E_nIf*nIf)) * (10^(h + h_E_If*If + h_E_nIf*nIf)) * R)
             
  )
  return(FR)
  
}

median(p_linkINF(post, nIf = 0,If = 0, 5))


xR = seq(from=0.1, to= 25, length= 100)


CM <- sapply(1:length(xR), function(i) p_linkINF(post, nIf = 0, If = 0, xR[i]))
nIfM <- sapply(1:length(xR), function(i) p_linkINF(post, nIf = 1, If = 0 , xR[i]))
IfM <- sapply(1:length(xR), function(i) p_linkINF(post, nIf = 0, If = 1, xR[i]))



LOS(nIfM-CM)# 63.4% that the curve of exposure is higher than the control
LOS(IfM-CM)# 83.72% that the curve of infected curve is higher than the control

LOS((IfM-CM) - (nIfM-CM))# 79% that the curve of infected curve is higher than the non infected

apply(Exp1,2, mean)



sdf <- data.frame(R= c(xR,xR,xR), 
                  Infection = c(rep("C",length(xR)),rep("I-",length(xR)), rep("I+",length(xR))),
                  Eaten = c(apply(CM, 2, median),apply(nIfM, 2, median),apply(IfM, 2, median)),
                  lc = c(t(apply(CM, 2, HPDI, prob = .68))[,1], t(apply(nIfM, 2, HPDI, prob = .68))[,1] , t(apply(IfM, 2, HPDI, prob = .68))[,1]),
                  uc = c(t(apply(CM, 2, HPDI, prob = .68))[,2], t(apply(nIfM, 2, HPDI, prob = .68))[,2] , t(apply(IfM, 2, HPDI, prob = .68))[,2]))


MP = ggplot(B2, aes(x = R, y = Eaten, group= Infection, fill = Infection)) + 
  geom_point(size = 1.5, shape = 21, position = position_jitterdodge(dodge.width=1.2), alpha = 0.2, show.legend = F) + 
  geom_line(data=sdf, aes(x=R, y=Eaten, colour = Infection), size = 1) + 
  geom_ribbon(data=sdf, aes(ymin = lc, ymax = uc, fill= Infection, group= Infection,), 
              alpha = 0.2, size = 0.01, show.legend = F) +
  scale_fill_manual(values=c("deepskyblue",
                             "black", "red"))+
  scale_color_manual(values=c("deepskyblue",
                              "black", "red")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ylab("Eaten prey \n (mean, 68% CI)") +
  xlab("Prey density") + ylim(c(0,13))

MP
### Results Males
LOS(10^(post$a + post$a_E_If) - 10^(post$a)) # 52.8% that the infected have higher attack rate.  there are no major difference in the attack rate
LOS(10^(post$h + post$h_E_If) - 10^(post$h)) # 20.95% that the infected have a higher handing time. no differences

LOS(10^(post$a + post$a_E_nIf) - 10^(post$a)) # 69.95% that the expose but not infected have higher attack rate.  there are no major difference in the attack rate
LOS(10^(post$h + post$h_E_nIf) - 10^(post$h)) # 56.5% that the expose but not infected have a higher handing time. no differences




###

library("ggpubr")
 theme_set(
    theme_bw() +
 theme(legend.position = "top")
 )



figure <- ggarrange(C, MP,
                    labels = c("A)", "B)"),
                    ncol = 2, nrow = 1)
figure


svg("plots/PlotB2 and B4.svg", width = 7, height = 5 )
figure 
graphics.off()
