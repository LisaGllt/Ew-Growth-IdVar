library(tidyverse)
library(readxl)
library(brms)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# 1. Data importation ----

df_data <- read_excel("0_Growth_EPX_DMX/data/Data_growth_EPX_DMX.xlsx", sheet="Growth")
head(df_data)

df_data$ID_f <- as.factor(df_data$ID)
df_data <- df_data%>%
  mutate(L=w^(1/3))%>%
  mutate(exposition=0+1*(exposition=="Trt"))
df_data$exposition <- as.factor(df_data$exposition)

df_data <- df_data[df_data$t<=56,] # data from t=0 to t=56

dead <- c(1,9,28,37,41,43,57,59) # Earthworms dead in the first 14 days
df_data <- df_data[!df_data$ID %in% dead,] # Removal of early deaths

ggplot(data=df_data, aes(x=t, y=w, color=exposition, group=ID_f))+
  geom_point(alpha=0.5)+
  geom_line(lwd=0.8, alpha=0.3)+
  labs(x="Temps (j)", y=expression(paste("Masse (mg)")))+
  theme_bw(11)+
  theme(legend.position = "bottom")

# 2. Draws ----

###############################################################################
n_draws <- 10000 # A CHANGER
n_rep <- 10
###############################################################################

individuals <- data.frame(ID = unique(df_data$ID)) # Individuals that can be drawn
ctrl_individuals <- individuals[individuals$ID <=30,]
trt_individuals <- individuals[individuals$ID > 30,]

set.seed(121212)
ctrl_individuals_draws <- t(replicate(n_draws, sample(ctrl_individuals, n_rep)))
colnames(ctrl_individuals_draws) <- c("ctrl_1", "ctrl_2", "ctrl_3")
trt_individuals_draws <- t(replicate(n_draws, sample(trt_individuals, n_rep)))
colnames(trt_individuals_draws) <- c("trt_1", "trt_2", "trt_3")

df_draws <- data.frame(ID_draws = seq(1:n_draws))

df_draws <- cbind(df_draws, ctrl_individuals_draws, trt_individuals_draws)

# 3. Model ----

# Model definition
bf.mod.VI <- bf(L ~ L0+a*t,
                L0~1 + (1||ID), # no expo at t0 
                a~0+exposition+(0+exposition||ID),
                center=T,
                nl=T)

priors <- prior(normal(2.25, 0.45), class = b, nlpar = L0, lb = 0) + # CV 20%
  prior(normal(2.25*0.2, 2.25*0.2), class = sd, nlpar = L0)+
  # Mean values of a for the two cohorts
  prior(normal(0.075, 0.015), class=b, coef=exposition0, nlpar = a) + # CV 20%
  prior(normal(0.075, 0.015), class=b, coef=exposition1, nlpar = a) + # CV 20%
  # Individual variations values of a for the two cohorts
  prior(normal(0.075*0.2, 0.075*0.2), class=sd, 
        coef=exposition0, group=ID, nlpar = a)+
  prior(normal(0.075*0.2, 0.075*0.2), class=sd, 
        coef=exposition1, group=ID, nlpar = a)+
  # Résidus
  prior(exponential(1), sigma)


# Mode of a list
Mode <- function(x, na.rm = FALSE) {
  if(na.rm){
    x = x[!is.na(x)]
  }
  
  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}

fit.model <- function(data){
  m.i <- brm(data = data, 
             bf.mod.VI, 
             backend = "rstan",
             prior = priors,
             sample_prior = "yes",
             seed = 42,
             chains = 4, cores = 4,
             iter = 6000, warmup = 1000,
             #threads = threading(4),
             control = list(adapt_delta = .95,
                            max_treedepth = 12))
}

# 4. Boucle ----

Rhat_a_ctrl <- c()
Rhat_a_trt <- c()
Rhat_asd_ctrl <- c()
Rhat_asd_trt <- c()
Rhat_L0 <- c()
Rhat_L0sd <- c()
L0_est <- c()
a_ctrl_est <- c()
a_trt_est <- c()
L0sd_est <- c()
asd_ctrl_est <- c()
asd_trt_est <- c()

sigma_est <- c()
Rhat_sigma <- c()


for(i in 1:n_draws){

  draws_i <- df_draws[i, 2:length(df_draws)]
  df_data_draws <- df_data[df_data$ID %in% draws_i,]
  
  m.i <- fit.model(df_data_draws)
  
  post_i <- as_draws_df(m.i)
  rhat_i <- rhat(m.i)
  
  L0_i <- Mode(post_i$b_L0_Intercept)
  a_ctrl_i <- Mode(post_i$b_a_exposition0)
  a_trt_i <- Mode(post_i$b_a_exposition1)
  Rhat_L0_i <- as.numeric(rhat_i[1])
  Rhat_a_ctrl_i <- as.numeric(rhat_i[2])
  Rhat_a_trt_i <- as.numeric(rhat_i[3])
  
  L0sd_i <- Mode(post_i$b_L0_Intercept)
  asd_ctrl_i <- Mode(post_i$b_a_exposition0)
  asd_trt_i <- Mode(post_i$b_a_exposition1)
  Rhat_L0sd_i <- as.numeric(rhat_i[4])
  Rhat_asd_ctrl_i <- as.numeric(rhat_i[5])
  Rhat_asd_trt_i <- as.numeric(rhat_i[6])
  
  sigma_i <- Mode(post_i$sigma)
  Rhat_sigma_i <- as.numeric(rhat_i[7])
  
  L0_est <- rbind(L0_est, L0_i)
  a_ctrl_est <- rbind(a_ctrl_est, a_ctrl_i)
  a_trt_est <- rbind(a_trt_est, a_trt_i)
  Rhat_L0 <- rbind(Rhat_L0, Rhat_L0_i)
  Rhat_a_ctrl <- rbind(Rhat_a_ctrl, Rhat_a_ctrl_i)
  Rhat_a_trt <- rbind(Rhat_a_trt, Rhat_a_trt_i)
  
  L0sd_est <- rbind(L0sd_est, L0sd_i)
  asd_ctrl_est <- rbind(asd_ctrl_est, asd_ctrl_i)
  asd_trt_est <- rbind(asd_trt_est, asd_trt_i)
  Rhat_L0sd <- rbind(Rhat_L0sd, Rhat_L0sd_i)
  Rhat_asd_ctrl <- rbind(Rhat_asd_ctrl, Rhat_asd_ctrl_i)
  Rhat_asd_trt <- rbind(Rhat_asd_trt, Rhat_asd_trt_i)
  
  sigma_est <- rbind(sigma_est, sigma_i)
  Rhat_sigma <- rbind(Rhat_sigma, Rhat_sigma_i)
  
}

df_draws <- cbind(df_draws, 
                  L0_est, a_ctrl_est, a_trt_est, Rhat_L0, Rhat_a_ctrl, Rhat_a_trt,
                  L0sd_est, asd_ctrl_est, asd_trt_est, Rhat_L0sd, Rhat_asd_ctrl, Rhat_asd_trt,
                  sigma_est, Rhat_sigma)

save(df_draws, file = "Bootstrap_VI_10rep.RData")

