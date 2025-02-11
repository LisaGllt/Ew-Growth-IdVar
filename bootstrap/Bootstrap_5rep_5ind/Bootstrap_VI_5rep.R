
library(readxl)
library(brms)
library(rethinking)
library(tidyverse)
library(tidybayes)
library(posterior)

# 1. Data importation ----

df_data <- read_excel("data/Data_growth_EPX_DMX.xlsx", sheet="Growth") |> 
  mutate(
    ID_f = as.factor(ID),
    exposition = as.factor(
      case_when(
        ID <= 30 ~ 0,
        ID > 30 ~ 1 
        )
    ),
    L = w^(1/3)
    )

df_data <- df_data[df_data$t<=56,] # data from t=0 to t=56

dead <- subset(df_data, t==14 & status == "D")$ID # Earthworms dead in the first 14 days
df_data <- df_data[!df_data$ID %in% dead,] # Removal of early deaths

# 2. Draws ----

###############################################################################
n_draws <- 500 # A CHANGER
n_rep <- 5
n_ind_cosm <- 5
###############################################################################

individuals <- data.frame(ID = unique(df_data$ID)) # Individuals that can be drawn
ctrl_individuals <- individuals[individuals$ID <=30,]
trt_individuals <- individuals[individuals$ID > 30,]

set.seed(121212)

# ctrl_individuals_draws <- t(replicate(n_draws, sample(ctrl_individuals, n_rep)))
# colnames(ctrl_individuals_draws) <- rep("Ctrl", n_rep)
# trt_individuals_draws <- t(replicate(n_draws, sample(trt_individuals, n_rep)))
# colnames(trt_individuals_draws) <- rep("Trt", n_rep)
# 
df_draws <- data.frame(ID_draws = seq(1:n_draws))
# 
# df_draws <- cbind(df_draws, ctrl_individuals_draws, trt_individuals_draws)

# Fonction pour effectuer le tirage aléatoire avec remise
rand_draws <- function(ids, n_rep, n_ind_cosm, n_draws) {
    replicate(
      n_draws, 
      replicate(
        n_rep, 
        sample(ids, n_ind_cosm, replace = TRUE), 
        simplify = FALSE
        ), 
      simplify = FALSE
      )
}

# Tirages aléatoires
ctrl_individuals_draws <- rand_draws(ctrl_individuals, n_rep, n_ind_cosm, n_draws)
trt_individuals_draws <- rand_draws(trt_individuals, n_rep, n_ind_cosm, n_draws)

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


fit.model <- function(data){
  m.i <- brm(data = data, 
             bf.mod.VI, 
             backend = "rstan",
             prior = priors,
             sample_prior = "yes",
             seed = 42,
             chains = 4, cores = 4,
             iter = 4000, warmup = 1000,
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
  
  print(i)
  
  # Data calculation
  
  ctrl_draws_i <- ctrl_individuals_draws[i]
  trt_draws_i <- trt_individuals_draws[i]
  
  df_data_i <- data.frame()
  
  # Mean of earthworms in one replicate
  for(rep in 1:n_rep){
    ID_rep <- c(ctrl_draws_i[[1]][[rep]], trt_draws_i[[1]][[rep]])
    df_data_rep_full <- subset(df_data, ID %in% ID_rep)
    df_data_rep  <- df_data_rep_full |> 
      aggregate(w ~ t + exposition, FUN = mean) |> 
      mutate(
        ID_draw = i,
        ID = rep,
        L = w^(1/3)
      )
    
    df_data_i <- rbind(df_data_i, df_data_rep)
  }
  
  if(i==1){
    m.i <- fit.model(df_data_i)
    m.1 <- m.i
  }else{
    m.i <- update(m.1, newdata = df_data_i)
  }
  
  #post_i <- as_draws_df(m.i)
  post_i <- posterior_summary(m.i)
  rhat_i <- brms::rhat(m.i)
  
  L0_i <- post_i[1, 1]
  a_ctrl_i <- post_i[2, 1]
  a_trt_i <- post_i[3, 1]
  Rhat_L0_i <- as.numeric(rhat_i[1])
  Rhat_a_ctrl_i <- as.numeric(rhat_i[2])
  Rhat_a_trt_i <- as.numeric(rhat_i[3])
  
  L0sd_i <- post_i[4, 1]
  asd_ctrl_i <- post_i[5, 1]
  asd_trt_i <- post_i[6, 1]
  Rhat_L0sd_i <- as.numeric(rhat_i[4])
  Rhat_asd_ctrl_i <- as.numeric(rhat_i[5])
  Rhat_asd_trt_i <- as.numeric(rhat_i[6])
  
  sigma_i <- post_i[7, 1]
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
  
  res <- cbind(df_draws[1:i,], 
                    L0_est, a_ctrl_est, a_trt_est, Rhat_L0, Rhat_a_ctrl, Rhat_a_trt,
                    L0sd_est, asd_ctrl_est, asd_trt_est, Rhat_L0sd, Rhat_asd_ctrl, Rhat_asd_trt,
                    sigma_est, Rhat_sigma)
  
  save(res, file = "Bootstrap_VI_5rep_seed121212.RData")
  
}



res_f <- cbind(df_draws, 
                  L0_est, a_ctrl_est, a_trt_est, Rhat_L0, Rhat_a_ctrl, Rhat_a_trt,
                  L0sd_est, asd_ctrl_est, asd_trt_est, Rhat_L0sd, Rhat_asd_ctrl, Rhat_asd_trt,
                  sigma_est, Rhat_sigma)

save(res_f, file = "Bootstrap_VI_5rep_seed121212_f.RData")

