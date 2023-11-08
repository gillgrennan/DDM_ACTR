# Simulation 4: multiple participants, vary dataset size


library(dplyr)
library(ggplot2)
library(rtdists)
library(purrr)
library(furrr)
library(tidyr)
library(truncdist)
library(cowplot)
library(grid)
library(rlang)
library(ggsci)

future::plan("multicore", workers = 6) # Set to desired number of cores

theme_paper <- theme_classic(base_size = 24) + 
  theme(axis.text = element_text(colour = "black"))


use_cached_results <- FALSE # Set to FALSE to rerun simulations (takes time!)


n_participants <- 25

set.seed(2021)

#Set the ACT-R parameters:

# Latency factor F
lf <- rtrunc(n_participants, spec = "norm", mean = 1, sd = .5, a = 0, b = Inf)
lf_range <- rtrunc(n_participants, spec = "norm", mean = .1, sd = .05, a = 0, b = Inf)
lf_lower <- lf - .5 * lf_range
lf_upper <- lf + .5 * lf_range

# Non-retrieval time t_er
t_er <- rtrunc(n_participants, spec = "norm", mean = .75, sd = .5, a = 0, b = Inf)

# Activation of correct answer
a_c_mu <- rtrunc(n_participants, spec = "norm", mean = -.5, sd = .5, a = -Inf, b = 0)
a_c_sd <- rep(1, n_participants)

# Activation of incorrect answer
a_f_mu <- rtrunc(n_participants, spec = "norm", mean = -1.5, sd = .5, a = -Inf, b = 0)
a_f_sd <- rep(1, n_participants)

dataset_sizes <- c(25, 50, 100, 250, 500, 1000, 2500, 5000)
dataset_sizes_short <- c(25, 50, 100)

ll_diffusion <- function(pars, rt, response)
{
  densities <- ddiffusion(rt, response=response,
                          a=pars[1], v=pars[2], t0=pars[3],z=pars[7],d = 0,
                          sz=pars[4],sv=pars[5],
                          st0=pars[6])
  if (any(densities == 0)) return(1e6)
  return(-sum(log(densities)))
}


  #s4_param_recov <- tibble()
  
  for(i in 1:length(dataset_sizes_short)) {
    
    # Generate dataset
    n_trials <- dataset_sizes[i]
    
    sim_actr <- tibble(participant = rep(1:n_participants, each = n_trials),
                       f = runif(n_participants * n_trials, min = rep(lf_lower, each = n_trials), max = rep(lf_upper, each = n_trials)),
                       a_c = rnorm(n_participants * n_trials, mean = rep(a_c_mu, each = n_trials), sd = rep(a_c_sd, each = n_trials)),
                       a_f = rnorm(n_participants * n_trials, mean = rep(a_f_mu, each = n_trials), sd = rep(a_f_sd, each = n_trials)),
                       t = rep(t_er, each = n_trials)
    ) %>%
      rowwise() %>%
      mutate(rt = f * exp(-max(a_c, a_f)) + t,
             response = ifelse(a_c > a_f, "upper", "lower")) %>%
      ungroup()
    
    head(sim_actr)
    prop.table(table(sim_actr$response))
    
    # Fit each participant separately
    Q1 <- quantile(sim_actr$rt, .25)
    Q3 <- quantile(sim_actr$rt, .75)
    IQR <- IQR(sim_actr$rt)
    
    sim_actr_no_outliers <- subset(sim_actr,sim_actr$rt> (Q1-1.5*IQR) & sim_actr$rt< (Q3 + 1.5*IQR))
    
    sim_actr_split <- split(sim_actr_no_outliers, sim_actr_no_outliers$participant)
    
    dat <- list()  	
    dat <- vector("list",length =n_participants) #change per subject #
    iterations <- list()
    
    for (d in 1:n_participants) {
      data <- sim_actr_split[[d]]
      start1 <- c(runif(1,0.5,1),runif(1,0,2),0.1,runif(3,0,0.5)) 
      z1 <- start1[1]*0.5
      start <- c(start1,z1)
      names(start) <- c("a","v","t0","sz","sv","st0","z")
      recov <- nlminb(start,ll_diffusion,lower=0, rt=data$rt,response=data$response)
      if (recov$par[2] == 0) {
        recov <- nlminb(start,ll_diffusion,lower=0, rt=data$rt,response=data$response)
      }
      if (recov$par[2] == 0) {
        recov <- nlminb(start,ll_diffusion,lower=0, rt=data$rt,response=data$response)
      }
      if (recov$par[2] == 0) {
        recov <- nlminb(start,ll_diffusion,lower=0, rt=data$rt,response=data$response)
      }
      dat[[d]] <- c(round(recov$par,5))
      iterations[[d]] <- recov$iterations

    }
      
      #mutate(dataset_size = dataset_sizes[i])
    
    s4_param_recov <- bind_rows(s4_param_recov, dat)
    
  }

## run lines if need to redo any dataset sizes
  # replace_count <- 75 
  # s4_param_recov[1:replace_count, ] <- s4_param_recov[(275 - replace_count + 1):275, ]
  # s4_param_recov <- s4_param_recov[1:(nrow(s4_param_recov)-replace_count), ]
  # 
  s4_param_recov$dataset_size = c(rep(dataset_sizes[1],each = n_participants),rep(dataset_sizes[2],each = n_participants),rep(dataset_sizes[3],each = n_participants),rep(dataset_sizes[4],each = n_participants),rep(dataset_sizes[5],each = n_participants),rep(dataset_sizes[6],each = n_participants),rep(dataset_sizes[7],each = n_participants),rep(dataset_sizes[8],each = n_participants))
  s4_param_recov$participant = c(1:25,1:25,1:25,1:25,1:25,1:25,1:25,1:25)

  # parameter recovery 
  f_recov <- numeric()
  t_recov <- numeric()
  a_c_recov <- numeric()
  a_f_recov <- numeric()
  ratio1 <- numeric()
  a_diff <- numeric()
  
  dat <- s4_param_recov

  # updated equations as of 10/16/2023 
    results <- dat
    f_recov <- results$a/2
    t_recov <- results$t0
    a_c_recov <- log(results$v)
    a_f_recov <- log(results$v)-(2*results$a*results$v)

s4_param_recov_long <- s4_param_recov %>%
  transmute(
    n = dataset_size,
    participant,
    A_c_recov = log(v),
    A_f_recov = log(v)-(2*a*v),
    t_er = t0,    
    F = a/2
  ) %>%
  pivot_longer(A_c_recov:F, names_to = "Parameter", values_to = "Recovered")

s4_param_orig_long <- tibble(n = rep(c(dataset_sizes, 100), each = n_participants),
                             participant = rep(1:n_participants, length(dataset_sizes)+1),
                             A_c_recov = rep(a_c_mu, length(dataset_sizes)+1),
                             A_f_recov = rep(a_f_mu, length(dataset_sizes)+1),
                             t_er = rep(t_er, length(dataset_sizes)+1),
                             F = rep(lf, length(dataset_sizes)+1)) %>%
  pivot_longer(A_c_recov:F, names_to = "Parameter", values_to = "Original")
s4_param_orig_long <- head(s4_param_orig_long, - 100)
s4_par_comp <- left_join(s4_param_recov_long, s4_param_orig_long, by = c("n", "participant", "Parameter"))

#Comparison of original ACT-R parameters to recovered parameter estimates from DDM:
ggplot(s4_par_comp, aes(x = Original, y = Recovered, colour = Parameter)) +
  facet_wrap(~ n, scales = "free", labeller = function(x) label_both(labels = x, sep = " = ")) +
  geom_abline(lty = 2) +
  geom_point() +
  scale_color_d3(palette = "category20") +
  labs(x = "Original parameter value",
       y = "Recovered parameter value",
       colour = "Parameter") +
  theme_paper


#Difference between the original ACT-R parameters and the recovered parameter estimates from the DDM, as a function of dataset size:
s4_par_comp_error <- s4_par_comp %>%
  mutate(Parameter = factor(Parameter, 
                            levels = c("A_c_recov","A_f_recov","F", "t_er"),
                            labels  = c(expression(A[c]),expression(A[f]), expression(bar(F)), expression(t[er]))))
s4_par_comp_error$ae <- abs(s4_par_comp_error$Original - s4_par_comp_error$Recovered)
s4_par_comp_error$ae[s4_par_comp_error$ae == Inf] <- 5

s4_par_comp_error_param <- s4_par_comp_error %>%
  group_by(n, Parameter) %>%
  #summarise(sdae = sd(ae)/sqrt(25))
  summarise(mae = mean(ae))
  
s4_par_comp_error_mean <- s4_par_comp_error %>%
  group_by(n, Parameter) %>%
  summarise(mae = mean(ae)) %>%
  group_by(n) %>%
  #summarise(mae = mean(mae))
  summarise(sdae = sd(mae)/sqrt(4))


### absolute error plot!! 
ggplot() +
  geom_jitter(data = s4_par_comp_error, aes(x = n, y = ae, colour = Parameter), height = 0, width = .025, alpha = .05) +
  geom_line(data = s4_par_comp_error_param, aes(x = n, y = mae, group = Parameter, colour = Parameter), lty = 2, size = 1.5) +
  geom_point(data = s4_par_comp_error_param, aes(x = n, y = mae, group = Parameter, colour = Parameter), size = rel(2), alpha = .9) +
  geom_line(data = s4_par_comp_error_mean, aes(x = n, y = mae), lty = 2, lwd = rel(2), size = 1.5) +
  geom_point(data = s4_par_comp_error_mean, aes(x = n, y = mae), size = rel(2)) +
  scale_x_log10(breaks = c(25, 50, 100, 250, 500, 1000, 2500, 5000),
                labels = c(25, 50, 100, 250, 500, "1K", "2.5K", "5K")) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_color_d3(palette = "category20") +
  coord_cartesian(ylim = c(0, 3.8)) +
  labs(x = "Number of trials",
       y = "Absolute error",
       colour = NULL) +
  theme_paper +
  theme(axis.text.x = element_text(face = c("plain", "plain", "plain", "plain", "plain", "plain", "plain", "plain")),
        legend.background = element_blank(),
        legend.direction = "horizontal",
        legend.position = "top",
        legend.justification = "right",
        legend.box.margin = unit(c(0, 0, -40, 0), "pt"))


#correlation between original and recovered parameter values:
s4_par_comp_cor <- s4_par_comp %>%
  mutate(parameter = factor(Parameter, 
                            levels = c("A_c_recov","A_f_recov","F", "t_er"),
                            labels  = c(expression(A[c]), expression(A[f]), expression(bar(F)), expression(t[er])))) %>%
  group_by(n, Parameter) %>%
  #summarise(r = abs(original-recovered))
  summarise(r = abs(cor(Original, Recovered, method = "spearman")))

s4_par_comp_cor_mean <- s4_par_comp_cor %>%
  group_by(n) %>%
  summarise(r = mean(r))

ggplot() +
  #geom_rect(aes(xmin = exp(log(250) - .25), xmax = exp(log(250) + .25), ymin = 0, ymax = 2.4), colour = NA, fill = "grey90") +
  geom_jitter(data = s4_par_comp_cor, aes(x = n, y = r, colour = Parameter), height = 0, width = .025, alpha = .05) +
  geom_line(data = s4_par_comp_cor, aes(x = n, y = r, group = Parameter, colour = Parameter), lty = 2, size = 1.5) +
  geom_point(data = s4_par_comp_cor, aes(x = n, y = r, group = Parameter, colour = Parameter), size = rel(2), alpha = .9) +
  geom_line(data = s4_par_comp_cor_mean, aes(x = n, y = r), lty = 2, lwd = rel(2), size = 1.5) +
  geom_point(data = s4_par_comp_cor_mean, aes(x = n, y = r), size = rel(2)) +
  scale_x_log10(breaks = c(25, 50, 100, 250, 500, 1000, 2500, 5000),
                labels = c(25, 50, 100, 250, 500, "1K", "2.5K", "5K")) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_color_d3(palette = "category20") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Number of trials",
       y = "Correlation",
       colour = NULL) +
  theme_paper +
  theme(axis.text.x = element_text(face = c("plain", "plain", "plain", "plain", "plain", "plain", "plain", "plain")),
        legend.background = element_blank(),
        legend.direction = "horizontal",
        legend.position = "top",
        legend.justification = "right",
        legend.box.margin = unit(c(0, 0, -40, 0), "pt"))


