---
# final equations for ms version 10/06/2023! woo  
  ---
  
  ## Setup
  
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
library(GGally)
library(scales)
library(ggh4x)
library(ggsci)

future::plan("multicore", workers = 6) # Set to desired number of cores

theme_paper <- theme_classic(base_size = 20) + 
  theme(axis.text = element_text(colour = "black"))


use_cached_results <- FALSE # Set to FALSE to rerun simulations (takes time!)

## Simulation 3: 25 participants, 100  trials (typical empirical task size)

## Generate data

set.seed(2021)

# Number of trials to simulate per participant:
n_trials <- 100

# Number of participants to simulate:
n_participants <- 25

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


#Distribution of parameters:
tibble(participant = 1:n_participants, lf_lower, lf_upper, t_er, a_c_mu, a_c_sd, a_f_mu, a_f_sd) %>%
  pivot_longer(-participant, "parameter") %>%
  ggplot(aes(x = parameter, y = value, colour = parameter)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(height = 0, width = .1) +
  labs(x = "Parameter",
       y = "Value",
       colour = "Parameter") +
  guides(colour = "none") +
  theme_paper


#Alternative visualization:
tibble(participant = 1:n_participants, lf_lower, lf_upper, t_er, a_c_mu, a_c_sd, a_f_mu, a_f_sd) %>%
  pivot_longer(-participant, "parameter") %>%
  mutate(parameter = factor(parameter, 
                            levels = c("a_c_mu", "a_f_mu", "a_c_sd", "a_f_sd", "lf_lower", "lf_upper", "t_er"),
                            labels  = c(expression(mu[c]), expression(mu[f]), expression(sigma[c]), expression(sigma[f]), expression(F[a]), expression(F[b]), expression(t[er])))) %>%
  ggplot(aes(y = parameter, x = value, colour = (parameter))) +
  facet_grid(parameter ~ ., scales = "free", switch = "y", labeller = labeller(parameter = label_parsed)) +
  geom_hline(aes(yintercept = parameter), colour = "grey90", lty = 3) +
  geom_jitter(height = .1, width = 0, size = .5) +
  labs(x = NULL,
       y = NULL) +
  guides(colour = "none") +
  scale_colour_viridis_d() +
  theme_paper +
  theme(strip.background = element_blank(),
        strip.text.y.left = element_text(angle = 0, size = rel(1.25)),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

#Generate the data (rt/response):
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

#write.csv(sim_actr,file = "/Users/gilliangrennan/Desktop/data.csv")

## outlier removal for rt -- making unreasonably large v or a 
Q1 <- quantile(sim_actr$rt, .25)
Q3 <- quantile(sim_actr$rt, .75)
IQR <- IQR(sim_actr$rt)

sim_actr_no_outliers <- subset(sim_actr,sim_actr$rt> (Q1-1.5*IQR) & sim_actr$rt< (Q3 + 1.5*IQR))

#RT distributions (figure):
ggplot(mutate(sim_actr_no_outliers, rt = ifelse(response == "upper", rt, -rt)), aes(x = rt, group = participant, colour = participant)) +
  geom_line(stat = "density", adjust = .5, n = 2^12, na.rm = F, alpha = .5) +
  geom_vline(aes(xintercept = 0), lty = 2, colour = "grey80") +
  scale_x_continuous(breaks = c(-15, 0, 15)) +
  scale_y_continuous(expand = c(0,0), breaks = NULL) +
  coord_cartesian(xlim = c(-20, 20)) +
  labs(x = "RT (s)",
       y = "Density") +
  guides(colour = "none") +
  theme_paper +
  theme(strip.background = element_blank(),
        strip.text = element_blank())

ll_diffusion <- function(pars, rt, response)
{
  densities <- ddiffusion(rt, response=response,
                          a=pars[1], v=pars[2], t0=pars[3],z=pars[7],d = 0,
                          sz=pars[4],sv=pars[5],
                          st0=pars[6])
  if (any(densities == 0)) return(1e6)
  return(-sum(log(densities)))
}

# run drift diffusion model
sim_actr_split <- split(sim_actr_no_outliers, sim_actr_no_outliers$participant) #split by participant

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
  dat[[d]] <- c(round(recov$par,5))
  iterations[[d]] <- recov$iterations
  
}

# parameter recovery 
f_recov <- numeric()
t_recov <- numeric()
a_c_recov <- numeric()
a_f_recov <- numeric()
ratio1 <- numeric()
a_diff <- numeric()

for (x in 1:n_participants) {
  results <- as.data.frame(t(dat[[x]]))
  t_recov <- append(t_recov,results$t0) #t0 = ter 
  af1 <- log(results$v)-(2*results$a*results$v) #af = ln(v)-2av
  #af1 <- log((results$v/(exp(2*results$a*results$v)))+1-results$v)
  ac1 <- log(results$v) #ac = ln(v)
  f_recov1 <- results$a/2 # f = a-z = 1/2(a)
  f_recov <- append(f_recov,f_recov1) 
  a_c_recov <- append(a_c_recov,ac1)
  a_f_recov <- append(a_f_recov,af1)
}

dat1 <- data.frame(matrix(unlist(dat),nrow=length(dat),byrow=TRUE))
dat1$participant <- seq(1,n_participants) 
colnames(dat1) <- c('a','v','t0','sz','st0','sv','z','participant')

s3_par_comp <- tibble(
  participant = rep(1:n_participants, 4),
  parameter = rep(c("F", "t_er", "A_correct", "A_error"), each = n_participants),
  original = c(lf, t_er, a_c_mu, a_f_mu),
  recovered = c(f_recov, t_recov, a_c_recov, a_f_recov)
)

s3_par_comp_plot <- s3_par_comp %>%
  filter(parameter != "A_correct_sd") %>%
  mutate(parameter = factor(parameter, 
                            levels = c("A_correct","A_error","F", "t_er"),
                            labels  = c(expression(A[c]), expression(A[f]), expression(bar(F)), expression(t[er])))) 
#filter(!parameter %in% c("F[a]", "F[b]")) %>%

p_par_comp <- ggplot(s3_par_comp_plot, aes(x = original, y = recovered)) +
  facet_wrap(~parameter, ncol = 4, scales = "free", labeller = labeller(parameter = label_parsed))+
  geom_abline(lty = 2, colour = "grey80") +
  geom_smooth(method = "lm", formula = y ~ x, alpha = .5, colour = "black", lwd = rel(.8)) +
  geom_point(aes(colour = parameter), size = 3, alpha = .5) +
  labs(x = "Original",
       y = "Recovered",
       colour = NULL) +
  #scale_x_continuous(limits = c(-2, max(s3_par_comp$recovered))) +
  scale_color_d3(palette = "category20") +
  guides(colour = "none") +
  theme_paper +
  theme_paper +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = rel(1.5)))

p_par_comp


### leaving out drift
s3_par_comp <- tibble(
  participant = rep(1:n_participants, 2),
  parameter = rep(c("F", "t_er"), each = n_participants),
  original = c(lf, t_er),
  recovered = c(f_recov, t_recov)
)

s3_par_comp_plot <- s3_par_comp %>%
  filter(parameter != "A_correct_sd") %>%
  mutate(parameter = factor(parameter, 
                            levels = c("F", "t_er"),
                            labels  = c(expression(bar(F)), expression(t[er])))) 
#filter(!parameter %in% c("F[a]", "F[b]")) %>%

p_par_comp <- ggplot(s3_par_comp_plot, aes(x = original, y = recovered)) +
  facet_wrap(~parameter, ncol = 2, scales = "free_x", labeller = labeller(parameter = label_parsed))+
  geom_abline(lty = 2, colour = "grey80") +
  geom_smooth(method = "lm", formula = y ~ x, alpha = .5, colour = "black", lwd = rel(.8)) +
  geom_point(aes(colour = parameter), alpha = .5) +
  labs(x = "Original",
       y = "Recovered",
       colour = NULL) +
  scale_x_continuous(limits = c(min(s3_par_comp$recovered), max(s3_par_comp$recovered))) +
  scale_colour_viridis_d() +
  guides(colour = "none") +
  theme_paper +
  theme_paper +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = rel(1)))

p_par_comp

# alternative visualization
ggplot(s3_par_comp_plot, aes(x = original, y = recovered)) +
  facet_wrap(~ parameter, ncol = 5, labeller = labeller(parameter = label_parsed))+
  geom_abline(lty = 2, colour = "grey80") +
  geom_smooth(method = "lm", formula = y ~ x, alpha = .5, colour = "black", lwd = rel(.8)) +
  geom_point(aes(colour = parameter), alpha = .5) +
  coord_fixed() +
  labs(x = "Original",
       y = "Recovered",
       colour = NULL) +
  scale_x_continuous(limits = c(min(s3_par_comp$recovered), max(s3_par_comp$recovered))) +
  scale_colour_viridis_d() +
  guides(colour = "none") +
  theme_paper +
  theme_paper +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = rel(1)))


library(scales)
library(ggh4x)


## another alternative visualization 
par_limits <- s3_par_comp_plot %>%
  group_by(parameter) %>%
  summarise(limits = list(c(min(c(original, recovered)), max(c(original, recovered)))))

p_par_comp <- ggplot(s3_par_comp_plot, aes(x = original, y = recovered)) +
  facet_wrap(~ parameter, ncol = 3, scales = "free", labeller = labeller(parameter = label_parsed)) +
  geom_abline(lty = 2, colour = "grey80") +
  geom_smooth(method = "lm", formula = y ~ x, alpha = .2, colour = "black", lwd = rel(.8)) +
  geom_point(aes(colour = parameter), alpha = .9) +
  labs(x = "Original",
       y = "Recovered",
       colour = NULL) +
  facetted_pos_scales(
    x = list(
      scale_x_continuous(limits = par_limits$limits[[1]]),
      scale_x_continuous(limits = par_limits$limits[[2]]),
      scale_x_continuous(limits = par_limits$limits[[3]])
    ),
    y = list(
      scale_y_continuous(limits = par_limits$limits[[1]]),
      scale_y_continuous(limits = par_limits$limits[[2]]),
      scale_y_continuous(limits = par_limits$limits[[3]])
    )) +
  scale_colour_viridis_d() +
  guides(colour = "none") +
  theme_paper +
  theme_paper +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = rel(1)),
        aspect.ratio = 1)

p_par_comp


########### 

#Compare fitted DDM density to actual data:
ddm_dat <- crossing(rt = seq(0, 20, by = .01),
                    response = c("upper", "lower"))

sim_ddm <- dat1 %>%
  split(.$participant) %>% 
  future_map_dfr(function (x) {
    bind_cols(ddm_dat, density =
                ddiffusion(rt = ddm_dat$rt,
                           response = ddm_dat$response, ## change to response = 'upper'?
                           a = x$a,
                           v = x$v,
                           t0 = x$t0,
                           z = x$z,
                           st0 = x$st0,
                           sv = x$sv,
                           sz = x$sz)) %>%
      mutate(participant = x$participant[1])
  }) %>%
  mutate(rt = ifelse(response == 'upper', rt, -rt),
         model = "DDM")

sim_actr_no_outliers %>%
  mutate(rt = ifelse(response == 'upper', rt, -rt),
         model = "Data") %>%
  ggplot(aes(x = rt, colour = model)) +
  facet_wrap(~ participant, ncol = 5) +
  geom_vline(xintercept = 0, lty = 2, colour = "grey80") +
  geom_histogram(aes(y = ..density..), binwidth = .5, colour = "black", fill = "white", size = .1) +
  geom_line(data = sim_ddm, aes(y = density)) +
  scale_x_continuous(limits = c(-20, 20), breaks = c(-15, 0, 15)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_colour_manual(values = c("#e66101", "#5e3c99")) +
  labs(x = "RT (s)",
       y = "Density",
       colour = NULL) +
  guides(colour = FALSE) +
  theme_paper +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank())

#Alternative plot for presentation:
draw_key_custom <- function(data, params, size) {
  if(data$colour == "#000000" && data$size == .1) { # ACT-R
    grobTree(
      linesGrob(
        c(.1, .1, .3, .3, .3, .5, .5, .5, .7, .7, .7, .9, .9),
        c(0, .5, .5, 0, .8, .8, 0, .65, .65, 0, .4, .4, 0)
      ),
      gp = gpar(
        col = data$colour %||% "grey20",
        fill = alpha(data$fill %||% "white", data$alpha),
        lwd = (data$size %||% 0.5) * .pt,
        lty = data$linetype %||% 1
      )
    )
  } 
  else if (data$colour == "#e66101") { # LBA
    grobTree(
      linesGrob(
        c(0, 1),
        c(.5, .5)
      ),
      gp = gpar(
        col = alpha(data$colour %||% "grey20", data$alpha),
        fill = alpha(data$fill %||% "white", data$alpha),
        lwd = (data$size %||% 0.5) * .pt,
        lty = data$linetype %||% 1
      )
    )
  }
  else {
    grobTree() # Don't draw
  }
}

### figure for mapping 1 model participant (cdf vs. actual RTs)
test <- 3
sim_actr %>%
  filter(participant == test) %>%
  mutate(rt = ifelse(response == "upper", rt, -rt),
         model = "ACT-R") %>%
  ggplot(aes(x = rt, colour = model)) +
  facet_wrap(~ participant, ncol = 1) +
  geom_vline(xintercept = 0, lty = 2, colour = "grey80") +
  geom_histogram(aes(y = ..density..), binwidth = .5, fill = "white", size = .1, key_glyph = draw_key_custom) +
  geom_line(data = filter(sim_ddm, participant == test), aes(y = density), key_glyph = draw_key_custom) +
  scale_x_continuous(limits = c(-20, 20), breaks = c(-15, 0, 15)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_colour_manual(values = c("#000000", "#e66101")) +
  labs(x = "RT (s)",
       y = "Density",
       colour = NULL) +
  theme_paper +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "top",
        legend.justification = "left",
        legend.direction = "vertical",
        legend.box.margin = unit(c(-20, 0, -40, -30), "pt"))


#Plot correlation between parameters:
library(GGally)

p_par_cor <- s3_par_comp_plot %>%
  select(-original) %>%
  pivot_wider(names_from = parameter, values_from = recovered) %>%
  select(participant, `A[c]`,`A[f]`,`bar(F)`, `t[er]`) %>%
  ggpairs(columns = 2:5,
          labeller = "label_parsed",
          lower = list(continuous = wrap("smooth_lm", se = FALSE, alpha = .5)),
          upper = list(continuous = wrap("cor", display_grid = TRUE, digits = 2, stars = TRUE)),
          axisLabels = "show",
          switch = "both",
          progress = FALSE) +
  theme_paper +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = rel(1)),
        strip.text.y.left = element_text(angle = 0),
        axis.text = element_blank())

p_par_cor

## correlation between original and recovered: 
s4_par_comp_cor <- s3_par_comp %>%
  mutate(parameter = factor(parameter, 
                            levels = c("A_correct","A_error","F", "t_er"),
                            labels  = c(expression(A[c]), expression(A[f]), expression(bar(F)), expression(t[er])))) %>%
  group_by(parameter) %>%
  #summarise(r = abs(original-recovered))
  summarise(r = abs(cor(original, recovered, method = "spearman")))





