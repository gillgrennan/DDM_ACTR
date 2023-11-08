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
library(scales)
library(ggh4x)

future::plan("multicore", workers = 6) # Set to desired number of cores

theme_paper <- theme_classic(base_size = 20) + 
  theme(axis.text = element_text(colour = "black"))

use_cached_results <- FALSE # Set to FALSE to rerun simulations (takes time!)

####### dont run above lines if you arent running ddm in r (recommend to run in fas-dm)

set.seed(2021)

### load in reaction time files in a obs format separated by participant: 
### ex: sim_actr has every single trial including: participant, rt, response

## dat1 will be your table with all parameters separated by participant -- can just load in log 
## only run this line if need to calculate densities for cdf plot 

ddm_dat <- crossing(rt = seq(0, 20, by = .01),
                    response = c("upper","lower"))

dat1 <- actr_recovery

dat1$participant <- seq(1:28)

### import trial type .log file from fasdm (i.e. placebo_neutral_invalid)
sim_ddm <- dat1 %>%
  split(.$dataset) %>% ### idk if we are splitting by actual participant here?
  future_map_dfr(function (x) {
    bind_cols(ddm_dat, density =
                ddiffusion(rt = ddm_dat$rt,
                           response = ddm_dat$response, 
                           a = x$a,
                           v = x$v,
                           t0 = x$t0,
                           z = x$szr,
                           st0 = x$st0, 
                           #st0 = 0,
                           #sv = 0,
                           sv = x$sv,
                           #sz = 0
                           sz = x$szr
                )) %>%
      mutate(participant = x$participant[1])
  }) %>%
  mutate(rt = ifelse(response == "upper", rt, -rt),
         model = "DDM")


### import in reaction times: 

n_trials = nrow(sub.001)
label = rep(1,each = n_trials)
sub.001$participant <- c(t(label))

n_trials = nrow(sub.002)
label = rep(2,each = n_trials)
sub.002$participant <- c(t(label))

n_trials = nrow(sub.003)
label = rep(3,each = n_trials)
sub.003$participant <- c(t(label))

n_trials = nrow(sub.004)
label = rep(4,each = n_trials)
sub.004$participant <- c(t(label))

n_trials = nrow(sub.005)
label = rep(5,each = n_trials)
sub.005$participant <- c(t(label))

n_trials = nrow(sub.006)
label = rep(6,each = n_trials)
sub.006$participant <- c(t(label))

n_trials = nrow(sub.007)
label = rep(7,each = n_trials)
sub.007$participant <- c(t(label))

n_trials = nrow(sub.008)
label = rep(8,each = n_trials)
sub.008$participant <- c(t(label))

n_trials = nrow(sub.009)
label = rep(9,each = n_trials)
sub.009$participant <- c(t(label))

n_trials = nrow(sub.010)
label = rep(10,each = n_trials)
sub.010$participant <- c(t(label))

n_trials = nrow(sub.011)
label = rep(11,each = n_trials)
sub.011$participant <- c(t(label))

n_trials = nrow(sub.012)
label = rep(12,each = n_trials)
sub.012$participant <- c(t(label))

n_trials = nrow(sub.013)
label = rep(13,each = n_trials)
sub.013$participant <- c(t(label))

n_trials = nrow(sub.014)
label = rep(14,each = n_trials)
sub.014$participant <- c(t(label))

n_trials = nrow(sub.015)
label = rep(15,each = n_trials)
sub.015$participant <- c(t(label))

n_trials = nrow(sub.016)
label = rep(16,each = n_trials)
sub.016$participant <- c(t(label))


n_trials = nrow(sub.017)
label = rep(17,each = n_trials)
sub.017$participant <- c(t(label))


n_trials = nrow(sub.018)
label = rep(18,each = n_trials)
sub.018$participant <- c(t(label))

n_trials = nrow(sub.019)
label = rep(19,each = n_trials)
sub.019$participant <- c(t(label))

n_trials = nrow(sub.020)
label = rep(20,each = n_trials)
sub.020$participant <- c(t(label))

n_trials = nrow(sub.021)
label = rep(21,each = n_trials)
sub.021$participant <- c(t(label))

n_trials = nrow(sub.022)
label = rep(22,each = n_trials)
sub.022$participant <- c(t(label))

n_trials = nrow(sub.023)
label = rep(23,each = n_trials)
sub.023$participant <- c(t(label))

n_trials = nrow(sub.024)
label = rep(24,each = n_trials)
sub.024$participant <- c(t(label))

n_trials = nrow(sub.025)
label = rep(25,each = n_trials)
sub.025$participant <- c(t(label))

n_trials = nrow(sub.026)
label = rep(26,each = n_trials)
sub.026$participant <- c(t(label))

n_trials = nrow(sub.027)
label = rep(27,each = n_trials)
sub.027$participant <- c(t(label))

n_trials = nrow(sub.028)
label = rep(28,each = n_trials)
sub.028$participant <- c(t(label))


sim_actr <- tibble(rt = c(sub.001$V1,sub.002$V1,sub.003$V1,sub.004$V1,sub.005$V1,sub.006$V1,sub.007$V1,sub.008$V1,sub.009$V1,sub.010$V1,sub.011$V1,sub.012$V1,sub.013$V1,sub.014$V1,sub.015$V1,sub.016$V1,sub.017$V1,sub.018$V1,sub.019$V1,sub.020$V1,sub.021$V1,sub.022$V1,sub.023$V1,sub.024$V1,sub.025$V1,sub.026$V1,sub.027$V1,sub.028$V1),
                   response = c(sub.001$V2,sub.002$V2,sub.003$V2,sub.004$V2,sub.005$V2,sub.006$V2,sub.007$V2,sub.008$V2,sub.009$V2,sub.010$V2,sub.011$V2,sub.012$V2,sub.013$V2,sub.014$V2,sub.015$V2,sub.016$V2,sub.017$V2,sub.018$V2,sub.019$V2,sub.020$V2,sub.021$V2,sub.022$V2,sub.023$V2,sub.024$V2,sub.025$V2,sub.026$V2,sub.027$V2,sub.028$V2),
                   participant = c(sub.001$participant,sub.002$participant,sub.003$participant,sub.004$participant,sub.005$participant,sub.006$participant,sub.007$participant,sub.008$participant,sub.009$participant,sub.010$participant,sub.011$participant,sub.012$participant,sub.013$participant,sub.014$participant,sub.015$participant,sub.016$participant,sub.017$participant,sub.018$participant,sub.019$participant,sub.020$participant,sub.021$participant,sub.022$participant,sub.023$participant,sub.024$participant,sub.025$participant,sub.026$participant,sub.027$participant,sub.028$participant)
) 


### graphing all cdf plots 
sim_actr %>%
  mutate(rt = ifelse(response == 1, rt, -rt),
         model = "Data") %>%
  ggplot(aes(x = rt, colour = model)) +
  facet_wrap(~ participant, ncol = 6) +
  geom_vline(xintercept = 0, lty = 2, colour = "grey80") +
  geom_histogram(aes(y = ..density..), binwidth = .1, colour = "black", fill = "white", size = .1) +
  geom_line(data = sim_ddm, aes(y = density)) +
  scale_x_continuous(limits = c(-20, 20), breaks = c(-15, 0, 15)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_colour_manual(values = c("#e66101", "#5e3c99")) +
  labs(x = "RT (s)",
       y = "Density",
       colour = NULL) +
  guides(colour = "none") +
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
        lwd = (data$size %||% 1) * .pt,
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
        lwd = (data$size %||% 1) * .pt,
        lty = data$linetype %||% 1
      )
    )
  }
  else {
    grobTree() # Don't draw
  }
}

### figure for mapping just 1 / a few participants at a time 
par1 <- 2 ## select which participant to graph
sim_actr %>%
  filter(participant == par1) %>%
  mutate(rt = ifelse(response == 1, rt, -rt),
         model = "ACT-R") %>%
  ggplot(aes(x = rt, colour = model)) +
  facet_wrap(~ participant, ncol = 1) +
  geom_vline(xintercept = 0, lty = 2, colour = "grey80") +
  geom_histogram(aes(y = ..density..), binwidth = .2, fill = "white", size = .15, key_glyph = draw_key_custom) +
  geom_line(data = filter(sim_ddm, participant == par1), aes(y = density), key_glyph = draw_key_custom) +
  scale_x_continuous(limits = c(-2, 2), breaks = c(-1, 0, 1)) +
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

### recovered act-r parameters 

# parameter recovery equations as of oct. 2023 
# parameter recovery 
results_actr <- dat1
results_actr$f <- dat1$a/2
results_actr$ter <- dat1$t0
results_actr$a_correct <- log(dat1$v)
results_actr$a_incorrect <- log(dat1$v)-(2*dat1$a*dat1$v)


statsf <- c(mean(results_actr$f), sd(results_actr$f), range(results_actr$f))
statster <- c(mean(results_actr$ter), sd(results_actr$ter), range(results_actr$ter))
statsa_correct <- c(mean(results_actr$a_correct), sd(results_actr$a_correct), range(results_actr$a_correct))

## graph ACTR estimated parameters (all 3 stimulus types): 

## Analyse ACT-R parameters per stimulus type
## load in all work-spaces and save ACT-R vars accordingly 
# vars should be results_actr_neutral, results_actr_positive, results_actr_negative
# for (i in 1:3)
#   
#   results_actr_positive   
# 
# param_infer_best <- tibble(
#   list = c(1,2,3),
#   parameter <- tibble(participant = c(results_actr_negative$participant),
#                       f_recov = c(results_actr_negative$f),
#                       t_recov = c(results_actr_negative$ter),
#                       u_recov = c(results_actr_negative$u)),
#   parameter <- tibble(participant = c(results_actr_negative$participant),
#                       f_recov = c(results_actr_negative$f),
#                       t_recov = c(results_actr_negative$ter),
#                       u_recov = c(results_actr_negative$u)),
#   parameter <- tibble(participant = c(results_actr_negative$participant),
#                       f_recov = c(results_actr_negative$f),
#                       t_recov = c(results_actr_negative$ter),
#                       u_recov = c(results_actr_negative$u))
#   
# )
# 
# 
# parameter <- tibble(participant = c(results_actr_positive$participant,results_actr_neutral$participant,results_actr_negative$participant),
#                     f_recov = c(results_actr_positive$f,results_actr_neutral$f,results_actr_negative$f),
#                     t_recov = c(results_actr_positive$ter,results_actr_neutral$ter,results_actr_negative$ter),
#                     u_recov = c(results_actr_positive$u,results_actr_neutral$u,results_actr_negative$u))
# 
# df1 <- tibble(
#   g = c(1, 2, 3),
#   data = list(
#     tibble(x = 1, y = 2),
#     tibble(x = 4:5, y = 6:7),
#     tibble(x = 10)
#   )
# )


# # group median 
# param_infer_summary <- param_infer_best %>%
#   group_by(parameter,list) %>%
#   summarise(median = median(value))
# 
# 
# s4_par_comp_error_param <- s4_par_comp_error %>%
#   group_by(n, parameter) %>%
#   summarise(mae = mean(ae))
# 
# s4_par_comp_error_mean <- s4_par_comp_error %>%
#   group_by(n, parameter) %>%
#   summarise(mae = mean(ae)) %>%
#   group_by(n) %>%
#   summarise(mae = mean(mae))
# 
# 
# ggplot(param_infer_plotdat, aes(x = list_jitter, y = value, group = participant, colour = parameter)) +
#   facet_wrap(~ parameter, ncol = 3, scales = "free_y", labeller = labeller(parameter = label_parsed))+
#   geom_line(alpha = .15) +
#   geom_point(alpha = .25) +
#   geom_line(data = param_infer_summary, aes(x = list, y = median, group = parameter), colour = "black", lty = 2) +
#   geom_point(data = param_infer_summary, aes(x = list, y = median, group = parameter), colour = "black", size = rel(2.5)) +
#   scale_x_continuous(breaks = c(1, 2, 3)) +
#   scale_colour_viridis_d() +
#   labs(x = "Session",
#        y = "Parameter value") +
#   guides(colour = FALSE) +
#   theme_paper +
#   theme(strip.background = element_blank(),
#         strip.text = element_text(size = rel(1)))
# 



#### ACT-R Simulation Compare 

### scatter plot of mean RT vs simulated mean RTs across participants 

# import ACT-R data: 
n_trials = nrow(`001`)
label = rep(1,each = n_trials)
`001`$participant <- c(t(label))

n_trials = nrow(`002`)
label = rep(2,each = n_trials)
`002`$participant <- c(t(label))

n_trials = nrow(`003`)
label = rep(3,each = n_trials)
`003`$participant <- c(t(label))

n_trials = nrow(`004`)
label = rep(4,each = n_trials)
`004`$participant <- c(t(label))

n_trials = nrow(`005`)
label = rep(5,each = n_trials)
`005`$participant <- c(t(label))

n_trials = nrow(`006`)
label = rep(6,each = n_trials)
`006`$participant <- c(t(label))

n_trials = nrow(`007`)
label = rep(7,each = n_trials)
`007`$participant <- c(t(label))

n_trials = nrow(`008`)
label = rep(8,each = n_trials)
`008`$participant <- c(t(label))

n_trials = nrow(`009`)
label = rep(9,each = n_trials)
`009`$participant <- c(t(label))

n_trials = nrow(`010`)
label = rep(10,each = n_trials)
`010`$participant <- c(t(label))

n_trials = nrow(`011`)
label = rep(11,each = n_trials)
`011`$participant <- c(t(label))

n_trials = nrow(`012`)
label = rep(12,each = n_trials)
`012`$participant <- c(t(label))

n_trials = nrow(`013`)
label = rep(13,each = n_trials)
`013`$participant <- c(t(label))

n_trials = nrow(`014`)
label = rep(14,each = n_trials)
`014`$participant <- c(t(label))

n_trials = nrow(`015`)
label = rep(15,each = n_trials)
`015`$participant <- c(t(label))

n_trials = nrow(`016`)
label = rep(16,each = n_trials)
`016`$participant <- c(t(label))

n_trials = nrow(`017`)
label = rep(17,each = n_trials)
`017`$participant <- c(t(label))

n_trials = nrow(`018`)
label = rep(18,each = n_trials)
`018`$participant <- c(t(label))

n_trials = nrow(`019`)
label = rep(19,each = n_trials)
`019`$participant <- c(t(label))

n_trials = nrow(`020`)
label = rep(20,each = n_trials)
`020`$participant <- c(t(label))

n_trials = nrow(`021`)
label = rep(21,each = n_trials)
`021`$participant <- c(t(label))

n_trials = nrow(`022`)
label = rep(22,each = n_trials)
`022`$participant <- c(t(label))

n_trials = nrow(`023`)
label = rep(23,each = n_trials)
`023`$participant <- c(t(label))

n_trials = nrow(`024`)
label = rep(24,each = n_trials)
`024`$participant <- c(t(label))

n_trials = nrow(`025`)
label = rep(25,each = n_trials)
`025`$participant <- c(t(label))

n_trials = nrow(`026`)
label = rep(26,each = n_trials)
`026`$participant <- c(t(label))

n_trials = nrow(`027`)
label = rep(27,each = n_trials)
`027`$participant <- c(t(label))

n_trials = nrow(`028`)
label = rep(28,each = n_trials)
`028`$participant <- c(t(label))

sim_new <- tibble(accuracy = c(`001`$accuracy,`002`$accuracy,`003`$accuracy,`004`$accuracy,`005`$accuracy,`006`$accuracy,`007`$accuracy,`008`$accuracy,`009`$accuracy,`010`$accuracy,`011`$accuracy,`012`$accuracy,`013`$accuracy,`014`$accuracy,`015`$accuracy,`016`$accuracy,`017`$accuracy,`018`$accuracy,`019`$accuracy,`020`$accuracy,`021`$accuracy,`022`$accuracy,`023`$accuracy,`024`$accuracy,`025`$accuracy,`026`$accuracy,`027`$accuracy,`028`$accuracy),
                  reaction_time = c(`001`$response_time,`002`$response_time,`003`$response_time,`004`$response_time,`005`$response_time,`006`$response_time,`007`$response_time,`008`$response_time,`009`$response_time,`010`$response_time,`011`$response_time,`012`$response_time,`013`$response_time,`014`$response_time,`015`$response_time,`016`$response_time,`017`$response_time,`018`$response_time,`019`$response_time,`020`$response_time,`021`$response_time,`022`$response_time,`023`$response_time,`024`$response_time,`025`$response_time,`026`$response_time,`027`$response_time,`028`$response_time),
                  participant = c(`001`$participant,`002`$participant,`003`$participant,`004`$participant,`005`$participant,`006`$participant,`007`$participant,`008`$participant,`009`$participant,`010`$participant,`011`$participant,`012`$participant,`013`$participant,`014`$participant,`015`$participant,`016`$participant,`017`$participant,`018`$participant,`019`$participant,`020`$participant,`021`$participant,`022`$participant,`023`$participant,`024`$participant,`025`$participant,`026`$participant,`027`$participant,`028`$participant)
)

orig_rt <- sim_actr %>%
  group_by(participant) %>%
  summarise(accuracy = mean(response),
            rt = mean(rt))

sim_rt <- sim_new %>%
  group_by(participant) %>%
  summarise(accuracy = mean(accuracy),
            rt = mean(reaction_time))

n_participants <- 28 

s3_par_comp <- tibble(
  participant = rep(1:n_participants, 2),
  parameter = rep(c("Accuracy", "RT"), each = n_participants),
  original = c(orig_rt$accuracy,orig_rt$rt),
  recovered = c(sim_rt$accuracy,sim_rt$rt)
)

s3_par_comp_plot <- s3_par_comp %>%
  mutate(parameter = factor(parameter, 
                            levels = c("Accuracy","RT"),
                            labels  = c(expression(Accuracy), expression(RT)))) 

library(ggsci)
#filter(!parameter %in% c("F[a]", "F[b]")) %>%
p_par_comp <- ggplot(s3_par_comp_plot, aes(x = original, y = recovered)) +
  facet_wrap(~parameter, ncol = 2, scales = "free", labeller = labeller(parameter = label_parsed))+
  geom_abline(lty = 2, colour = "grey80") +
  geom_smooth(method = "lm", formula = y ~ x, alpha = .5, colour = "black", lwd = rel(.8)) +
  geom_point(aes(colour = parameter), size = 3, alpha = .5) +
  labs(x = "Original",
       y = "Recovered",
       colour = NULL,
       size = rel(2)) +
  scale_x_continuous(limits = c(min(s3_par_comp$original), max(s3_par_comp$recovered))) +
  scale_color_d3(palette = "category20") +
  guides(colour = "none") +
  theme_paper +
  theme_paper +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = rel(1.5)))

p_par_comp