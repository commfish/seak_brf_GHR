# black rockfish selectivity at age
# ben.williams@alaska.gov

# load ----
source(here::here("R/helper.r"))

# data ----

# note if the lengths are not rounded there are occasionaly more precise estimates
# these will not allow the model to converge
# combine sport and comm age data

plus_group <- 30

# Female only data for sport and commercial fisheries in CSEO
# combined both fleets as there is not enough data currently to evaluate them separately

read_csv(here::here("data/br_bio.csv"), guess = 50000) %>% 
  rename_all(tolower) %>% 
  dplyr::select(year, Area = g_management_area_code, 
                length = length_millimeters, age, sex = sex_code,
                weight = weight_kilograms) %>% 
  mutate(length = round(length / 10),
         fish = "comm") %>% 
  bind_rows(read_csv(here::here("data/sport_brf_bio_se.csv"), guess = 50000) %>%
              rename_all(tolower) %>%
              dplyr::select(year, Area = area,
                            length, age, sex , weight = wt_kg) %>%
              mutate(length = round(length / 10),
                     fish = "sport",
                     sex = case_when(sex=="F" ~ 2,
                                     sex=="M" ~ 1))) -> brf

# 
brf %>% 
  mutate(age = ifelse(age>plus_group, plus_group, age)) %>% 
  filter(Area %in% c("CSEO"), !is.na(age), !is.na(length), 
         sex %in% 2, fish %in% c("comm", "sport")) -> dat

brf %>%
  filter(Area %in% c("CSEO"), !is.na(age), !is.na(length), 
         sex %in% 1:2, fish %in% c("sport", "comm")) %>% 
  group_by(year, sex) %>% 
  tally() %>% 
  pivot_wider(names_from = sex, values_from = n) %>% 
  dplyr::select(male = `1`, female = `2`) -> ratio_dat


# estimate length/weight relationship
# a * length^b

lw <- unname(lm(log(weight) ~ log(length), data = filter(brf, sex==2))$coef)
lw_int = (lw[1])
lw_slope = (lw[2])

# maturity at age from Kodiak study 
# A50 - 10.47 years may not apply...

mat_int = -7.521637
mat_slope = 0.717806

# von Bertalanffy model----

# TMB model
setwd(here::here("tmb"))

# compile("vonb.cpp")
dyn.load(dynlib("vonb"))

# inputs 
data <- list(age = dat$age,
            length = dat$length)

# starting parameters
params <- list(logLinf = log(55), 
               logkappa = log(0.15), 
               t0 = -0.5, 
               logSigma = 0.001)

# build model
model <- MakeADFun(data = data, 
                   parameters = params, 
                   DLL="vonb")

# optimize the model
fit <- nlminb(model$par, 
              model$fn, 
              model$gr)

rep <- sdreport(model)
(f_summary <- summary(rep, select = "all", p.value = T))

# female vonb output 
vb_fit <- model$report()$fit
(Linf <- model$report()$Linf)
Linf_se = (f_summary[1,2])
(kappa <- model$report()$kapp)
kappa_se = (f_summary[2,2])
(t0 <- model$report()$t0)
t0_se = (f_summary[3,2])
(vb_Sigma <- model$report()$Sigma)

dat %>% 
  mutate(fit = Linf * (1 - exp(-kappa * (age - t0)))) %>% 
  ggplot(aes(age, length)) + 
  geom_point() + 
  geom_line(aes(y = fit)) +
  expand_limits(y = 0)

# sex ratio -----
# compile("ratio.cpp")
dyn.load(dynlib("ratio"))

data <- list(male = ratio_dat$male, 
             female = ratio_dat$female,
             model = 2)
params <- list(p = 0.5, tau = 0.8)

L = c(0, # p
      0) # tau

U = c(1, # p
      Inf) # tau

# build model
model <- MakeADFun(data = data, 
                   parameters = params, 
                   DLL="ratio")

# optimize the model
fit <- nlminb(model$par, 
              model$fn, 
              model$gr,
              lower = L,
              upper = U)


rep <- sdreport(model)
(f_summary <- summary(rep, select = "all", p.value = T))
(ratio <- unname(rep$par.fixed[1]))

data.frame(p = model$report()$p,
          p_est = as.list(rep, what = "Estimate")$`p`,
          p_std = as.list(rep, what = "Std")$`p`) %>% 
  mutate(ul = p_est + (1.96 * p_std),
         ll = p_est - (1.96 * p_std)) %>% 
  ggplot(aes(x = p)) +
  geom_histogram(aes(x = p, y = ..density..),
                 alpha = 0.5) +
  geom_vline(aes(xintercept = p_est)) +
  geom_vline(aes(xintercept = ll), lty = 3) +
  geom_vline(aes(xintercept = ul), lty = 3) +
  ggtitle(label = "Sex ratio")


# potential M
brf %>% 
  group_by(sex) %>% 
  summarise(max_age = max(age, na.rm = T))

brf %>% 
  filter(sex==2) %>% 
  ggplot(aes(age)) + 
  geom_histogram(bins = 100)

# max age-based estimate of M
4.899 * 51 ^-0.916
4.899 * 45 ^-0.916
4.899 * 40 ^-0.916
4.899 * 35 ^-0.916

#von B parameter estimates of M
4.118 * kappa^ 0.73 * Linf^-0.33

# looks like M is likely to range somewhere between 0.13 & 0.22
# that is higher than the 0.123 that has been used for AK previously

# clean data for selectivity ----
# get a weighted proportion at age
# using all data lumped together

dat %>% 
  dplyr::select(age, year) %>% 
  mutate(total_n = n()) %>% 
  group_by(year) %>% 
  mutate(annual_n = n()) %>% 
  group_by(age, year) %>% 
  mutate(n = n()) %>% 
  ungroup %>% 
  mutate(prop = range01(n / annual_n) * annual_n) %>% 
  group_by(age) %>% 
  summarise(prop = sum(prop) / mean(total_n)) %>% 
  mutate(prop = range01(prop)) %>% 
  left_join(data.frame(age = 0:plus_group), .) %>% 
  mutate(prop = replace_na(prop, 0)) -> select_dat


# Selectivity model ----
setwd(here::here("tmb"))

compile("select.cpp")
dyn.load(dynlib("select"))

data = list(ages = select_dat$age,
            paaC = select_dat$prop,
            mu_M = 0.15,
            sd_M = 0.03,
            mu_F = 0.1,
            sd_F = 0.05)

params = list(logM = log(0.15),	
              logF = log(0.10),
              logmu = log(8),		    
              logupsilon = log(1.5),
              logsigR = 0.02)

map = list()


# parameter bounds

L <- c(logM = log(.04), 
       logF = log(.02), 
       logmu = log(2), 
       logupsilon = log(.07), 
       logsigR = log(0.0001))

U <- c(logM = log(.4), 
       logF = log(.4),
       logmu = log(15), 
       logupsilon = log(5), 
       logsigR = log(10))

# build model
model <- MakeADFun(data = data, 
                   parameters = params, 
                   DLL = "select",
                   lower = L,
                   upper = U, 
                   map = map)


# optimize the model
fit <- nlminb(model$par, 
              model$fn, 
              model$gr)
rep <- sdreport(model)
(ss = summary(rep, select = "all", p.value = T))


propC <- model$report()$propC
(M <- model$report()$M)
(F <- model$report()$F)
Fa <- model$report()$Fa
saC <- model$report()$saC
(mu <- model$report()$mu)
(upsilon <- model$report()$upsilon)
Ca <- model$report()$Ca 
Ua <- model$report()$Va
Na <- model$report()$Na

data.frame(saf = saC,
           M = M,
           F = F) %>%
  write_csv(here::here("data/select.csv"))

tibble(age = 0:plus_group) %>% 
  mutate(length = Linf * (1 - exp(-kappa * (age - t0))),
         weight = exp(lw_int + lw_slope * log(length)),
         mature = 1 / (1 + exp(mat_int + mat_slope * age)) * exp(mat_int + mat_slope * age),
         unfished = Ua * weight * mature,
         fished = Na * weight * mature,
         Ca = Ca * weight,
         expN = Na * saC,
         tot_bio = Na * weight,
         exp_bio = expN * weight,
         Fa = Fa,
         Sa = saC) -> report

report %>% 
  dplyr::select(age, Sa, mature) %>% 
  pivot_longer(-age, names_to = "measure", "value") %>% 
  ggplot(aes(age, value, color = measure)) + 
  geom_line() +
  scale_color_viridis_d(name = "", end = 0.75) +
  theme(legend.position = c(0.8, 0.2))

report %>% 
  ggplot(aes(age, propC)) + 
  geom_point() +
  geom_bar(aes(y = select_dat$prop), stat = "identity", alpha = 0.3) + 
  geom_line(aes(y = saC), lty = 3) 

report %>% 
  ggplot(aes(age, unfished)) + 
  geom_line() +
  geom_line(aes(y = fished), col = 4)

report %>% 
  summarise(spr = sum(fished) / sum(unfished))


# exploring shape of input and output

data.frame(M = seq(0.05, 0.3, by = 0.01)) %>% 
  ggplot(aes(M)) + 
  stat_function(fun=dnorm, args = list(mean=data$mu_M, sd=data$sd_M), 
                geom="density", fill = "#1B9E77", alpha = 0.3) + 
  stat_function(fun=dnorm, args = list(mean=exp(-1.4), sd=.1), 
                geom="density", fill = "#7570B3", alpha = 0.3) 

data.frame(F = seq(0.0, 0.3, by = 0.01)) %>% 
  ggplot(aes(F)) + 
  stat_function(fun=dnorm, args = list(mean=data$mu_F, sd=data$sd_F), 
                geom="density", fill = "#1B9E77", alpha = 0.3) + 
  stat_function(fun=dnorm, args = list(mean=exp(ss[2,1]), sd=ss[2,2]), 
                geom="density", fill = "#7570B3", alpha = 0.3) 

# additional model exploration
library(adnuts)
test <- sample_tmb(obj=model, init=NULL)
post <- extract_samples(test)
sp <- extract_sampler_params(test)

launch_shinytmb(test)

library(tmbstan)
fit_mcmc <- tmbstan(model, chains = 4, iter = 4000)
pairs(fit_mcmc, pars=names(model$par))
traceplot(fit_mcmc, pars=names(model$par), inc_warmup=FALSE)
