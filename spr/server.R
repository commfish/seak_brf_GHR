#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)

# data ----
read.csv(here::here("data/select.csv")) -> slx

theme_set(theme_light() +
              theme(
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  strip.background = element_rect(fill = NA, colour = NA),
                  strip.text.x = element_text(colour = "black"),
                  strip.text.y = element_text(colour = "black"),
                  panel.border = element_rect(fill = NA),
                  legend.key.size = unit(0.9, "lines"),
                  legend.key = element_rect(colour = NA, fill = NA),
                  legend.background = element_rect(colour = NA, fill = NA)
              ))

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
    mutate(age = ifelse(age>30, 30, age)) %>% 
    filter(Area %in% c("CSEO"), !is.na(age), !is.na(length), 
           sex %in% 2, fish %in% c("comm", "sport")) -> dat

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
    left_join(data.frame(age = 0:max(dat$age)), .) %>% 
    mutate(prop = replace_na(prop, 0)) -> select_dat

# functions ----
# scale data to 1
range01 <- function(x){
    (x - min(x)) / (max(x) - min(x))
}


# shiny ----
# Define server logic required to draw a histogram
shinyServer(function(input, output) {
    
    
    # fig 1 ----
    output$distPlot <- renderPlot({
        
        tibble(area = 1000:2500, 
               density = input$density, 
               M = input$M, 
               avg_wt = input$wt, 
               ratio = input$ratio,
               wt_se = input$wt_se, 
               density_se = input$density_se) %>% 
            mutate(wt_low = avg_wt * 0.001 - 2 * wt_se * 0.001,
                   wt_high = avg_wt * 0.001 + 2 * wt_se * 0.001,
                   density_low = density - 2 * density_se,
                   density_high = density + 2 * density_se,
                   biomass = area * density * M * avg_wt * ratio * 0.001, 
                   cols = ifelse(area == 1781,  "1", "0"),
                   sized = ifelse(area == 1781,  3, 1),
                   ll = area * density_low * M * wt_low * ratio,
                   ul = area * density_high * M * wt_high * ratio,
                   bio = 1781 * density * M * avg_wt * ratio * 0.001,
                   bll = 1781 * density * M * wt_low * ratio,
                   bul = 1781 * density * M * wt_high * ratio) -> out
        
        out %>% 
            ggplot(aes(area, biomass, group = 1)) + 
                geom_point(aes(size = sized, color = cols)) +
                geom_ribbon(aes(ymin = ll, ymax = ul), alpha = 0.3) +
                scale_y_continuous(labels = scales::comma) +
                geom_hline(yintercept = 100000  * 0.0004535924, lty = 3) +
                expand_limits(y = c(0)) +
                scale_color_grey(start = 0.7, end = 0.2) +
                guides(color = FALSE, size = FALSE) +
                geom_vline(xintercept = 1781, lty = 3) +
                ylab("biomass t") +
            ggtitle(paste0("female spawning biomass = ", round(out$bio), " t (",round(out$bll),
                           "-",round(out$bul),")"))
        
    })
    
    
    # fig 2 ----
    output$sprPlot <- renderPlot({
        
        Na = Ua = 1781 * input$density * input$ratio
        M = input$M
        F = input$Fc
        ww = input$wt + input$wt_se
        # maturity at age from Kodiak study - may not apply...
        mat_int = -7.521637
        mat_slope = 0.717806
        
        
            # von b
        vonb <- read_csv(here::here("output/vonb.csv"))
        Linf <- vonb$value[vonb$param=="Linf"]
        kappa <- vonb$value[vonb$param=="kappa"]
        t0 <- vonb$value[vonb$param=="t0"]
        
        # weight 
        lw_int = -10.65245
        lw_slope = 2.915867
        
        Sa = slx$saf

        
        for(i in 2:31){
            if(i<31) Na[i] <- exp(-(M + Sa[i-1] * F)) * Na[i-1]
            if(i==31) Na[i] <- (Na[i-1] * exp(-(M + F * Sa[i-1]))) / (1 - exp(-(M + F * Sa[i])))
            if(i<31) Ua[i] <- exp(-M) * Ua[i-1]
            if(i==31) Ua[i] <- (Ua[i-1] * exp(-M)) / (1 - exp(-M))
        }
        
        if(F==0) Na = Ua
        
        Ca = NA
        for(i in 1:31){
            Ca[i] = Na[i] * (1.0 - exp(-(M + F * Sa[i]))) * F * Sa[i] / (M + F * Sa[i]);
        }
        data.frame(age = 0:30, Na, Ua) %>% 
            mutate(RLS = -1.58 + 2.58 * (1-exp(-0.247 * age)),
                   length = Linf * (1 - exp(-kappa * (age - t0))),
                   weight = exp(lw_int + lw_slope * log(length)),
                   mature = exp(mat_int + mat_slope * age) / (1 + exp(mat_int + mat_slope * age)),
                   weight = ifelse(weight>ww, ww, weight), 
                   unfished = Ua * weight * mature * 0.001 * RLS,
                   fished = Na * weight * mature * 0.001 * RLS,
                   Ca = round(Ca * weight * 0.001),
                   expN = Na * Sa,
                   tot_bio = Na * weight,
                   exp_bio = expN * weight,
                   Fa = F * Sa) -> report
        
        report %>% 
            summarise(spr = round(sum(fished) / sum(unfished), digits = 2)) -> spr
        
        report %>% 
            summarise(wt = sum(Ca)) %>% pull() %>% round() -> wt
        
        report %>% 
            ggplot(aes(age, unfished)) + 
            geom_line() +
            geom_line(aes(y = fished), col = 4) + 
            ggtitle(paste0("spr = ", spr, "; catch = ", wt, " t")) + 
            scale_y_continuous(labels = scales::comma) 
    })
    
    # fig 3 ----
    output$sprPlot2 <- renderPlot({
        
        Na = Ua = 1781 * input$density * input$ratio
        M = input$M
        F = input$Fc
        ww = input$wt + input$wt_se
        # maturity at age from Kodiak study - may not apply...
        mat_int = -7.521637
        mat_slope = 0.717806
        
        # von b
        # von b
        vonb <- read_csv(here::here("output/vonb.csv"))
        Linf <- vonb$value[vonb$param=="Linf"]
        kappa <- vonb$value[vonb$param=="kappa"]
        t0 <- vonb$value[vonb$param=="t0"]
        
        # weight 
        lw_int = -10.65245
        lw_slope = 2.915867
        
        Sa = slx$saf
        
        
        for(i in 2:31){
            if(i<31) Na[i] <- exp(-(M + Sa[i-1] * F)) * Na[i-1]
            if(i==31) Na[i] <- (Na[i-1] * exp(-(M + F * Sa[i-1]))) / (1 - exp(-(M + F * Sa[i])))
            if(i<31) Ua[i] <- exp(-M) * Ua[i-1]
            if(i==31) Ua[i] <- (Ua[i-1] * exp(-M)) / (1 - exp(-M))
        }
        if(F==0) Na = Ua
        
        Ca = NA
        for(i in 1:31){
            Ca[i] = Na[i] * (1.0 - exp(-(M + F * Sa[i]))) * F * Sa[i] / (M + F * Sa[i]);
        }
        data.frame(age = 0:30, Na, Ua) %>% 
            mutate(RLS = -1.58 + 2.58 * (1-exp(-0.247 * age)),
                   length = Linf * (1 - exp(-kappa * (age - t0))),
                   weight = exp(lw_int + lw_slope * log(length)),
                   mature = exp(mat_int + mat_slope * age) / (1 + exp(mat_int + mat_slope * age)),
                   weight = ifelse(weight>ww, ww, weight), 
                   unfished = Ua * weight * mature * 0.001 * RLS,
                   fished = Na * weight * mature * 0.001 * RLS,
                   Ca = round(Ca * weight * 0.001),
                   expN = Na * Sa,
                   tot_bio = Na * weight,
                   exp_bio = expN * weight,
                   Fa = F * Sa) -> report
        
        report %>% 
            summarise(spr = round(sum(fished) / sum(unfished), digits = 2)) -> spr
        
        report %>% 
            summarise(wt = sum(Ca)) %>% pull() %>% round() -> wt
        report %>% 
            ggplot(aes(age, Ca / max(Ca))) + 
            # ggtitle(paste0("spr = ", spr, "; catch = ", wt, " t")) + 
            scale_y_continuous(labels = scales::comma) +
            geom_bar(stat = "identity", alpha = 0.3, fill = "#1B9E77") +
            geom_bar(aes(y = select_dat$prop), stat = "identity", alpha = 0.5, fill = "#7570B3") 
    })
    
    
})

