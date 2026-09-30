pacman::p_load(EpiEstim, 
               slider,
               rio, 
               tidyverse)

linelist <- rio::import(here::here("data", "Ebola_North_kivu_10th_outbreak_2018", 
                       "master_linelist_current - deidentified fo 2026.xlsx"))

ditch_weeks <- map(seq(1,isoweek(min(linelist$onset_calc, na.rm = TRUE))-1,1), 
                   ~paste0(2018, "_W",  str_pad(.x,
                                                side = "left", pad = "0", width = 2))
) %>% 
  append(map(seq(isoweek(max(linelist$onset_calc, na.rm = TRUE))+1, 53), 
             ~paste0(2020, "_W", str_pad(.x, 
                                         side = "left", 
                                         pad = "0", 
                                         width = 2)))
  )

inc_pre <- linelist %>% 
  mutate(date_onset = case_when(is.na(onset_calc) & community_death == 1 & 
                                  !is.na(date_death) ~ date_death - ddays(8), 
                                .default = onset_calc), 
         year = isoyear(date_onset))

years <- unique(inc_pre$year)
years <- years[!is.na(years)]

years_weeks <- map(years, 
         function(x){
           paste0(x, "_W",  str_pad(seq(1,isoweek(paste0(x, "-12-28")), 1),
                                    side = "left", pad = "0", width = 2))
         }
         ) %>% 
  unlist() %>% 
  setdiff(ditch_weeks) %>% 
  as_tibble()

inc <- inc_pre %>% 
  mutate(
         date_onset = ymd(onset_calc), 
         iso_week = paste0(isoyear(onset_calc), "_W", str_pad(isoweek(onset_calc), side = "left", pad ="0", width =  2)), 
         isolation_binary = case_when(community_death == 1|delay_admission > 3 ~ "no", 
                                      delay_admission <= 3 ~ "yes")
         ) %>% 
  group_by(iso_week) %>% 
  summarise(count = n(), 
            n_community_death = sum(community_death, na.rm = TRUE)/sum(!is.na(community_death)), 
            n_isolated_early = sum(isolation_binary == "yes", na.rm = TRUE)/sum(!is.na(isolation_binary)), 
            n_contact_registered = sum(contact_registered == "yes", na.rm = TRUE)/
              sum(contact_registered %in% c("no", "yes"), na.rm = TRUE)) %>% 
  ungroup() %>% 
  right_join(years_weeks, by = c("iso_week" = "value")) %>% 
  mutate(count = case_when(is.na(count) ~ 0, 
                           .default = count), 
         iso_week = factor(iso_week, levels = sort(.$iso_week))) %>% 
  arrange(iso_week) %>% 
  mutate(n = row_number(), 
         community_death_slider = slide_dbl(n_community_death, mean, .before = 2), 
         isolation_slider = slide_dbl(n_isolated_early, mean, .before = 2), 
         contact_slider = slide_dbl(n_contact_registered, mean, .before = 2))

cases_d_i <- inc$count

window_width <- 2
t_start <- seq(2, 104 - window_width + 1, by = window_width)
t_end   <- t_start + window_width - 1

r_est_i1 <- wallinga_teunis(
  incid = cases_d_i,
  method = "parametric_si", 
  config = list(mean_si = 15.3/7,
                std_si = 7.3/7,
                method = "parametric_si",
                n_sim = 100,
                t_start = t_start,
                t_end = t_end, 
                seed = 3684
  )
)

r_est_i <- as.data.frame(r_est_i1[["R"]])[, c("t_start","Mean(R)","Std(R)", 
                                             "Quantile.0.025(R)", "Quantile.0.975(R)")] %>% 
  left_join(inc %>% 
              select(n, iso_week), 
            by = c("t_start" = "n")) %>% 
  mutate(shrink = case_when(`Mean(R)` < 1 ~ "1", 
                            .default = "2"))

r <- ggplot(r_est_i)+
  geom_ribbon(aes(x = t_start, 
                  ymin = `Quantile.0.025(R)`, 
                  ymax = `Quantile.0.975(R)`
                  ), 
              fill = "lightblue")+
  geom_line(aes(x = t_start, y = `Mean(R)`, group = 1))+
  geom_hline(yintercept = 1, linetype = "dashed")+
  coord_cartesian(xlim = c(0, 104))+
  theme_minimal()

underwater <- r_est_i %>% 
  filter(shrink == "1") %>% 
  mutate(rect_start = t_start -1,
         rect_end = t_start + 3, 
         lag = lag(t_start)+2, 
         lead = lead(t_start), 
         left = case_when(rect_start > lag | is.na(lag) ~ rect_start), 
         right = case_when(rect_end < lead ~ rect_end)
      #   overlap = case_when(lead == rect_end ~ TRUE), 
      #   lag_overlap = lag(overlap), 
      #   lead_overlap = lead(overlap)
      ) %>% 
  select(rect_start, t_start, rect_end, lag, lead, left, right) %>% 
  fill(left, .direction = "down") %>%
  fill(right, .direction = "up") %>% 
  mutate(right = case_when(is.na(right) ~ left + 3, 
                           .default = right)) %>% 
  distinct(left, right)

death <- ggplot()+
  geom_rect(data = underwater, aes(xmin = left, xmax = right, ymin = 0, ymax = 1), 
            fill = "grey", alpha = 0.3)+
  geom_line(data = inc, aes(x = n, y = n_community_death), color = "black")+
  coord_cartesian(xlim = c(0, 104), ylim = c(0,1))+
  theme_minimal()

isolation <- ggplot()+
  geom_rect(data = underwater, aes(xmin = left, xmax = right, ymin = 0, ymax = 1), 
            fill = "grey", alpha = 0.3)+
  geom_line(data = inc, 
            aes(x = n, y = n_isolated_early), color = "cadetblue")+
  coord_cartesian(xlim = c(0, 104), ylim = c(0,1))+
  theme_minimal()

known_c <- ggplot(inc)+
  geom_rect(data = underwater, aes(xmin = left, xmax = right, ymin = 0, ymax = 1), 
            fill = "grey", alpha = 0.3)+
  geom_line(data = inc, 
            aes(x = n, y = n_contact_registered), color = "blueviolet")+
  coord_cartesian(xlim = c(0, 104), ylim = c(0,1))+
  theme_minimal()

cowplot::plot_grid(r, death, isolation, known_c, 
                   ncol = 1, align = "hv", axis = "tblr")

