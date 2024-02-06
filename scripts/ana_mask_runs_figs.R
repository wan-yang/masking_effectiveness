# Model code used in Yang & Shaman "Reconciling the efficacy and effectiveness of masking on epidemic outcomes"
# author: Wan Yang
# date: March 2023
# note: this code is for research use only; may not be entirely optimized nor neatly organized. 

# script to analyze the main runs and plot the figures

# TO USE - FIRST SET THE WORKING DIRECTORY TO SOURCE FILE LOCATION
dir_data = '../data/'
dir_code = '../scripts/'
dir_res = '../results/'
dir_fig = '../figures/'

#######################################################
# load the results
# main analysis
load("../results/res_main_diff.bysn.RData")  # population-level effect
load("../results/res_main_diff.byeff.bysn.RData") # individual level effect
load("../results/res_main_daily.RData") # simulated epidemic curves

# counterfactual - mask at home
load("../results/res_maskathome_diff.bysn.RData")
# sensitivity analysis - assuming different weights for nonphysical contact
load("../results/res_nonphywt_diff.bysn.RData")
# sensitivity analysis - assuming R0 is 2x the value used in the main & counterfactual analyses
load("../results/res_R0x2_diff.bysn.RData")
load("../results/res_R0x2_diff.byeff.bysn.RData") # individual level effect

#######################################################
# format for plotting
sz.t = 10
theme.t = theme(plot.title = element_text(v=0, size = sz.t, margin=margin(0,0,3,0)), 
                strip.placement = "outside", strip.text = element_text(size = sz.t-1, margin=margin(1.5,0,1.5,0), hjust = 0),
                axis.title = element_text(size =sz.t, margin=margin(0,0.2,0,0)), 
                axis.title.x = element_text(margin = margin(t = -5, r = 0, b = 0, l = 0)),
                axis.text.y = element_text(size=sz.t -1, margin=margin(0,0.2,0,0)), 
                axis.text.x = element_text(size=sz.t -1,angle = 35, hjust =1, vjust = 1.3),
                plot.margin=unit(c(c(.3, 1, .1, .5)), units="line"), # top, right, bottom, left
                legend.title = element_text(size=sz.t -1), legend.text=element_text(size=sz.t-1, vjust = 1),
                legend.margin=margin(0,0,0,0),
                legend.box.margin=margin(5,-10,0,-6),
                legend.key.size = unit(.2, 'cm'), #change legend key size
                legend.key.height = unit(.7, 'cm'), #change legend key height
                legend.key.width = unit(.2, 'cm')) #change legend key width)

theme.t2 = theme(plot.title = element_text(v=0, size = sz.t, margin=margin(0,0,3,0)), 
                 strip.placement = "outside", strip.text = element_text(size = sz.t-1, margin=margin(1.5,0,1.5,0), hjust = 0),
                 axis.title = element_text(size =sz.t, margin=margin(0,0.2,0,0)), 
                 axis.title.x = element_text(margin = margin(t = -5, r = 0, b = 0, l = 0)),
                 axis.text.y = element_text(size=sz.t -1, margin=margin(0,0.2,0,0)), 
                 axis.text.x = element_text(size=sz.t -1,angle = 0, hjust =.5, vjust = 1.3),
                 plot.margin=unit(c(.3, .1, .1, .5), units="line"), # top, right, bottom, left
                 legend.title = element_text(size=sz.t -1), legend.text=element_text(size=sz.t-1, vjust = .5),
                 legend.margin=margin(0,0,0,0),
                 legend.box.margin=margin(t = 5, r = 2, b = 0, l = -10),
                 legend.key.size = unit(.2, 'cm'), #change legend key size
                 legend.key.height = unit(.5, 'cm'), #change legend key height
                 legend.key.width = unit(.2, 'cm')) #change legend key width)

theme.t3 = theme(plot.title = element_text(v=0, size = sz.t, margin=margin(0,0,3,0)), 
                 panel.spacing = unit(0, "lines"),
                 strip.placement = "outside", strip.text = element_text(size = sz.t-1, margin=margin(1.5,0,1.5,0), hjust = 0),
                 axis.title = element_text(size =sz.t, margin=margin(0,0.2,0,0)), 
                 axis.title.x = element_text(margin = margin(t = -5, r = 0, b = 0, l = 0)),
                 axis.text.y = element_text(size=sz.t -1, margin=margin(0,0.2,0,0)), 
                 axis.text.x = element_blank(), # element_text(size=sz.t -1,angle = 30, hjust =.5, vjust = 1.3)
                 plot.margin=unit(c(.1, -0.5, 0.1, 0.5), units="line"), # top, right, bottom, left
                 legend.title = element_text(size=sz.t -1), legend.text=element_text(size=sz.t-1, vjust = .5),
                 legend.margin=margin(0,0,0,0),
                 legend.box.margin=margin(t = 2, r = 10, b = 0, l = -6),
                 legend.key.size = unit(.2, 'cm'), #change legend key size
                 legend.key.height = unit(.5, 'cm'), #change legend key height
                 legend.key.width = unit(.2, 'cm')) #change legend key width)
theme.t4b = theme(plot.title = element_text(v=0, size = sz.t, margin=margin(5,0,3,0)), 
                  strip.placement = "outside", strip.text = element_text(size = sz.t-1, margin=margin(1.5,0,1.5,0), hjust = 0),
                  axis.title = element_text(size =sz.t, margin=margin(0,0.2,0,0)), 
                  axis.title.x = element_text(margin = margin(t = -5, r = 0, b = 0, l = 0)),
                  axis.text.y = element_text(size=sz.t -1, margin=margin(0,0.2,0,0)), 
                  axis.text.x = element_text(size=sz.t -1,angle = 30, hjust =1, vjust = 1.3),
                  plot.margin=unit(c(.3, .5, .5, .5), units="line"), # top, right, bottom, left
                  legend.position = 'null') #change legend key width)
#######################################################

#######################################################
# labeling etc. 
eff.cut = 25 # typical upper bound of observed effectiveness

timeframes = res_main_diff.bysn$timeframe %>% unique
names(timeframes) = factor(timeframes, levels = c("cum2date", "this.season"), labels = c('cumulative', 'this season'))
mea.metrics = res_main_diff.bysn$mea.metric %>% unique
seasons.t = unique(res_main_diff.bysn$season) %>% sort
seasons = unique(res_main_diff.bysn$season) %>% sort
seasons.labs = gsub('2019-2020','March-June 2020', seasons)

frac.masking_vec = seq(0, 1, by = .1)  # fm
cut.mask.start_vec = c(0, .01/100, .02/100, .05/100, .1/100, .2/100, .5/100, 1/100, 2/100) # cm
cut.mask.start_vec.lab = c('always','0.01%','0.02%','0.05%','0.1%','0.2%','0.5%', '1%', '2%')

age.grps = read.csv(paste0(dir_code,'age.grps8.csv'))

############################################################
# FIGURE 1
# Fig 1 time series + combined % reduction 
mea.t = 'infections'
diff.ctc.redn.bymask.t = F
diff.ctc.bymask.t = F
n.tot = 1e5
metric.t = 'perc.redn'; mask.grp.t = 'all'; age.grp.t = 'all'

cut.mask.start_vec.t = c(0, .1/100, .5/100, 1/100)
cut.mask.start_vec.lab.t = c('always wear masks when outside the home', 'wear masks when prevalence >0.1%', 
                             'wear masks when prevalence >0.5%', 'wear masks when prevalence >1%')
it = 1
cut.mask.start.t = cut.mask.start_vec.t[it]
tda = res_main_daily %>% 
  filter(diff.ctc.redn.bymask == diff.ctc.redn.bymask.t & diff.ctc.bymask == diff.ctc.bymask.t &
           measure == mea.t & mask.grp %in% mask.grps.t)  # mask.grp == mask.grp.t &
# need to add no mask           
tda = tda %>% filter((cut.mask.start == cut.mask.start.t & ave.frac.masking %in% seq(.2, 1, by=.2)) | ave.frac.masking == 0) %>% 
  melt(., id.vars = c("diff.ctc.redn.bymask", "diff.ctc.bymask","ave.frac.masking","cut.mask.start",
                      "measure", "mask.grp","date")) %>%
  # need to normalize the numbers
  mutate(n.ppl = case_when(grepl('no mask', mask.grp) ~ n.tot * (1 - ave.frac.masking),
                           grepl('use mask', mask.grp) ~ n.tot * ave.frac.masking,
                           grepl('all', mask.grp) ~ n.tot)) %>%
  mutate(perc = value / n.ppl * 100) %>%
  dcast(., diff.ctc.redn.bymask + diff.ctc.bymask + ave.frac.masking + cut.mask.start + measure + mask.grp + date ~ variable, value.var = 'perc') %>%
  mutate(perc.mask = factor(ave.frac.masking * 100))
dates.t = tda %>% .$date %>% unique %>% as.Date %>% sort
dates.t = seq(dates.t[1], tail(dates.t, 1), by = '3 months')

pp1a = ggplot(tda, aes(x = date, y = median, color = perc.mask)) +
  geom_line() +  # data = tda %>% filter(mask.grp == 'no mask (assigned)')
  geom_ribbon(aes(x = date, ymin = ci50lwr, ymax = ci50upr, fill = perc.mask), alpha = .3, color = NA)  +
  # facet_rep_wrap(~cut.mask.start, repeat.tick.labels = T, ncol = 2) +
  scale_x_date(breaks = dates.t, labels = format(dates.t,'%b\n%Y'), limits = as.Date(c('2020/3/1','2024/12/31')),expand=c(.01,.01)) +
  labs(y = paste('%', mea.t, ' per day'), title = paste0('(A) Simulated daily infection rates: ', cut.mask.start_vec.lab.t[it]), 
       x = '', color = '%', fill = '%') + 
  theme_minimal() + theme.t2

it = 4
cut.mask.start.t = cut.mask.start_vec.t[it]
tda = res_main_daily %>% 
  filter(diff.ctc.redn.bymask == diff.ctc.redn.bymask.t & diff.ctc.bymask == diff.ctc.bymask.t &
           measure == mea.t & mask.grp %in% mask.grps.t)  # mask.grp == mask.grp.t &
# need to add no mask           
tda = tda %>% filter((cut.mask.start == cut.mask.start.t & ave.frac.masking %in% seq(.2, 1, by=.2)) | ave.frac.masking == 0) %>% 
  melt(., id.vars = c("diff.ctc.redn.bymask", "diff.ctc.bymask","ave.frac.masking","cut.mask.start",
                      "measure", "mask.grp","date")) %>%
  # need to normalize the numbers
  mutate(n.ppl = case_when(grepl('no mask', mask.grp) ~ n.tot * (1 - ave.frac.masking),
                           grepl('use mask', mask.grp) ~ n.tot * ave.frac.masking,
                           grepl('all', mask.grp) ~ n.tot)) %>%
  mutate(perc = value / n.ppl * 100) %>%
  dcast(., diff.ctc.redn.bymask + diff.ctc.bymask + ave.frac.masking + cut.mask.start + measure + mask.grp + date ~ variable, value.var = 'perc') %>%
  mutate(perc.mask = factor(ave.frac.masking * 100))

pp1b = ggplot(tda, aes(x = date, y = median, color = perc.mask)) +
  geom_line() +  # data = tda %>% filter(mask.grp == 'no mask (assigned)')
  geom_ribbon(aes(x = date, ymin = ci50lwr, ymax = ci50upr, fill = perc.mask), alpha = .3, color = NA)  +
  # facet_rep_wrap(~cut.mask.start, repeat.tick.labels = T, ncol = 2) +
  scale_x_date(breaks = dates.t, labels = format(dates.t,'%b\n%Y'), limits = as.Date(c('2020/3/1','2024/12/31')),expand=c(.01,.01)) +
  labs(y = paste('%', mea.t, ' per day'), 
       title = paste0('(B) Simulated daily infection rates: ', cut.mask.start_vec.lab.t[it]), 
       x = '', color = '%', fill = '%') + 
  theme_minimal() + theme.t2

# add mobility and seasonality
da.p.ctc.redn = read.csv(paste0(dir_data,'da.p.ctc.redn.csv')) %>% data.table()
da.p.ctc.redn$date = da.p.ctc.redn$date %>% as.Date 
dates.t = da.p.ctc.redn %>% .$date %>% unique %>% as.Date %>% sort
dates.t = seq(dates.t[1], tail(dates.t, 1), by = '3 months')

date.start = as.Date('2020/3/1'); num.yr.sim = 5
date.end = as.Date(paste0(2020+num.yr.sim-1,'/12/31')) # as.Date('2119/12/31')
tab_dates = data.table(date = seq(date.start, date.end, by = 'day'))
tab_dates$tt = 1:nrow(tab_dates) # used as a time stamp to save space
tab_dates$year = tab_dates$date %>% year() # would cause problem: %>% MMWRweek() %>% .$MMWRyear # 
# tab_dates$week = tab_dates$date %>% MMWRweek() %>% .$MMWRweek
tab_dates$EndOfYr = grepl('-12-31', tab_dates$date)
tab_dates$MidYr = grepl('-06-30', tab_dates$date) # tally by season July 1 to June 30 the next year
tab_dates$DayOfYr = tab_dates$date %>% lubridate::yday()
tab_dates$DayOfWk = tab_dates$date %>% wday()
# add season
tab_dates = tab_dates %>% mutate(season = case_when(date <= as.Date(paste0(year,'-06-30')) ~ paste0(year-1,'-', year),
                                                    T ~ paste0(year,'-', year+1)))
setDT(tab_dates)
tx.sn.amp = .25
tx.sn = data.table(date = da.p.ctc.redn$date, rel.sn = 1 + tx.sn.amp * cos((tab_dates$DayOfYr %% 365)/365 * 2 * pi)) # seasonal risk

pp1c = ggplot(da.p.ctc.redn, aes(x = date, y = value)) +
  geom_line(color = 'blue') +  # data = tda %>% filter(mask.grp == 'no mask (assigned)')
  geom_line(data = tx.sn, aes(x=date, y = rel.sn), color = 'red') +
  scale_x_date(breaks = dates.t, labels = format(dates.t,'%b\n%Y'), limits = as.Date(c('2020/3/1','2024/12/31')),expand=c(.01,.01)) +
  scale_y_continuous(name = "Mobility relative to pre-COVID-19 level", # Features of the first axis
                     sec.axis = sec_axis(trans=~., name="Relative seasonal infection risk") # Add a second axis and specify its features
  ) +
  labs(title = '(C) Changes in population mobility and seasonal infection risk', 
       x = '') + 
  theme_minimal() + theme.t2 + theme( axis.line.y.left = element_line(color = "blue"), 
                                      axis.ticks.y.left = element_line(color = "blue"),
                                      axis.text.y.left = element_text(color = 'blue'),
                                      axis.title.y.left = element_text(colour = 'blue'),
                                      axis.line.y.right = element_line(color = "red"), 
                                      axis.ticks.y.right = element_line(color = "red"),
                                      axis.text.y.right = element_text(color = 'red'),
                                      axis.title.y.right = element_text(colour = 'red'))



pdf(paste0(dir_fig, 'Fig1.pdf'), width = 10, height = 2.5*3)
grid.arrange(grobs = list(pp1a, pp1b, pp1c),
             nrow = 3, ncol = 1)
dev.off()
############################################################



############################################################
# FIGURE 2
ifig.start = 1
timeframe.t = 'cum2date'; mea.metric.t = 'total'; mea.t = 'infections'
tda = res_main_diff.bysn %>% 
  filter(diff.ctc.bymask == diff.ctc.bymask.t & diff.ctc.redn.bymask == diff.ctc.redn.bymask.t) %>%
  mutate(cut.mask.start1 = factor(cut.mask.start1, levels = cut.mask.start_vec, labels = cut.mask.start_vec.lab),
         ave.frac.masking1 = factor(ave.frac.masking1, levels = frac.masking_vec)
  ) %>%
  mutate(measure = case_when(grepl('AR', variable) ~ 'infections',
                             grepl('Hosp', variable) ~ 'hospitalizations'),
         mea.metric = case_when(grepl('tot', variable) ~ 'total',
                                grepl('1plus', variable) ~ 'individuals affected >=1 time',
                                grepl('2plus', variable) ~ 'individuals affected >=2 times',
                                grepl('3plus', variable) ~ 'individuals affected >=3 times')) %>%
  mutate(season = factor(season, levels = seasons, labels = seasons.labs))%>% 
  filter(mask.grp == mask.grp.t & age.grp == age.grp.t & 
           metric == metric.t & timeframe == timeframe.t & 
           measure == mea.t & mea.metric  %in% mea.metric.t)
layer = tda %>% filter(round(med,0) < eff.cut) # above a certain cut off
# significantly diff
d.sig = tda %>% filter(ci95lwr > 0)
pp2a = ggplot(tda, 
              aes(x = cut.mask.start1, y = ave.frac.masking1)) +  # use the median - the mean can be very skewed 
  geom_tile(aes(fill = med)) + 
  geom_tile(data=layer, alpha = 0.0, color = "black", size = 1, linejoin = "round") +
  geom_tile(data=layer, alpha = 1, aes(fill = med)) +
  geom_text(aes(x = cut.mask.start1, y = ave.frac.masking1, label = med %>% round(., 0)), size = 3) +
  geom_text(data=d.sig, aes(x = cut.mask.start1, y = ave.frac.masking1, label = med %>% round(., 0)), size = 3, fontface = "bold") +
  facet_rep_wrap(~ season, repeat.tick.labels = T, ncol = length(seasons.t), labeller = label_wrap_gen(multi_line=FALSE)) +
  scale_fill_distiller(palette = "RdPu", direction = 1, na.value = 'lightgrey', limits = c(0, 100)) +
  labs(x = 'time to start masking', y = 'fraction pop masking', fill = '% redn', 
       title = paste0('(',LETTERS[ifig.start],') Percent reduction in total ', mea.t, ': ', ifelse(timeframe.t=='cum2date','tallies from March 2020 through this season','tallies for this season'))) +
  theme_minimal() + theme.t


timeframe.t = 'cum2date'; mea.metric.t = 'individuals affected >=1 time'; 
tda = res_main_diff.bysn %>% 
  filter(diff.ctc.bymask == diff.ctc.bymask.t & diff.ctc.redn.bymask == diff.ctc.redn.bymask.t) %>%
  mutate(cut.mask.start1 = factor(cut.mask.start1, levels = cut.mask.start_vec, labels = cut.mask.start_vec.lab),
         ave.frac.masking1 = factor(ave.frac.masking1, levels = frac.masking_vec)
  ) %>%
  mutate(measure = case_when(grepl('AR', variable) ~ 'infections',
                             grepl('Hosp', variable) ~ 'hospitalizations'),
         mea.metric = case_when(grepl('tot', variable) ~ 'total',
                                grepl('1plus', variable) ~ 'individuals affected >=1 time',
                                grepl('2plus', variable) ~ 'individuals affected >=2 times',
                                grepl('3plus', variable) ~ 'individuals affected >=3 times')) %>%
  mutate(season = factor(season, levels = seasons, labels = seasons.labs))%>% 
  filter(mask.grp == mask.grp.t & age.grp == age.grp.t & 
           metric == metric.t & timeframe == timeframe.t & 
           measure == mea.t & mea.metric  %in% mea.metric.t)
layer = tda %>% filter(round(med,0) < eff.cut) # above a certain cut off
# significantly diff
d.sig = tda %>% filter(ci95lwr > 0)
pp2b = ggplot(tda, 
              aes(x = cut.mask.start1, y = ave.frac.masking1)) +  # use the median - the mean can be very skewed 
  geom_tile(aes(fill = med)) + 
  geom_tile(data=layer, alpha = 0.0, color = "black", size = 1, linejoin = "round") +
  geom_tile(data=layer, alpha = 1, aes(fill = med)) +
  geom_text(aes(x = cut.mask.start1, y = ave.frac.masking1, label = med %>% round(., 0)), size = 3) +
  geom_text(data=d.sig, aes(x = cut.mask.start1, y = ave.frac.masking1, label = med %>% round(., 0)), size = 3, fontface = "bold") +
  facet_rep_wrap(~ season, repeat.tick.labels = T, ncol = length(seasons.t), labeller = label_wrap_gen(multi_line=FALSE)) +
  scale_fill_distiller(palette = "RdPu", direction = 1, na.value = 'lightgrey', limits = c(0, 100)) +
  labs(x = 'time to start masking', y = 'fraction pop masking', fill = '% redn', 
       title = paste0('(',LETTERS[ifig.start+1],') Percent reduction in individuals ever-infected: ', 
                      # gsub('affected', ifelse(mea.t=='infections','infected','hospitalized'),mea.metric.t),': ', 
                      ifelse(timeframe.t=='cum2date','tallies from March 2020 through this season','tallies for this season'))) +
  theme_minimal() + theme.t


timeframe.t = 'this.season'; mea.metric.t = 'total'; 
tda = res_main_diff.bysn %>% 
  filter(diff.ctc.bymask == diff.ctc.bymask.t & diff.ctc.redn.bymask == diff.ctc.redn.bymask.t) %>%
  mutate(cut.mask.start1 = factor(cut.mask.start1, levels = cut.mask.start_vec, labels = cut.mask.start_vec.lab),
         ave.frac.masking1 = factor(ave.frac.masking1, levels = frac.masking_vec)
  ) %>%
  mutate(measure = case_when(grepl('AR', variable) ~ 'infections',
                             grepl('Hosp', variable) ~ 'hospitalizations'),
         mea.metric = case_when(grepl('tot', variable) ~ 'total',
                                grepl('1plus', variable) ~ 'individuals affected >=1 time',
                                grepl('2plus', variable) ~ 'individuals affected >=2 times',
                                grepl('3plus', variable) ~ 'individuals affected >=3 times')) %>%
  mutate(season = factor(season, levels = seasons, labels = seasons.labs))%>% 
  filter(mask.grp == mask.grp.t & age.grp == age.grp.t & 
           metric == metric.t & timeframe == timeframe.t & 
           measure == mea.t & mea.metric  %in% mea.metric.t)
layer = tda %>% filter(round(med,0) < eff.cut) # above a certain cut off
# significantly diff
d.sig = tda %>% filter(ci95lwr > 0)
pp2c = ggplot(tda, 
              aes(x = cut.mask.start1, y = ave.frac.masking1)) +  # use the median - the mean can be very skewed 
  geom_tile(aes(fill = med)) + 
  geom_tile(data=layer, alpha = 0.0, color = "black", size = 1, linejoin = "round") +
  geom_tile(data=layer, alpha = 1, aes(fill = med)) +
  geom_text(aes(x = cut.mask.start1, y = ave.frac.masking1, label = med %>% round(., 0)), size = 3) +
  geom_text(data=d.sig, aes(x = cut.mask.start1, y = ave.frac.masking1, label = med %>% round(., 0)), size = 3, fontface = "bold") +
  facet_rep_wrap(~ season, repeat.tick.labels = T, ncol = length(seasons.t), labeller = label_wrap_gen(multi_line=FALSE)) +
  scale_fill_distiller(palette = "RdPu", direction = 1, na.value = 'lightgrey', limits = c(0, 100)) +
  labs(x = 'time to start masking', y = 'fraction pop masking', fill = '% redn', 
       title = paste0('(',LETTERS[ifig.start+2],') Percent reduction in total ', mea.t, ': ', ifelse(timeframe.t=='cum2date','tallies from March 2020 through this season','tallies for this season'))) +
  theme_minimal() + theme.t

pdf(paste0(dir_fig, 'Fig2_redn_main.pdf'), width = 10, height = 2.5*3)
grid.arrange(grobs = list(pp2a, pp2b, pp2c),
             nrow = 3, ncol = 1)
dev.off()
############################################################

############################################################
# FIGURE 3
# Fig 3 - plot example distribution of contact
source(paste0(dir_code,'R0_by_ctc.R'))
ctc = rec %>% dplyr::select(idx, age, num.HHctc, num.nonHHctc) %>% 
  {eval(parse(text = paste('dplyr::mutate(., age.grp = case_when(', str_age_recode, ', T ~ NA_character_))', sep='')))} %>%
  mutate(age.grp = factor(age.grp, levels = age.grps$label, labels = paste('Age',age.grps$label))) %>% mutate(age = NULL)
ctc2 = ctc %>% melt(., id.vars = c('idx','age.grp')) %>%
  mutate(ctc_type = factor(variable, levels = c('num.HHctc','num.nonHHctc'), labels = c('HH','nonHH'))) %>% mutate(variable = NULL)


theme.t3 = theme(plot.title = element_text(v=0, size = sz.t, margin=margin(0,0,3,0)), 
                 strip.placement = "outside", strip.text = element_text(size = sz.t-1, margin=margin(1.5,0,1.5,0), hjust = 0),
                 axis.title = element_text(size =sz.t, margin=margin(0,0.2,0,0)), 
                 axis.title.x = element_text(margin = margin(t = 2, r = 0, b = 0, l = 0)),
                 axis.text.y = element_text(size=sz.t -1, margin=margin(0,0.2,0,0)), 
                 axis.text.x = element_text(size=sz.t -1,angle = 0, margin = margin(t = -1, r = 0, b = 0, l = 0)),
                 plot.margin=unit(c(c(.3, 1, 1, .5)), units="line"), # top, right, bottom, left
                 legend.title = element_text(size=sz.t -1), legend.text=element_text(size=sz.t-1, vjust = 1),
                 legend.margin=margin(0,0,0,0),
                 legend.box.margin=margin(5,-10,0,-6),
                 legend.key.size = unit(.2, 'cm'), #change legend key size
                 legend.key.height = unit(.3, 'cm'), #change legend key height
                 legend.key.width = unit(.2, 'cm')) #change legend key width)

pp3a = ggplot(ctc2, aes(x = value, fill = ctc_type, color = ctc_type)) +
  geom_histogram(aes(y = ..density..), position = 'identity', alpha = .3, binwidth = .2, size = .3) +
  labs(title = '(A) Example contact distributions: Household (HH) vs. non-household (nonHH) contacts', 
       x = 'weighted cumulative contact duration (days per day)') + 
  facet_rep_wrap(~age.grp, repeat.tick.labels = T, ncol = nrow(age.grps),
                 labeller = label_wrap_gen(multi_line=FALSE))+
  theme_minimal() + theme.t3



timeframe.t = 'cum2date'; mea.metric.t = 'total'; mea.t = 'infections'
tda = res_maskathome_diff.bysn %>% 
  filter(diff.ctc.bymask == diff.ctc.bymask.t & diff.ctc.redn.bymask == diff.ctc.redn.bymask.t) %>%
  mutate(cut.mask.start1 = factor(cut.mask.start1, levels = cut.mask.start_vec, labels = cut.mask.start_vec.lab),
         ave.frac.masking1 = factor(ave.frac.masking1, levels = frac.masking_vec)
  ) %>%
  mutate(measure = case_when(grepl('AR', variable) ~ 'infections',
                             grepl('Hosp', variable) ~ 'hospitalizations'),
         mea.metric = case_when(grepl('tot', variable) ~ 'total',
                                grepl('1plus', variable) ~ 'individuals affected >=1 time',
                                grepl('2plus', variable) ~ 'individuals affected >=2 times',
                                grepl('3plus', variable) ~ 'individuals affected >=3 times')) %>%
  mutate(season = factor(season, levels = seasons, labels = seasons.labs))%>% 
  filter(mask.grp == mask.grp.t & age.grp == age.grp.t & 
           metric == metric.t & timeframe == timeframe.t & 
           measure == mea.t & mea.metric  %in% mea.metric.t)
layer = tda %>% filter(round(med,0) < eff.cut) # above a certain cut off

# ctc2 %>% group_by(age.grp, ctc_type) %>%
#  summarise(mn = mean(value), med = median(value))

# significantly diff
d.sig = tda %>% filter(ci95lwr > 0)
pp3b = ggplot(tda, 
              aes(x = cut.mask.start1, y = ave.frac.masking1)) +  # use the median - the mean can be very skewed 
  geom_tile(aes(fill = med)) + 
  geom_tile(data=layer, alpha = 0.0, color = "black", size = 1, linejoin = "round") +
  geom_tile(data=layer, alpha = 1, aes(fill = med)) +
  geom_text(aes(x = cut.mask.start1, y = ave.frac.masking1, label = med %>% round(., 0)), size = 3) +
  geom_text(data=d.sig, aes(x = cut.mask.start1, y = ave.frac.masking1, label = med %>% round(., 0)), size = 3, fontface = "bold") +
  facet_rep_wrap(~ season, repeat.tick.labels = T, ncol = length(seasons.t), labeller = label_wrap_gen(multi_line=FALSE)) +
  scale_fill_distiller(palette = "RdPu", direction = 1, na.value = 'lightgrey', limits = c(0, 100)) +
  labs(x = 'time to start masking', y = 'fraction pop masking', fill = '% redn', 
       title = paste0('(B) Assume people use masks at home: percent reduction in total infections, ', ifelse(timeframe.t=='cum2date','tallies from March 2020 through this season','tallies for this season'))) +
  theme_minimal() + theme.t

pdf(paste0(dir_fig, 'Fig3_ctc_maskathome.pdf'), width = 10, height = 5)
grid.arrange(grobs = list(pp3a, pp3b),
             nrow = 2, ncol = 1)
dev.off()
############################################################



############################################################
# FIGURE 4 - individual level effect, and reduction by mask efficacy group
parms = read.csv(paste0(dir_code, 'parm.bounds.csv'))
parm.t = parms %>% filter(parm == "ave.p.redn.mask.byage" & age <=85)
parm.t$default %>% mean

diff.ctc.redn.bymask.t = F # basic baseline settings
diff.ctc.bymask.t = F
# for the diff measures, b/c repeated outcomes take time to observed, only plot outcomes for the last 3 seasons
sn.t = tail(seasons, 1)
mea.t = 'infections'
timeframe.t = 'cum2date'; mask.grp.t = 'use mask'; age.grp.t = 'all'
mask.grps.t = c("mask.lowE","mask.midE","mask.highE")  # , 'use mask', 'all'
mask.grps.labs.t = c('<50%', '50-90%', '>90%') #  c('low', 'medium', 'high')
age.grps.labs = paste('Age',age.grps$label)
age.grps.t = age.grps.labs # [c(3,5,8)]
ave.frac.masking.t = 0.5; cut.mask.start.t = '0.1%'

RES = res_main_diff.byeff.bysn %>% 
  filter(diff.ctc.bymask == diff.ctc.bymask.t & diff.ctc.redn.bymask == diff.ctc.redn.bymask.t) %>%
  mutate(cut.mask.start = factor(cut.mask.start, levels = cut.mask.start_vec, labels = cut.mask.start_vec.lab),
         ave.frac.masking = factor(ave.frac.masking, levels = frac.masking_vec)
  ) %>%
  mutate(measure = case_when(grepl('AR', variable) ~ 'infections',
                             grepl('Hosp', variable) ~ 'hospitalizations'),
         mea.metric = case_when(grepl('tot', variable) ~ 'total',
                                grepl('1plus', variable) & grepl('AR', variable) ~ 'individuals infected >=1 time',
                                grepl('2plus', variable) & grepl('AR', variable) ~ 'individuals infected >=2 times',
                                grepl('3plus', variable) & grepl('AR', variable) ~ 'individuals infected >=3 times',
                                grepl('1plus', variable) & grepl('Hosp', variable) ~ 'individuals hospitalized >=1 time',
                                grepl('2plus', variable) & grepl('Hosp', variable) ~ 'individuals hospitalized >=2 times',
                                grepl('3plus', variable) & grepl('Hosp', variable) ~ 'individuals hospitalized >=3 times')) %>%
  mutate(season = factor(season, levels = seasons, labels = seasons.labs))
# plot it
mea.metrics = RES$mea.metric %>% unique
measures = RES$measure %>% unique

mea.metrics.t = paste('individuals',ifelse(mea.t=='infections','infected','hospitalized'),c('>=1 time','>=2 times','>=3 times'))

tda = RES %>% filter(metric == metric.t & timeframe == timeframe.t & 
                       measure == mea.t & mea.metric  %in% mea.metrics.t &
                       season == sn.t & mask.grp == mask.grp.t & 
                       age.grp == age.grp.t & ave.frac.masking !=1)
layer = tda %>% filter(round(med,0) < eff.cut) # above a certain cut off
d.sig = tda %>% filter(ci95lwr > 0)
pp4a = ggplot(tda, 
             aes(x = cut.mask.start, y = ave.frac.masking)) +  # use the median - the mean can be very skewed 
  geom_tile(aes(fill = med)) + 
  # geom_tile(data=layer, alpha = 0.0, color = "black", size = 1, linejoin = "round") +
  # geom_tile(data=layer, alpha = 1, aes(fill = med)) +
  geom_text(aes(x = cut.mask.start, y = ave.frac.masking, label = med %>% round(., 0)), size = 3) +
  geom_text(data=d.sig, aes(x = cut.mask.start, y = ave.frac.masking, label = med %>% round(., 0)), size = 3, fontface = "bold") +
  facet_rep_wrap(~ mea.metric, repeat.tick.labels = T, ncol = length(mea.metrics.t), labeller = label_wrap_gen(multi_line=FALSE, width = 50)) +
  scale_fill_distiller(palette = "RdPu", direction = 1, na.value = 'lightgrey', limits = c(0, 100)) +
  labs(x = 'time to start masking', y = 'fraction pop masking', fill = '% redn', 
       title = paste0('(A) All age groups and all mask efficacy groups combined, percent reduction in')) +
  theme_minimal() + theme.t

# dose response
tda2 = RES %>% filter(metric == metric.t & timeframe == timeframe.t & season == sn.t & 
                        measure == mea.t & mea.metric == mea.metric.t &
                        age.grp != 'all' & 
                        ave.frac.masking == ave.frac.masking.t & cut.mask.start == cut.mask.start.t) %>% 
  filter(mask.grp %in% mask.grps.t) %>%
  mutate(mask.grp = factor(mask.grp, levels = mask.grps.t, labels = mask.grps.labs.t)) %>%
  mutate(age.grp = factor(age.grp, levels = age.grps$label, labels = age.grps.labs)) %>%
  filter(age.grp %in% age.grps.t)
layer = tda2 %>% filter(round(med,0) < eff.cut) # above a certain cut off
d.sig = tda2 %>% filter(ci95lwr > 0)
ymax.t = tda2$ci95upr %>% max %>% ceiling()
pp4b = ggplot(tda2, aes(x = mask.grp, color = mask.grp)) +
  geom_boxplot(aes(
    lower = ci50lwr, 
    upper = ci50upr, 
    middle = med, 
    ymin = ci95lwr, 
    ymax = ci95upr,
    group= mask.grp),
    stat = "identity"
  ) + 
  facet_rep_wrap(~age.grp, ncol = 7, labeller = label_wrap_gen(multi_line=FALSE)) + # repeat.tick.labels = 'x', 
  labs(x = 'mask efficacy', y = '% reduction', color = 'efficacy', 
       title = paste0('(B) Precent reduction in total ', mea.t, ': by age group, by mask efficacy\n(example: ',ave.frac.masking.t*100,'% people use mask when prevalence is >', cut.mask.start.t,')')) +
  lims(y = c(0, ymax.t)) + 
  theme_minimal() + theme.t4b

pdf(paste0(dir_fig, 'Fig4.pdf'), width = length(mea.metrics.t)*2.5, height = 5)
grid.arrange(grobs = list(pp4a, pp4b),
             nrow = 2, ncol = 1)
dev.off()
############################################################


# SUPPLEMENTAL FIGURES
################################################
## Fig S1 - sensitivity analysis re weights of nonphysical contact
wt_nonphys_vec = c(.1, .2, .5, 1); # assume non-phyical contact has 20% of the chance resulting infection
metric.t = 'perc.redn'; mask.grp.t = 'all'; age.grp.t = 'all'
timeframes = res_nonphywt_diff.bysn$timeframe %>% unique
names(timeframes) = factor(timeframes, levels = c("cum2date", "this.season"), labels = c('cumulative', 'this season'))
mea.metrics = res_nonphywt_diff.bysn$mea.metric %>% unique
seasons.t = unique(res_nonphywt_diff.bysn$season) %>% sort
seasons = unique(res_nonphywt_diff.bysn$season) %>% sort
seasons.labs = gsub('2019-2020','March-June 2020', seasons)
diff.ctc.bymask.t = F
diff.ctc.redn.bymask.t = F
timeframe.t = 'cum2date'; mea.metric.t = 'total'; mea.t = 'infections'

RES = res_nonphywt_diff.bysn %>% 
  filter(diff.ctc.bymask == diff.ctc.bymask.t & diff.ctc.redn.bymask == diff.ctc.redn.bymask.t) %>%
  mutate(cut.mask.start1 = factor(cut.mask.start1, levels = cut.mask.start_vec, labels = cut.mask.start_vec.lab),
         ave.frac.masking1 = factor(ave.frac.masking1, levels = frac.masking_vec, labels = paste0(frac.masking_vec*100,'% use masks')),
         wt_nonphys = factor(wt_nonphys, levels = wt_nonphys_vec)
  ) %>%
  mutate(measure = case_when(grepl('AR', variable) ~ 'infections',
                             grepl('Hosp', variable) ~ 'hospitalizations'),
         mea.metric = case_when(grepl('tot', variable) ~ 'total',
                                grepl('1plus', variable) ~ 'individuals affected >=1 time',
                                grepl('2plus', variable) ~ 'individuals affected >=2 times',
                                grepl('3plus', variable) ~ 'individuals affected >=3 times')) %>%
  mutate(season = factor(season, levels = seasons, labels = seasons.labs))

tda = RES %>% filter(mask.grp == mask.grp.t & age.grp == age.grp.t & 
                       metric == metric.t & timeframe == timeframe.t & 
                       measure == mea.t & mea.metric  %in% mea.metric.t)
layer = tda %>% filter(round(med,0) < eff.cut) # above a certain cut off
# significantly diff
d.sig = tda %>% filter(ci95lwr > 0)

pp.s1 = ggplot(tda, 
            aes(x = cut.mask.start1, y = wt_nonphys)) +  # use the median - the mean can be very skewed 
  geom_tile(aes(fill = med)) + 
  geom_tile(data=layer, alpha = 0.0, color = "black", size = 1, linejoin = "round") +
  geom_tile(data=layer, alpha = 1, aes(fill = med)) +
  geom_text(aes(x = cut.mask.start1, y = wt_nonphys, label = med %>% round(., 0)), size = 3) +
  geom_text(data=d.sig, aes(x = cut.mask.start1, y = wt_nonphys, label = med %>% round(., 0)), size = 3, fontface = "bold") +
  facet_rep_wrap(~ave.frac.masking1 + season, repeat.tick.labels = T, ncol = length(seasons.t), labeller = label_wrap_gen(multi_line=FALSE, width = 50)) +
  scale_fill_distiller(palette = "RdPu", direction = 1, na.value = 'lightgrey', limits = c(0, 100)) +
  labs(x = 'time to start masking', y = 'relative probability of transmission:\nnonphysical vs. physical contact', fill = '% redn', 
       title = paste0('Percent reduction in ', ifelse(mea.metric.t=='total',
                                                      paste('total', mea.t),
                                                      gsub('affected', ifelse(mea.t=='infections','infected','hospitalized'),mea.metrics[im])
       ),  ': ', ifelse(timeframe.t=='cum2date','tallies from March 2020 through this season','tallies for this season'))) +
  theme_minimal() + theme.t

pdf(paste0(dir_fig, 'FigS1_nonphywt.pdf'), width = length(seasons.t) * 2.5, height = 5)
print(pp.s1)
dev.off()
################################################

################################################
## Fig S2 - contact reduction
ifig.start = 1
timeframe.t = 'cum2date'; mea.metric.t = 'total'; mea.t = 'infections'
age.grp.t = 'all'; mask.grp.t = 'all'
tda = res_main_diff.bysn %>% 
  filter(diff.ctc.bymask == diff.ctc.bymask.t & diff.ctc.redn.bymask == T) %>%
  mutate(cut.mask.start1 = factor(cut.mask.start1, levels = cut.mask.start_vec, labels = cut.mask.start_vec.lab),
         ave.frac.masking1 = factor(ave.frac.masking1, levels = frac.masking_vec)
  ) %>%
  mutate(measure = case_when(grepl('AR', variable) ~ 'infections',
                             grepl('Hosp', variable) ~ 'hospitalizations'),
         mea.metric = case_when(grepl('tot', variable) ~ 'total',
                                grepl('1plus', variable) ~ 'individuals affected >=1 time',
                                grepl('2plus', variable) ~ 'individuals affected >=2 times',
                                grepl('3plus', variable) ~ 'individuals affected >=3 times')) %>%
  mutate(season = factor(season, levels = seasons, labels = seasons.labs)) %>% 
  filter(mask.grp == mask.grp.t & age.grp == age.grp.t & 
           metric == metric.t & timeframe == timeframe.t & 
           measure == mea.t & mea.metric  %in% mea.metric.t)
layer = tda %>% filter(round(med,0) < eff.cut) # above a certain cut off
# significantly diff
d.sig = tda %>% filter(ci95lwr > 0)
pp.s2a = ggplot(tda, 
              aes(x = cut.mask.start1, y = ave.frac.masking1)) +  # use the median - the mean can be very skewed 
  geom_tile(aes(fill = med)) + 
  geom_tile(data=layer, alpha = 0.0, color = "black", size = 1, linejoin = "round") +
  geom_tile(data=layer, alpha = 1, aes(fill = med)) +
  geom_text(aes(x = cut.mask.start1, y = ave.frac.masking1, label = med %>% round(., 0)), size = 3) +
  geom_text(data=d.sig, aes(x = cut.mask.start1, y = ave.frac.masking1, label = med %>% round(., 0)), size = 3, fontface = "bold") +
  facet_rep_wrap(~ season, repeat.tick.labels = T, ncol = length(seasons.t), labeller = label_wrap_gen(multi_line=FALSE)) +
  scale_fill_distiller(palette = "RdPu", direction = 1, na.value = 'lightgrey', limits = c(0, 100)) +
  labs(x = 'time to start masking', y = 'fraction pop masking', fill = '% redn', 
       title = paste0('(',LETTERS[ifig.start],') Percent reduction in total ', mea.t, ': ', ifelse(timeframe.t=='cum2date','tallies from March 2020 through this season','tallies for this season'))) +
  theme_minimal() + theme.t


timeframe.t = 'cum2date'; mea.metric.t = 'individuals affected >=1 time'; 
tda = res_main_diff.bysn %>% 
  filter(diff.ctc.bymask == diff.ctc.bymask.t & diff.ctc.redn.bymask == T) %>%
  mutate(cut.mask.start1 = factor(cut.mask.start1, levels = cut.mask.start_vec, labels = cut.mask.start_vec.lab),
         ave.frac.masking1 = factor(ave.frac.masking1, levels = frac.masking_vec)
  ) %>%
  mutate(measure = case_when(grepl('AR', variable) ~ 'infections',
                             grepl('Hosp', variable) ~ 'hospitalizations'),
         mea.metric = case_when(grepl('tot', variable) ~ 'total',
                                grepl('1plus', variable) ~ 'individuals affected >=1 time',
                                grepl('2plus', variable) ~ 'individuals affected >=2 times',
                                grepl('3plus', variable) ~ 'individuals affected >=3 times')) %>%
  mutate(season = factor(season, levels = seasons, labels = seasons.labs))%>% 
  filter(mask.grp == mask.grp.t & age.grp == age.grp.t & 
           metric == metric.t & timeframe == timeframe.t & 
           measure == mea.t & mea.metric  %in% mea.metric.t)
layer = tda %>% filter(round(med,0) < eff.cut) # above a certain cut off
# significantly diff
d.sig = tda %>% filter(ci95lwr > 0)
pp.s2b = ggplot(tda, 
              aes(x = cut.mask.start1, y = ave.frac.masking1)) +  # use the median - the mean can be very skewed 
  geom_tile(aes(fill = med)) + 
  geom_tile(data=layer, alpha = 0.0, color = "black", size = 1, linejoin = "round") +
  geom_tile(data=layer, alpha = 1, aes(fill = med)) +
  geom_text(aes(x = cut.mask.start1, y = ave.frac.masking1, label = med %>% round(., 0)), size = 3) +
  geom_text(data=d.sig, aes(x = cut.mask.start1, y = ave.frac.masking1, label = med %>% round(., 0)), size = 3, fontface = "bold") +
  facet_rep_wrap(~ season, repeat.tick.labels = T, ncol = length(seasons.t), labeller = label_wrap_gen(multi_line=FALSE)) +
  scale_fill_distiller(palette = "RdPu", direction = 1, na.value = 'lightgrey', limits = c(0, 100)) +
  labs(x = 'time to start masking', y = 'fraction pop masking', fill = '% redn', 
       title = paste0('(',LETTERS[ifig.start+1],') Percent reduction in individuals ever-infected: ', 
                      # gsub('affected',ifelse(mea.t=='infections','infected','hospitalized'),mea.metric.t) ,': ', 
                      ifelse(timeframe.t=='cum2date','tallies from March 2020 through this season','tallies for this season'))) +
  theme_minimal() + theme.t

pdf(paste0(dir_fig, 'FigS2_ctc.redn.pdf'), width = 10, height = 2.5*2)
grid.arrange(grobs = list(pp.s2a, pp.s2b),
             nrow = 2, ncol = 1)
dev.off()
################################################


################################################
# Fig S3 and Fig S4 - higher R0 (x2 the value used in the main analysis)
metric.t = 'perc.redn'; mask.grp.t = 'all'; age.grp.t = 'all'
timeframes = res_R0x2_diff.bysn$timeframe %>% unique
names(timeframes) = factor(timeframes, levels = c("cum2date", "this.season"), labels = c('cumulative', 'this season'))
mea.metrics = res_R0x2_diff.bysn$mea.metric %>% unique
seasons.t = unique(res_R0x2_diff.bysn$season) %>% sort
seasons = unique(res_R0x2_diff.bysn$season) %>% sort
seasons.labs = gsub('2019-2020','March-June 2020', seasons)

timeframe.t = 'cum2date'; mea.metric.t = 'total'; mea.t = 'infections'
tda = res_R0x2_diff.bysn %>% 
  filter(diff.ctc.bymask == diff.ctc.bymask.t & diff.ctc.redn.bymask == diff.ctc.redn.bymask.t) %>%
  mutate(cut.mask.start1 = factor(cut.mask.start1, levels = cut.mask.start_vec, labels = cut.mask.start_vec.lab),
         ave.frac.masking1 = factor(ave.frac.masking1, levels = frac.masking_vec)
  ) %>%
  mutate(measure = case_when(grepl('AR', variable) ~ 'infections',
                             grepl('Hosp', variable) ~ 'hospitalizations'),
         mea.metric = case_when(grepl('tot', variable) ~ 'total',
                                grepl('1plus', variable) ~ 'individuals affected >=1 time',
                                grepl('2plus', variable) ~ 'individuals affected >=2 times',
                                grepl('3plus', variable) ~ 'individuals affected >=3 times')) %>%
  mutate(season = factor(season, levels = seasons, labels = seasons.labs))%>% 
  filter(mask.grp == mask.grp.t & age.grp == age.grp.t & 
           metric == metric.t & timeframe == timeframe.t & 
           measure == mea.t & mea.metric  %in% mea.metric.t)
layer = tda %>% filter(round(med,0) < eff.cut) # above a certain cut off
# significantly diff
d.sig = tda %>% filter(round(ci95lwr,0) > 0)
pp.s3a = ggplot(tda, 
              aes(x = cut.mask.start1, y = ave.frac.masking1)) +  # use the median - the mean can be very skewed 
  geom_tile(aes(fill = med)) + 
  geom_tile(data=layer, alpha = 0.0, color = "black", size = 1, linejoin = "round") +
  geom_tile(data=layer, alpha = 1, aes(fill = med)) +
  geom_text(aes(x = cut.mask.start1, y = ave.frac.masking1, label = med %>% round(., 0)), size = 3) +
  geom_text(data=d.sig, aes(x = cut.mask.start1, y = ave.frac.masking1, label = med %>% round(., 0)), size = 3, fontface = "bold") +
  facet_rep_wrap(~ season, repeat.tick.labels = T, ncol = length(seasons.t), labeller = label_wrap_gen(multi_line=FALSE)) +
  scale_fill_distiller(palette = "RdPu", direction = 1, na.value = 'lightgrey', limits = c(0, 100)) +
  labs(x = 'time to start masking', y = 'fraction pop masking', fill = '% redn', 
       title = paste0('(',LETTERS[1],') Percent reduction in total ', mea.t, ': ', ifelse(timeframe.t=='cum2date','tallies from March 2020 through this season','tallies for this season'))) +
  theme_minimal() + theme.t


timeframe.t = 'cum2date'; mea.metric.t = 'individuals affected >=1 time'; 
tda = res_R0x2_diff.bysn %>% 
  filter(diff.ctc.bymask == diff.ctc.bymask.t & diff.ctc.redn.bymask == diff.ctc.redn.bymask.t) %>%
  mutate(cut.mask.start1 = factor(cut.mask.start1, levels = cut.mask.start_vec, labels = cut.mask.start_vec.lab),
         ave.frac.masking1 = factor(ave.frac.masking1, levels = frac.masking_vec)
  ) %>%
  mutate(measure = case_when(grepl('AR', variable) ~ 'infections',
                             grepl('Hosp', variable) ~ 'hospitalizations'),
         mea.metric = case_when(grepl('tot', variable) ~ 'total',
                                grepl('1plus', variable) ~ 'individuals affected >=1 time',
                                grepl('2plus', variable) ~ 'individuals affected >=2 times',
                                grepl('3plus', variable) ~ 'individuals affected >=3 times')) %>%
  mutate(season = factor(season, levels = seasons, labels = seasons.labs))%>% 
  filter(mask.grp == mask.grp.t & age.grp == age.grp.t & 
           metric == metric.t & timeframe == timeframe.t & 
           measure == mea.t & mea.metric  %in% mea.metric.t)
layer = tda %>% filter(round(med,0) < eff.cut) # above a certain cut off
# significantly diff
d.sig = tda %>% filter(round(ci95lwr,0) > 0)
pp.s3b = ggplot(tda, 
              aes(x = cut.mask.start1, y = ave.frac.masking1)) +  # use the median - the mean can be very skewed 
  geom_tile(aes(fill = med)) + 
  geom_tile(data=layer, alpha = 0.0, color = "black", size = 1, linejoin = "round") +
  geom_tile(data=layer, alpha = 1, aes(fill = med)) +
  geom_text(aes(x = cut.mask.start1, y = ave.frac.masking1, label = med %>% round(., 0)), size = 3) +
  geom_text(data=d.sig, aes(x = cut.mask.start1, y = ave.frac.masking1, label = med %>% round(., 0)), size = 3, fontface = "bold") +
  facet_rep_wrap(~ season, repeat.tick.labels = T, ncol = length(seasons.t), labeller = label_wrap_gen(multi_line=FALSE)) +
  scale_fill_distiller(palette = "RdPu", direction = 1, na.value = 'lightgrey', limits = c(0, 100)) +
  labs(x = 'time to start masking', y = 'fraction pop masking', fill = '% redn', 
       title = paste0('(',LETTERS[2],') Percent reduction in individuals ever-infected: ', 
                      # gsub('affected', ifelse(mea.t=='infections','infected','hospitalized'),mea.metric.t),': ', 
                      ifelse(timeframe.t=='cum2date','tallies from March 2020 through this season','tallies for this season'))) +
  theme_minimal() + theme.t


timeframe.t = 'this.season'; mea.metric.t = 'total'; 
tda = res_R0x2_diff.bysn %>% 
  filter(diff.ctc.bymask == diff.ctc.bymask.t & diff.ctc.redn.bymask == diff.ctc.redn.bymask.t) %>%
  mutate(cut.mask.start1 = factor(cut.mask.start1, levels = cut.mask.start_vec, labels = cut.mask.start_vec.lab),
         ave.frac.masking1 = factor(ave.frac.masking1, levels = frac.masking_vec)
  ) %>%
  mutate(measure = case_when(grepl('AR', variable) ~ 'infections',
                             grepl('Hosp', variable) ~ 'hospitalizations'),
         mea.metric = case_when(grepl('tot', variable) ~ 'total',
                                grepl('1plus', variable) ~ 'individuals affected >=1 time',
                                grepl('2plus', variable) ~ 'individuals affected >=2 times',
                                grepl('3plus', variable) ~ 'individuals affected >=3 times')) %>%
  mutate(season = factor(season, levels = seasons, labels = seasons.labs))%>% 
  filter(mask.grp == mask.grp.t & age.grp == age.grp.t & 
           metric == metric.t & timeframe == timeframe.t & 
           measure == mea.t & mea.metric  %in% mea.metric.t)
layer = tda %>% filter(round(med,0) < eff.cut) # above a certain cut off
# significantly diff
d.sig = tda %>% filter(round(ci95lwr,0) > 0)
pp.s3c = ggplot(tda, 
              aes(x = cut.mask.start1, y = ave.frac.masking1)) +  # use the median - the mean can be very skewed 
  geom_tile(aes(fill = med)) + 
  geom_tile(data=layer, alpha = 0.0, color = "black", size = 1, linejoin = "round") +
  geom_tile(data=layer, alpha = 1, aes(fill = med)) +
  geom_text(aes(x = cut.mask.start1, y = ave.frac.masking1, label = med %>% round(., 0)), size = 3) +
  geom_text(data=d.sig, aes(x = cut.mask.start1, y = ave.frac.masking1, label = med %>% round(., 0)), size = 3, fontface = "bold") +
  facet_rep_wrap(~ season, repeat.tick.labels = T, ncol = length(seasons.t), labeller = label_wrap_gen(multi_line=FALSE)) +
  scale_fill_distiller(palette = "RdPu", direction = 1, na.value = 'lightgrey', limits = c(0, 100)) +
  labs(x = 'time to start masking', y = 'fraction pop masking', fill = '% redn', 
       title = paste0('(',LETTERS[3],') Percent reduction in total ', mea.t, ': ', ifelse(timeframe.t=='cum2date','tallies from March 2020 through this season','tallies for this season'))) +
  theme_minimal() + theme.t

pdf(paste0(dir_fig, 'FigS3_redn_R0x2.pdf'), width = 10, height = 2.5*3)
grid.arrange(grobs = list(pp.s3a, pp.s3b, pp.s3c),
             nrow = 3, ncol = 1)
dev.off()

## Fig S4
diff.ctc.redn.bymask.t = F # basic baseline settings
diff.ctc.bymask.t = F
# for the diff measures, b/c repeated outcomes take time to observed, only plot outcomes for the last 3 seasons
sn.t = tail(seasons, 1)
mea.t = 'infections'
timeframe.t = 'cum2date'; mask.grp.t = 'use mask'; age.grp.t = 'all'
mask.grps.t = c("mask.lowE","mask.midE","mask.highE")  # , 'use mask', 'all'
mask.grps.labs.t = c('<50%', '50-90%', '>90%') #  c('low', 'medium', 'high')
age.grps.labs = paste('Age',age.grps$label)
age.grps.t = age.grps.labs # [c(3,5,8)]
ave.frac.masking.t = 0.5; cut.mask.start.t = '0.1%'

RES = res_R0x2_diff.byeff.bysn %>% 
  filter(diff.ctc.bymask == diff.ctc.bymask.t & diff.ctc.redn.bymask == diff.ctc.redn.bymask.t) %>%
  mutate(cut.mask.start = factor(cut.mask.start, levels = cut.mask.start_vec, labels = cut.mask.start_vec.lab),
         ave.frac.masking = factor(ave.frac.masking, levels = frac.masking_vec)
  ) %>%
  mutate(measure = case_when(grepl('AR', variable) ~ 'infections',
                             grepl('Hosp', variable) ~ 'hospitalizations'),
         mea.metric = case_when(grepl('tot', variable) ~ 'total',
                                grepl('1plus', variable) & grepl('AR', variable) ~ 'individuals infected >=1 time',
                                grepl('2plus', variable) & grepl('AR', variable) ~ 'individuals infected >=2 times',
                                grepl('3plus', variable) & grepl('AR', variable) ~ 'individuals infected >=3 times',
                                grepl('1plus', variable) & grepl('Hosp', variable) ~ 'individuals hospitalized >=1 time',
                                grepl('2plus', variable) & grepl('Hosp', variable) ~ 'individuals hospitalized >=2 times',
                                grepl('3plus', variable) & grepl('Hosp', variable) ~ 'individuals hospitalized >=3 times')) %>%
  mutate(season = factor(season, levels = seasons, labels = seasons.labs))

mea.metrics.t = paste('individuals',ifelse(mea.t=='infections','infected','hospitalized'),c('>=1 time','>=2 times','>=3 times'))

tda = RES %>% filter(metric == metric.t & timeframe == timeframe.t & 
                       measure == mea.t & mea.metric  %in% mea.metrics.t &
                       season == sn.t & mask.grp == mask.grp.t & 
                       age.grp == age.grp.t & ave.frac.masking !=1)
layer = tda %>% filter(round(med,0) < eff.cut) # above a certain cut off
d.sig = tda %>% filter(ci95lwr > 0)
pp.s4a = ggplot(tda, 
              aes(x = cut.mask.start, y = ave.frac.masking)) +  # use the median - the mean can be very skewed 
  geom_tile(aes(fill = med)) + 
  # geom_tile(data=layer, alpha = 0.0, color = "black", size = 1, linejoin = "round") +
  # geom_tile(data=layer, alpha = 1, aes(fill = med)) +
  geom_text(aes(x = cut.mask.start, y = ave.frac.masking, label = med %>% round(., 0)), size = 3) +
  geom_text(data=d.sig, aes(x = cut.mask.start, y = ave.frac.masking, label = med %>% round(., 0)), size = 3, fontface = "bold") +
  facet_rep_wrap(~ mea.metric, repeat.tick.labels = T, ncol = length(mea.metrics.t), labeller = label_wrap_gen(multi_line=FALSE, width = 50)) +
  scale_fill_distiller(palette = "RdPu", direction = 1, na.value = 'lightgrey', limits = c(0, 100)) +
  labs(x = 'time to start masking', y = 'fraction pop masking', fill = '% redn', 
       title = paste0('(A) All age groups and all mask efficacy groups combined, percent reduction in')) +
  theme_minimal() + theme.t

# dose response
tda2 = RES %>% filter(metric == metric.t & timeframe == timeframe.t & season == sn.t & 
                        measure == mea.t & mea.metric == mea.metric.t &
                        age.grp != 'all' & 
                        ave.frac.masking == ave.frac.masking.t & cut.mask.start == cut.mask.start.t) %>% 
  filter(mask.grp %in% mask.grps.t) %>%
  mutate(mask.grp = factor(mask.grp, levels = mask.grps.t, labels = mask.grps.labs.t)) %>%
  mutate(age.grp = factor(age.grp, levels = age.grps$label, labels = age.grps.labs)) %>%
  filter(age.grp %in% age.grps.t)
layer = tda2 %>% filter(round(med,0) < eff.cut) # above a certain cut off
d.sig = tda2 %>% filter(ci95lwr > 0)
ymax.t = tda2$ci95upr %>% max %>% ceiling()
pp.s4b = ggplot(tda2, aes(x = mask.grp, color = mask.grp)) +
  geom_boxplot(aes(
    lower = ci50lwr, 
    upper = ci50upr, 
    middle = med, 
    ymin = ci95lwr, 
    ymax = ci95upr,
    group= mask.grp),
    stat = "identity"
  ) + 
  facet_rep_wrap(~age.grp, ncol = 7, labeller = label_wrap_gen(multi_line=FALSE)) + # repeat.tick.labels = 'x', 
  labs(x = 'mask efficacy', y = '% reduction', color = 'efficacy', 
       title = paste0('(B) Precent reduction in total ', mea.t, ': by age group, by mask efficacy\n(example: ',ave.frac.masking.t*100,'% people use mask when prevalence is >', cut.mask.start.t,')')) +
  lims(y = c(0, ymax.t)) + 
  theme_minimal() + theme.t4b

pdf(paste0(dir_fig, 'FigS4.pdf'), width = length(mea.metrics.t)*2.5, height = 5)
grid.arrange(grobs = list(pp.s4a, pp.s4b),
             nrow = 2, ncol = 1)
dev.off()

################################################


################################################
## Fig S5 - hospitalizations
mea.t = 'hospitalizations'
diff.ctc.redn.bymask.t = F
diff.ctc.bymask.t = F
n.tot = 1e5
mask.grp.t = 'all'; age.grp.t = 'all'
ifig.start = 1
timeframe.t = 'cum2date'; mea.metric.t = 'total'; 
pp.s5 = list(3)
tda = res_main_diff.bysn %>% 
  filter(diff.ctc.bymask == diff.ctc.bymask.t & diff.ctc.redn.bymask == diff.ctc.redn.bymask.t) %>%
  mutate(cut.mask.start1 = factor(cut.mask.start1, levels = cut.mask.start_vec, labels = cut.mask.start_vec.lab),
         ave.frac.masking1 = factor(ave.frac.masking1, levels = frac.masking_vec)
  ) %>%
  mutate(measure = case_when(grepl('AR', variable) ~ 'infections',
                             grepl('Hosp', variable) ~ 'hospitalizations'),
         mea.metric = case_when(grepl('tot', variable) ~ 'total',
                                grepl('1plus', variable) ~ 'individuals affected >=1 time',
                                grepl('2plus', variable) ~ 'individuals affected >=2 times',
                                grepl('3plus', variable) ~ 'individuals affected >=3 times')) %>%
  mutate(season = factor(season, levels = seasons, labels = seasons.labs))%>% 
  filter(mask.grp == mask.grp.t & age.grp == age.grp.t & 
           metric == metric.t & timeframe == timeframe.t & 
           measure == mea.t & mea.metric  %in% mea.metric.t)
layer = tda %>% filter(round(med,0) < eff.cut) # above a certain cut off
# significantly diff
d.sig = tda %>% filter(ci95lwr > 0)
pp.s5[[1]] = ggplot(tda, 
              aes(x = cut.mask.start1, y = ave.frac.masking1)) +  # use the median - the mean can be very skewed 
  geom_tile(aes(fill = med)) + 
  geom_tile(data=layer, alpha = 0.0, color = "black", size = 1, linejoin = "round") +
  geom_tile(data=layer, alpha = 1, aes(fill = med)) +
  geom_text(aes(x = cut.mask.start1, y = ave.frac.masking1, label = med %>% round(., 0)), size = 3) +
  geom_text(data=d.sig, aes(x = cut.mask.start1, y = ave.frac.masking1, label = med %>% round(., 0)), size = 3, fontface = "bold") +
  facet_rep_wrap(~ season, repeat.tick.labels = T, ncol = length(seasons.t), labeller = label_wrap_gen(multi_line=FALSE)) +
  scale_fill_distiller(palette = "RdPu", direction = 1, na.value = 'lightgrey', limits = c(0, 100)) +
  labs(x = 'time to start masking', y = 'fraction pop masking', fill = '% redn', 
       title = paste0('(',LETTERS[1],') Percent reduction in total ', mea.t,': ', ifelse(timeframe.t=='cum2date','tallies from March 2020 through this season','tallies for this season'))) +
  theme_minimal() + theme.t


timeframe.t = 'cum2date'; # mea.metric.t = 'individuals affected >=1 time'; 
cnt = 1
mea.metric_vec.t = c('individuals affected >=1 time','individuals affected >=2 times','individuals affected >=3 times')
for(mea.metric.t in mea.metric_vec.t){
  cnt = cnt +1
  tda = res_main_diff.bysn %>% 
    filter(diff.ctc.bymask == diff.ctc.bymask.t & diff.ctc.redn.bymask == diff.ctc.redn.bymask.t) %>%
    mutate(cut.mask.start1 = factor(cut.mask.start1, levels = cut.mask.start_vec, labels = cut.mask.start_vec.lab),
           ave.frac.masking1 = factor(ave.frac.masking1, levels = frac.masking_vec)
    ) %>%
    mutate(measure = case_when(grepl('AR', variable) ~ 'infections',
                               grepl('Hosp', variable) ~ 'hospitalizations'),
           mea.metric = case_when(grepl('tot', variable) ~ 'total',
                                  grepl('1plus', variable) ~ 'individuals affected >=1 time',
                                  grepl('2plus', variable) ~ 'individuals affected >=2 times',
                                  grepl('3plus', variable) ~ 'individuals affected >=3 times')) %>%
    mutate(season = factor(season, levels = seasons, labels = seasons.labs))%>% 
    filter(mask.grp == mask.grp.t & age.grp == age.grp.t & 
             metric == metric.t & timeframe == timeframe.t & 
             measure == mea.t & mea.metric  %in% mea.metric.t)
  
  if(mea.metric.t %in% mea.metric_vec.t[2:3]){
    tda = tda %>%
      mutate(med = case_when(season == 'March-June 2020' ~ NA_real_, 
                             T ~ med))
  }
  
  layer = tda %>% filter(round(med,0) < eff.cut) # above a certain cut off
  d.sig = tda %>% filter(ci95lwr > 0)
  pp.s5[[cnt]] = ggplot(tda, 
               aes(x = cut.mask.start1, y = ave.frac.masking1)) +  # use the median - the mean can be very skewed 
    geom_tile(aes(fill = med)) + 
    geom_tile(data=layer, alpha = 0.0, color = "black", size = 1, linejoin = "round") +
    geom_tile(data=layer, alpha = 1, aes(fill = med)) +
    geom_text(aes(x = cut.mask.start1, y = ave.frac.masking1, label = med %>% round(., 0)), size = 3) +
    geom_text(data=d.sig, aes(x = cut.mask.start1, y = ave.frac.masking1, label = med %>% round(., 0)), size = 3, fontface = "bold") +
    facet_rep_wrap(~ season, repeat.tick.labels = T, ncol = length(seasons.t), labeller = label_wrap_gen(multi_line=FALSE)) +
    scale_fill_distiller(palette = "RdPu", direction = 1, na.value = 'lightgrey', limits = c(0, 100)) +
    labs(x = 'time to start masking', y = 'fraction pop masking', fill = '% redn', 
         title = paste0('(',LETTERS[cnt],') Percent reduction in', gsub('affected',ifelse(mea.t=='infections','infected','hospitalized'),mea.metric.t),': ', ifelse(timeframe.t=='cum2date','tallies from March 2020 through this season','tallies for this season'))) +
    theme_minimal() + theme.t
  
}

pdf(paste0(dir_fig, 'FigS5_hospitalizations.pdf'), width = 10, height = 7.5)
grid.arrange(grobs = list(pp.s5[[1]],pp.s5[[2]],pp.s5[[3]]),
             nrow = 3, ncol = 1)
dev.off()
################################################

################################################
# Fig S6 (vaccination)
load("../data/da.vx.daily.RData")  # vaccination settings: NYC vax data used in the simulation
tda = da.vx.daily %>% left_join(x = ., y = age.grps %>% dplyr::select(age.start, label), 
                                by = c('age.start')) %>%
  mutate(age.grp = factor(label, levels = age.grps$label, labels = paste('Age',age.grps$label)))
dates.t = tda %>% .$date %>% unique %>% as.Date %>% sort
dates.t = seq(dates.t[1], tail(dates.t, 1), by = '6 months')

pp.s6 = ggplot(tda, aes(x = date, y = value)) +
   geom_line() +
  scale_x_date(breaks = dates.t, labels = format(dates.t,'%b\n%Y'), limits = as.Date(c('2020/3/1','2024/12/31')),expand=c(.01,.01)) +
  labs(y = '% population vaccinated per day', 
       x = '') +
  facet_rep_wrap(~age.grp, ncol=2, repeat.tick.labels = T) + theme_minimal() + theme.t2


pdf(paste0(dir_fig, 'FigS6_vacc.pdf'), width = 8, height = 8)
print(pp.s6)
dev.off()
################################################

