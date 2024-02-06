# Model code used in Yang & Shaman "Reconciling the efficacy and effectiveness of masking on epidemic outcomes"
# author: Wan Yang
# date: March 2023
# note: this code is for research use only; may not be entirely optimized nor neatly organized. 

# script to load parameter settings based on external data sources/estimates


# set simulation time
# more detailed date/time setting to much with the calendar (for leap years)
if(!exists('date.start'))
  date.start = as.Date('2020/3/1')
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
ndays = nrow(tab_dates)
ndays / 365

if(!exists('MASKATHOME'))
  MASKATHOME = F # set this to TRUE only when testing a very specific scenario that ppl also mask at home

if(!exists('parm.type.t'))
  parm.type.t = 'default' # which parm set to use


if(MASKATHOME){
  print('WARNING: MASKATHOME SET TO TRUE')
}


# to reduce model stochasticity and considering high chance of epidemic and introduction
days.epi.seeding = seq(40 * 7, length.out = 21) # make sure at least 1 introduction during these days

# load parameters
parm.bounds = read.csv(paste0(dir_code, 'parm.bounds.csv')) %>% data.table()
age.ranges = parm.bounds %>% filter(!is.na(age)) %>% .$age %>% unique %>% sort


# non age specific parms
parm_vec.t = c('frac.ini.inf0', # initial infection rate, 
               'seed0', # baseline seeding 
               'VE',
               'Tei', 'Tir', 
               # 'prob.tx.per.ctc0',  # Assign this using different scripts per specific model simulation settings
               'tx.sn.amp', 'vx.sn.amp', # seasonality
               'ave.HHctc.dur.per.day', 'ave.nonHHctc.dur.per.day', 
               'p.ve.v.i', 'p.inf2ve', 'p.ve.v.dis','p.inf.v.dis'  # imm protection related
               )
for(parm.t in parm_vec.t){
  tmp = fn_get.parmx(parm.t, type = parm.type.t) %>% dplyr::select(-parm, -age) %>% as.numeric()
  assign(parm.t, tmp)
}


# age-group specific parms
parm_vec.t = c('ave.num.nonHHctc.byage','sd.num.nonHHctc.byage','ave.num.HHctc.byage','sd.num.HHctc.byage', # contact rates
               'ave.p.redn.mask.byage','sd.p.redn.mask.byage', # mask effectiveness
               'severity.by.age', 'suscept.by.age', 'adj.ctc.redn.mask', 'adj.ctc.redn.nomask'
               )

for(parm.t in parm_vec.t){
  tmp = fn_get.parmx(parm.t, type = parm.type.t) %>% dplyr::select(-parm) %>% setnames(.,'default','value')
  assign(parm.t, tmp)
}

if(diff.ctc.redn.bymask){
  # use diff redn rate
  adj.ctc.redn.bymask = data.table(masking = 0:1, value = c(adj.ctc.redn.nomask$value, adj.ctc.redn.mask$value))
  
} else {
  adj.ctc.redn.bymask = data.table(masking = 0:1, value = 1)
  
}

# convert 1-year age group to coarser age groups
# if not specific for each age, convert it to that first
if(nrow(ave.num.nonHHctc.byage) > nrow(age.grps)){ # by age group
  parm_vec.t = c('ave.num.nonHHctc.byage','sd.num.nonHHctc.byage','ave.num.HHctc.byage','sd.num.HHctc.byage', # contact rates
                 'ave.p.redn.mask.byage','sd.p.redn.mask.byage' # mask effectiveness
  )
  for(parm.t in parm_vec.t){
    parm.grp.t = get(parm.t)
    if(nrow(parm.grp.t) == nrow(age.grps))
      next
    tmp = fn_map.age2coarser(parm.grp.t, 
                     age.grps)
    assign(parm.t, tmp)
  }
} 


# load VE estimates
# for tuning the vaccine-induced immunity against infection
# function form: VEadj.wt ~ 1 - 1 / (1+exp(-k * (day - tm.imm/2)))
est_parmVimmLoss = read.csv(paste0(dir_data, 'est_parmVimmLoss.csv')) %>% data.table()
# wt: for wildtype, alpha
variant.t = 'omicron'  # get set it to other variant per assumption
tm.imm.t0 = est_parmVimmLoss[variant == variant.t]$tm.imm; # during of vaccine-induced protection against SYMPTOMATIC infection
# use shorter time for any infection
tm.imm.t = tm.imm.t0 # / 2 # 180 # 
p.imm.wane.max.t =  est_parmVimmLoss[variant == variant.t]$p.imm.wane.max; # maximal level of immunity loss (=1 or lower)
k.t =  est_parmVimmLoss[variant == variant.t]$k; 

# based on these estimates, assume p.inf2ve times duration of protection from prior infection
# based on these estimates, assume p.ve.v.dis times duration of protection against severe disease from vaccination
# based on these estimates, assume p.inf.v.dis times duration of protection against severe disease from vaccination
# b/c uncertainty in these assumptions, put the ratios in a file to allow testing of multiple settings
Trs = tm.imm.t0 * p.inf2ve # immunity period
Tvs = tm.imm.t0 * p.ve.v.i  # duration of vaccine protection against infection
Tvd = tm.imm.t0 * p.ve.v.dis  # duration of vaccine protection against severe disease
Tid = tm.imm.t0 * p.inf.v.dis  # duration of protection against severe disease from prior infection

print(paste('tm.imm prior inf v new infection:', round(Trs,0)))
print(paste('tm.imm ve v infection:', round(Tvs,0)))
print(paste('tm.imm ve v severe dis:', round(Tvd,0)))
print(paste('tm.imm prior inf v severe:', round(Tid,0)))

parm.immLossInf = list(Dur = Trs, # mean immunity period
                       p.imm.wane.max = p.imm.wane.max.t, # parm for the logistic function, maximal level of immune waving within the mean imm period
                       k = k.t,  # parm for the logistic function 
                       p.imm.wane.midpt = .5)
# compute prob table 
prob.table = fn_calProbTable(parm.immLossInf, reverse.prob = F)
parm.immLossInf$prob.table = prob.table

parm.immLossVax = list(Dur = Tvs, # mean immunity period
                       p.imm.wane.max = p.imm.wane.max.t, # parm for the logistic function, maximal level of immune waving within the mean imm period
                       k = k.t,  # parm for the logistic function 
                       p.imm.wane.midpt = .5)
# compute prob table 
prob.table = fn_calProbTable(parm.immLossVax, reverse.prob = F)
parm.immLossVax$prob.table = prob.table

# parms for immune waning against severity, for both prior infection and vaccination
parm.immLossVseverity = list(Dur.vax = Tvd, # mean duration of vaccine protection
                             vax.p.max = p.imm.wane.max.t, # parm for the logistic function, maximal level of immune waving within the mean imm period
                             vax.k = k.t,  # parm for the logistic function 
                             vax.p.midpt = .5, 
                             Dur.imm = Tid, # mean duration of protection from prior infection
                             imm.p.max = p.imm.wane.max.t, # parm for the logistic function, maximal level of immune waving within the mean imm period
                             imm.k = k.t,  # parm for the logistic function 
                             imm.p.midpt = .5)

# parms for the fn_getVaxProbLogistic function
parm.getVaxProb = list(t2vax = 365, # mean duration of subsequent vaccinations
                       p.max = 1, # parm for the logistic function, maximal level of immune waving within the mean imm period
                       k = .05,  # parm for the logistic function 
                       p.midpt = 1,
                       vax.interval = 365, # period b/w two series of vaccination
                       tm.min.inf2vx = 14 # policy: cannot get vx w/in 90 days of recent infection
                       # but not all infected know they are infected due to e.g., asymptomatic infection
                       # so just limited to the recent 2 weeks (asked before vax)
)


# load vax coverage
load(paste0(dir_data, "vx_input_agegrps",n.age.grp.t,".RData"))
# fill in the vx rate for days w/o data, if needed
if(tail(tab_dates$date,1) > tail(vx.daily$date,1)){
  da.vx.daily = NULL
  for(ia in 1:nrow(age.grps)){
    vx.t = vx.daily %>% filter(age.start == age.grps$age.start[ia]) 
    if(nrow(vx.t) < 1)
      next
    vx.t = vx.t %>% ungroup %>% dplyr::select(date, value) %>% 
      full_join(x = ., y = tab_dates, by = 'date')
    
    vx_by_day = vx.t %>% # dplyr::filter(date > as.Date('2020/4/1')) %>% # exclude pre-pandemic
      reshape2::dcast(., DayOfYr ~ year, value.var = 'value') # %>% dplyr::select(-`2020`) # exclude year 2020
    d0 = ifelse(age.grps$age.start[ia]<5, '2022/3/1', '2021/9/1') # younger age groups had later vax program
    tmp = vx.t %>% dplyr::filter(date >= as.Date(d0)) %>% 
      reshape2::dcast(., DayOfYr ~ year, value.var = 'value') %>% 
      dplyr::select(-DayOfYr) %>% apply(., 1, mean, na.rm=T)
    if(length(tmp) < nrow(vx_by_day)){
      # no day 366
      tmp = c(tmp, tail(tmp,1))
    }
    vx_by_day$mn = tmp
    
    # replace later part of 2021 & year 2022 with the long term cycle
    # vx_by_day[271:366,'2021'] = NA_real_
    # vx_by_day[,'2022'] = NA_real_
    
    # fill in for na's
    j0 = which(names(vx_by_day)=='2022')
    for(j in j0:(ncol(vx_by_day)-1)){
      for(i in 1:nrow(vx_by_day)){
        vx_by_day[i,j] = case_when(is.na(vx_by_day[i,j]) ~ vx_by_day$mn[i],
                                   T ~ vx_by_day[i,j])
      }
    }
    
    
    
    vx_by_day$mn = NULL
    vx.t = vx_by_day %>% reshape2::melt(., id.vars = 'DayOfYr') %>%
      suppressMessages() %>%
      dplyr::mutate(date = as.Date(DayOfYr-1, origin = paste0(variable, '-1-1'))) %>%
      dplyr::mutate(year = format(date,'%Y')) %>%
      filter(variable == year) %>%
      dplyr::select(date, value) %>%
      filter(!is.na(value))
    
    vx.t = vx.t %>% right_join(x = ., y = tab_dates, by = 'date') %>%
      dplyr::arrange(., date)
    
    da.vx.daily = rbind(da.vx.daily,
                        data.table(vx.t, age.start = age.grps$age.start[ia],
                                   age.end = age.grps$age.end[ia]))
    rm(vx.t, vx_by_day)
    
  }
  
} else{
  da.vx.daily = vx.daily
}
  
# ggplot(da.vx.daily, aes(x = date, y = value)) +
#   geom_line() +
#   facet_rep_wrap(~age.start) + theme_minimal()
setDT(da.vx.daily)
da.vx.weekly = NULL
# easier just aggregate from the daily estimates
da.vx.weekly = da.vx.daily %>%
  dplyr::mutate(year = MMWRweek(date)$MMWRyear, week = MMWRweek(date)$MMWRweek) %>%
  group_by(year, week, age.start, age.end) %>%
  summarise(value = sum(value, na.rm = T)) %>%
  dplyr::mutate(date = MMWRweek2Date(year, week, 1) + 6) %>% # shift to end of week
  ungroup %>% dplyr::select(-year, week) %>% 
  right_join(x = ., y = tab_dates %>% filter(DayOfWk == 7), by = 'date') %>%
  dplyr::arrange(., age.start, date)

# adjust for ve
da.vx.daily = da.vx.daily %>%
  dplyr::mutate(value = (value * VE / 100)) # not these are %
da.vx.weekly = da.vx.weekly %>%
  dplyr::mutate(value = (value * VE / 100))

tt.vax.start = tab_dates %>% filter (date == (da.vx.daily$date %>% min)) %>% .$tt

# save a copy
save(da.vx.daily, file = paste0(dir_data, 'da.vx.daily.RData'))

rm(vx.daily, vx.daily.1st, vx.daily.booster1, vx.daily.booster1plus,
   vx.weekly, vx.weekly.1st, vx.weekly.booster1, vx.weekly.booster1plus)

# load contact pattern
load(paste0(dir_data, 'contact_input_agegrps',n.age.grp.t,'_maxcut0.95.RData'))  # set upper bound at the 95% quintile
da.proc.t = 'raw' # based on raw data, assume that those w/o specific type of contact was b/c it was not recorded, rather than having 0 contact

if(!exists('cont_freq.t'))
  cont_freq.t = 'all'; # including all contacts, regardless of freq

if(!exists('wt_nonphys.t'))
  wt_nonphys.t = .2; # assume non-phyical contact has 20% of the chance resulting infection

mea.t = 'cont.dur.wt' # use weighted contact duraction based on the assumed wt_nonphys
method.t = 'mle' # fitted using MLE
ctc.fit = res.cnt.fit %>% filter(da.proc == da.proc.t & method == method.t & cont_freq == cont_freq.t &
                                     wt_nonphys == wt_nonphys.t & measure == mea.t & 
                                   is.na(hh_size) # in case hh_size specfic model is also run
                                 )

rm(res.cnt.fit, res.hhsz.fit, res.stat)

# load mobility data to represent reduction in contact rate
da.mob = read.csv(paste0(dir_data, 'da_mobility_nyc_daily.csv')) %>% # data.table() %>%
  dplyr::mutate(date = date %>% as.Date())
# dates w/o data, fill in based on historical pattern
da.mob = da.mob %>% right_join(x = ., y = tab_dates, by = 'date')
da.mob.t = da.mob %>% dplyr::select(date, year, DayOfYr, mob.bus) %>% setnames(., 'mob.bus','mob')
  # dplyr::mutate(mob = stats::filter(mob.bus, filter=rep(1/3,3)))  # use moving average instead, as some days had very irregular mobility
# incomplete data, use the latest for the missing week
# but if the latest week is the summer, due to hot weather mob might be low
# use the same day of year instead
mob_by_day = da.mob.t %>% # dplyr::filter(date > as.Date('2020/4/1')) %>% # exclude pre-pandemic
  reshape2::dcast(., DayOfYr ~ year, value.var = 'mob') # %>% dplyr::select(-`2020`) # exclude year 2020

if(any(is.na(da.mob.t$mob))){
  mob_by_day$mx = da.mob.t %>% dplyr::filter(date > as.Date('2020/4/1')) %>% 
    reshape2::dcast(., DayOfYr ~ year, value.var = 'mob') %>% 
    dplyr::select(-DayOfYr) %>% apply(., 1, max, na.rm=T)
  
  
  # set years after July 2021 to 1 
  j0 = which(names(mob_by_day)=='2021')
  for(j in j0:(ncol(mob_by_day)-1)){
    if(j == j0){
      mob_by_day[181:nrow(mob_by_day),j] = 1
    } else {
      mob_by_day[,j] = 1
    }
  }
  # to avoid a large jump, gradually increase in 2021
  # j21 = which(names(mob_by_day)=='2021')
  # mob_by_day[1:90,j0] = seq(max(tail(mob_by_day$`2020`,3)),1, length.out = 90)
  mob_by_day[180+(1:90),j0] = seq(mean(tail(mob_by_day$`2021`[1:180],3)),1, length.out = 90)
  
  mob_by_day$mx = NULL
  
  da.mob.t = mob_by_day %>% reshape2::melt(., id.vars = 'DayOfYr') %>%
    suppressMessages() %>%
    dplyr::mutate(date = as.Date(DayOfYr-1, origin = paste0(variable, '-1-1'))) %>%
    dplyr::mutate(year = format(date,'%Y')) %>%
    filter(variable == year) %>%
    dplyr::select(date, value) %>%
    filter(!is.na(value))
} else {
  da.mob.t = da.mob.t %>% setnames(., 'mob', 'value') %>%
    dplyr::select(date, value)
}



da.p.ctc.redn = da.mob.t %>% right_join(x = ., y = tab_dates, by = 'date') %>%
  dplyr::arrange(., date)
rm(da.mob, da.mob.t, mob_by_day)

# save a copy
write.csv(da.p.ctc.redn, paste0(dir_data, 'da.p.ctc.redn.csv'), row.names = F)

if(MASKATHOME){
  print('WARNING: MASKATHOME SET TO TRUE')
}


