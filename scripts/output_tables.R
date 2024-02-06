# Model code used in Yang & Shaman "Reconciling the efficacy and effectiveness of masking on epidemic outcomes"
# author: Wan Yang
# date: March 2023
# note: this code is for research use only; may not be entirely optimized nor neatly organized. 


## OUPPUT TABLES

# TO USE - FIRST SET THE WORKING DIRECTORY TO SOURCE FILE LOCATION
dir_data = '../data/'
dir_code = '../scripts/'
dir_res = '../results/'
dir_table = '../tables/'

#######################################################
# load the results
# main analysis
load("../results/res_main_diff.bysn.RData")  # population-level effect

### use the script load_parms.R to get the param settings

num_pop = 100000 # population size
num_ens = 1 # number of ensemble runs
num.yr.sim = 5 # number of years to run the model

n.age.grp.t = 8
date.start = as.Date('2020/3/1')

agg.by.cohort = F # whether to compute and return tallies by birth cohort (no need if run for a short time)
agg.by.sn = T # whether to compute and return tallies by season  (no need if run for a short time)
wt_nonphys.main = .2 # main setting for wt of nonphysical contact

diff.ctc.redn.bymask = F # whether assume diff level of redn in contact rate by masking status
diff.ctc.bymask = F # whether to assume diff (e.g. lower) contact rate for those prefer masking
source(paste0(dir_code, 'get_main.runs.settings.R'))

inclVax = inclVax.t # whether to include vaccination
UseContactDuration = UseContactDuration.t #  whether to use contact duration directly or use number of contacts
UseVxData = UseVxData.t # whether to use vax data for the vaccination module
UseAgeSpecSuscept = UseAgeSpecSuscept.t # whether to assume diff susceptibility by age. set it to FALSE to reduce uncertainty
# diff.ctc.redn.bymask = F # whether assume diff level of redn in contact rate by masking status
# diff.ctc.bymask = F # whether to assume diff (e.g. lower) contact rate for those prefer masking
MASKATHOME = MASKATHOME.t # set this to TRUE only when testing a very specific scenario that ppl also mask at home


if(inclVax){
  status_vec = c('S','E','I','R','V')  # all possible status here
} else {
  status_vec = c('S','E','I','R')  # all possible status here
}

source(paste0(dir_code,'loadPackages.R'))
source(paste0(dir_code, 'Fns_util_AMBmasking.R'))
source(paste0(dir_code, 'Fn_ABMmasking.R'))
options(dplyr.summarise.inform = FALSE)
options(datatable.verbose = F)


# read in age group setting
age.grps = read.csv(paste0(dir_code,'age.grps',n.age.grp.t,'.csv')) %>% data.table()
str_age_recode = age.grps[,c('code','value')] %>% apply(1, paste, collapse = ' ~ ') %>% paste(collapse = ', ')

birth.rate = 1/age.grps[nrow(age.grps)]$age.end/365

mask.E.cut = c(.5, .9) # cut point for stratification of mask effectiveness
method.calHHrisk = 'like4like' # set method to calculate HH risk (like4like is probably more accurate than adjHHmask)

tm1 = Sys.time()

# load parameters for the simulation
source(paste0(dir_code, 'load_parms.R'))

################################################
# Table S1
tab1a = data.table(`Initital infection rate (% population)` = paste0(frac.ini.inf0 * 100,'%'), # initial infection rate, 
                   `Daily seeding rate` = seed0 * num_pop, # baseline seeding 
                   `Vacccination Effectiveness (VE)` = paste0(VE*100,'%'),
                   `Latency period` = paste(Tei,'days'), 
                   `Infectious period` = paste(Tir,'days'), 
                   `Immune duration of prior infection against infection` = paste(round(Trs,0),'days'),
                   `Immune duration of vaccination against infection` = paste(round(Tvs,0),'days'),
                   `Immune duration of prior infection against severe disease` = paste(round(Tid,0),'days'),
                   `Immune duration of vaccination against severe disease` = paste(round(Tvd,0),'days')) %>% t

tab1a = data.table(Parameter = c('Initial infection rate (% population)',
                                 'Daily seeding rate',
                                 'Vaccination Effectiveness (VE)',
                                 'Latency period',
                                 'Infectious period',
                                 'Immune duration of prior infection against infection',
                                 'Immune duration of prior infection against severe disease',
                                 'Immune duration of vaccination against infection',
                                 'Immune duration of vaccination against severe disease'),
                   Age = 'NA',
                   Value = c(paste0(frac.ini.inf0 * 100,'%'),
                             seed0 * num_pop,
                             paste0(VE*100,'%'),
                             paste(Tei,'days'),
                             paste(Tir,'days'), 
                             paste(round(Trs,0),'days'),
                             paste(round(Tid,0),'days'),
                             paste(round(Tvs,0),'days'),
                             paste(round(Tvd,0),'days')
                   ))

# ave.p.redn.mask.byage = c(0, rep(.5, 4), rep(.5, 8), rep(.5, 5), rep(.5, 7), rep(.7, 20), rep(.7, 20), rep(.8, 36)) # c(0, .5, .5, .5, .7, .7, .8) # c(0, .5, .5, .5, .7, .7, .8, .8)
# sd.p.redn.mask.byage = rep(.25, 101)

mask.eff = data.table(Parameter = 'Mask efficacy (age specific)', Age = age.grps$label, Value = c('no masks', rep('50% (SD = 25%)',4),rep('70% (SD = 25%)',2), '80% (SD = 25%)'))
hosp.rate = data.table(Parameter = 'Hospitalization rate (age specific)', Age = c('<20', '20'), Value = c(.5/1e3, 'exp(-7.5 + 0.075Age)'))
tab1 = rbind(tab1a, mask.eff, hosp.rate)
write.csv(tab1, paste0(dir_table,'TableS1_parameters.csv'), row.names = F)


# Table S2 - estimates
sn.t = "2023-2024"; tm.t = 'cum2date'; mask.grp.t = 'all'; age.grp.t = 'all'; 
vars.t = c('totAR', 'totHospRate'); vars_labs.t = c('Infections','Hospitalizations')
cut.mask.start_vec.lab.t = c('Always', paste0('Prevalence >',c('0.01%','0.02%','0.05%','0.1%','0.2%','0.5%', '1%', '2%')))

tab2 = res_main_diff.bysn %>% filter(season == sn.t & timeframe == tm.t & mask.grp == mask.grp.t & age.grp == age.grp.t &
                                       variable %in% vars.t & metric == 'perc.redn' & !diff.ctc.bymask & !diff.ctc.redn.bymask) %>%
  rename(., cut.mask.start = cut.mask.start1, 'Fraction of population masking' = ave.frac.masking1) %>%
  mutate(cut.mask.start = factor(cut.mask.start, levels = cut.mask.start_vec, labels = cut.mask.start_vec.lab.t))
tab2$value = tab2 %>% dplyr::select(med, ci95lwr, ci95upr) %>% apply(., 1, fn_format, rounddigt=1)
tab2 = tab2 %>% dcast(., variable + `Fraction of population masking` ~ cut.mask.start, value.var = 'value') %>%
  mutate(variable = factor(variable, levels = vars.t, labels = vars_labs.t))
write.csv(tab2, paste0(dir_table,'TableS2_full_est.eff.csv'), row.names = F)

# subset
sn.t = "2023-2024"; tm.t = 'cum2date'; mask.grp.t = 'all'; age.grp.t = 'all'; 
vars.t = c('totAR', 'totHospRate'); vars_labs.t = c('Infections','Hospitalizations')

cut.mask.start_vec.lab.t = c('Always', paste0('Prevalence >',c('0.01%','0.02%','0.05%','0.1%','0.2%','0.5%', '1%', '2%')))
cut.mask_vec.t = cut.mask.start_vec.lab.t[c(1,2,3,5,6,8,9)]
tab2 = res_main_diff.bysn %>% filter(season == sn.t & timeframe == tm.t & mask.grp == mask.grp.t & age.grp == age.grp.t &
                                       variable %in% vars.t & metric == 'perc.redn' & !diff.ctc.bymask & !diff.ctc.redn.bymask) %>%
  rename(., cut.mask.start = cut.mask.start1, 'Fraction of population masking' = ave.frac.masking1) %>%
  mutate(cut.mask.start = factor(cut.mask.start, levels = cut.mask.start_vec, labels = cut.mask.start_vec.lab.t)) %>%
  filter(cut.mask.start %in% cut.mask_vec.t)
tab2$value = tab2 %>% dplyr::select(med, ci95lwr, ci95upr) %>% apply(., 1, fn_format, rounddigt=1)
tab2 = tab2 %>% dcast(., variable + `Fraction of population masking` ~ cut.mask.start, value.var = 'value') %>%
  mutate(variable = factor(variable, levels = vars.t, labels = vars_labs.t))
write.csv(tab2, paste0(dir_table,'TableS2_est.eff.csv'), row.names = F)

# Table S3 - ctc distributions
# load contact pattern
n.age.grp.t = 8
load(paste0(dir_data, 'contact_input_agegrps',n.age.grp.t,'_maxcut0.95.RData'))  # set upper bound at the 95% quintile
da.proc.t = 'raw' # based on raw data, assume that those w/o specific type of contact was b/c it was not recorded, rather than having 0 contact
mea.t = 'cont.dur.wt' # use weighted contact duraction based on the assumed wt_nonphys
method.t = 'mle' # fitted using MLE
ctc.fit = res.cnt.fit %>% filter(da.proc == da.proc.t & method == method.t & cont_freq == cont_freq.t &
                                   measure == mea.t & 
                                   is.na(hh_size) # in case hh_size specfic model is also run
)
ctc.fit = ctc.fit %>% 
  mutate(distribution = case_when(prob.dist == 'gamma' ~ paste0('Gamma (rate=', round(rate,digits=2),', shape=',round(shape,2),')'),
                                  prob.dist == 'lnorm' ~ paste0('Log-normal (meanlog=',round(meanlog,2),', sdlog=',round(sdlog,2),')'),
                                  prob.dist == 'weibull' ~ paste0('Weibull (scale=',round(scale,2), ', shape=',round(shape,2),')'),
                                  prob.dist == 'norm' ~ paste0('Normal (mean=',round(mean,2),', sd=',round(sd,2),')'),
                                  prob.dist == 'unif' ~ paste0('Uniform (min=',round(min, 2), ', max=',round(max,2),')'),
                                  prob.dist == 'exp'~ paste0('Exponential (rate=',round(rate,2),')')
  )) %>%
  mutate(age.grp = factor(age.grp, levels = age.grps$label))
tabS3 = ctc.fit %>% dcast(., wt_nonphys + age.grp ~ cnt_home, value.var = 'distribution') 
colnames(tabS3) = c('Relative infection risk:\nnonphysical vs physical contact','Age group (year)','Non-household contact','Household contact')
write.xlsx(tabS3, paste0(dir_table,'TableS3.xlsx'))
