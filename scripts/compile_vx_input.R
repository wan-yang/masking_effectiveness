# Model code used in Yang & Shaman "Reconciling the efficacy and effectiveness of masking on epidemic outcomes"
# author: Wan Yang
# date: March 2023
# note: this code is for research use only; may not be entirely optimized nor neatly organized. 

n.age.grp.t = 8

# TO USE - FIRST SET THE WORKING DIRECTORY TO SOURCE FILE LOCATION
dir_data = '../data/'
dir_code = '../scripts/'
dir_res = '../results/'
dir_plot = '../plot_vx/' # plot the distributions if you'd like

if(! file.exists(dir_plot))
  dir.create(dir_plot, recursive = T)


source(paste0(dir_code, 'loadPackages.R'))
  
age.grps = read.csv(paste0(dir_code, 'age.grps',n.age.grp.t ,'.csv')) %>% data.table()

url.t = 'https://raw.githubusercontent.com/nychealth/covid-vaccine-data/main/people/coverage-by-demo-allages.csv'
vx.tot = read_csv(url(url.t)) %>% data.table() %>% janitor::clean_names()

url.t = 'https://raw.githubusercontent.com/nychealth/covid-vaccine-data/main/people/trends-byage.csv'
vx.trend = read_csv(url(url.t)) %>% data.table() %>% janitor::clean_names()


tda.vx = vx.tot %>% filter(group == 'Age') %>% dplyr::select(subgroup, perc_fully)
tmp = tda.vx$subgroup %>% sub("'",'',.) %>% str_split(., '-') %>% unlist %>% matrix(., nrow = nrow(tda.vx), ncol = 2, byrow = T)
# last row: 85+
tmp[nrow(tmp),] = c(85, 100)
tda.vx$age.start = tmp[,1]
tda.vx$age.end = tmp[,2]

n.dose_vec = c('1st', 'booster1', 'booster1plus') # 1st = primary series; booster1 = 3 dose; booster1plus = subsequent boosters
v.1st = 'perc_fully'
v.booster1 = 'perc_additional'
v.booster1plus = 'perc_bivalent_additional'
for(n.dose in n.dose_vec){
  v.t = get(paste0('v.',n.dose))
  
  # tda.trend = vx.trend %>% dplyr::select(date, grep('perc_fully', names(vx.trend)))
  tda.trend = vx.trend %>% dplyr::select(date, grep(v.t, names(vx.trend)))
  # covert cumulative to daily or weekly and add the lag (?)
  # no lag here b/c we are counting the first 2 doses together
  tda.trend.daily = vx.trend %>% dplyr::select(grep(v.t, names(vx.trend))) 
  tda.trend.daily = rbind(tda.trend.daily[1,],
                          tda.trend.daily[-1,] - tda.trend.daily[-nrow(tda.trend.daily),])
  tda.trend.daily$date = tda.trend$date
  # tda.trend.daily2 = tda.trend.daily %>% melt(., id.vars = 'date') %>% filter(variable != 'city_perc_fully')
  tda.trend.daily2 = tda.trend.daily %>% melt(., id.vars = 'date') %>% filter(! grepl('city', variable))
  lab.t = unique(tda.trend.daily2$variable)
  lab.a0 = lab.t %>% gsub(v.t, '', .) %>% str_split(., '_') %>% lapply(., function(x){x[2]}) %>% unlist %>% gsub('up','',.)
  lab.a1 = lab.t %>% gsub(v.t, '', .)  %>% str_split(., '_') %>% lapply(., function(x){x[3]}) %>% unlist %>% replace_na(.,85)
  tda.age.grps = data.table(age.start = lab.a0 %>% as.numeric(),
                            age.end = lab.a1 %>% as.numeric())
  tda.trend.daily2 = tda.trend.daily2 %>% 
    mutate(age.start = factor(variable, levels = lab.t, labels = lab.a0),
           age.end = factor(variable, levels = lab.t, labels = lab.a1)) %>% 
    dplyr::select(-variable)
  # match with the age groups in the model
  # for each age group in the model, search in the data, which ones are relevant
  tda.trend.daily3 = NULL
  for(ia in 1:nrow(age.grps)){
    ia0.t = age.grps[ia]$age.start
    ia1.t = age.grps[ia]$age.end
    if(ia1.t <1)
      next
    
    grp.t = c()
    for(ja in 1:nrow(tda.age.grps)){
      ja0.t = tda.age.grps[ja]$age.start
      ja1.t = tda.age.grps[ja]$age.end
      
      if((ja0.t:ja1.t) %in% (ia0.t:ia1.t) %>% sum() > length((ja0.t:ja1.t))/2){
        # at least half the group in the dataset is in this age group
        grp.t = c(grp.t, ja)
      }
    }
    
    if(length(grp.t) < 1)
      next
    
    tda.age.grps[grp.t]
    
    # assemble all relevant groups
    tda.t = tda.trend.daily2 %>% filter(age.start %in% tda.age.grps[grp.t]$age.start)
    
    if(length(grp.t)==1){
      # if only 1 grp is relevant
      tda.t = tda.t %>% mutate(age.start = ia0.t, age.end = ia1.t) # re-assign age group
      
    } else{
      # multiple groups, take the wt'ed average
      wts.t = (tda.age.grps[grp.t]$age.end - tda.age.grps[grp.t]$age.start + 1) 
      wts.t = wts.t/sum(wts.t)
      tda.t = tda.t %>% 
        mutate(wt = factor(age.start, levels = tda.age.grps[grp.t]$age.start, labels = wts.t)) %>%
        mutate(wt = wt %>% as.character() %>% as.numeric())
      
      # tda.t = tda.t %>% dcast(., date ~ age.start, value.var = 'value') %>%
      #   setnames(., old = c('date', tda.age.grps[grp.t]$age.start),
      #            new = c('date',paste0('gr',1:length(grp.t))))
      tda.t = tda.t %>% group_by(date) %>% 
        summarise(value = weighted.mean(value, w = wt)) %>%
        mutate(age.start = ia0.t, age.end = ia1.t)
    }
    tda.trend.daily3 = rbind(tda.trend.daily3,
                             tda.t)
  }
  
  
  tmp = tda.trend.daily3$date %>% MMWRweek()
  
  tda.trend.daily3 = tda.trend.daily3 %>% 
    mutate(year = tmp$MMWRyear, week = tmp$MMWRweek) 
  tda.trend.weekly = tda.trend.daily3 %>%
    group_by(year, week, age.start, age.end) %>%
    summarise(value = sum(value)) %>% 
    mutate(date = MMWRweek2Date(year, week, 1))
  
  assign(paste0('vx.daily.',n.dose), tda.trend.daily3)
  assign(paste0('vx.weekly.',n.dose), tda.trend.weekly)
  
  rm(tda.trend.daily, tda.trend.daily2, tda.trend.daily3, tda.trend.weekly)
}


# plot
tda = rbind(data.table(dose = 'primary', vx.weekly.1st),
            data.table(dose = '3rd', vx.weekly.booster1),
            data.table(dose = '4th', vx.weekly.booster1plus))
tda = tda %>% # mutate(date = MMWRweek2Date(year, week, 1)) %>%
  mutate(dose = factor(dose, levels = c('primary', '3rd', '4th')))
tda$age.grp = tda[,c('age.start','age.end')] %>% apply(1, paste, collapse='-')

pdf(paste0(dir_plot,'fig_vx.by.age.group.pdf'), width = 8, height = 12)
ggplot(tda, aes(x = date, y = value, color = dose)) +
  geom_line() +
  facet_rep_wrap(~age.grp, repeat.tick.labels = T, ncol=2) +
  theme_minimal()
dev.off()

# result: almost no overlap of the different doses
# could combine to simplify
vx.daily = rbind(vx.daily.1st, vx.daily.booster1, vx.daily.booster1plus) %>%
  group_by(date, year, week, age.start, age.end) %>%
  summarise(value = sum(value))

vx.weekly = rbind(vx.weekly.1st, vx.weekly.booster1, vx.weekly.booster1plus) %>%
  group_by(date, year, week, age.start, age.end) %>%
  summarise(value = sum(value))


tda = rbind(data.table(dose = 'primary', vx.weekly.1st),
            data.table(dose = '3rd', vx.weekly.booster1),
            data.table(dose = '4th', vx.weekly.booster1plus),
            data.table(dose = 'combined', vx.weekly))
tda = tda %>% # mutate(date = MMWRweek2Date(year, week, 1)) %>%
  mutate(dose = factor(dose, levels = c('primary', '3rd', '4th','combined')))
tda$age.grp = tda[,c('age.start','age.end')] %>% apply(1, paste, collapse='-')

pdf(paste0(dir_plot,'fig2_vx.by.age.group.pdf'), width = 8, height = 12)
ggplot(tda, aes(x = date, y = value, color = dose)) +
  geom_line() +
  facet_rep_wrap(~age.grp, repeat.tick.labels = T, ncol=2) +
  theme_minimal()
dev.off()

vars = ls()
vars = vars[grepl('vx.daily|vx.weekly', vars, perl = T)]
vars = vars %>% paste(., collapse = ', ')
eval(parse(text = paste0('save(', vars, ', file="', dir_data,'vx_input_agegrps',n.age.grp.t ,'.RData")')))


