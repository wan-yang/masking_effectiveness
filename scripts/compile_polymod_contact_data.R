# Model code used in Yang & Shaman "Reconciling the efficacy and effectiveness of masking on epidemic outcomes"
# author: Wan Yang
# date: March 2023
# note: this code is for research use only; may not be entirely optimized nor neatly organized. 

# code used to compile POLYMOD contact data for parameterization

# POLYMOD social contact data, from
# https://zenodo.org/record/3874557#.ZAiucy3MyF0


# TO USE - FIRST SET THE WORKING DIRECTORY TO SOURCE FILE LOCATION
dir_data = '../data/'
dir_code = '../scripts/'
dir_plot = '../plot_polymod/' # plot the distributions if you'd like

if(! file.exists(dir_plot))
  dir.create(dir_plot, recursive = T)

source(paste0(dir_code, 'loadPackages.R'))
library(fitdistrplus)
library(logspline)

if(F){
  # creating the age groups 
  age.grps = data.table(age.start = c(0, 1, 5, 15, 25, 45, 65),
                        age.end = c(0, 4, 14, 24, 44, 64, 80))
  write.csv(age.grps, paste0(dir_code, 'age.grps7.csv'), row.names = F)
  
  age.grps = data.table(age.start = c(0, 1, 5, 13, 18, 25, 45, 65),
                        age.end = c(0, 4, 12, 17, 24, 44, 64, 80))
  write.csv(age.grps, paste0(dir_code, 'age.grps8.csv'), row.names = F)
}

n.age.grp.t = 8 # number of age groups

max.cut = .95 # upper quantile for the contact data; tried .975, looked very high/extreme

# standardize number of contact
ave.HHctc.dur.per.day = 12/24 # average contact duration among HH members per day
ave.nonHHctc.dur.per.day = .5/24

age.grps = read.csv(paste0(dir_code,'age.grps',n.age.grp.t,'.csv'))
# put it in code
age.grps$value = age.grps[,c('age.start','age.end')] %>% apply(1, paste, collapse = '-') %>% sQuote(q='') %>% 
  gsub(pattern = '0-0', x = ., replacement = '<1') %>% 
  gsub(pattern = '65-80', x = ., replacement = '65+')
age.grps$label = age.grps[,c('age.start','age.end')] %>% apply(1, paste, collapse = '-') %>% 
  gsub(pattern = '0-0', x = ., replacement = '<1') %>% 
  gsub(pattern = '65-80', x = ., replacement = '65+')
age.grps$lower = paste("`part_age` >=", age.grps$age.start)
age.grps$upper = paste("`part_age` <=", age.grps$age.end)
age.grps$code = age.grps[, c('lower', 'upper')] %>% apply(1, paste, collapse = ' & ')
str_age_recode = age.grps[,c('code','value')] %>% apply(1, paste, collapse = ' ~ ') %>% paste(collapse = ', ')

# THESE ARE DATA FROM THE POLYMOD STUDY
code.part.common = read_xls(paste0(dir_data, '2008_Mossong_POLYMOD_dictionary.xls'), sheet = 1)
code.part.extra  = read_xls(paste0(dir_data, '2008_Mossong_POLYMOD_dictionary.xls'), sheet = 2)
code.cont.common = read_xls(paste0(dir_data, '2008_Mossong_POLYMOD_dictionary.xls'), sheet = 3)
code.hh.common = read_xls(paste0(dir_data, '2008_Mossong_POLYMOD_dictionary.xls'), sheet = 4)
code.hh.extra = read_xls(paste0(dir_data, '2008_Mossong_POLYMOD_dictionary.xls'), sheet = 5)

sday = read.csv(paste0(dir_data, '2008_Mossong_POLYMOD_sday.csv'))  # day of week survey done
# 0 = Sun; 6 = Sat

hh.common = read.csv(paste0(dir_data, '2008_Mossong_POLYMOD_hh_common.csv'))
hh.extra = read.csv(paste0(dir_data, '2008_Mossong_POLYMOD_hh_extra.csv'))

cont.common = read.csv(paste0(dir_data, '2008_Mossong_POLYMOD_contact_common.csv')) # individual contact info
part.common = read.csv(paste0(dir_data, '2008_Mossong_POLYMOD_participant_common.csv'))
part.extra = read.csv(paste0(dir_data, '2008_Mossong_POLYMOD_participant_extra.csv')) 

prob_bounds = c(.5, .25, .75, .025, .975)
prob_names = c('median', 'ci50lwr', 'ci50upr', 'ci95lwr', 'ci95upr')

idx = which(code.cont.common$`Data repository name` == 'duration_multi')
cont.duration = code.cont.common[idx:nrow(code.cont.common), ncol(code.cont.common)]
cont.duration$grp = 1:nrow(cont.duration)
cont.duration$lwr = c(1, 5, 15, 60, 60*4)
cont.duration$upr.nonhh = c(5, 15, 60, 60*4, 60 * 8) # upper bound for nonHH contacts
cont.duration$upr.hh = c(5, 15, 60, 60*4, 60 * 16)  # upper bound for HH contacts
cont.duration$mean.nonhh = cont.duration[, c('lwr', 'upr.nonhh')] %>% apply(., 1, mean)
cont.duration$mean.hh = cont.duration[, c('lwr', 'upr.hh')] %>% apply(., 1, mean)
# for outside contact
cont.duration$code.nonhh = paste("`cnt_home` == F & `duration_multi` == ", cont.duration$grp)
cont.duration$code.hh = paste("`cnt_home` == T & `duration_multi` == ", cont.duration$grp)
str_dur.nonhh_recode = cont.duration[,c('code.nonhh','mean.nonhh')] %>% apply(1, paste, collapse = ' ~ ') %>% paste(collapse = ', ')
str_dur.hh_recode = cont.duration[,c('code.hh','mean.hh')] %>% apply(1, paste, collapse = ' ~ ') %>% paste(collapse = ', ')
str_dur_recode = paste(str_dur.nonhh_recode, str_dur.hh_recode, sep = ', ')

# check missing data for cont duration
tmp = cont.common %>% filter(is.na(duration_multi) | is.na(phys_contact)) 
tmp = cont.common %>% filter(is.na(duration_multi)) # 1423 missing data
nrow(tmp)
hist(tmp$cnt_age_exact)

# note the weighted cont.dur account for number of contacts, duration of contact, and nature of contact

prob.dists = c('unif','exp', 'norm', 'lnorm', 'weibull', 'gamma') # list of prob distr to test the data on
dur0 = .1 # to add to 0
method_vec = c('mle') # different estimation method
qme.probs = seq(.1, .9, by = .1)
# for plotting
ncol.t = nrow(age.grps)
nrow.t = 2
cnt.max = nrow.t * ncol.t
# test and save diff weighting
cont_freq = c('daily', 'weekly_plus','monthly_plus', 'all')
cont_daily = 1
cont_weekly_plus = 1:2
cont_monthly_plus = 1:3
cont_all = c(1:5, NA)
wts_nonphy = c(.1, .2, .5, 1) 

res.cnt.fit = NULL
# fill in 0 contacts
tag.t = 'fill0'
for(method.t in method_vec){
  for(freq.t in cont_freq){
    print(freq.t)
    cont_freq.t = get(paste0('cont_',freq.t))
    for(iwt in 1:length(wts_nonphy)){
      wt.nonphy = wts_nonphy[iwt]
      
      wts = c(1, wt.nonphy) # relative weights for physical vs. non-physical contact
      wt.na = mean(wts)
      idx = which(code.cont.common$`Data repository name` == 'phys_contact')
      cont.phy = code.cont.common[idx+0:1, ncol(code.cont.common)] 
      cont.phy$grp = 1:nrow(cont.phy)
      cont.phy$value = wts
      cont.phy$code = paste("`phys_contact` == ", cont.phy$grp)
      str_wt_recode = cont.phy[,c('code','value')] %>% apply(1, paste, collapse = ' ~ ') %>% paste(collapse = ', ')
      
      # exclude missing data
      cont.common.t = cont.common %>% filter(!is.na(duration_multi)) %>%
        # filter based on contact freqency
        filter(frequency_multi %in% cont_freq.t) %>% 
        # fill in those with 0 contacts?
        full_join(., y = rbind(data.table(part_id = part.common$part_id, cnt_home = F),
                               data.table(part_id = part.common$part_id, cnt_home = T)),
                  by = c('part_id', 'cnt_home')) %>%
        replace_na(., list(cont.dur = 0))  %>%
        left_join(., part.common, by = 'part_id') %>%  
        left_join(., hh.common, by = 'hh_id') %>%
        # day of week
        left_join(., sday %>% dplyr::select(part_id, dayofweek), by = 'part_id') %>%
        # assign age group
        {eval(parse(text = paste('dplyr::mutate(., age.grp = case_when(', str_age_recode, ', T ~ NA_character_))', sep='')))} %>%
        mutate(age.grp = factor(age.grp, levels = age.grps$label)) %>%
        # set numerical duration of contact
        {eval(parse(text = paste('dplyr::mutate(., cont.dur = case_when(', str_dur_recode, ', T ~ NA_real_))', sep='')))} %>%
        # weight physical vs non-physical contact
        {eval(parse(text = paste('dplyr::mutate(., dur.wt = case_when(', str_wt_recode, ', T ~ wt.na))', sep='')))} %>%
        # convert to hour
        mutate(cont.dur = cont.dur / 60) %>%
        # weighted duration
        mutate(wt.cont.dur = cont.dur * dur.wt)
      
      # da.cont.all = cont.common.t %>% 
      #   # sum duration across all contacts
      #   group_by(part_id) %>% # per participant
      #   summarise(age.grp = first(age.grp),
      #             hh_size = first(hh_size), # household size including participants
      #             n.cont = n(),
      #             cont.dur = sum(cont.dur, na.rm = T), 
      #             wt.cont.dur = sum(wt.cont.dur, na.rm = T)) 
      
      da.cont = cont.common.t %>% # filter(cnt_home ==T) %>%  
        # sum duration across all contacts
        group_by(cnt_home, part_id) %>%
        summarise(age.grp = first(age.grp),
                  hh_size = first(hh_size),
                  dayofweek = first(dayofweek),
                  cont.dur = sum(cont.dur, na.rm = T), 
                  wt.cont.dur = sum(wt.cont.dur, na.rm = T)) %>%
        mutate(cont.dur.per.hhctc = case_when(cnt_home == F ~ NA_real_,
                                              cnt_home == T ~ cont.dur / pmax(1, (hh_size - 1))), # in case div 0
               wt.cont.dur.per.hhctc = case_when(cnt_home == F ~ NA_real_,
                                                 cnt_home == T ~ wt.cont.dur / pmax(1, (hh_size - 1)))) # exclude self
      any(is.na(da.cont$wt.cont.dur))
      any(is.na(da.cont$cont.dur))
      
      # fit the data to a certain distribution? 
      fit.cont = NULL
      
      # plot it
      cnt = 0; pp=list(); pp.wt = list()
      pdf(paste0(dir_plot, 'fig_polymod_',tag.t,'_maxcut',max.cut,'_',method.t,'_',freq.t,'_wt',wt.nonphy,'.pdf'), width = ncol.t*2, height = nrow.t*2.5)
      for(cnt_home.t in c(TRUE, FALSE)){
        for(age.grp.t in age.grps$label){
          tda = da.cont %>% filter(age.grp == age.grp.t & cnt_home == cnt_home.t) 
          # for some ppl, hhctc.dur == 0 but hh_size >>0 -> likely b/c hhctc was not recorded, rather than no hhctc at all
          # exclude those cases
          if(cnt_home.t){
            tda = tda  %>% filter(!(cont.dur == 0 & hh_size >= 2))
          } else{
            # for nonhh contact, exclude 0 likely due to weekends
            tda = tda  %>% filter(!(cont.dur == 0 & dayofweek %in% c(0,6,NA)))
          }
          
          
          # exclude the long tails
          dur.max = tda$cont.dur %>% quantile(., prob = max.cut)
          wt.dur.max = tda$wt.cont.dur %>% quantile(., prob = max.cut)
          
          
          # descdist(tmp$wt.cont.dur, discrete = F)
          fits = fits.wt = list(length(prob.dists))
          bic.fits = bic.fits.wt = numeric(length(prob.dists))
          for(ip in 1:length(prob.dists)){
            fits[[ip]]=fitdist(data = tda %>% filter(cont.dur<dur.max) %>% .$cont.dur+dur0, 
                               distr = prob.dists[ip], method = method.t) # , probs =qme.probs
            fits.wt[[ip]]=fitdist(data = tda %>% filter(wt.cont.dur<wt.dur.max) %>% .$wt.cont.dur+dur0, 
                                  distr = prob.dists[ip], method = method.t)
            
            bic.fits[ip] = fits[[ip]]$bic
            bic.fits.wt[ip] = fits.wt[[ip]]$bic
          }
          # get the best fit dist based on bic
          # all possible parameters
          prob.parm.names = lapply(fits, function(x){names(x$estimate)}) %>% unlist() %>% c() %>% unique %>% sort
          prob.parms0 = data.table(1)[,`:=`(c(prob.parm.names),NA)][,V1:=NULL][.0]
          ibest = which.min(bic.fits)
          ibest.wt = which.min(bic.fits.wt)
          prob.parms = rbind(prob.parms0, fits[[ibest]]$estimate %>% t %>% as.data.table(), fill=T)
          prob.parms.wt = rbind(prob.parms0, fits.wt[[ibest.wt]]$estimate %>% t %>% as.data.table(), fill=T)
          prob.parms = prob.parms %>% 
            mutate(method = method.t, cnt_home = cnt_home.t, age.grp = age.grp.t, hh_size = NA, n.smp = nrow(tda), max.cut = dur.max, # add upper bound
                   measure = 'cont.dur', prob.dist = prob.dists[ibest]) %>%
            setcolorder(., neworder = c('cnt_home', 'age.grp', 'measure', 'prob.dist'))
          prob.parms.wt = prob.parms.wt %>% 
            mutate(method = method.t, cnt_home = cnt_home.t, age.grp = age.grp.t, hh_size = NA, n.smp = nrow(tda), max.cut = wt.dur.max, # add upper bound
                   measure = 'cont.dur.wt', prob.dist = prob.dists[ibest.wt]) %>%
            setcolorder(., neworder = c('cnt_home', 'age.grp', 'measure', 'prob.dist'))
          
          fit.cont = rbind(fit.cont, 
                           prob.parms,
                           prob.parms.wt)
          
          # plotting
          cnt = cnt + 1
          pp[[cnt]] = denscomp(fits[[ibest]], plotstyle = "ggplot") + 
            labs(title = paste0('age: ', age.grp.t, '; ', ifelse(cnt_home.t, 'hhctc','nonhh'),'\nbest.dist: ', prob.dists[ibest]), x = 'cont dur (hours)') +
            theme_minimal() + theme(legend.position = 'none')
          pp.wt[[cnt]] = denscomp(fits.wt[[ibest.wt]], plotstyle = "ggplot") + 
            labs(title = paste0('age: ', age.grp.t, '; ',ifelse(cnt_home.t, 'hhctc','nonhh'),'\nbest.dist: ', prob.dists[ibest.wt]), x = 'weighted cont dur (hours)') +
            theme_minimal() + theme(legend.position = 'none')
          
          if(cnt == cnt.max){
            code.t = 'p = ggarrange('
            code.wt.t = 'p.wt = ggarrange('
            for(i in 1:cnt){
              if(i < cnt){
                code.t = paste0(code.t, 'pp[[',i,']],')
                code.wt.t = paste0(code.wt.t, 'pp.wt[[',i,']],')
              } else {
                code.t = paste0(code.t, 'pp[[',i,']], ', 'ncol = ncol.t, nrow= nrow.t)')
                code.wt.t = paste0(code.wt.t, 'pp.wt[[',i,']], ', 'ncol = ncol.t, nrow=nrow.t)')
              }
              
            }
            eval(parse(text = code.t))
            eval(parse(text = code.wt.t))
            print(p)
            print(p.wt)
            # reset counter 
            cnt = 0; pp = list(); pp.wt = list()
          }
        } # age
        
      } # cnt_home
      if(cnt >= 1){
        code.t = 'p = ggarrange('
        code.wt.t = 'p.wt = ggarrange('
        for(i in 1:cnt){
          if(i < cnt){
            code.t = paste0(code.t, 'pp[[',i,']],')
            code.wt.t = paste0(code.wt.t, 'pp.wt[[',i,']],')
          } else {
            code.t = paste0(code.t, 'pp[[',i,']], ', 'ncol = ncol.t, nrow=nrow.t)')
            code.wt.t = paste0(code.wt.t, 'pp.wt[[',i,']], ', 'ncol = ncol.t, nrow=nrow.t)')
          }
          
        }
        eval(parse(text = code.t))
        eval(parse(text = code.wt.t))
        print(p)
        print(p.wt)
        # reset counter 
        cnt = 0; pp = list(); pp.wt = list()
      }
      dev.off()
      
      # for hh cont further stratify by hh size
      hh_size.max = da.cont$hh_size %>% max(., na.rm = T)
      if(wt.nonphy == min(wts_nonphy)){
        # par(mfrow = c(4, 3), mar = c(2, 2, .5, .5), mgp = c(1, .3, 0), tck=-.02, cex = .8)
        cnt = 0; pp=list()
        # plot it
        pdf(paste0(dir_plot, 'fig_polymod_', tag.t,'_maxcut',max.cut,'_',method.t,'_',freq.t,'_hhcnt_byhhsz_byage.pdf'), width = ncol.t * 2.5, height = nrow.t * 2.5)
        for(age.grp.t in age.grps$label){
          for(hh_size.t in 2:hh_size.max){ 
            # note: hh_size includes self
            tda = da.cont %>% filter(age.grp == age.grp.t & cnt_home == T & hh_size == hh_size.t) 
            # for some ppl, hhctc.dur == 0 but hh_size >>0 -> likely b/c hhctc was not recorded, rather than no hhctc at all
            # exclude those cases
            if(cnt_home.t){
              tda = tda  %>% filter(!(cont.dur == 0 & hh_size >= 2))
            }
            
            # exclude the long tails
            dur.max = tda$cont.dur %>% quantile(., prob = max.cut)
            wt.dur.max = tda$wt.cont.dur %>% quantile(., prob = max.cut)
            
            if(nrow(tda) < 50)
              next
            
            # descdist(tmp$wt.cont.dur, discrete = F)
            fits = fits.wt = list(length(prob.dists))
            bic.fits = bic.fits.wt = numeric(length(prob.dists))
            for(ip in 1:length(prob.dists)){
              fits[[ip]]=fitdist(data = tda %>% filter(cont.dur<dur.max) %>% .$cont.dur+dur0, 
                                 distr = prob.dists[ip], method = method.t) # , probs =qme.probs
              fits.wt[[ip]]=fitdist(data = tda %>% filter(wt.cont.dur<wt.dur.max) %>% .$wt.cont.dur+dur0, 
                                    distr = prob.dists[ip], method = method.t)
              
              bic.fits[ip] = fits[[ip]]$bic
              bic.fits.wt[ip] = fits.wt[[ip]]$bic
            }
            # get the best fit dist based on bic
            # all possible parameters
            prob.parm.names = lapply(fits, function(x){names(x$estimate)}) %>% unlist() %>% c() %>% unique %>% sort
            prob.parms0 = data.table(1)[,`:=`(c(prob.parm.names),NA)][,V1:=NULL][.0]
            ibest = which.min(bic.fits)
            ibest.wt = which.min(bic.fits.wt)
            prob.parms = rbind(prob.parms0, fits[[ibest]]$estimate %>% t %>% as.data.table(), fill=T)
            prob.parms.wt = rbind(prob.parms0, fits.wt[[ibest.wt]]$estimate %>% t %>% as.data.table(), fill=T)
            prob.parms = prob.parms %>% 
              mutate(method = method.t, cnt_home = cnt_home.t, age.grp = age.grp.t, hh_size = hh_size.t, n.smp = nrow(tda), max.cut = dur.max, # add upper bound
                     measure = 'cont.dur', prob.dist = prob.dists[ibest]) %>%
              setcolorder(., neworder = c('cnt_home', 'age.grp', 'measure', 'prob.dist'))
            prob.parms.wt = prob.parms.wt %>% 
              mutate(method = method.t, cnt_home = cnt_home.t, age.grp = age.grp.t, hh_size = hh_size.t, n.smp = nrow(tda), max.cut = wt.dur.max, # add upper bound
                     measure = 'cont.dur.wt', prob.dist = prob.dists[ibest.wt]) %>%
              setcolorder(., neworder = c('cnt_home', 'age.grp', 'measure', 'prob.dist'))
            
            fit.cont = rbind(fit.cont, 
                             prob.parms,
                             prob.parms.wt)
            # plot it
            # denscomp(fits[[ibest]], legendtext = prob.dists[ibest],
            #          main = paste0('hh_sz=',hh_size.t, '; age: ', age.grp.t), xlab = 'cont dur (hours)')
            
            cnt = cnt + 1
            pp[[cnt]] = denscomp(fits[[ibest]], plotstyle = "ggplot") + 
              labs(title = paste0('age: ', age.grp.t, '; hh_sz=',hh_size.t,'\nbest.dist: ', prob.dists[ibest]), x = 'cont dur (hours)') +
              theme_minimal() + theme(legend.position = 'none')
            pp.wt[[cnt]] = denscomp(fits.wt[[ibest.wt]], plotstyle = "ggplot") + 
              labs(title = paste0('age: ', age.grp.t, '; hh_sz=',hh_size.t,'\nbest.dist: ', prob.dists[ibest.wt]), x = 'weighted cont dur (hours)') +
              theme_minimal() + theme(legend.position = 'none')
            
            if(cnt == cnt.max){
              code.t = 'p = ggarrange('
              code.wt.t = 'p.wt = ggarrange('
              for(i in 1:cnt){
                if(i < cnt){
                  code.t = paste0(code.t, 'pp[[',i,']],')
                  code.wt.t = paste0(code.wt.t, 'pp.wt[[',i,']],')
                } else {
                  code.t = paste0(code.t, 'pp[[',i,']], ', 'ncol = ncol.t, nrow= nrow.t)')
                  code.wt.t = paste0(code.wt.t, 'pp.wt[[',i,']], ', 'ncol = ncol.t, nrow=nrow.t)')
                }
                
              }
              eval(parse(text = code.t))
              eval(parse(text = code.wt.t))
              print(p)
              print(p.wt)
              # reset counter 
              cnt = 0; pp = list(); pp.wt = list()
            }
            
          } # hh_size
        } # age
        # last batch
        if(cnt >= 1){
          code.t = 'p = ggarrange('
          code.wt.t = 'p.wt = ggarrange('
          for(i in 1:cnt){
            if(i < cnt){
              code.t = paste0(code.t, 'pp[[',i,']],')
              code.wt.t = paste0(code.wt.t, 'pp.wt[[',i,']],')
            } else {
              code.t = paste0(code.t, 'pp[[',i,']], ', 'ncol = ncol.t, nrow=nrow.t)')
              code.wt.t = paste0(code.wt.t, 'pp.wt[[',i,']], ', 'ncol = ncol.t, nrow=nrow.t)')
            }
            
          }
          eval(parse(text = code.t))
          eval(parse(text = code.wt.t))
          print(p)
          print(p.wt)
          # reset counter 
          cnt = 0; pp = list(); pp.wt = list()
        }
        dev.off()
      } # by hh_size
      
      fit.cont = fit.cont %>%
        arrange(., cnt_home, age.grp) %>% # hh_size, 
        mutate(da.proc = tag.t, cont_freq = freq.t, wt_phys = wts[1], wt_nonphys = wts[2], wt_na = wt.na) %>%
        setcolorder(., neworder = c('da.proc', 'method','cont_freq', 'wt_phys','wt_nonphys','wt_na'))
      
      
      res.cnt.fit = rbind(res.cnt.fit, fit.cont)
    
    }
  }
} # method


method.t = 'mle'
freq.t = 'all'
iwt = 1
# do not fill in 0
tag.t = 'raw'
for(method.t in method_vec){
  for(freq.t in cont_freq){
    print(freq.t)
    cont_freq.t = get(paste0('cont_',freq.t))
    for(iwt in 1:length(wts_nonphy)){
      wt.nonphy = wts_nonphy[iwt]
      
      wts = c(1, wt.nonphy) # relative weights for physical vs. non-physical contact
      wt.na = mean(wts)
      idx = which(code.cont.common$`Data repository name` == 'phys_contact')
      cont.phy = code.cont.common[idx+0:1, ncol(code.cont.common)] 
      cont.phy$grp = 1:nrow(cont.phy)
      cont.phy$value = wts
      cont.phy$code = paste("`phys_contact` == ", cont.phy$grp)
      str_wt_recode = cont.phy[,c('code','value')] %>% apply(1, paste, collapse = ' ~ ') %>% paste(collapse = ', ')
      
      # exclude missing data
      cont.common.t = cont.common %>% filter(!is.na(duration_multi)) %>%
        # filter based on contact freqency
        filter(frequency_multi %in% cont_freq.t) %>% 
        # do not fill in those with 0 contacts?
        # full_join(., y = rbind(data.table(part_id = part.common$part_id, cnt_home = F),
        #                        data.table(part_id = part.common$part_id, cnt_home = T)), 
        #           by = c('part_id', 'cnt_home')) %>% 
        # replace_na(., list(cont.dur = 0))  %>% 
        left_join(., part.common, by = 'part_id') %>%  
        left_join(., hh.common, by = 'hh_id') %>%
        # day of week
        left_join(., sday %>% dplyr::select(part_id, dayofweek), by = 'part_id') %>%
        # assign age group
        {eval(parse(text = paste('dplyr::mutate(., age.grp = case_when(', str_age_recode, ', T ~ NA_character_))', sep='')))} %>%
        mutate(age.grp = factor(age.grp, levels = age.grps$label)) %>%
        # set numerical duration of contact
        {eval(parse(text = paste('dplyr::mutate(., cont.dur = case_when(', str_dur_recode, ', T ~ NA_real_))', sep='')))} %>%
        # weight physical vs non-physical contact
        {eval(parse(text = paste('dplyr::mutate(., dur.wt = case_when(', str_wt_recode, ', T ~ wt.na))', sep='')))} %>%
        # convert to hour
        mutate(cont.dur = cont.dur / 60) %>%
        # weighted duration
        mutate(wt.cont.dur = cont.dur * dur.wt)
      
      # da.cont.all = cont.common.t %>% 
      #   # sum duration across all contacts
      #   group_by(part_id) %>% # per participant
      #   summarise(age.grp = first(age.grp),
      #             hh_size = first(hh_size), # household size including participants
      #             n.cont = n(),
      #             cont.dur = sum(cont.dur, na.rm = T), 
      #             wt.cont.dur = sum(wt.cont.dur, na.rm = T)) 
      
      da.cont = cont.common.t %>% # filter(cnt_home ==T) %>%  
        # sum duration across all contacts
        group_by(cnt_home, part_id) %>%
        summarise(age.grp = first(age.grp),
                  hh_size = first(hh_size),
                  dayofweek = first(dayofweek),
                  cont.dur = sum(cont.dur, na.rm = T), 
                  wt.cont.dur = sum(wt.cont.dur, na.rm = T)) %>%
        mutate(cont.dur.per.hhctc = case_when(cnt_home == F ~ NA_real_,
                                              cnt_home == T ~ cont.dur / pmax(1, (hh_size - 1))), # in case div 0
               wt.cont.dur.per.hhctc = case_when(cnt_home == F ~ NA_real_,
                                                 cnt_home == T ~ wt.cont.dur / pmax(1, (hh_size - 1)))) # exclude self
      any(is.na(da.cont$wt.cont.dur))
      any(is.na(da.cont$cont.dur))
      
      # fit the data to a certain distribution? 
      fit.cont = NULL
      
      # plot it
      cnt = 0; pp=list(); pp.wt = list()
      pdf(paste0(dir_plot, 'fig_polymod_',tag.t,'_maxcut',max.cut,'_',method.t,'_',freq.t,'_wt',wt.nonphy,'.pdf'), width = ncol.t*2, height = nrow.t*2.5)
      for(cnt_home.t in c(TRUE, FALSE)){
        for(age.grp.t in age.grps$label){
          tda = da.cont %>% filter(age.grp == age.grp.t & cnt_home == cnt_home.t) 
          # for some ppl, hhctc.dur == 0 but hh_size >>0 -> likely b/c hhctc was not recorded, rather than no hhctc at all
          # exclude those cases
          if(cnt_home.t){
            tda = tda  %>% filter(!(cont.dur == 0 & hh_size >= 2))
          } else{
            # for nonhh contact, exclude 0 likely due to weekends
            # tda = tda  %>% filter(!(cont.dur == 0 & dayofweek %in% c(0,6,NA)))
          }
          
          
          # exclude the long tails
          dur.max = tda$cont.dur %>% quantile(., prob = max.cut)
          wt.dur.max = tda$wt.cont.dur %>% quantile(., prob = max.cut)
          
          
          # descdist(tmp$wt.cont.dur, discrete = F)
          fits = fits.wt = list(length(prob.dists))
          bic.fits = bic.fits.wt = numeric(length(prob.dists))
          for(ip in 1:length(prob.dists)){
            fits[[ip]]=fitdist(data = tda %>% filter(cont.dur<dur.max) %>% .$cont.dur+dur0, 
                               distr = prob.dists[ip], method = method.t) # , probs =qme.probs
            fits.wt[[ip]]=fitdist(data = tda %>% filter(wt.cont.dur<wt.dur.max) %>% .$wt.cont.dur+dur0, 
                                  distr = prob.dists[ip], method = method.t)
            
            bic.fits[ip] = fits[[ip]]$bic
            bic.fits.wt[ip] = fits.wt[[ip]]$bic
          }
          # get the best fit dist based on bic
          # all possible parameters
          prob.parm.names = lapply(fits, function(x){names(x$estimate)}) %>% unlist() %>% c() %>% unique %>% sort
          prob.parms0 = data.table(1)[,`:=`(c(prob.parm.names),NA)][,V1:=NULL][.0]
          ibest = which.min(bic.fits)
          ibest.wt = which.min(bic.fits.wt)
          prob.parms = rbind(prob.parms0, fits[[ibest]]$estimate %>% t %>% as.data.table(), fill=T)
          prob.parms.wt = rbind(prob.parms0, fits.wt[[ibest.wt]]$estimate %>% t %>% as.data.table(), fill=T)
          prob.parms = prob.parms %>% 
            mutate(method = method.t, cnt_home = cnt_home.t, age.grp = age.grp.t, hh_size = NA, n.smp = nrow(tda), max.cut = dur.max, # add upper bound
                   measure = 'cont.dur', prob.dist = prob.dists[ibest]) %>%
            setcolorder(., neworder = c('cnt_home', 'age.grp', 'measure', 'prob.dist'))
          prob.parms.wt = prob.parms.wt %>% 
            mutate(method = method.t, cnt_home = cnt_home.t, age.grp = age.grp.t, hh_size = NA, n.smp = nrow(tda), max.cut = wt.dur.max, # add upper bound
                   measure = 'cont.dur.wt', prob.dist = prob.dists[ibest.wt]) %>%
            setcolorder(., neworder = c('cnt_home', 'age.grp', 'measure', 'prob.dist'))
          
          fit.cont = rbind(fit.cont, 
                           prob.parms,
                           prob.parms.wt)
          
          # plotting
          cnt = cnt + 1
          pp[[cnt]] = denscomp(fits[[ibest]], plotstyle = "ggplot") + 
            labs(title = paste0('age: ', age.grp.t, '; ', ifelse(cnt_home.t, 'hhctc','nonhh'),'\nbest.dist: ', prob.dists[ibest]), x = 'cont dur (hours)') +
            theme_minimal() + theme(legend.position = 'none')
          pp.wt[[cnt]] = denscomp(fits.wt[[ibest.wt]], plotstyle = "ggplot") + 
            labs(title = paste0('age: ', age.grp.t, '; ',ifelse(cnt_home.t, 'hhctc','nonhh'),'\nbest.dist: ', prob.dists[ibest.wt]), x = 'weighted cont dur (hours)') +
            theme_minimal() + theme(legend.position = 'none')
          
          if(cnt == cnt.max){
            code.t = 'p = ggarrange('
            code.wt.t = 'p.wt = ggarrange('
            for(i in 1:cnt){
              if(i < cnt){
                code.t = paste0(code.t, 'pp[[',i,']],')
                code.wt.t = paste0(code.wt.t, 'pp.wt[[',i,']],')
              } else {
                code.t = paste0(code.t, 'pp[[',i,']], ', 'ncol = ncol.t, nrow= nrow.t)')
                code.wt.t = paste0(code.wt.t, 'pp.wt[[',i,']], ', 'ncol = ncol.t, nrow=nrow.t)')
              }
              
            }
            eval(parse(text = code.t))
            eval(parse(text = code.wt.t))
            print(p)
            print(p.wt)
            # reset counter 
            cnt = 0; pp = list(); pp.wt = list()
          }
        } # age
        
      } # cnt_home
      if(cnt >= 1){
        code.t = 'p = ggarrange('
        code.wt.t = 'p.wt = ggarrange('
        for(i in 1:cnt){
          if(i < cnt){
            code.t = paste0(code.t, 'pp[[',i,']],')
            code.wt.t = paste0(code.wt.t, 'pp.wt[[',i,']],')
          } else {
            code.t = paste0(code.t, 'pp[[',i,']], ', 'ncol = ncol.t, nrow=nrow.t)')
            code.wt.t = paste0(code.wt.t, 'pp.wt[[',i,']], ', 'ncol = ncol.t, nrow=nrow.t)')
          }
          
        }
        eval(parse(text = code.t))
        eval(parse(text = code.wt.t))
        print(p)
        print(p.wt)
        # reset counter 
        cnt = 0; pp = list(); pp.wt = list()
      }
      dev.off()
      
      # for hh cont further stratify by hh size
      hh_size.max = da.cont$hh_size %>% max(., na.rm = T)
      if(wt.nonphy == min(wts_nonphy)){
        # par(mfrow = c(4, 3), mar = c(2, 2, .5, .5), mgp = c(1, .3, 0), tck=-.02, cex = .8)
        cnt = 0; pp=list()
        # plot it
        pdf(paste0(dir_plot, 'fig_polymod_', tag.t,'_maxcut',max.cut,'_',method.t,'_',freq.t,'_hhcnt_byhhsz_byage.pdf'), width = ncol.t * 2.5, height = nrow.t * 2.5)
        for(age.grp.t in age.grps$label){
          for(hh_size.t in 2:hh_size.max){ 
            # note: hh_size includes self
            tda = da.cont %>% filter(age.grp == age.grp.t & cnt_home == T & hh_size == hh_size.t) 
            # for some ppl, hhctc.dur == 0 but hh_size >>0 -> likely b/c hhctc was not recorded, rather than no hhctc at all
            # exclude those cases
            if(cnt_home.t){
              tda = tda  %>% filter(!(cont.dur == 0 & hh_size >= 2))
            }
            
            # exclude the long tails
            dur.max = tda$cont.dur %>% quantile(., prob = max.cut)
            wt.dur.max = tda$wt.cont.dur %>% quantile(., prob = max.cut)
            
            if(nrow(tda) < 50)
              next
            
            # descdist(tmp$wt.cont.dur, discrete = F)
            fits = fits.wt = list(length(prob.dists))
            bic.fits = bic.fits.wt = numeric(length(prob.dists))
            for(ip in 1:length(prob.dists)){
              fits[[ip]]=fitdist(data = tda %>% filter(cont.dur<dur.max) %>% .$cont.dur+dur0, 
                                 distr = prob.dists[ip], method = method.t) # , probs =qme.probs
              fits.wt[[ip]]=fitdist(data = tda %>% filter(wt.cont.dur<wt.dur.max) %>% .$wt.cont.dur+dur0, 
                                    distr = prob.dists[ip], method = method.t)
              
              bic.fits[ip] = fits[[ip]]$bic
              bic.fits.wt[ip] = fits.wt[[ip]]$bic
            }
            # get the best fit dist based on bic
            # all possible parameters
            prob.parm.names = lapply(fits, function(x){names(x$estimate)}) %>% unlist() %>% c() %>% unique %>% sort
            prob.parms0 = data.table(1)[,`:=`(c(prob.parm.names),NA)][,V1:=NULL][.0]
            ibest = which.min(bic.fits)
            ibest.wt = which.min(bic.fits.wt)
            prob.parms = rbind(prob.parms0, fits[[ibest]]$estimate %>% t %>% as.data.table(), fill=T)
            prob.parms.wt = rbind(prob.parms0, fits.wt[[ibest.wt]]$estimate %>% t %>% as.data.table(), fill=T)
            prob.parms = prob.parms %>% 
              mutate(method = method.t, cnt_home = cnt_home.t, age.grp = age.grp.t, hh_size = hh_size.t, n.smp = nrow(tda), max.cut = dur.max, # add upper bound
                     measure = 'cont.dur', prob.dist = prob.dists[ibest]) %>%
              setcolorder(., neworder = c('cnt_home', 'age.grp', 'measure', 'prob.dist'))
            prob.parms.wt = prob.parms.wt %>% 
              mutate(method = method.t, cnt_home = cnt_home.t, age.grp = age.grp.t, hh_size = hh_size.t, n.smp = nrow(tda), max.cut = wt.dur.max, # add upper bound
                     measure = 'cont.dur.wt', prob.dist = prob.dists[ibest.wt]) %>%
              setcolorder(., neworder = c('cnt_home', 'age.grp', 'measure', 'prob.dist'))
            
            fit.cont = rbind(fit.cont, 
                             prob.parms,
                             prob.parms.wt)
            # plot it
            # denscomp(fits[[ibest]], legendtext = prob.dists[ibest],
            #          main = paste0('hh_sz=',hh_size.t, '; age: ', age.grp.t), xlab = 'cont dur (hours)')
            
            cnt = cnt + 1
            pp[[cnt]] = denscomp(fits[[ibest]], plotstyle = "ggplot") + 
              labs(title = paste0('age: ', age.grp.t, '; hh_sz=',hh_size.t,'\nbest.dist: ', prob.dists[ibest]), x = 'cont dur (hours)') +
              theme_minimal() + theme(legend.position = 'none')
            pp.wt[[cnt]] = denscomp(fits.wt[[ibest.wt]], plotstyle = "ggplot") + 
              labs(title = paste0('age: ', age.grp.t, '; hh_sz=',hh_size.t,'\nbest.dist: ', prob.dists[ibest.wt]), x = 'weighted cont dur (hours)') +
              theme_minimal() + theme(legend.position = 'none')
            
            if(cnt == cnt.max){
              code.t = 'p = ggarrange('
              code.wt.t = 'p.wt = ggarrange('
              for(i in 1:cnt){
                if(i < cnt){
                  code.t = paste0(code.t, 'pp[[',i,']],')
                  code.wt.t = paste0(code.wt.t, 'pp.wt[[',i,']],')
                } else {
                  code.t = paste0(code.t, 'pp[[',i,']], ', 'ncol = ncol.t, nrow= nrow.t)')
                  code.wt.t = paste0(code.wt.t, 'pp.wt[[',i,']], ', 'ncol = ncol.t, nrow=nrow.t)')
                }
                
              }
              eval(parse(text = code.t))
              eval(parse(text = code.wt.t))
              print(p)
              print(p.wt)
              # reset counter 
              cnt = 0; pp = list(); pp.wt = list()
            }
            
          } # hh_size
        } # age
        # last batch
        if(cnt >= 1){
          code.t = 'p = ggarrange('
          code.wt.t = 'p.wt = ggarrange('
          for(i in 1:cnt){
            if(i < cnt){
              code.t = paste0(code.t, 'pp[[',i,']],')
              code.wt.t = paste0(code.wt.t, 'pp.wt[[',i,']],')
            } else {
              code.t = paste0(code.t, 'pp[[',i,']], ', 'ncol = ncol.t, nrow=nrow.t)')
              code.wt.t = paste0(code.wt.t, 'pp.wt[[',i,']], ', 'ncol = ncol.t, nrow=nrow.t)')
            }
            
          }
          eval(parse(text = code.t))
          eval(parse(text = code.wt.t))
          print(p)
          print(p.wt)
          # reset counter 
          cnt = 0; pp = list(); pp.wt = list()
        }
        dev.off()
      } # by hh_size
      fit.cont = fit.cont %>%
        arrange(., cnt_home, age.grp) %>% # hh_size, 
        mutate(da.proc = tag.t, cont_freq = freq.t, wt_phys = wts[1], wt_nonphys = wts[2], wt_na = wt.na) %>%
        setcolorder(., neworder = c('da.proc', 'method','cont_freq', 'wt_phys','wt_nonphys','wt_na'))
      
      
      res.cnt.fit = rbind(res.cnt.fit, fit.cont)
      
    }
  }
} # method


# fit hh size by age group
res.hhsz.fit = NULL
for(method.t in method_vec){
  # exclude missing data
  hh.common.t = part.common %>% 
    left_join(., hh.common, by = 'hh_id') %>%
    # assign age group
    {eval(parse(text = paste('dplyr::mutate(., age.grp = case_when(', str_age_recode, ', T ~ NA_character_))', sep='')))} 
  
  cnt = 0; pp = list()
  pdf(paste0(dir_plot,'fig_hhsz.byage_',method.t,'_maxcut',max.cut,'.pdf'), width = 10, height = 6)
  for(age.grp.t in age.grps$label){
    # note: hh_size includes self
    tda = hh.common.t %>% filter(age.grp == age.grp.t) 
    # exclude the long tails
    hh.max = tda$hh_size %>% quantile(., prob = max.cut)
    
    if(nrow(tda) < 50)
      next
    
    # descdist(tmp$wt.cont.dur, discrete = F)
    fits = list(length(prob.dists))
    # fits.wt = list(length(prob.dists))
    bic.fits = numeric(length(prob.dists))
    # bic.fits.wt = numeric(length(prob.dists))
    for(ip in 1:length(prob.dists)){
      fits[[ip]]=fitdist(tda %>% filter(hh_size<=hh.max) %>% .$hh_size, prob.dists[ip], method = method.t)
      # fits.wt[[ip]]=fitdist(tda$wt.cont.dur+ dur0, prob.dists[ip])
      
      bic.fits[ip] = fits[[ip]]$bic
      # bic.fits.wt[ip] = fits.wt[[ip]]$bic
    }
    # get the best fit dist based on bic
    # all possible parameters
    prob.parm.names = lapply(fits, function(x){names(x$estimate)}) %>% unlist() %>% c() %>% unique %>% sort
    prob.parms0 = data.table(1)[,`:=`(c(prob.parm.names),NA)][,V1:=NULL][.0]
    ibest = which.min(bic.fits)
    # ibest.wt = which.min(bic.fits.wt)
    prob.parms = rbind(prob.parms0, fits[[ibest]]$estimate %>% t %>% as.data.table(), fill=T)
    # prob.parms.wt = rbind(prob.parms0, fits.wt[[ibest.wt]]$estimate %>% t %>% as.data.table(), fill=T)
    prob.parms = prob.parms %>% 
      mutate(method = method.t, age.grp = age.grp.t, n.smp = nrow(tda), max.cut = hh.max, # add upper bound
             measure = 'hh_size', prob.dist = prob.dists[ibest]) %>%
      setcolorder(., neworder = c('age.grp', 'measure', 'prob.dist'))
    
    res.hhsz.fit = rbind(res.hhsz.fit, 
                     prob.parms)  # prob.parms.wt
    
    # plot it
    # denscomp(fits[[ibest]], legendtext = prob.dists[ibest],
    #          main = paste0('hh_sz=',hh_size.t, '; age: ', age.grp.t), xlab = 'cont dur (hours)')
    
    cnt = cnt + 1
    pp[[cnt]] = denscomp(fits[[ibest]], plotstyle = "ggplot") + 
      labs(title = paste0('age: ', age.grp.t, '\nbest.dist: ', prob.dists[ibest]), x = 'hh size') +
      theme_minimal() + theme(legend.position = 'none')
    
    if(cnt == nrow(age.grps)){
      code.t = 'p = ggarrange('
      for(i in 1:cnt){
        if(i < cnt){
          code.t = paste0(code.t, 'pp[[',i,']],')
        } else {
          code.t = paste0(code.t, 'pp[[',i,']], ', 'ncol = 4, nrow=2)')
        }
        
      }
      eval(parse(text = code.t))
      print(p)
      # reset counter 
      cnt = 0; pp = list()
    }
  } # age
  dev.off()
} # method


res.stat = NULL
for(freq.t in cont_freq){
  print(freq.t)
  cont_freq.t = get(paste0('cont_',freq.t))
  for(iwt in 1:length(wts_nonphy)){
    wt.nonphy = wts_nonphy[iwt]
    
    wts = c(1, wt.nonphy) # relative weights for physical vs. non-physical contact
    wt.na = mean(wts)
    idx = which(code.cont.common$`Data repository name` == 'phys_contact')
    cont.phy = code.cont.common[idx+0:1, ncol(code.cont.common)] 
    cont.phy$grp = 1:nrow(cont.phy)
    cont.phy$value = wts
    cont.phy$code = paste("`phys_contact` == ", cont.phy$grp)
    str_wt_recode = cont.phy[,c('code','value')] %>% apply(1, paste, collapse = ' ~ ') %>% paste(collapse = ', ')
    
    # exclude missing data
    cont.common.t = cont.common %>% filter(!is.na(duration_multi)) %>%
      # filter based on contact freqency
      filter(frequency_multi %in% cont_freq.t) %>% 
      # need to fill in those with 0 contacts 
      # full_join(., y = rbind(data.table(part_id = part.common$part_id, cnt_home = F),
      #                        data.table(part_id = part.common$part_id, cnt_home = T)), 
      #           by = c('part_id', 'cnt_home')) %>% 
      # replace_na(., list(cont.dur = 0))  %>% 
      left_join(., part.common, by = 'part_id') %>%  
      left_join(., hh.common, by = 'hh_id') %>%
      # assign age group
      {eval(parse(text = paste('dplyr::mutate(., age.grp = case_when(', str_age_recode, ', T ~ NA_character_))', sep='')))} %>%
      mutate(age.grp = factor(age.grp, levels = age.grps$label)) %>%
      # set numerical duration of contact
      {eval(parse(text = paste('dplyr::mutate(., cont.dur = case_when(', str_dur_recode, ', T ~ NA_real_))', sep='')))} %>%
      # weight physical vs non-physical contact
      {eval(parse(text = paste('dplyr::mutate(., dur.wt = case_when(', str_wt_recode, ', T ~ wt.na))', sep='')))} %>%
      # convert to hour
      mutate(cont.dur = cont.dur / 60) %>%
      # weighted duration
      mutate(wt.cont.dur = cont.dur * dur.wt)
    
    # da.cont.all = cont.common.t %>% 
    #   # sum duration across all contacts
    #   group_by(part_id) %>% # per participant
    #   summarise(age.grp = first(age.grp),
    #             hh_size = first(hh_size), # household size including participants
    #             n.cont = n(),
    #             cont.dur = sum(cont.dur, na.rm = T), 
    #             wt.cont.dur = sum(wt.cont.dur, na.rm = T)) 
    
    da.cont = cont.common.t %>% # filter(cnt_home ==T) %>%  
      # sum duration across all contacts
      group_by(cnt_home, part_id) %>%
      summarise(age.grp = first(age.grp),
                hh_size = first(hh_size),
                cont.dur = sum(cont.dur, na.rm = T), 
                wt.cont.dur = sum(wt.cont.dur, na.rm = T)) %>%
      mutate(cont.dur.per.hhctc = case_when(cnt_home == F ~ NA_real_,
                                            cnt_home == T ~ cont.dur / pmax(1, (hh_size - 1))), # in case div 0
             wt.cont.dur.per.hhctc = case_when(cnt_home == F ~ NA_real_,
                                               cnt_home == T ~ wt.cont.dur / pmax(1, (hh_size - 1)))) # exclude self
    any(is.na(da.cont$wt.cont.dur))
    any(is.na(da.cont$cont.dur))
    
    stat.cont = da.cont %>%
      group_by(cnt_home, age.grp) %>% # hh_size, 
      summarise(cont.dur.mn = mean(cont.dur, na.rm = T), 
                wt.cont.dur.mn = mean(wt.cont.dur, na.rm = T),
                cont.dur.per.hhctc.mn = mean(cont.dur.per.hhctc, na.rm = T),
                wt.cont.dur.per.hhctc.mn = mean(wt.cont.dur.per.hhctc, na.rm = T),
                cont.dur.median = median(cont.dur, na.rm = T), 
                wt.cont.dur.median = median(wt.cont.dur, na.rm = T),
                cont.dur.per.hhctc.median = median(cont.dur.per.hhctc, na.rm = T),
                wt.cont.dur.per.hhctc.median = median(wt.cont.dur.per.hhctc, na.rm = T),
                cont.dur.sd = sd(cont.dur, na.rm = T), 
                wt.cont.dur.sd = sd(wt.cont.dur, na.rm = T),
                cont.dur.per.hhctc.sd = sd(cont.dur.per.hhctc, na.rm = T),
                wt.cont.dur.per.hhctc.sd = sd(cont.dur.per.hhctc, na.rm = T)
      ) %>%
      filter(!is.na(cnt_home)) %>%
      filter(!is.na(age.grp)) %>%
      arrange(., cnt_home, age.grp) %>% # hh_size, 
      mutate(cont_freq = freq.t, wt_phys = wts[1], wt_nonphys = wts[2], wt_na = wt.na) %>%
      setcolorder(., neworder = c('cont_freq', 'wt_phys','wt_nonphys','wt_na'))
    
    
    # assign(paste0('sheet', iwt), stat.cont)
    res.stat = rbind(res.stat, stat.cont)
    # res.cnt.fit = rbind(res.cnt.fit, fit.cont)
  }
}


save(res.cnt.fit, res.hhsz.fit, res.stat, file = paste0(dir_data, 'contact_input_agegrps',n.age.grp.t,'_maxcut',max.cut,'.RData'))

