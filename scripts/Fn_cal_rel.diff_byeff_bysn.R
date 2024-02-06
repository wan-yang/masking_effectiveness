# Model code used in Yang & Shaman "Reconciling the efficacy and effectiveness of masking on epidemic outcomes"
# author: Wan Yang
# date: March 2023
# note: this code is for research use only; may not be entirely optimized nor neatly organized. 

# compute relative diff vs. baseline, by mask effectiveness group, by season

{
  fn_stat = function(x){
    round.git = 3
    return(list(mn = mean(x, na.rm = T) %>% round(., round.git),
                med = median(x, na.rm = T) %>% round(., round.git),
                sd = sd(x, na.rm = T) %>% round(., round.git),
                ci50lwr = quantile(x, .25, na.rm = T) %>% round(., round.git),
                ci50upr = quantile(x, .75, na.rm = T) %>% round(., round.git),
                ci95lwr = quantile(x, .025, na.rm = T) %>% round(., round.git),
                ci95upr = quantile(x, .975, na.rm = T) %>% round(., round.git)))
  }
  
}

tmp = try(load(paste0(dir_out, 'out_diff.byeff.bysn.',
                      ifelse(diff.ctc.bymask.t,'CTCbymask','CTCo'),'_',
                      ifelse(diff.ctc.redn.bymask.t,'MOBbymask','MOBo'),
                      'fm',fm.no,
                      '.RData')))

if(class(tmp) == 'try-error'){ # not there yet
  RES = NULL
  for(fm.t in fm.no){  # 2:length(frac.masking_vec)
    for(cm.t in 1:length(cut.mask.start_vec)){ # length(cut.mask.start_vec)
      
      print(paste('fm.t =', fm.t, 'cm.t =', cm.t))
      if(fm.t == 1 & cm.t == 1)
        next
      
      
      fsce = paste0(dir_res, 'res_',
                    ifelse(inclVax.t,'vx','novx'),'_',
                    ifelse(UseAgeSpecSuscept.t,'Sbyage','So'),'_',
                    ifelse(diff.ctc.bymask.t,'CTCbymask','CTCo'),'_',
                    ifelse(diff.ctc.redn.bymask.t,'MOBbymask','MOBo'),
                    '_fm', fm.t, '_cm', cm.t)  # , '_ens', ens_no, '.RData'
      
      
      res.cp = NULL
      for(ens.t in 1:num_ens){
        
        fname.t = paste0(fsce, '_ens', ens.t, '.RData')
        
        tmp = try(load(fname.t), silent = T)
        if(any(class(tmp) == 'try-error'))
          next
        
        # outs = c('res.daily','res.wkly','res.yrly','res.byyr.stat', 'res.bycohort.byage.stat','res.stat','ctc.stat','parm.list')
        outs = c('res.bysn.stat')
        for(oo in outs){
          eval(parse(text = paste0('out.t = res$',oo)))
          # compute the diff
          # res.daily: daily tallies
          # res.wkly: % of population in each dis states at the end of each week
          # res.yrly: yearly tallies
          # res.byyr.stat: both yearly, and cumulative stats
          # res.snly:sesason tallies
          # res.bysn.stat: both season, and cumulative stats
          # res.bycohort.byage.stat: stats by birth cohort
          # res.stat: cumulative, and across all ens
          # ctc.stat: summary stats for checking ctc rate by age group, and masking assignment
          # parm.list: save the stats seperately to reduce space
          stat.type.t = 'cum2date'
          mask.grp.t = 'all' # b/c the baseline has no mask group
          age.grp.t = '25-44'
          # compute for each group, vs the ref group, level of reduction
          ref.t = 'no mask'
          tda0 = out.t %>% filter(mask.grp == ref.t) %>% 
            reshape2::melt(., id.vars = c('season', 'ens', 'timeframe', 'mask.grp','age.grp')) %>% 
            mutate(mask.grp = NULL) %>% 
            rename(., base = value) 
          tda1 = out.t %>% filter(mask.grp != ref.t) %>% 
            reshape2::melt(., id.vars = c('season', 'ens', 'timeframe', 'mask.grp','age.grp')) %>%
            rename(., sce = value)
          tda = tda1 %>% full_join(x = ., 
                                   y = tda0, 
                                   by = c('season', 'ens', 'timeframe', 'age.grp', 'variable')) %>% 
            group_by(season, ens, timeframe, mask.grp, age.grp, variable) %>%
            summarise(perc.redn = case_when(base == 0 & sce == 0 ~ 0,
                                            base == 0 & sce >0 ~ ((base - sce)/(base+.1)) * 100, # NA_real_, 
                                            T ~ ((base - sce)/base) * 100) %>% round(., 2),
                      # ratio
                      ratio = case_when(base == 0 & sce == 0 ~ 1,
                                        base == 0 & sce >0 ~ (sce/(base+.1)), # NA_real_, 
                                        T ~ (sce/base)) %>% round(., 3)) %>%
            ungroup() %>%
            mutate(ens = ens.t)
          
          
          
          res.cp = rbind(res.cp, tda)
          
        }
      } # ens
      
      if(is.null(res.cp))
        next
      
      # summarize across the ensemble
      setDT(res.cp)
      ens.stat = rbind(res.cp[, fn_stat(perc.redn), 
                              by = c("season",'timeframe','mask.grp','age.grp', 'variable')] %>%
                         mutate(metric = 'perc.redn'),
                       res.cp[, fn_stat(ratio), 
                              by = c("season",'timeframe','mask.grp','age.grp', 'variable')] %>%
                         mutate(metric = 'ratio')) %>%
        mutate(# diff.ctc.bymask = diff.ctc.bymask.t,
          # diff.ctc.redn.bymask = diff.ctc.redn.bymask.t,
          ave.frac.masking = frac.masking_vec[fm.t],
          cut.mask.start = cut.mask.start_vec[cm.t])
      
      RES = rbind(RES, ens.stat)
    } # cut.mask.start
  } # frac.masking
  common.parms = list(inclVax = inclVax.t, # whether to include vaccination
                      UseContactDuration = UseContactDuration.t, #  whether to use contact duration directly or use number of contacts
                      UseVxData = UseVxData.t, # whether to use vax data for the vaccination module
                      UseAgeSpecSuscept = UseAgeSpecSuscept.t, # whether to assume diff susceptibility by age. set it to FALSE to reduce uncertainty
                      diff.ctc.redn.bymask = diff.ctc.redn.bymask.t, # whether assume diff level of redn in contact rate by masking status
                      diff.ctc.bymask = diff.ctc.bymask.t, # whether to assume diff (e.g. lower) contact rate for those prefer masking
                      MASKATHOME = MASKATHOME.t, # set this to TRUE only when testing a very specific scenario that ppl also mask at home
                      parm.type = parm.type.t, # which parm set to use
                      cont_freq = cont_freq.t, # for the contact data, diff contact frequency
                      wt_nonphys = wt_nonphys.t  # assume non-phyical contact has 20% of the chance resulting infection
  )
  
  save(RES, common.parms, file = paste0(dir_out, 'out_diff.byeff.bysn.',
                                        ifelse(diff.ctc.bymask.t,'CTCbymask','CTCo'),'_',
                                        ifelse(diff.ctc.redn.bymask.t,'MOBbymask','MOBo'),
                                        'fm',fm.no,
                                        '.RData'))
}

