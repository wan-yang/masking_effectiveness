# Model code used in Yang & Shaman "Reconciling the efficacy and effectiveness of masking on epidemic outcomes"
# author: Wan Yang
# date: March 2023
# note: this code is for research use only; may not be entirely optimized nor neatly organized. 

# compute relative diff vs. baseline, by season

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

tmp = try(load(paste0(dir_out, 'out_diff.bysn.',
                      ifelse(vsBaseCTCoMOBo,'vsBase0','vsBase1'),'_',
                      ifelse(diff.ctc.bymask.t,'CTCbymask','CTCo'),'_',
                      ifelse(diff.ctc.redn.bymask.t,'MOBbymask','MOBo'),
                      '.RData')))

if(class(tmp) == 'try-error'){ # not there yet
  RES = NULL
  
  idx.wt = 2 # main setting
  if(!exists('wt_nonphys.t'))
    idx.wt = 1:length(wt_nonphys_vec)
  
  for(wt.t in idx.wt){
    wt_nonphys.t = wt_nonphys_vec[wt.t]
    for(fm.t in 2:length(frac.masking_vec)){  # length(frac.masking_vec)
      for(cm.t in 1:length(cut.mask.start_vec)){ # length(cut.mask.start_vec)
        
        print(paste('wt.t =', wt.t, 'fm.t =', fm.t, 'cm.t =', cm.t))
        if(fm.t == 1 & cm.t == 1)
          next
        
        fm0 = 1; cm0 = 1 # baseline: no mask
        if(vsBaseCTCoMOBo){ # use the same baseline regardless of diff.ctc.bymask.t / diff.ctc.redn.bymask.t
          fbase = paste0(dir_res, 'res_',
                         ifelse(inclVax.t,'vx','novx'),'_',
                         ifelse(UseAgeSpecSuscept.t,'Sbyage','So'),'_CTCo_MOBo',
                         # ifelse(diff.ctc.bymask.t,'CTCbymask','CTCo'),'_',
                         # ifelse(diff.ctc.redn.bymask.t,'MOBbymask','MOBo'),
                         ifelse(wt_nonphys.t == wt_nonphys.main, '', paste0('_wt', wt.t)),
                         '_fm', fm0, '_cm', cm0)
        } else {
          fbase = paste0(dir_res, 'res_',
                         ifelse(inclVax.t,'vx','novx'),'_',
                         ifelse(UseAgeSpecSuscept.t,'Sbyage','So'),'_',
                         ifelse(diff.ctc.bymask.t,'CTCbymask','CTCo'),'_',
                         ifelse(diff.ctc.redn.bymask.t,'MOBbymask','MOBo'),
                         ifelse(wt_nonphys.t == wt_nonphys.main, '', paste0('_wt', wt.t)),
                         '_fm', fm0, '_cm', cm0)  # , '_ens', ens_no, '.RData'
        }
        
        
        fm1 = fm.t; cm1 = cm.t # a different scenario to compare 
        fsce = paste0(dir_res, 'res_',
                      ifelse(inclVax.t,'vx','novx'),'_',
                      ifelse(UseAgeSpecSuscept.t,'Sbyage','So'),'_',
                      ifelse(diff.ctc.bymask.t,'CTCbymask','CTCo'),'_',
                      ifelse(diff.ctc.redn.bymask.t,'MOBbymask','MOBo'),
                      ifelse(wt_nonphys.t == wt_nonphys.main, '', paste0('_wt', wt.t)),
                      '_fm', fm1, '_cm', cm1)  # , '_ens', ens_no, '.RData'
        
        
        res.cp.bysn = NULL
        for(ens.t in 1:num_ens){
          fname0 = paste0(fbase, '_ens', ens.t, '.RData')
          fname1 = paste0(fsce, '_ens', ens.t, '.RData')
          
          tmp = try(load(fname1), silent = T)
          if(any(class(tmp) == 'try-error'))
            next
          
          res1 = res; rm(res)
          
          
          load(fname0); res0 = res; rm(res)
          
          # outs = c('res.daily','res.wkly','res.yrly','res.byyr.stat', 'res.bycohort.byage.stat','res.stat','ctc.stat','parm.list')
          outs = c('res.bysn.stat')
          for(oo in outs){
            for(otype in c('0', '1')){
              eval(parse(text = paste0('tmp = res',otype,'$',oo)))
              assign(paste0(oo,otype), tmp)
            }
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
            
            cp.bysn = res.bysn.stat0 %>% filter(# timeframe == stat.type.t & 
              mask.grp == mask.grp.t) %>% 
              reshape2::melt(., id.vars = c("season",'ens','timeframe','mask.grp','age.grp')) %>%
              rename(., base = value) %>% 
              left_join(., 
                        y = res.bysn.stat1 %>% filter(# timeframe == stat.type.t & 
                          mask.grp == mask.grp.t) %>% 
                          reshape2::melt(., id.vars = c("season",'ens','timeframe','mask.grp','age.grp')) %>%
                          rename(., sce = value),
                        by = c("season",'ens','timeframe','mask.grp','age.grp', 'variable')
              ) %>% 
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
            
            res.cp.bysn = rbind(res.cp.bysn, cp.bysn)
            
          }
        } # ens
        
        if(is.null(res.cp.bysn))
          next
        
        # summarize across the ensemble
        setDT(res.cp.bysn)
        ens.stat = rbind(res.cp.bysn[, fn_stat(perc.redn), 
                                     by = c("season",'timeframe','mask.grp','age.grp', 'variable')] %>%
                           mutate(metric = 'perc.redn'),
                         res.cp.bysn[, fn_stat(ratio), 
                                     by = c("season",'timeframe','mask.grp','age.grp', 'variable')] %>%
                           mutate(metric = 'ratio')) %>%
          mutate(# diff.ctc.bymask = diff.ctc.bymask.t,
            # diff.ctc.redn.bymask = diff.ctc.redn.bymask.t,
            ave.frac.masking0 = frac.masking_vec[fm0],
            ave.frac.masking1 = frac.masking_vec[fm1],
            cut.mask.start0 = cut.mask.start_vec[cm0],
            cut.mask.start1 = cut.mask.start_vec[cm1])
        
        if(wt_nonphys.t != wt_nonphys.main)
          ens.stat$wt_nonphys = wt_nonphys.t
        
        RES = rbind(RES, ens.stat)
      } # cut.mask.start
    } # frac.masking
    if(wt_nonphys.t == wt_nonphys.main){
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
    } else {
      common.parms = list(inclVax = inclVax.t, # whether to include vaccination
                          UseContactDuration = UseContactDuration.t, #  whether to use contact duration directly or use number of contacts
                          UseVxData = UseVxData.t, # whether to use vax data for the vaccination module
                          UseAgeSpecSuscept = UseAgeSpecSuscept.t, # whether to assume diff susceptibility by age. set it to FALSE to reduce uncertainty
                          diff.ctc.redn.bymask = diff.ctc.redn.bymask.t, # whether assume diff level of redn in contact rate by masking status
                          diff.ctc.bymask = diff.ctc.bymask.t, # whether to assume diff (e.g. lower) contact rate for those prefer masking
                          MASKATHOME = MASKATHOME.t, # set this to TRUE only when testing a very specific scenario that ppl also mask at home
                          parm.type = parm.type.t, # which parm set to use
                          cont_freq = cont_freq.t # for the contact data, diff contact frequency
      )
    }
  }
  
  
  save(RES, common.parms, file = paste0(dir_out, 'out_diff.bysn.',
                                        ifelse(vsBaseCTCoMOBo,'vsBase0','vsBase1'),'_',
                                        ifelse(diff.ctc.bymask.t,'CTCbymask','CTCo'),'_',
                                        ifelse(diff.ctc.redn.bymask.t,'MOBbymask','MOBo'),
                                        '.RData'))
}


