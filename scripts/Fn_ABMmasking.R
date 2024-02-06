# Model code used in Yang & Shaman "Reconciling the efficacy and effectiveness of masking on epidemic outcomes"
# author: Wan Yang
# date: March 2023
# note: this code is for research use only; may not be entirely optimized nor neatly organized. 

# Agent-based model (ABM) for modeling masking impact

# could run multiple models (model ensemble) simultaneously
# why the ensemble - potentially used in conjunction with the filter


Fn_ABMmasking = function(model_inputs, model_settings, agg.by.cohort = F, agg.by.sn = T){ # 
  
  # get the model inputs from model_inputs
  curEnv=environment()
  list2env(model_inputs, envir = curEnv)
  list2env(model_settings, envir = curEnv)
  
  # set parms accordingly
  parms.t = list(num_ens = num_ens,
                 num_pop = num_pop, # population size
                 year.cur = tab_dates[1]$year, # year simulation starts
                 UseContactDuration = UseContactDuration, #  whether to use contact duration directly or use number of contacts
                 frac.ini.inf = frac.ini.inf0, # initial fraction of population being infectious 
                 frac.masking = c(0, rep(ave.frac.masking, nrow(age.grps)-1)), # initial fraction of population adopt masking
                 diff.ctc.bymask = diff.ctc.bymask,
                 ave.p.redn.mask = ave.p.redn.mask.byage$value, # average effectiveness, mean
                 sd.p.redn.mask = sd.p.redn.mask.byage$value, # effectiveness, standard deviation (sd)
                 age.grps = age.grps,
                 ave.num.nonHHctc = ave.num.nonHHctc.byage$value, # average number of non hh contact for diff age groups, mean
                 sd.num.nonHHctc = sd.num.nonHHctc.byage$value,  # number of non hh contact for diff age groups, sd
                 ave.num.HHctc = ave.num.HHctc.byage$value, # average number of hh contact (exclude self)
                 sd.num.HHctc = sd.num.HHctc.byage$value, # number of hh contact (exclude self), sd
                 birth.rate =birth.rate)  
  
  ndays = tab_dates %>% nrow()
  
  # record summary results
  res.daily = NULL
  res.wkly = NULL
  res.yrly = NULL
  res.snly = NULL
  res.byyr.stat = NULL  # to record stats by year
  res.bysn.stat = NULL  # to record stats by season 
  res.bycohort.byage.stat = NULL # to record birth cohort and age
  
  rec = NULL # ppl of the entire ensemble
  num_pop.mask = NULL
  num_pop.nomask = NULL
  num_births.per.day = NULL
  
  # initialize the system for each ens member and combine
  tmp = fn_initalizeABM_ens(parms.t)
  rec = tmp$rec %>% data.table()
  n.byage_ens.t = tmp$n.byage_ens
  # save a record for comparison with next year
  rec.last = rec %>% 
    dplyr::select(-status, -tm.last.event, -tm.last.inf, -tm.last.hosp,-tm.last.vax, -num.nonHHctc, -num.HHctc, -masking, p.redn.mask) %>% data.table()
  
  rec.last.sn = rec.last  # for recording by season
  
  R0ini.t = tmp$R0ini %>% round(2)
  R0ini.no.mask.t = tmp$R0ini.no.mask  %>% round(2)
  year.cur = parms.t$year.cur
  
  num_pop.bymask = rec %>% group_by(ens, masking) %>% summarise(num_pop.bymask = n()) %>% ungroup
  # num_pop.mask = tmp %>% dplyr::filter(masking == 1) %>% arrange(., by = 'ens') # %>% .$n
  # num_pop.nomask = tmp %>% dplyr::filter(masking == 0) %>% arrange(., by = 'ens') #  %>% .$n
  
  num_births.per.day = with(parms.t, birth.rate * num_pop)
  
  
  newi.all = newi.nomask = newi.mask = newi.nomask.itt = newi.mask.itt = matrix(0, ndays, num_ens)
  newh.all = newh.nomask = newh.mask = newh.nomask.itt = newh.mask.itt = matrix(0, ndays, num_ens)
  newv.all = newv.nomask = newv.mask = matrix(0, ndays, num_ens)  # track num of vaccination
  
  res.yrly = res.wkly = NULL  # for quality check 
  
  for(tt in 1:(ndays)){
    
    # update year, used to set birth year
    if(tab_dates[tt]$DayOfYr==1)
      parms.t$year.cur = tab_dates[tt]$year # parms.t$year.cur + 1 # used to set birth year
    
    # random seeding - do it before getting the current status or at the end of the time step
    # to reduce model stochastic, do it by week - see below
    
    
    # compute population prevalence at time step
    Itot.t = rec %>% 
      dplyr::filter(status == 'I') %>% group_by(ens) %>% 
      summarise(AVEprev = n() / num_pop) %>% # do not use: + seed0
      full_join(., y = data.table(ens = 1:num_ens), by = 'ens') %>% 
      replace_na(., list(AVEprev = 0)) # !!! set NA to 0, so as not to introduce difference due to seeding
    
    Istra.t = rec %>% dplyr::filter(status == 'I') %>%
      group_by(ens, masking) %>% 
      summarise(HHprev = n()) %>% # may have missing entries, if no one has status == 'I'
      left_join(., y = num_pop.bymask, by = c('ens', 'masking')) %>%
      dplyr::mutate(HHprev = HHprev / num_pop.bymask) %>% # may have missing entries, if no one has status == 'I'
      full_join(., y = data.table(ens = 1:num_ens), by = 'ens') %>% 
      replace_na(., list(HHprev = 0)) # also set NA to 0, so as not to introduce difference due to seeding
    
    
    print(paste(tab_dates[tt]$year, tab_dates[tt]$DayOfYr, mean(Itot.t$AVEprev) %>% round(5)))
    # update for each type of individuals
    
    # because upstream update will affect those updated later, get all the diff type of ppl at this time step first
    # for the susceptibles
    ppl.t = rec[status == 'S']; nppl.t = ppl.t %>% nrow #  %>% group_by(ens) %>% summarise(n = n())
    S.ppl.t = ppl.t; S.nppl.t = nppl.t
    # for the exposed
    ppl.t = rec[status == 'E'];  nppl.t = ppl.t %>% nrow #   %>% group_by(ens) %>% summarise(n = n())
    E.ppl.t = ppl.t; E.nppl.t = nppl.t
    # for the infectious
    ppl.t = rec[status == 'I'];  nppl.t = ppl.t %>% nrow #  %>% group_by(ens) %>% summarise(n = n())
    I.ppl.t = ppl.t; I.nppl.t = nppl.t
    # for the recovered
    ppl.t = rec[status == 'R'];  nppl.t = ppl.t %>% nrow #  %>% group_by(ens) %>% summarise(n = n())
    R.ppl.t = ppl.t; R.nppl.t = nppl.t
    # for the vaccinated and still having protection against infection
    ppl.t = rec[status == 'V'];  nppl.t = ppl.t %>% nrow #  %>% group_by(ens) %>% summarise(n = n())
    V.ppl.t = ppl.t; V.nppl.t = nppl.t
    
    # NOW DO THE UPDATES
    # for the susceptible
    ppl.t = S.ppl.t; nppl.t = S.nppl.t
    ppl.t = ppl.t %>% left_join(., y = Itot.t, by = c('ens')) %>%
      left_join(., y = Istra.t, by = c('ens', 'masking')) %>%
      replace_na(., list(AVEprev = 0, HHprev = 0))
    
    # for babies (no masking), assign them a lower prev, if masking is on (?) given the lower background force of infection
    # can't do this when no one masks (baseline)
    ave.prev.t = (Itot.t  %>% .$AVEprev %>% mean())
    ave.m.prev.t = Istra.t %>% filter(masking == 1) %>% .$HHprev %>% mean()
    if(is.na(ave.m.prev.t)){  # no one masks in the population
      HHprev.babies.t = ave.prev.t * .6 # reduce it by 40%
    } else {
      HHprev.babies.t = pmin(ave.m.prev.t, ave.prev.t * .6)
    }
    ppl.t = ppl.t %>% 
      dplyr::mutate(HHprev = case_when(age == 0 ~ HHprev.babies.t,
                                T ~ HHprev)) %>%
      data.table()
    
    
    # set mask and risk avoidance (cnt redn) for this time
    # note this will also be used to set adj.ctc.redn in the fn_foi function
    for(iens in 1:num_ens){
      if((Itot.t %>% filter(ens==iens) %>% .$AVEprev) < cut.mask.start |
         ave.frac.masking == 0 # need to do this for scenario that no one mask too, otherwise, systematically lower the prev 
         ){
        # not activated
        ppl.t[ens==iens]$masking = 0
        ppl.t[ens==iens]$HHprev = ppl.t[ens==iens]$AVEprev
      }
    }
    
    # add household (HH) specific risk based on masking
    # ppl.t = ppl.t %>% mutate(HHprev = case_when(masking == 0 ~ prev.nomask.t,
    #                                             masking == 1 ~ prev.mask.t))
    
    if(nppl.t > 0){
      # attach numbers to the data frame for matching
      tx.sn.t = 1 + tx.sn.amp * cos((tab_dates[tt]$DayOfYr %% 365)/365 * 2 * pi) # seasonal risk
      tt.t = tt
      p.ctc.redn.t = da.p.ctc.redn %>% filter(tt == tt.t) %>% .$value
      # print(paste('tx redn:', round(p.ctc.redn.t, 2)))
      ppl.t = ppl.t %>% 
        dplyr::mutate(p.t = fn_foi(ppl.t, 
                                   prob.tx.per.ctc0 = prob.tx.per.ctc0, # baseline probability of infection per infectious contact
                                   method.calHHrisk = method.calHHrisk,
                                   suscept.by.age = suscept.by.age, UseAgeSpecSuscept = UseAgeSpecSuscept, # allow age-specific susceptibility
                                   p.ctc.redn.t = p.ctc.redn.t, # default, no change
                                   adj.ctc.redn.bymask = adj.ctc.redn.bymask, # could tie to mask behavior
                                   UseContactDuration = UseContactDuration, # use contact duration from polymod data
                                   MASKATHOME = MASKATHOME
        ) * tx.sn.t,  # probability of infection this time step
        u = runif(nppl.t))
      idx.t = ppl.t %>% filter(u < p.t) %>% .$idx # use the ppl idx to match individuals
      
      if(length(idx.t) > 0){ # S -> E 
        rec[idx %in% idx.t]$status = 'E' # update status
        rec[idx %in% idx.t]$tm.last.event = tt # record time of event
        rec[idx %in% idx.t]$tm.last.inf = tt # record time of new infection
        rec = rec %>% dplyr::mutate(num.inf = case_when(idx %in% idx.t ~ num.inf + 1,
                                                 T ~ num.inf)) %>% data.table()
        
        # count new infection, without a delay (just for counting)
        # ppl.t %>% filter(idx %in% idx.t) %>% group_by(ens) %>% summarise(n = n())
        
        
        if(F){ # error: this will count number of rows
          newi.all[tt, ] = ppl.t %>% filter(idx %in% idx.t) %>% 
            full_join(., y = data.table(ens = 1:num_ens), by = 'ens') %>% 
            group_by(ens) %>%
            summarise(n = n()) %>%
            replace_na(., list(n = 0)) %>%
            dplyr::arrange(., ens) %>% .$n
        }
        newi.all[tt, ] = ppl.t %>% filter(idx %in% idx.t) %>% 
          group_by(ens) %>%
          summarise(n = n()) %>%
          full_join(., y = data.table(ens = 1:num_ens), by = 'ens') %>% 
          replace_na(., list(n = 0)) %>%
          dplyr::arrange(., ens) %>% .$n
        # based whether mask is in use at the time
        newi.mask[tt, ] = ppl.t %>% filter(idx %in% idx.t & masking == 1) %>% 
          group_by(ens) %>%
          summarise(n = n()) %>%
          full_join(., y = data.table(ens = 1:num_ens), by = 'ens') %>% 
          replace_na(., list(n = 0)) %>%
          dplyr::arrange(., ens) %>% .$n
        # count new infection among those use mask
        newi.nomask[tt, ] = ppl.t %>% filter(idx %in% idx.t & masking == 0) %>% 
          group_by(ens) %>%
          summarise(n = n()) %>%
          full_join(., y = data.table(ens = 1:num_ens), by = 'ens') %>% 
          replace_na(., list(n = 0)) %>%
          dplyr::arrange(., ens) %>% .$n
        # based on assigned 'treatment'
        newi.mask.itt[tt, ] = rec %>% filter(idx %in% idx.t & masking == 1) %>% 
          group_by(ens) %>%
          summarise(n = n()) %>%
          full_join(., y = data.table(ens = 1:num_ens), by = 'ens') %>% 
          replace_na(., list(n = 0)) %>%
          dplyr::arrange(., ens) %>% .$n
        # count new infection among those use mask
        newi.nomask.itt[tt, ] = rec %>% filter(idx %in% idx.t & masking == 0) %>% 
          group_by(ens) %>%
          summarise(n = n()) %>%
          full_join(., y = data.table(ens = 1:num_ens), by = 'ens') %>% 
          replace_na(., list(n = 0)) %>%
          dplyr::arrange(., ens) %>% .$n
        
        # compute severe case 
        ppl2.t = ppl.t %>% filter(idx %in% idx.t); 
        nppl2.t = ppl2.t %>% nrow() # group_by(ens) %>% summarise(n = n())
        ppl2.t = ppl2.t %>%
          dplyr::mutate(p.t = fn_severity(tt, ppl.t = ., severity.by.age, parm.immLossVseverity), 
                 u = runif(nppl2.t))
        idx2.t = ppl2.t %>% filter(u < p.t) %>% .$idx # severe case
        
        # record it
        if(length(idx2.t) > 0){
          rec[idx %in% idx2.t]$num.hosp = rec[idx %in% idx2.t]$num.hosp  + 1
          rec[idx %in% idx2.t]$tm.last.hosp = tt # record time of new hospitalization
          
          # count new hospitalization, without a delay (just for counting)
          newh.all[tt, ] = ppl2.t %>% filter(idx %in% idx2.t) %>% 
            group_by(ens) %>%
            summarise(n = n()) %>%
            full_join(., y = data.table(ens = 1:num_ens), by = 'ens') %>% 
            replace_na(., list(n = 0)) %>%
            dplyr::arrange(., ens) %>% .$n
          newh.mask[tt, ] = ppl2.t %>% filter(idx %in% idx2.t & masking == 1) %>% 
            group_by(ens) %>%
            summarise(n = n()) %>%
            full_join(., y = data.table(ens = 1:num_ens), by = 'ens') %>% 
            replace_na(., list(n = 0)) %>%
            dplyr::arrange(., ens) %>% .$n
          # count new infection among those use mask
          newh.nomask[tt, ] = ppl2.t %>% filter(idx %in% idx2.t & masking == 0) %>% 
            group_by(ens) %>%
            summarise(n = n()) %>%
            full_join(., y = data.table(ens = 1:num_ens), by = 'ens') %>% 
            replace_na(., list(n = 0)) %>%
            dplyr::arrange(., ens) %>% .$n
          # based on assigned 'treatment'
          newh.mask.itt[tt, ] = rec %>% filter(idx %in% idx2.t & masking == 1) %>% 
            group_by(ens) %>%
            summarise(n = n()) %>%
            full_join(., y = data.table(ens = 1:num_ens), by = 'ens') %>% 
            replace_na(., list(n = 0)) %>%
            dplyr::arrange(., ens) %>% .$n
          newh.nomask.itt[tt, ] = rec %>% filter(idx %in% idx2.t & masking == 0) %>% 
            group_by(ens) %>%
            summarise(n = n()) %>%
            full_join(., y = data.table(ens = 1:num_ens), by = 'ens') %>% 
            replace_na(., list(n = 0)) %>%
            dplyr::arrange(., ens) %>% .$n
        }
        
      }
    }
    rm(ppl.t, nppl.t, idx.t, ppl2.t, nppl2.t, idx2.t)
    
    # for the exposed
    # ppl.t = rec[status == 'E'];  nppl.t = ppl.t %>% nrow
    ppl.t = E.ppl.t; nppl.t = E.nppl.t
    if(nppl.t > 0){
      # attach numbers to the data frame for matching
      ppl.t = ppl.t %>% 
        dplyr::mutate(p.t = 1 - exp(-1/Tei * (tt - ppl.t$tm.last.event)),  # probability of moving to being infectious
               u = runif(nppl.t))
      idx.t = ppl.t %>% filter(u < p.t) %>% .$idx # use the ppl idx to match individuals
      
      if(length(idx.t) > 0){ # E -> I
        rec[idx %in% idx.t]$status = 'I' # update status
        rec[idx %in% idx.t]$tm.last.event = tt # record time of event
        
      } 
    }
    rm(ppl.t, nppl.t, idx.t)
    
    # for the infectious
    ppl.t = I.ppl.t; nppl.t = I.nppl.t
    if(nppl.t > 0){
      ppl.t = ppl.t %>% 
        dplyr::mutate(p.t = 1 - exp(-1/Tir * (tt - ppl.t$tm.last.event)),  # probability of moving to being recovered
               u = runif(nppl.t))
      idx.t = ppl.t %>% filter(u < p.t) %>% .$idx # use the ppl idx to match individuals
      
      if(length(idx.t) > 0){ # I -> R
        rec[idx %in% idx.t]$status = 'R' # update status
        rec[idx %in% idx.t]$tm.last.event = tt # record time of event
      } 
    }
    rm(ppl.t, nppl.t, idx.t)
    
    # for the recovered
    ppl.t = R.ppl.t; nppl.t = R.nppl.t
    if(nppl.t > 0){
      ppl.t = ppl.t %>% 
        dplyr::mutate(p.t = fn_immLossLogistic(tt, ppl.t, parm.immLossInf),  # probability of moving to being susceptible
               u = runif(nppl.t))
      idx.t = ppl.t %>% dplyr::filter(u < p.t) %>% .$idx # use the ppl idx to match individuals
      
      if(length(idx.t) > 0){ # R -> S
        rec[idx %in% idx.t]$status = 'S' # update status - back to susceptible
        rec[idx %in% idx.t]$tm.last.event = tt # record time of event
      }
    }
    rm(ppl.t, nppl.t, idx.t)
    
    # now do vaccine waning protection
    # for the vax'ed and protected
    ppl.t = V.ppl.t; nppl.t = V.nppl.t
    if(nppl.t > 0){
      ppl.t = ppl.t %>% 
        dplyr::mutate(p.t = fn_immLossLogistic(tt, ppl.t, parm.immLossVax),  # probability of moving to being susceptible
               u = runif(nppl.t))
      idx.t = ppl.t %>% dplyr::filter(u < p.t) %>% .$idx # use the ppl idx to match individuals
      if(length(idx.t) > 0){ # V -> S
        rec[idx %in% idx.t]$status = 'S' # update status - back to susceptible
        rec[idx %in% idx.t]$tm.last.event = tt # record time of event
        # rec[idx %in% idx.t]$num.vax = rec[idx %in% idx.t]$num.vax - .5 # allow partial protection against severe disease to remain
        # no need to update number of vax, time since last vaccination will help gauge the level of protection
      }
    }
    rm(ppl.t, nppl.t, idx.t)
    
    # NOW DO VACCINATION
    # n.vax.t = num_vax[tt]  # number of effective vaccination for this day
    # sample n.vax.t ppl from those S, E, or R
    # if an 'S' person is selected for vax, update status to V, 
    #   update tm.vax, update tm.last.event, ++num.vax 
    # if an 'E' person is selected for vax, keep status as E (i.e. course of infection contintue), 
    #   keep tm.last.event, update tm.vax, num.vax = +.5 (allow vax to have partial effect on disease severity)
    # if an 'R' person is selected for vax, update status to V, 
    #   update tm.last.event, update tm.vax, ++num.vax 
    # vax of 'I' not allowed
    if(inclVax & 
       tt >= tt.vax.start &
       tab_dates[tt]$DayOfWk == 7  # do it weekly
    ){  # n.vax.t > 0
      
      if(tt == tt.vax.start)
        print('start vaccination')
      
      # reason for "over-vaccination": individuals are checked for eligibility (asked to get vaccinated) many times before they get vaccinated
      # and assume all those times are independent
      # need a diff function for likelihood of vaccination, after the eligibility, prob to get vax go down again
      # could be periodic 
      ppl.t = rec %>% dplyr::filter(# ((tt - tm.last.vax) > t2revax.min | is.na(tm.last.vax)) & # also restrict timing 
        status %in% c('S', 'E', 'R'));  # pool of eligible ppl
      nppl.t = ppl.t %>% nrow
      
      if(UseVxData){
        # based on daily or weekly vx data/trends
        tt.t = tt
        da.vx.t = da.vx.weekly %>% filter(tt == tt.t)
        all.idx.t = fn_doVaccination(ppl.t, tt, da.vx.t, n.byage_ens.t)
        
      } else {
        # simulated vax, w/o data
        # set eligibility and vax.rate for this time step
        vax.rate.t = ini.vax.rate[pmax(0, pmin(tt -tt.vax.start+1, length(ini.vax.rate)))]
        if(tt %in% (tt.vax.start + 0:365)){  # 0:365 limit to the first three months(?) so don't over vaccinate
          # first year, primary series
          vax.coverage.byage.t = vax.coverage.byage.1st %>%
            dplyr::mutate(value = value * vax.rate.t)
          vx.sn.t = 1 # no seasonal vx for the first year
        } else {
          # booster
          vax.coverage.byage.t = vax.coverage.byage.boost 
          vx.sn.t = 1 + vx.sn.amp * cos(((tab_dates[tt]$DayOfYr + 30*2) %% 365)/365 * 2 * pi) # seasonal vx, assume vx occurs ~2 months prior
          
        }
        
        vax.elig.byage.t = vax.elig.byage
        ppl.t = ppl.t %>% dplyr::mutate(p.t = fn_getVaxProbLogistic(tt, ppl.t,vax.elig.byage.t, # set eligibility by age
                                                             vax.coverage.byage.t,
                                                             parm.getVaxProb) * 
                                   vx.sn.t, # add the same seasonality, b/c vx is often tied to timing of outbreak
                                 # include further eligibility for vax 
                                 u = runif(nppl.t)) 
        
        
      }
      
      
      # for those Susceptible or recovered - combine them
      if(UseVxData){
        idx.t = ppl.t %>% filter(idx %in% all.idx.t & status %in% c('S', 'R')) %>% .$idx
      } else {
        idx.t = ppl.t %>% dplyr::filter(u < p.t & status %in% c('S', 'R')) %>% .$idx # use the ppl idx to match individuals
      }
      
      if(length(idx.t) > 0){
        # if an 'S' person is selected for vax, update status to V, 
        #   update tm.vax, update tm.last.event, ++num.vax
        rec[idx %in% idx.t]$status = 'V' # update status 
        rec[idx %in% idx.t]$tm.last.event = tt # record time of event
        rec[idx %in% idx.t]$tm.last.vax = tt
        rec[idx %in% idx.t]$num.vax = rec[idx %in% idx.t]$num.vax + 1
      }
      # for those newly infected (E)
      if(UseVxData){
        idx.t = ppl.t %>% filter(idx %in% all.idx.t & status == 'E') %>% .$idx
      } else {
        idx.t = ppl.t %>% dplyr::filter(u < p.t & status == 'E') %>% .$idx # use the ppl idx to match individuals
      }
      if(length(idx.t) > 0){
        # if an 'E' person is selected for vax, keep status as E (i.e. course of infection contintue), 
        #   keep tm.last.event, update tm.vax, num.vax = +.5 (allow vax to have partial effect on disease severity)
        # do not update status to E or update tm.last.event for the infection to run its course
        rec[idx %in% idx.t]$tm.last.vax = tt
        rec[idx %in% idx.t]$num.vax = rec[idx %in% idx.t]$num.vax + .5 
      }
      
      
      # count vaccination, no do it when use simulated vx prob?
      if(UseVxData){
        idx.t = all.idx.t
      } else {
        idx.t = ppl.t %>% dplyr::filter(u < p.t) %>% .$idx # use the ppl idx to match individuals
        
      }
      newv.all[tt, ] = ppl.t %>% filter(idx %in% idx.t) %>% 
        group_by(ens) %>%
        summarise(n = n()) %>%
        full_join(., y = data.table(ens = 1:num_ens), by = 'ens') %>% 
        replace_na(., list(n = 0)) %>%
        dplyr::arrange(., ens) %>% .$n
      newv.mask[tt, ] = ppl.t %>% filter(idx %in% idx.t & masking == 1) %>% 
        group_by(ens) %>%
        summarise(n = n()) %>%
        full_join(., y = data.table(ens = 1:num_ens), by = 'ens') %>% 
        replace_na(., list(n = 0)) %>%
        dplyr::arrange(., ens) %>% .$n
      # count new infection among those use mask
      newv.nomask[tt, ] = ppl.t %>% filter(idx %in% idx.t & masking == 0) %>% 
        group_by(ens) %>%
        summarise(n = n()) %>%
        full_join(., y = data.table(ens = 1:num_ens), by = 'ens') %>% 
        replace_na(., list(n = 0)) %>%
        dplyr::arrange(., ens) %>% .$n
      
    } # end vaccination
    rm(all.idx.t, idx.t, ppl.t, nppl.t)
    
    # tally at the end of the year (do it before updating demographics)
    # if(tt %% 365 == 0 & tt !=0){
    if(tab_dates[tt]$EndOfYr & tt !=0){ 
      
      # idx.t = tab_dates %>% filter(year == tab_dates[tt]$year) %>% .$tt. # dunno why, somehow did not work
      yr.t = tab_dates[tt]$year
      idx.t = tab_dates %>% dplyr::filter(year == yr.t) %>% .$tt
      # iy = tab_dates[tt]$year #  iy + 1
      res.yrly = rbind(res.yrly,
                       data.table(year = tab_dates[tt]$year, 
                                  ens = 1:num_ens,
                                  newi.all = newi.all[idx.t, ,drop=F] %>% colSums, 
                                  newi.nomask = newi.nomask[idx.t, ,drop=F] %>% colSums, 
                                  newi.mask = newi.mask[idx.t, ,drop=F] %>% colSums,
                                  newi.nomask.itt = newi.nomask.itt[idx.t, ,drop=F] %>% colSums, 
                                  newi.mask.itt = newi.mask.itt[idx.t, ,drop=F] %>% colSums,
                                  newh.all = newh.all[idx.t, ,drop=F] %>% colSums, 
                                  newh.nomask = newh.nomask[idx.t, ,drop=F] %>% colSums, 
                                  newh.mask = newh.mask[idx.t, ,drop=F] %>% colSums,
                                  newh.nomask.itt = newh.nomask.itt[idx.t, ,drop=F] %>% colSums, 
                                  newh.mask.itt = newh.mask.itt[idx.t, ,drop=F] %>% colSums,
                                  newv.all = newv.all[idx.t, ,drop=F] %>% colSums, 
                                  newv.nomask = newv.nomask[idx.t, ,drop=F] %>% colSums, 
                                  newv.mask = newv.mask[idx.t, ,drop=F] %>% colSums))
      
      # get detailed stats for this year
      tmp = fn_getStatYrly(rec, rec.last, idx.t)
      tmp.cum = tmp$stat.cum
      tmp.this = tmp$stat.this
      
      res.byyr.stat = rbind(res.byyr.stat, 
                            data.table(year = tab_dates[tt]$year,
                                                      tmp.cum),
                            data.table(year = tab_dates[tt]$year,
                                       tmp.this))
      rm(tmp, tmp.this, tmp.cum)
      
      # by birth cohort and age, as mask behavior change over time, hard to classify
      # only do it when agg.by.cohort == T
      if(agg.by.cohort){
        tmp = fn_getStatYrlyByCohort(rec, rec.last, idx.t)
        tmp.cum = tmp$stat.cum
        tmp.this = tmp$stat.this
        
        res.bycohort.byage.stat = rbind(res.bycohort.byage.stat, 
                                        data.table(year = tab_dates[tt]$year,
                                                   tmp.cum),
                                        data.table(year = tab_dates[tt]$year,
                                                   tmp.this))
        rm(tmp, tmp.this, tmp.cum)
      }
      
      
      # save a record for comparison with next year
      rec.last = rec %>% dplyr::select(-status, -tm.last.event, -tm.last.inf, -tm.last.hosp,-tm.last.vax, -num.nonHHctc, -num.HHctc) %>% data.table()
      
    }
    
    # tally by season 
    if(tab_dates[tt]$MidYr & tt !=0 & 
       agg.by.sn  # only do it when agg.by.sn == T
       ){ 
      
      season.t = tab_dates[tt]$season
      idx.t = tab_dates %>% dplyr::filter(season == season.t) %>% .$tt
      
      res.snly = rbind(res.snly,
                       data.table(season = tab_dates[tt]$season, 
                                  ens = 1:num_ens,
                                  newi.all = newi.all[idx.t, ,drop=F] %>% colSums, 
                                  newi.nomask = newi.nomask[idx.t, ,drop=F] %>% colSums, 
                                  newi.mask = newi.mask[idx.t, ,drop=F] %>% colSums,
                                  newi.nomask.itt = newi.nomask.itt[idx.t, ,drop=F] %>% colSums, 
                                  newi.mask.itt = newi.mask.itt[idx.t, ,drop=F] %>% colSums,
                                  newh.all = newh.all[idx.t, ,drop=F] %>% colSums, 
                                  newh.nomask = newh.nomask[idx.t, ,drop=F] %>% colSums, 
                                  newh.mask = newh.mask[idx.t, ,drop=F] %>% colSums,
                                  newh.nomask.itt = newh.nomask.itt[idx.t, ,drop=F] %>% colSums, 
                                  newh.mask.itt = newh.mask.itt[idx.t, ,drop=F] %>% colSums,
                                  newv.all = newv.all[idx.t, ,drop=F] %>% colSums, 
                                  newv.nomask = newv.nomask[idx.t, ,drop=F] %>% colSums, 
                                  newv.mask = newv.mask[idx.t, ,drop=F] %>% colSums))
      
      # get detailed stats for this season
      tmp = fn_getStatYrly(rec, rec.last.sn, idx.t)
      tmp.cum = tmp$stat.cum
      tmp.this = tmp$stat.this
      tmp.this$timeframe = 'this.season'
      
      res.bysn.stat = rbind(res.bysn.stat, 
                            data.table(season = tab_dates[tt]$season,
                                       tmp.cum),
                            data.table(season = tab_dates[tt]$season,
                                       tmp.this))
      rm(tmp, tmp.this, tmp.cum)
      
      
      # save a record for comparison with next season
      rec.last.sn = rec %>% dplyr::select(-status, -tm.last.event, -tm.last.inf, -tm.last.hosp,-tm.last.vax, -num.nonHHctc, -num.HHctc) %>% data.table()
      
    }
    
    
    # add birth and death - important for long term simulation (severity, etc, so wouldn't over count # reinfection)
    # also do random seeding here
    if(tt %% 7 == 0 & tt !=0 & tt!=ndays){
      # also record state variables
      tmp = rec %>% group_by(ens, status) %>% summarise(n = n()) %>%
        full_join(., y = data.table(ens = rep(1:num_ens, e=length(status_vec)), status = rep(status_vec, num_ens)), by = c('ens','status')) %>%
        replace_na(., list(n = 0)) %>%
        dplyr::mutate(perc = n / num_pop * 100, week.ending = tab_dates[tt]$date) %>%
        dplyr::mutate(n = NULL) %>%
        setcolorder(., c("ens", "week.ending", "status", "perc"))
      
      res.wkly = rbind(res.wkly, 
                       tmp)
      
      # do it by week 
      # b/c number births per day is not a round number, first sample to get the number 
      num_births.per.day.t = rpois(num_ens, num_births.per.day * 7) # sample from a Possion dist with the average
      # idx.t = sample(1:num_pop, num_births.per.day.t, replace = F) # sample these people: die and reborn
      
      if(tab_dates[tt]$year == tab_dates[1]$year){
        # if it is the first year of simulation
        # b/c the newborns have been generated at the start of the simulation for this year
        # to maintain the same cohort size, rather than replacing the new deaths with newborns
        # randomly sample from the <1 yr and replace them to account for the differential exposure for this birth cohort
        # as a result, no deaths the initial year
        ppl.t = rec %>% filter(age < 1) %>% # sample from the initial newborns
          group_by(ens) %>%
          dplyr::mutate(p.t = 1/n()) # random sampling
        
      } else {
        # if not the first year, 
        # sample from the oldest age group and replace with newborns
        ppl.t = rec %>% filter(age >= tail(age.grps$age.start,1)) # sample from the oldest age group
        # compute probability based on age
        ppl.t = ppl.t %>% dplyr::mutate(p.t = fn_getProbDeath(ppl.t)) 
        
      }
      
      # sample by ens
      idx.newborn = NULL
      for(ii in 1:num_ens){
        
        if(num_births.per.day.t[ii]<1)
          next
        
        if(nrow(ppl.t) <= num_births.per.day.t[ii]){
          idx.newborn.t = ppl.t$idx
          print('Birth: small smp size: all get replace!')
        } else {
          idx.newborn.t = sample(x = ppl.t %>% dplyr::filter(ens == ii) %>% .$idx, size = num_births.per.day.t[ii], 
                                 prob = ppl.t %>% dplyr::filter(ens == ii) %>% .$p.t,
                                 replace = F) # sample these people: die and reborn
        }
        
        idx.newborn = c(idx.newborn, idx.newborn.t)
      }
      
      # re-initiate the newborns
      # ppl.t = rec[idx %in% idx.t]
      if(length(idx.newborn) > 0)
        rec = fn_initalizeABMnewborns(parms.t, rec, idx.newborn)  %>% data.table()
      # caution - this would affect counting of total # infection based on final rec
      # if running the modeling for many years, should tally yearly to minimize the changes due to population turn over 
      
      # do random seeding
      num_seed.per.week.t = rpois(num_ens, seed0 * 7 * num_pop) # sample from a Poisson dist with the average
      # to reduce model stochasticity and considering high chance of epidemic and introduction
      if(tab_dates[tt]$DayOfYr %in% days.epi.seeding){
        num_seed.per.week.t = num_seed.per.week.t + 1
        print('epi week seeding')
      }
      
      ens2seed = which(num_seed.per.week.t > 0) # only do it for the ones with >1 
      if(length(ens2seed) > 0){
        for(iens in ens2seed){
          # should sample from those susceptile
          smp.seed = rec %>% filter(ens == iens & status == 'S') %>% .$ens.idx
          idx.seed = sample(smp.seed, size = num_seed.per.week.t[iens], replace = F)
          rec[ens == iens & ens.idx %in% idx.seed]$status = 'E'  # exposed
          rec[ens == iens & ens.idx %in% idx.seed]$tm.last.event = tt # record time of event
          rec[ens == iens & ens.idx %in% idx.seed]$tm.last.inf = tt # record time of new infection
          rec[ens == iens & ens.idx %in% idx.seed]$num.inf = rec[ens == iens & ens.idx %in% idx.seed]$num.inf + 1  # record number times infected
        }
      }
    }
    
    # population turnover, done yearly
    if(tab_dates[tt]$EndOfYr & tt !=0 & tt!=ndays){
      print('update pop age') 
      rec = fn_pop.aging(parms.t, rec) %>% data.table()
    }
    
    
  } # end all times
  
  # save results
  colnames(newi.all) = colnames(newi.nomask) = colnames(newi.mask) = colnames(newi.nomask.itt) = colnames(newi.mask.itt) = paste0('ens', 1:num_ens)
  colnames(newh.all) = colnames(newh.nomask) = colnames(newh.mask) = colnames(newh.nomask.itt) = colnames(newh.mask.itt) = paste0('ens', 1:num_ens)
  colnames(newv.all) = colnames(newv.nomask) = colnames(newv.mask) = paste0('ens', 1:num_ens)
  res.daily = rbind(data.table(measure = 'infections', mask.grp = 'all', date = tab_dates$date, newi.all), 
                    data.table(measure = 'infections', mask.grp = 'no mask (actual)', date = tab_dates$date, newi.nomask),
                    data.table(measure = 'infections', mask.grp = 'use mask (actual)', date = tab_dates$date, newi.mask),
                    data.table(measure = 'infections', mask.grp = 'no mask (assigned)', date = tab_dates$date, newi.nomask.itt),
                    data.table(measure = 'infections', mask.grp = 'use mask (assigned)', date = tab_dates$date, newi.mask.itt),
                    data.table(measure = 'hospitalizations', mask.grp = 'all', date = tab_dates$date, newh.all), 
                    data.table(measure = 'hospitalizations', mask.grp = 'no mask (actual)', date = tab_dates$date, newh.nomask),
                    data.table(measure = 'hospitalizations', mask.grp = 'use mask (actual)', date = tab_dates$date, newh.mask),
                    data.table(measure = 'hospitalizations', mask.grp = 'no mask (assigned)', date = tab_dates$date, newh.nomask.itt),
                    data.table(measure = 'hospitalizations', mask.grp = 'use mask (assigned)', date = tab_dates$date, newh.mask.itt),
                    data.table(measure = 'vaccinations', mask.grp = 'all', date = tab_dates$date, newv.all), 
                    data.table(measure = 'vaccinations', mask.grp = 'no mask', date = tab_dates$date, newv.nomask),
                    data.table(measure = 'vaccinations', mask.grp = 'use mask', date = tab_dates$date, newv.mask))
  
  
  # get summary stats across all years
  res.stat = fn_getStat(rec)
  
  # got summary stats for the pop
  ctc.stat = fn_getPopStat(rec)
  
  parm.list = model_settings 
  parm.list$status_vec = NULL
  parm.list$str_age_recode = NULL
  parm.list$R0ini = R0ini.t
  parm.list$R0ini.no.mask = R0ini.no.mask.t
  
  return(list(
    res.daily = res.daily, # daily tallies
    res.wkly = res.wkly, # % of population in each dis states at the end of each week
    res.yrly = res.yrly,  # yearly tallies
    res.byyr.stat = res.byyr.stat, # both yearly, and cumulative stats
    # by season
    res.snly = res.snly,  # season tallies
    res.bysn.stat = res.bysn.stat, # both season, and cumulative stats
    res.bycohort.byage.stat = res.bycohort.byage.stat, # stats by birth cohort
    res.stat = res.stat, # cumulative, and across all ens
    ctc.stat = ctc.stat, # summary stats for checking ctc rate by age group, and masking assignment
    # save the stats seperately to reduce space
    parm.list = parm.list
  ))
}

