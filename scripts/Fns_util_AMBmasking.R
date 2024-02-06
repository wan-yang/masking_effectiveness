# Model code used in Yang & Shaman "Reconciling the efficacy and effectiveness of masking on epidemic outcomes"
# author: Wan Yang
# date: March 2023
# note: this code is for research use only; may not be entirely optimized nor neatly organized. 

# Functions for running the ABM 
# all functions used to run the masking ABM

# function to get the corresponding param value
fn_get.parmx1 = function(parm.name, type = 'default'){
  if(type == 'default'){
    parmx = parm.bounds[parm == parm.name]$default 
  } else if(type == 'smp'){
    # random sampling from lower, upper bound
    parmx = runif(1, parm.bounds[parm == parm.name]$lwr, parm.bounds[parm == parm.name]$upr)
  }
  parmx
}
fn_get.parmx = function(parm.name, type = 'default', nsmp = 100){
  if(type == 'default'){
    parmx = parm.bounds %>% dplyr::filter(parm %in% parm.name) %>% dplyr::select(parm, age, default) %>% dplyr::arrange(., parm, age) 
  } else if(type == 'smp'){
    # random sampling from lower, upper bound
    lwr = parm.bounds %>% dplyr::filter(parm %in% parm.name) %>% dplyr::select(parm, age, lwr) %>% dplyr::arrange(., parm, age)
    upr = parm.bounds %>% dplyr::filter(parm %in% parm.name) %>% dplyr::select(parm, age, upr) %>% dplyr::arrange(., parm, age) 
    parmx = cbind(parm.bounds %>% dplyr::filter(parm %in% parm.name) %>% dplyr::select(parm, age) %>% dplyr::arrange(., parm, age),
                  tgp::lhs(nsmp, 
                           rect = cbind(lwr$lwr, upr$upr)) %>% t)
  }
  parmx
  
}

# function to map age grouped variable to each 1-y age
fn_map.age2finer = function(parm.grp, age.grps, age.ranges = 0:100){
  if(is.null(dim(parm.grp))){
    # the variable is just a vector, add age
    parm.grp = data.table(age.grps, value = parm.grp)
  }
  parmx = data.table(age = age.ranges)
  tmp = age.grps
  tmp$value = parm.grp$value
  tmp[nrow(tmp)]$age.end = 200
  tmp$lower = paste("`age` >=", tmp$age.start)
  tmp$upper = paste("`age` <=", tmp$age.end)
  tmp$code = tmp[, c('lower', 'upper')] %>% apply(1, paste, collapse = ' & ')
  str_age_recode = tmp[,c('code','value')] %>% apply(1, paste, collapse = ' ~ ') %>% paste(collapse = ', ')
  parmx = parmx %>% {eval(parse(text = paste('dplyr::mutate(., value = case_when(', str_age_recode,'))', sep='')))}
  
  parmx
}

# function to map 1-y age specific variable to age grouped values
fn_map.age2coarser = function(parm.grp, age.grps, age.ranges = 0:100){
  if(is.null(dim(parm.grp))){
    # the variable is just a vector, add age
    parm.grp = data.table(age = age.ranges, value = parm.grp)
  }
  parmx = data.table(age.grps %>% dplyr::select(age.start, age.end), value = 0)
  for(ia in 1:nrow(age.grps)){
    parmx[ia]$value = parm.grp %>% dplyr::filter(age >= age.grps$age.start[ia] & age <= age.grps$age.end[ia]) %>% .$value %>% mean
  }
  parmx
}

# function to sample duration of contact from a specific prob distribution, based on polymod data
fn_getContDurSmp = function(age.grp.t, # age group
                            nsmp.t, # number sample
                            cnt_home.t, # whether it is hh contact
                            ctc.fit # data table storing the distributions and parameters
){
  tmp = ctc.fit %>% filter(age.grp == age.grps[age.grp.t]$label & cnt_home ==cnt_home.t)
  tmp.parm = tmp %>% 
    dplyr::select(-da.proc, -method, -cont_freq, -wt_phys, -wt_nonphys, -wt_na, -cnt_home, -age.grp, -measure, -prob.dist,-hh_size, -n.smp, -max.cut) %>%
    dplyr::mutate(id = 1) %>% melt(.,id.vars = 'id') %>% dplyr::mutate(id = NULL) %>% filter(!is.na(value))
  
  ctc.max.t = tmp$max.cut / 24 # div by 24 to convert to per day
  tmp.parm$code = tmp.parm %>% apply(.,1,paste,collapse='=')
  code.parm = tmp.parm$code %>% paste(., collapse=', ')
  code.t = paste0('r',tmp$prob.dist,'(',nsmp.t,',',code.parm,') / 24') # div by 24 to convert to per day
  # dist.t = get(paste0('r',tmp$prob.dist))
  if(cnt_home.t){
    ctc.min.t = 0 # possible to have 0 hh ctc
  } else {
    ctc.min.t = .01 # in days, non hh ctc, slightly >0
  }
  ctc.t = eval(parse(text=code.t)) %>% round(.,4)
  ctc.t = ctc.t %>% pmax(., ctc.min.t) 
  if(!is.null(ctc.max.t))
    ctc.t = ctc.t %>% pmin(., ctc.max.t)  # make sure it is >0
  
  ctc.t
}

# do the vaccination component per weekly or daily data
fn_doVaccination = function(ppl.t, tt, da.vx.t, n.byage_ens.t){
  idx.t = NULL
  ppl.t = ppl.t %>% 
    dplyr::mutate(p.t = fn_getVaxRelProb(tt, # current time
                                                ppl.t, # information for these individuals
                                                parm.getVaxProb)) %>%
    data.table
  
  for(ia in 1:nrow(age.grps)){
    n.vx.t = da.vx.t %>% filter(age.start == age.grps$age.start[ia]) %>% .$value
    
    if(length(n.vx.t) < 1)
      next  # no vax for this age group
    
    # do it for each ens
    for(iens in 1:num_ens){
      n.t = (n.vx.t * n.byage_ens.t[ia, iens]) %>% round(., 0)
      if(n.t < 1)
        next
      
      tmp.ppl.t = ppl.t %>% 
        filter(ens == iens & age >= age.grps$age.start[ia] & age <= age.grps$age.end[ia])
      
      if(nrow(tmp.ppl.t) <= n.t){
        smp.t = tmp.ppl.t$idx
        print("doVaccination: small smp size, all get vax'ed!!!")
      } else {
        smp.t = sample(tmp.ppl.t$idx, 
                       size = n.t, prob = tmp.ppl.t$p.t, replace = F)
      }
      
      
      idx.t = append(idx.t, smp.t)
    }
  }
  
  idx.t
}

fn_format = function(x, rounddigt = 2){
  paste0(round(x[1], digits = rounddigt), ' (', round(x[2], digits = rounddigt), ', ', round(x[3], digits = rounddigt),')')
}
fn_foi0 = function(indv.info.t, # information for this individual
                   prev.pop.t, # current population level prevalence
                   # prob.tx.per.ctc0 = .5, # baseline probability of infection per infectious contact
                   prob.tx.per.ctc0 = prob.tx.per.ctc0, # baseline probability of infection per infectious contact
                   p.redn.mask = .9, # assumed % reduction in risk of prob of infection per contact
                   ave.HHctc.dur.per.day = 12/24, # average contact duration among HH members per day
                   ave.nonHHctc.dur.per.day = .5/24, # average contact duration among nonHH members per day
                   dt = 1 # time step of simulation
){ # compute force of infection
  # prev.pop: population level infection prevalence, updated based on this time step
  # frac.tm.hh: fraction of time spent at home (with HH contacts) per day, set to 1/2 default
  # ave.HHctc.dur.per.day: average contact duration among HH members per day, set to 12 hours, or 12/24
  # ave.nonHHctc.dur.per.day: average contact duration among nonHH members per day, set to .5 hour, or .5/24
  # dt: time step, set to 1 by default
  # prob
  prob.tx.per.ctc.out = prob.tx.per.ctc0 * (1 - p.redn.mask * indv.info.t$masking) # prob of tx per contact depending on whether indv using mask
  prev.hh.t = prev.pop.t * (1 - p.redn.mask * indv.info.t$masking) # assume prevalence of infection w/i household scaled with maksing practice
  # masking = 1 (use mask): prev.hh =prev.pop * (1 - p.redn.mask); 
  # masking = 0 (no mask): prev.hh = prev.pop 
  1 - exp(-(prob.tx.per.ctc0 * indv.info.t$num.HHctc * ave.HHctc.dur.per.day * dt * prev.hh.t + # HHctc, assume no mask
              prob.tx.per.ctc.out * indv.info.t$num.nonHHctc * ave.nonHHctc.dur.per.day * dt * prev.pop.t # nonHHctc
  ))
}

# allow computation of multiple individuals

# ensemble version of the model
fn_foi = function(ppl.t, # information for these individuals
                  # prev.pop.t, # current population level prevalence
                  # prob.tx.per.ctc0 = .5, # baseline probability of infection per infectious contact
                  prob.tx.per.ctc0 = prob.tx.per.ctc0, # baseline probability of infection per infectious contact
                  ave.HHctc.dur.per.day = 12/24, # average contact duration among HH members per day
                  ave.nonHHctc.dur.per.day = .5/24, # average contact duration among nonHH members per day
                  dt = 1, # time step of simulation
                  method.calHHrisk = 'like4like', #  options: 'adjHHmask', 'like4like'
                  suscept.by.age = 1, # age specific susceptibility
                  # adjust for contact reduction
                  p.ctc.redn.t = 1, # default, no change
                  adj.ctc.redn.bymask = data.table(masking = 0:1, value = 1), # default, no difference
                  UseAgeSpecSuscept = F, # whether age specific susceptibility is considered
                  UseContactDuration = F, # whether contact duration is used directly
                  MASKATHOME = F # unless to test the scenario, in general won't be the case
){ # compute force of infection
  # ppl.t: data table storing all the info for related people
  # prev.pop: population level infection prevalence, updated based on this time step
  # frac.tm.hh: fraction of time spent at home (with HH contacts) per day, set to 1/2 default
  # ave.HHctc.dur.per.day: average contact duration among HH members per day, set to 12 hours, or 12/24
  # ave.nonHHctc.dur.per.day: average contact duration among nonHH members per day, set to .5 hour, or .5/24
  # dt: time step, set to 1 by default
  # could further account for maternal immunity (so lower risk of infection for newborns)
  # prob
  
  if(UseContactDuration){
    # if contact duration is used directly
    # set these parameters to 1, as it is already accounted for
    ave.HHctc.dur.per.day = 1 # average contact duration among HH members per day
    ave.nonHHctc.dur.per.day = 1
  }
  
  prob.tx.per.ctc.out = prob.tx.per.ctc0 * (1 - ppl.t$p.redn.mask * ppl.t$masking) # prob of tx per contact depending on whether indv using mask
  # prev.hh.t = prev.pop.t * (1 - ppl.t$p.redn.mask * ppl.t$masking) # assume prevalence of infection w/i household scaled with maksing practice
  if(method.calHHrisk == 'adjHHmask'){
    # method 'adjHHmask', adjust HH prevalence based on masking practice of the HH 
    # assume prevalence of infection w/i household scaled with masking practice of the HH (multiple HHctc's, assume multiplicative risk for the HH)
    # - it reduces the prev w/i HH using mask; but dose not increase the risk of HH not using mask
    # as a result, it would likely UNDERESTIMATE the overall epidemic intensity
    prev.hh.t = ppl.t$AVEprev * (1 - (ppl.t$p.redn.mask * ppl.t$masking)^ppl.t$num.HHctc) 
    
  } else if (method.calHHrisk == 'like4like'){
    # method 'like4like', use the population level prevalence for each stratum (mask vs no mask) to represent each type of HH 
    # - it could reduce the prev w/i HH using mask and increase the risk of HH not using mask
    # might flip depending on stratum specific susceptibility, e.g. if most those who don't use mask are immune
    prev.hh.t = ppl.t$HHprev
  }
  p.suscept = 1 # default - all equally susceptible
  if(UseAgeSpecSuscept){
    p.suscept = suscept.by.age[pmin(ppl.t$age, 100) + 1, ]$value
  }
  # further risk/ctc reduction tied to masking behavior 
  p.ctc.redn.bymask = adj.ctc.redn.bymask[ppl.t$masking+1, ]$value * p.ctc.redn.t %>% pmin(., 1) # make sure it does not exceed 1, as this is reduction
  # masking = 1 (use mask): prev.hh =prev.pop * (1 - p.redn.mask); 
  # masking = 0 (no mask): prev.hh = prev.pop
  if(MASKATHOME){ 
    foi.t = 1 - exp(-(ppl.t$num.HHctc * ave.HHctc.dur.per.day * dt * prev.hh.t + # HHctc, assume also WITH mask
                        ppl.t$num.nonHHctc * ave.nonHHctc.dur.per.day * dt * ppl.t$AVEprev * # nonHHctc
                        p.ctc.redn.bymask # reduction from contact reduction, only applied to nonHHctc
    ) * prob.tx.per.ctc.out * p.suscept # further adjust for susceptibility by age
    )
  } else {
    foi.t = 1 - exp(-(prob.tx.per.ctc0 * ppl.t$num.HHctc * ave.HHctc.dur.per.day * dt * prev.hh.t + # HHctc, assume no mask
                        prob.tx.per.ctc.out * ppl.t$num.nonHHctc * ave.nonHHctc.dur.per.day * dt * ppl.t$AVEprev * # nonHHctc
                        p.ctc.redn.bymask # reduction from contact reduction, only applied to nonHHctc
    ) * p.suscept # further adjust for susceptibility by age
    )
  }
  
  
  return(foi.t)
}

# ensemble version of the model
fn_foi_ens = function(ppl.t, # information for these individuals
                      # prev.pop.t, # current population level prevalence
                      # prob.tx.per.ctc0 = .5, # baseline probability of infection per infectious contact
                      prob.tx.per.ctc0 = prob.tx.per.ctc0, # baseline probability of infection per infectious contact
                      # p.redn.mask = .9, # assumed % reduction in risk of prob of infection per contact
                      # if pass in number contacts, rather than the total converted contact duration
                      ave.HHctc.dur.per.day = 12/24, # average contact duration among HH members per day
                      ave.nonHHctc.dur.per.day = .5/24, # average contact duration among nonHH members per day
                      dt = 1, # time step of simulation
                      method.calHHrisk = 'like4like', #  options: 'adjHHmask', 'like4like'
                      suscept.by.age = 1, # age specific susceptibility
                      # adjust for contact reduction
                      p.ctc.redn.t = 1, # default, no change
                      adj.ctc.redn.bymask = data.table(masking = 0:1, value = 1), # default, no difference
                      UseAgeSpecSuscept = F, # whether age specific susceptibility is considered
                      UseContactDuration = F # whether contact duration is used directly
){ # compute force of infection
  # ppl.t: data table storing all the info for related people
  # prev.pop: population level infection prevalence, updated based on this time step
  # frac.tm.hh: fraction of time spent at home (with HH contacts) per day, set to 1/2 default
  # ave.HHctc.dur.per.day: average contact duration among HH members per day, set to 12 hours, or 12/24
  # ave.nonHHctc.dur.per.day: average contact duration among nonHH members per day, set to .5 hour, or .5/24
  # dt: time step, set to 1 by default
  # could further account for maternal immunity (so lower risk of infection for newborns)
  # prob
  
  if(UseContactDuration){
    # if contact duration is used directly
    # set these parameters to 1, as it is already accounted for
    ave.HHctc.dur.per.day = 1 # average contact duration among HH members per day
    ave.nonHHctc.dur.per.day = 1
  }
  
  prob.tx.per.ctc.out = prob.tx.per.ctc0 * (1 - ppl.t$p.redn.mask * ppl.t$masking) # prob of tx per contact depending on whether indv using mask
  # prev.hh.t = prev.pop.t * (1 - ppl.t$p.redn.mask * ppl.t$masking) # assume prevalence of infection w/i household scaled with maksing practice
  if(method.calHHrisk == 'adjHHmask'){
    # method 'adjHHmask', adjust HH prevalence based on masking practice of the HH 
    # assume prevalence of infection w/i household scaled with masking practice of the HH (multiple HHctc's, assume multiplicative risk for the HH)
    # - it reduces the prev w/i HH using mask; but dose not increase the risk of HH not using mask
    # as a result, it would likely UNDERESTIMATE the overall epidemic intensity
    prev.hh.t = ppl.t$AVEprev * (1 - (ppl.t$p.redn.mask * ppl.t$masking)^ppl.t$num.HHctc) 
    
  } else if (method.calHHrisk == 'like4like'){
    # method 'like4like', use the population level prevalence for each stratum (mask vs no mask) to represent each type of HH 
    # - it could reduce the prev w/i HH using mask and increase the risk of HH not using mask
    # might flip depending on stratum specific susceptibility, e.g. if most those who don't use mask are immune
    prev.hh.t = ppl.t$HHprev
  }
  p.suscept = 1 # default - all equally susceptible
  if(UseAgeSpecSuscept){
    p.suscept = suscept.by.age[pmin(ppl.t$age, 100) + 1, ]$value
  }
  # further risk/ctc reduction tied to masking behavior 
  p.ctc.redn.bymask = adj.ctc.redn.bymask[ppl.t$masking+1, ]$value * p.ctc.redn.t %>% pmin(., 1) # make sure it does not exceed 1, as this is reduction
  # masking = 1 (use mask): prev.hh =prev.pop * (1 - p.redn.mask); 
  # masking = 0 (no mask): prev.hh = prev.pop 
  foi.t = 1 - exp(-(prob.tx.per.ctc0 * ppl.t$num.HHctc * ave.HHctc.dur.per.day * dt * prev.hh.t + # HHctc, assume no mask
                      prob.tx.per.ctc.out * ppl.t$num.nonHHctc * ave.nonHHctc.dur.per.day * dt * ppl.t$AVEprev * # nonHHctc
                      p.ctc.redn.bymask # reduction from contact reduction, only applied to nonHHctc
  ) * p.suscept # further adjust for susceptibility by age
  )
  
  return(foi.t)
}


fn_immLossExp = function(tt, # current time
                         ppl.t, # information for these individuals
                         Trs # mean immunity period
){
  # this function computes probability of an recovred individual returning to being susceptible per an exponential function
  1 - exp(-1/Trs * (tt - ppl.t$tm.last.event))
}

# to convert the cumulative cdf to pdf 
# tested and confirmed: results from fn_calProbTable MATCHE THE CUMULATIVE!
fn_calProbTable = function(parm.tm2event, reverse.prob = F){
  
  prob.table = with(parm.tm2event, {
    # look up table for the prob:
    cdf.t = p.imm.wane.max / (1+exp(-k*(1:(Dur) - Dur * p.imm.wane.midpt))) # this cumulative
    pdf.t = c(0, cdf.t[-1] - cdf.t[-length(cdf.t)]) # density
    pdf.t = pdf.t / sum(pdf.t)
    tmp = numeric(Dur)
    for(i in 1:Dur){
      tmp[i] = 1/(prod(1-tmp[1:i])) * pdf.t[i] # adjust for those already moved, with recursive update
    }
    if(reverse.prob)
      tmp = rev(tmp)
    data.table(tm2now = 1:Dur, prob = tmp)
  })
  
  prob.table
}

fn_immLossLogistic = function(tt, # current time
                              ppl.t, # information for these individuals
                              parm.immLoss # related parms
){
  # this function computes probability of a recovered individual returning to being susceptible per a logistic function
  p = with(parm.immLoss, {
    prob.table[pmin(tt - ppl.t$tm.last.event+1, Dur),]$prob
  })
  
  p
}


# note: it is a diff situation here, cp fn_immLossLogistic
# b/c not changes in status is made, so the p.redn should be cumulative over time
fn_severity = function(tt, # current time
                       ppl.t, # information for these individuals
                       severity.by.age,  # severity measure for each age
                       parm.immLossVseverity 
){
  # to compute severity / severe cases based on immune history (# infections, vaccination, etc)
  p = with(parm.immLossVseverity, {
    p.redn.vax = (vax.p.max / (1+exp(vax.k*(tt - ppl.t$tm.last.vax - Dur.vax * vax.p.midpt)))) %>% replace_na(., 0)
    p.redn.inf = (imm.p.max / (1+exp(imm.k*(tt - ppl.t$tm.last.inf - Dur.imm * imm.p.midpt)))) %>% replace_na(., 0)
    severity.by.age[pmin(ppl.t$age, 100) + 1, ]$value * # baseline severity by age
      (1 - (ppl.t$num.vax >=1) * p.redn.vax) *  # from vaccine protection
      (1 - (ppl.t$num.inf >=1) * p.redn.inf) # from prior infection
  })
  
  p
}


if(F){ # FOR TESTING
  imm.p.max = .95 # initial VE against severe disease
  imm.k = .025  # similar to estimated using VE data
  Dur.imm = 365; imm.p.midpt = .6
  plot(imm.p.max / (1+exp(imm.k*(1:365 - Dur.imm * imm.p.midpt))))
}



fn_getVaxProbLogistic = function(tt, # current time
                                 ppl.t, # information for these individuals
                                 vax.elig.byage, # set eligibility by age
                                 vax.coverage.byage, # use coverage to represent willingness, could change over time given availability etc
                                 parm.getVaxProb
){
  # this function computes probability of an individual getting vaccinated 
  # based on time since last vaccine dose
  # And additional eligibility criteria
  # reason for "over-vaccination": individuals are checked for eligibility (asked to get vaccinated) many times before they get vaccinated
  # and assume all those times are independent
  # need a diff function for likelihood of vaccination, after the eligibility, prob to get vax go down again
  # could be periodic 
  

  # use sort of a periodic spike for this probability, except for the first year(?)
  # use a high frequency sinusoidal function? 
  p.getVax = with(as.list(parm.getVaxProb),{
    # if 1st time, tm.last.vax would be NA's
    # p1 = p1 %>% replace_na(., 1)  # set to tt.vax.start ?
    # not good
    # p1b = p.max / (1+exp(-k*((tt - (ppl.t$tm.last.vax %>% replace_na(., tt.vax.start - 250))) %% vax.interval - t2vax * p.midpt))) # this cumulative up to this day (%%: allow it to cycle back next vax period)
    # p1a = p.max / (1+exp(-k*((tt - (ppl.t$tm.last.vax %>% replace_na(., tt.vax.start - 250)) -1) %% vax.interval - t2vax * p.midpt))) # cumulative up to 1 day ago
    
    p1b = p.max / (1+exp(-k*((tt - ppl.t$tm.last.vax) %% vax.interval - t2vax * p.midpt))) # this cumulative up to this day (%%: allow it to cycle back next vax period)
    p1a = p.max / (1+exp(-k*((tt - ppl.t$tm.last.vax -1) %% vax.interval - t2vax * p.midpt))) # cumulative up to 1 day ago
    p1 = p1b - p1a # density for today
    # if 1st time, tm.last.vax would be NA's
    # p1 = p1 %>% replace_na(., 1)  # set to tt.vax.start ?
    p1 = p1 %>% replace_na(., 1/60) # suppose one gets asked 60 times for the first vaccination
    
    # additional criteria
    # e.g., vax.eligible could be based on roll out of mass vax
    # vax.willingness could be based on vax coverage data
    # p2 = p1 * ppl.t$vax.eligible * ppl.t$vax.willingness
    p2 = p1 * vax.elig.byage[pmin(ppl.t$age, 100) + 1, ]$value * vax.coverage.byage[pmin(ppl.t$age, 100) + 1, ]$value
    
    p2
  })
  
  p.getVax
}

# to get the overall vax coverage, example #ppl getting vax at time t per relative probability instead
fn_getVaxRelProb = function(tt, # current time
                            ppl.t, # information for these individuals
                            # vax.elig.byage, # set eligibility by age
                            # vax.coverage.byage, # use coverage to represent willingness, could change over time given availability etc
                            parm.getVaxProb
){
  # this function computes probability of an individual getting vaccinated 
  # based on time since last vaccine dose
  # And additional eligibility criteria
  # reason for "over-vaccination": individuals are checked for eligibility (asked to get vaccinated) many times before they get vaccinated
  # and assume all those times are independent
  # need a diff function for likelihood of vaccination, after the eligibility, prob to get vax go down again
  # could be periodic 
  
  # use sort of a periodic spike for this probability, except for the first year(?)
  # use a high frequency sinusoidal function? 
  p.getVax = with(as.list(parm.getVaxProb),{
    
    p1 = p.max / (1+exp(-k*((tt - ppl.t$tm.last.vax) - t2vax * p.midpt))) # this relative only based on last vax
    
    p1 = p1 %>% replace_na(., 1) 
    
    # additional criteria
    p2 = (tt - ppl.t$tm.last.inf) > tm.min.inf2vx # * vax.coverage.byage[pmin(ppl.t$age, 100) + 1, ]$value
    p2 = p2 %>% replace_na(., 1)
    # but not all infected know they are infected due to e.g., asymptomatic infection
    
    p1 * p2
  })
  
  p.getVax
}




# initialization - to initialize the system, based on diff assumptions/interactions

fn_initalizeABM = function(parms){
  
  rec = with(as.list(parms),
             {
               # to store individual info
               rec = data.table(idx = 1:num_pop, 
                                age = 0, # to record age of the individual (for population turnover, severity, etc.)
                                birth.year = 0, # for analysis of birth cohort
                                status = 'S', # to record disease status: S, E, I, R, V
                                tm.last.event = NA_integer_,  # to record the time of last event
                                tm.last.inf = NA_integer_,  # to record the time of last infection
                                tm.last.hosp = NA_integer_,  # to record the time of last hospitalization
                                tm.last.vax = NA_integer_, # to record the time of last vaccination
                                # time when an individual of age x become eligible to get vax?
                                num.inf = 0, # to record number of infection (for adjusting severity, prob of reinfection etc.)
                                num.hosp = 0, # to record number of hospitalization 
                                num.vax = 0, # to record number of vax (for adjusting severity, prob of reinfection etc.)
                                num.nonHHctc = 0, # to record number of contacts outside of the household
                                num.HHctc = 0, # to number of HH/close contact
                                masking = 0, # record masking status
                                p.redn.mask = 0.5 # allow difference in mask effectiveness
               )
               
               # initialize the population
               n.agegrp = nrow(age.grps)
               p.age = age.grps$age.end - age.grps$age.start + 1
               p.age = p.age / sum(p.age)
               n.byage = rmultinom(1, size = num_pop, prob = p.age) # number of ppl in each age group
               # idx.list = list(n.agegrp)
               smp.t = 1:num_pop # initial pool
               for(ia in 1:n.agegrp){
                 
                 if(n.byage[ia] < 1)
                   next
                 
                 idx.t = sample(smp.t, n.byage[ia], replace = F)
                 
                 if(length(idx.t)>0){
                   # assign age etc
                   rec[idx %in% idx.t]$age = sample(age.grps$age.start[ia]:age.grps$age.end[ia], n.byage[ia], replace = T)
                   rec[idx %in% idx.t]$birth.year = year.cur - rec[idx %in% idx.t]$age
                   
                   # assign contact
                   if(!UseContactDuration){ # set value based on number of contacts
                     non.hhctc.t = rnorm(n.byage[ia], mean = ave.num.nonHHctc[ia], sd = sd.num.nonHHctc[ia]) %>% round(., 0) %>% pmax(., 1) # at least 1
                     rec.t[ens.idx %in% idx.t]$num.nonHHctc = 
                     rec.t[ens.idx %in% idx.t]$num.HHctc = rnorm(n.byage[ia], mean = ave.num.HHctc[ia], sd = sd.num.HHctc[ia]) %>% round(., 0) %>% pmax(., ifelse(ia <=3, 1, 0)) # make sure there is at least a guardian for children
                     
                   } else {
                     # set value based on duration of contacts
                     non.hhctc.t = fn_getContDurSmp(age.grp.t = ia, nsmp.t = n.byage[ia], cnt_home.t = F, ctc.fit)
                     rec.t[ens.idx %in% idx.t]$num.nonHHctc = non.hhctc.t
                     rec.t[ens.idx %in% idx.t]$num.HHctc = fn_getContDurSmp(age.grp.t = ia, nsmp.t = n.byage[ia], cnt_home.t = T, ctc.fit)
                     
                   }
                   # assign masking practice
                   if(diff.ctc.bymask){ # assume diff ctc rate for ppl with diff masking preference
                     # scale the prob of masking with the inverse of the ctc rate
                     m.wts = (1/(non.hhctc.t + .001))  %>% pmin(., max(non.hhctc.t))  # avoid extreme values # ^2
                     m.wts = m.wts/sum(m.wts)
                     idx.mask = sample(idx.t, size = round(frac.masking[ia] * n.byage[ia], 0), prob = m.wts, replace = F)
                   } else {
                     # random 
                     idx.mask = sample(idx.t, size = round(frac.masking[ia] * n.byage[ia], 0), replace = F)
                   }
                   
                   if(length(idx.mask) >0){
                     rec[idx.mask]$masking = 1  # use mask
                   }
                  
                   # mask effectiveness
                   rec[idx.t]$p.redn.mask = rnorm(n.byage[ia], mean = ave.p.redn.mask[ia], sd = sd.p.redn.mask[ia]) %>% round(., 3) %>% pmin(., .95) %>% pmax(., 0.05) # [.05 - .95]
                   
                   # update sampling pool
                   smp.t = smp.t[! smp.t %in% idx.t] # exclude those have been picked 
                 }
                 
               }
               
               # random seeding
               idx.inf = sample(1:num_pop, size = round(frac.ini.inf * num_pop, 0), replace = F)
               if(length(idx.inf)){
                 rec$status[idx.inf] = 'I'  # exposed
                 rec$tm.last.event[idx.inf] = 0  # exposed
                 rec$tm.last.event[idx.inf]$tm.last.inf = 0 # record time of new infection
                 rec$num.inf[idx.inf] = rec$num.inf[idx.inf] + 1  # record number times infected
                 
               }
               
               mn.num.HHctc = sum(ave.num.HHctc * n.byage/num_pop)
               mn.num.nonHHctc = sum(ave.num.nonHHctc * n.byage/num_pop)
               mn.frac.masking = sum(frac.masking * n.byage/num_pop)
               mn.p.redn.mask = sum(ave.p.redn.mask * n.byage/num_pop)
               R0ini = (prob.tx.per.ctc0 * mn.num.HHctc * ave.HHctc.dur.per.day * Tir + # assume mask, HHctc
                          prob.tx.per.ctc0 * (1 - mn.p.redn.mask) * mn.num.nonHHctc * ave.nonHHctc.dur.per.day * Tir) * mn.frac.masking +  # mask, nonHHctc
                 (prob.tx.per.ctc0 * mn.num.HHctc  * ave.HHctc.dur.per.day * Tir + # HHctc, assume no mask
                    prob.tx.per.ctc0  * mn.num.nonHHctc * ave.nonHHctc.dur.per.day * Tir) * (1 - mn.frac.masking)
               
               R0ini.no.mask = (prob.tx.per.ctc0 * Tir * (mn.num.HHctc * ave.HHctc.dur.per.day + # HHctc, assume no mask
                                                            mn.num.nonHHctc * ave.nonHHctc.dur.per.day))
               
               
               return(list(rec = rec, R0ini = R0ini, R0ini.no.mask = R0ini.no.mask))
             })
  
  return(rec)
  
  if(F){
    # quality check
    ia = 7
    ppl.t = rec %>% filter(age >= age.grps$age.start[ia] & age <= age.grps$age.end[ia])
    summary(ppl.t$num.nonHHctc); ave.num.nonHHctc[ia]
    summary(ppl.t$num.HHctc); ave.num.HHctc[ia]
    summary(ppl.t$p.redn.mask); ave.p.redn.mask[ia]
  }
  
}

fn_initalizeABM_ens = function(parms){
  
  rec = with(as.list(parms),
             {
               # to store individual info
               rec = data.table(ens =rep(1:num_ens, e=num_pop),
                                ens.idx = rep(1:num_pop, num_ens), 
                                idx = 1:(num_ens * num_pop), # global idx
                                age = NA_integer_, # to record age of the individual (for population turnover, severity, etc.)
                                birth.year = NA_integer_, # for analysis of birth cohort
                                status = 'S', # to record disease status: S, E, I, R, V
                                tm.last.event = NA_integer_,  # to record the time of last event
                                tm.last.inf = NA_integer_,  # to record the time of last infection
                                tm.last.hosp = NA_integer_,  # to record the time of last hospitalization
                                tm.last.vax = NA_integer_, # to record the time of last vaccination
                                # time when an individual of age x become eligible to get vax?
                                num.inf = 0, # to record number of infection (for adjusting severity, prob of reinfection etc.)
                                num.hosp = 0, # to record number of hospitalization 
                                num.vax = 0, # to record number of vax (for adjusting severity, prob of reinfection etc.)
                                num.nonHHctc = 0, # to record number (or duration) of contacts outside of the household
                                num.HHctc = 0, # to number (or duration) of HH/close contact
                                masking = 0, # record masking status
                                p.redn.mask = 0.5 # allow difference in mask effectiveness
               )
               
               # initialize the population
               n.agegrp = nrow(age.grps)
               p.age = age.grps$age.end - age.grps$age.start + 1
               p.age = p.age / sum(p.age)
               n.byage_ens = rmultinom(num_ens, size = num_pop, prob = p.age) # number of ppl in each age group
               # idx.list = list(n.agegrp)
               
               
               
               R0ini = NULL
               R0ini.no.mask = NULL
               
               # need to do it one by one
               for(ii in 1:num_ens){
                 rec.t = rec[(ii-1)*num_pop+1:num_pop, ]
                 smp.t = 1:num_pop # initial pool
                 n.byage = n.byage_ens[,ii]
                 for(ia in 1:n.agegrp){
                   
                   if(n.byage[ia] < 1)
                     next
                   
                   idx.t = sample(smp.t, n.byage[ia], replace = F)
                   
                   if(length(idx.t) > 0){
                     # assign age etc
                     rec.t[ens.idx %in% idx.t]$age = sample(age.grps$age.start[ia]:age.grps$age.end[ia], n.byage[ia], replace = T)
                     rec.t[ens.idx %in% idx.t]$birth.year = year.cur - rec.t[ens.idx %in% idx.t]$age
                     
                     # assign contact
                     if(!UseContactDuration){ # set value based on number of contacts
                       non.hhctc.t = rnorm(n.byage[ia], mean = ave.num.nonHHctc[ia], sd = sd.num.nonHHctc[ia]) %>% round(., 0) %>% pmax(., 1) # at least 1
                       rec.t[ens.idx %in% idx.t]$num.nonHHctc = non.hhctc.t
                       rec.t[ens.idx %in% idx.t]$num.HHctc = rnorm(n.byage[ia], mean = ave.num.HHctc[ia], sd = sd.num.HHctc[ia]) %>% round(., 0) %>% pmax(., ifelse(ia <=3, 1, 0)) # make sure there is at least a guardian for children
                       
                     } else {
                       # set value based on duration of contacts
                       non.hhctc.t = fn_getContDurSmp(age.grp.t = ia, nsmp.t = n.byage[ia], cnt_home.t = F, ctc.fit)
                       rec.t[ens.idx %in% idx.t]$num.nonHHctc = non.hhctc.t
                       rec.t[ens.idx %in% idx.t]$num.HHctc = fn_getContDurSmp(age.grp.t = ia, nsmp.t = n.byage[ia], cnt_home.t = T, ctc.fit)
                       
                     }
                     # assign masking practice
                     if(diff.ctc.bymask){ # assume diff ctc rate for ppl with diff masking preference
                       # scale the prob of masking with the inverse of the ctc rate
                       m.wts = (1/(non.hhctc.t + .001)) %>% pmin(., max(non.hhctc.t))  # avoid extreme values # ^2
                       m.wts = m.wts/sum(m.wts)
                       idx.mask = sample(idx.t, size = round(frac.masking[ia] * n.byage[ia], 0), prob = m.wts, replace = F)
                     } else {
                       # random 
                       idx.mask = sample(idx.t, size = round(frac.masking[ia] * n.byage[ia], 0), replace = F)
                     }
                     
                     if(length(idx.mask) > 0){
                       rec.t[ens.idx %in% idx.mask]$masking = 1  # use mask
                     }
                     
                     # mask effectiveness
                     rec.t[ens.idx %in% idx.t]$p.redn.mask = rnorm(n.byage[ia], mean = ave.p.redn.mask[ia], sd = sd.p.redn.mask[ia]) %>% round(., 3) %>% pmin(., .95) %>% pmax(., 0.05) # [.05 - .95]
                     
                     # update sampling pool
                     smp.t = smp.t[! smp.t %in% idx.t] # exclude those have been picked 
                   }
                   
                 }
                 # random seeding
                 idx.inf = sample(1:num_pop, size = round(frac.ini.inf * num_pop, 0), replace = F)
                 if(length(idx.inf) > 0){
                   rec.t[ens.idx %in% idx.inf]$status = 'I'  # exposed
                   rec.t[ens.idx %in% idx.inf]$tm.last.event = 0  # exposed
                   rec.t[ens.idx %in% idx.inf]$tm.last.inf = 0 # record time of new infection
                   rec.t[ens.idx %in% idx.inf]$num.inf = rec.t[ens.idx %in% idx.inf]$num.inf + 1  # record number times infected
                 }
                 
                 mn.num.HHctc = rec.t$num.HHctc %>% mean
                 mn.num.nonHHctc = rec.t$num.nonHHctc %>% mean
                 mn.frac.masking = rec.t$masking %>% mean
                 mn.p.redn.mask = rec.t$p.redn.mask %>% mean
                 if(!UseContactDuration){
                   
                   R0ini.t = (prob.tx.per.ctc0 * mn.num.HHctc * ave.HHctc.dur.per.day * Tir + # assume mask, HHctc
                                prob.tx.per.ctc0 * (1 - mn.p.redn.mask) * mn.num.nonHHctc * ave.nonHHctc.dur.per.day * Tir) * mn.frac.masking +  # mask, nonHHctc
                     (prob.tx.per.ctc0 * mn.num.HHctc  * ave.HHctc.dur.per.day * Tir + # HHctc, assume no mask
                        prob.tx.per.ctc0  * mn.num.nonHHctc * ave.nonHHctc.dur.per.day * Tir) * (1 - mn.frac.masking)
                   
                   R0ini.no.mask.t = (prob.tx.per.ctc0 * Tir * (mn.num.HHctc * ave.HHctc.dur.per.day + # HHctc, assume no mask
                                                                  mn.num.nonHHctc * ave.nonHHctc.dur.per.day))
                   
                 } else {
                   R0ini.t = prob.tx.per.ctc0 * Tir * ((mn.num.HHctc + # assume mask, HHctc
                                                          (1 - mn.p.redn.mask) * mn.num.nonHHctc) * mn.frac.masking +  # mask, nonHHctc
                                                         (mn.num.HHctc + # HHctc, assume no mask
                                                            mn.num.nonHHctc) * (1 - mn.frac.masking))
                   
                   R0ini.no.mask.t = (prob.tx.per.ctc0 * Tir * (mn.num.HHctc + # HHctc, assume no mask
                                                                  mn.num.nonHHctc))
                 }
                 
                 rec[(ii-1)*num_pop+1:num_pop,] =  rec.t
                 R0ini = c(R0ini, R0ini.t)
                 R0ini.no.mask = c(R0ini.no.mask, R0ini.no.mask.t)
               }
               
               
               return(list(rec = rec, R0ini = R0ini, R0ini.no.mask = R0ini.no.mask, n.byage_ens = n.byage_ens))
             })
  
  return(rec)
  
  if(F){
    # quality check
    ia = 7
    ppl.t = rec %>% filter(age >= age.grps$age.start[ia] & age <= age.grps$age.end[ia])
    summary(ppl.t$num.nonHHctc); ave.num.nonHHctc[ia]
    summary(ppl.t$num.HHctc); ave.num.HHctc[ia]
    summary(ppl.t$p.redn.mask); ave.p.redn.mask[ia]
    
    rec1 = rec$rec %>% {eval(parse(text = paste('dplyr::mutate(., age.grp = case_when(', str_age_recode, ', T ~ NA_character_))', sep='')))} 
    rec1 %>% group_by(ens, age.grp, masking) %>% # by masking
      summarise(num.nonHHctc.mn = num.nonHHctc %>% mean,num.nonHHctc.median = num.nonHHctc %>% median, num.nonHHctc.lwr = num.nonHHctc %>% quantile(., prob = .05), num.nonHHctc.upr = num.nonHHctc %>% quantile(., prob = .95),
                num.HHctc.mn = num.HHctc %>% mean,num.HHctc.median = num.HHctc %>% median) %>% view
    tmp = rec$rec %>% group_by(age, masking) %>%
      summarise(num.nonHHctc = mean(num.nonHHctc))
  }
  
}

# this works for ensemble too
fn_initalizeABMnewborns = function(parms, rec, idx.newborn){
  
  rec = with(as.list(parms),
             {
               rec[idx %in% idx.newborn]$age = 0 %>% as.integer()
               rec[idx %in% idx.newborn]$birth.year = year.cur
               rec[idx %in% idx.newborn]$status = 'S'
               rec[idx %in% idx.newborn]$tm.last.event = NA_integer_
               rec[idx %in% idx.newborn]$tm.last.inf = NA_integer_
               rec[idx %in% idx.newborn]$tm.last.hosp = NA_integer_
               rec[idx %in% idx.newborn]$tm.last.vax = NA_integer_
               rec[idx %in% idx.newborn]$num.inf = 0
               rec[idx %in% idx.newborn]$num.hosp = 0
               rec[idx %in% idx.newborn]$num.vax = 0
               rec[idx %in% idx.newborn]$masking = 0  # no mask for newborns
               # mask effectiveness
               rec[idx %in% idx.newborn]$p.redn.mask = 0
               
               # assign contact
               if(!UseContactDuration){ # set value based on number of contacts
                 rec[idx %in% idx.newborn]$num.nonHHctc = rnorm(length(idx.newborn), mean = ave.num.nonHHctc[1], sd = sd.num.nonHHctc[1]) %>% round(., 0) %>% pmax(., 1) # at least 1
                 rec[idx %in% idx.newborn]$num.HHctc = rnorm(length(idx.newborn), mean = ave.num.HHctc[1], sd = sd.num.HHctc[1]) %>% round(., 0) %>% pmax(., 1) # make sure there is at least a guardian for children
                 
                 } else {
                 # set value based on duration of contacts
                 rec[idx %in% idx.newborn]$num.nonHHctc = fn_getContDurSmp(age.grp.t = 1, nsmp.t = length(idx.newborn), cnt_home.t = F, ctc.fit)
                 rec[idx %in% idx.newborn]$num.HHctc = fn_getContDurSmp(age.grp.t = 1, nsmp.t = length(idx.newborn), cnt_home.t = T, ctc.fit)
                 
               }
               
               rec
             })
  
  return(rec)
}


# ok for the ensemble too
fn_pop.aging = function(parms, rec, NoMaskSwitch=T){
  # update age etc
  # for those change age group, also update contact info
  # NoMaskSwitch=T: same masking behavior over lifetime
  rec = with(as.list(parms),{
    rec$age = rec$age + 1 # everyone get older by 1 year
    for(ia in 2:nrow(age.grps)){
      idx.t = rec %>% filter(age == age.grps$age.start[ia]) %>% .$idx # new folks who just entering a new age group
      n.t = length(idx.t)
      
      if(n.t < 1){
        print(paste('!!! pop.aging: no one aged', age.grps$age.start[ia]))
        next
      }
        
      # update info
      # assign contact
      # rec[idx %in% idx.t]$num.nonHHctc = rnorm(n.t, mean = ave.num.nonHHctc[ia], sd = sd.num.nonHHctc[ia]) %>% round(., 0) %>% pmax(., 1) # at least 1
      # rec[idx %in% idx.t]$num.HHctc = rnorm(n.t, mean = ave.num.HHctc[ia], sd = sd.num.HHctc[ia]) %>% round(., 0) %>% pmax(., ifelse(ia <=3, 1, 0)) # make sure there is at least a guardian for children
      
      # to keep the same ctc pattern among those prefer masking, set ctc using the same ranking 
      if(diff.ctc.bymask){
        # rank ctc 
        if(ia > 2){ # i.e. when masking preference is already assigned, and thus need to match with that
          non.hhctc.last = rec[idx %in% idx.t]$num.nonHHctc  # ctc before the update
          ctc.order = non.hhctc.last %>% order()
        }
      }
      
      if(!UseContactDuration){ # set value based on number of contacts
        non.hhctc.t = rnorm(n.t, mean = ave.num.nonHHctc[ia], sd = sd.num.nonHHctc[ia]) %>% round(., 0) %>% pmax(., 1) # at least 1
        if(diff.ctc.bymask & ia > 2){
          non.hhctc.t = non.hhctc.t %>% sort # first sort it
          non.hhctc.t = non.hhctc.t[ctc.order]  # use the same order 
          # tmp = cbind(old = non.hhctc.last, new = non.hhctc.t)
        }
        rec[idx %in% idx.t]$num.nonHHctc = non.hhctc.t
        rec[idx %in% idx.t]$num.HHctc = rnorm(n.t, mean = ave.num.HHctc[ia], sd = sd.num.HHctc[ia]) %>% round(., 0) %>% pmax(., ifelse(ia <=3, 1, 0)) # make sure there is at least a guardian for children
        
      } else {
        # set value based on duration of contacts
        non.hhctc.t = fn_getContDurSmp(age.grp.t = ia, nsmp.t = n.t, cnt_home.t = F, ctc.fit)
        if(diff.ctc.bymask & ia > 2){
          non.hhctc.t = non.hhctc.t %>% sort # first sort it
          non.hhctc.t = non.hhctc.t[ctc.order]  # use the same order 
        }
        rec[idx %in% idx.t]$num.nonHHctc = non.hhctc.t 
        rec[idx %in% idx.t]$num.HHctc = fn_getContDurSmp(age.grp.t = ia, nsmp.t = n.t, cnt_home.t = T, ctc.fit)
        
      }
      
      # assign masking practice
      # to keep it consistent over a lifetime, only do it when entering grp2 (when masking is possible)
      if(!NoMaskSwitch | ia ==2){
        rec[idx %in% idx.t]$masking = 0 # set all to 0 first
        
        if(diff.ctc.bymask){ # assume diff ctc rate for ppl with diff masking preference
          # scale the prob of masking with the inverse of the ctc rate
          m.wts = (1/(non.hhctc.t + .001)) %>% pmin(., max(non.hhctc.t))  # avoid extreme values
          m.wts = m.wts/sum(m.wts)
          idx.mask = sample(idx.t, size = round(frac.masking[ia] * length(idx.t), 0), prob = m.wts, replace = F)
        } else {
          # random 
          idx.mask = sample(idx.t, size = round(frac.masking[ia] * length(idx.t), 0), replace = F)
        }
        
        if(length(idx.mask) > 0)
          rec[idx %in% idx.mask]$masking = 1  # use mask
      }
      
      # mask effectiveness
      rec[idx %in% idx.t]$p.redn.mask = rnorm(n.t, mean = ave.p.redn.mask[ia], sd = sd.p.redn.mask[ia]) %>% round(., 3) %>% pmin(., .95) %>% pmax(., 0.05) # [.05 - .95]
    } # end this age group
    rec
  })
  
  rec
}

fn_getProbDeath = function(ppl.t, life.expectancy = 75, age.death0 = 64){
  # compute probability of death based on age
  1 - exp(-2.358 + 6.556e-02 * ppl.t$age - 4.582e-04 * ppl.t$age^2) # based on fit to the life table
}

# get stats for this year
# need to keep a record of last year
fn_getStatYrly = function(rec, rec.last, idx.t){
  # record detailed stats for each year
  # NOTE: the tallies here are based on all those 'alive' (i.e., will miss those pass and get replace during population renewal)
  # for the total including the deceased, use the res.daily and res.yrly records
  
  rec1 = rec %>% dplyr::select(-status, -tm.last.event, -tm.last.vax) %>% # , -num.nonHHctc, -num.HHctc
    left_join(., y = rec.last %>% dplyr::select(-age), 
              by = c('ens', 'ens.idx', 'idx', 'birth.year'), suffix = c('','.t0')) %>%
    dplyr::mutate(num.inf.t = case_when(!is.na(num.inf.t0) ~ num.inf - num.inf.t0,
           is.na(num.inf.t0) & (tm.last.inf %in% idx.t) ~ 1, # in case the individual get replaced
           T ~ 0),
           num.hosp.t = case_when(!is.na(num.hosp.t0) ~ num.hosp - num.hosp.t0,
                                  is.na(num.hosp.t0) & (tm.last.hosp %in% idx.t) ~ 1, # in case the individual age replaced
                                  T ~ 0)
             ) %>% # record whether a new infection/hosp occurs this year
    dplyr::mutate(mask.grp = case_when(masking == 1 & (p.redn.mask <mask.E.cut[1]) ~ 'mask.lowE',
                                masking == 1 & (p.redn.mask >=mask.E.cut[1] & p.redn.mask < mask.E.cut[2]) ~ 'mask.midE',
                                masking == 1 & (p.redn.mask >=mask.E.cut[2]) ~ 'mask.highE',
                                T ~ 'no mask')) %>%
    {eval(parse(text = paste('dplyr::mutate(., age.grp = case_when(', str_age_recode, ', T ~ NA_character_))', sep='')))} %>%
    data.table()
  
  base.hosp = 1e5;  # hospitalization per 100K
  base.inf = 100; # infection, in %
  
  stat.cum = rbind(rec1 %>% # combine all, for each ensemble member
                group_by(ens) %>%
                summarise(totAR = (sum(num.inf) / n() * base.inf)  %>% round(., 3), 
                          indvAR1plus = (sum(num.inf>=1) / n() * base.inf) %>% round(., 3),
                          indvAR2plus = (sum(num.inf>=2) / n() * base.inf) %>% round(., 3),
                          indvAR3plus = (sum(num.inf>=3) / n() * base.inf) %>% round(., 3),
                          # for severe outcomes
                          totHospRate = (sum(num.hosp) / n() * base.hosp) %>% round(., 1), 
                          indvHospRate1plus = (sum(num.hosp>=1) / n() * base.hosp) %>% round(., 1),
                          indvHospRate2plus = (sum(num.hosp>=2) / n() * base.hosp) %>% round(., 1),
                          indvHospRate3plus = (sum(num.hosp>=3) / n() * base.hosp) %>% round(., 1)
                ) %>%
                dplyr::mutate(mask.grp = 'all',
                       age.grp = 'all') %>%
                ungroup(),
              rec1 %>% group_by(ens, masking) %>% # all age groups combined, by masking
                summarise(totAR = (sum(num.inf) / n() * base.inf)  %>% round(., 3), 
                          indvAR1plus = (sum(num.inf>=1) / n() * base.inf) %>% round(., 3),
                          indvAR2plus = (sum(num.inf>=2) / n() * base.inf) %>% round(., 3),
                          indvAR3plus = (sum(num.inf>=3) / n() * base.inf) %>% round(., 3),
                          # for severe outcomes
                          totHospRate = (sum(num.hosp) / n() * base.hosp) %>% round(., 1), 
                          indvHospRate1plus = (sum(num.hosp>=1) / n() * base.hosp) %>% round(., 1),
                          indvHospRate2plus = (sum(num.hosp>=2) / n() * base.hosp) %>% round(., 1),
                          indvHospRate3plus = (sum(num.hosp>=3) / n() * base.hosp) %>% round(., 1)
                ) %>%
                dplyr::mutate(mask.grp = factor(masking, levels = c(0, 1), labels = c('no mask', 'use mask')),
                       age.grp = 'all', masking = NULL) %>%
                ungroup(),
              rec1 %>% group_by(ens, age.grp) %>% # all mask groups combined, by age
                summarise(totAR = (sum(num.inf) / n() * base.inf)  %>% round(., 3), 
                          indvAR1plus = (sum(num.inf>=1) / n() * base.inf) %>% round(., 3),
                          indvAR2plus = (sum(num.inf>=2) / n() * base.inf) %>% round(., 3),
                          indvAR3plus = (sum(num.inf>=3) / n() * base.inf) %>% round(., 3),
                          # for severe outcomes
                          totHospRate = (sum(num.hosp) / n() * base.hosp) %>% round(., 1), 
                          indvHospRate1plus = (sum(num.hosp>=1) / n() * base.hosp) %>% round(., 1),
                          indvHospRate2plus = (sum(num.hosp>=2) / n() * base.hosp) %>% round(., 1),
                          indvHospRate3plus = (sum(num.hosp>=3) / n() * base.hosp) %>% round(., 1)
                ) %>%
                dplyr::mutate(mask.grp = 'all', masking = NULL) %>%
                ungroup(),
              rec1 %>% group_by(ens, age.grp, mask.grp) %>% # by age group, mask group, 
                summarise(totAR = (sum(num.inf) / n() * base.inf)  %>% round(., 3),  
                          indvAR1plus = (sum(num.inf>=1) / n() * base.inf) %>% round(., 3),
                          indvAR2plus = (sum(num.inf>=2) / n() * base.inf) %>% round(., 3),
                          indvAR3plus = (sum(num.inf>=3) / n() * base.inf) %>% round(., 3),
                          # for severe outcomes
                          totHospRate = (sum(num.hosp) / n() * base.hosp) %>% round(., 1), 
                          indvHospRate1plus = (sum(num.hosp>=1) / n() * base.hosp) %>% round(., 1),
                          indvHospRate2plus = (sum(num.hosp>=2) / n() * base.hosp) %>% round(., 1),
                          indvHospRate3plus = (sum(num.hosp>=3) / n() * base.hosp) %>% round(., 1)
                ) %>% ungroup()) %>% 
    # somehow did not work on the cluster...
    # dplyr::mutate(mask.grp = factor(mask.grp, levels = c('all', 'no mask', 'use mask', paste0('mask.', c('lowE', 'midE', 'highE')))),
    #               age.grp = factor(age.grp, levels = c('all', age.grps$label))) %>%
    dplyr::mutate(timeframe = 'cum2date') %>% 
    data.table::setcolorder(., c('ens', 'timeframe', 'mask.grp', 'age.grp', 'totAR', 'indvAR1plus', 'indvAR2plus', 'indvAR3plus', 
                     'totHospRate', 'indvHospRate1plus', 'indvHospRate2plus', 'indvHospRate3plus')) # %>%
    # dplyr::arrange(., ens, age.grp, mask.grp) 
  
  
  # count infections this year, based on timing of the lastest infection
  stat.this = rbind(rec1 %>% # combine all, for each ensemble member
                 group_by(ens) %>%
                 summarise(totAR = (sum(num.inf.t) / n() * base.inf)  %>% round(., 3), 
                           indvAR1plus = (sum(num.inf.t>=1) / n() * base.inf) %>% round(., 3),
                           indvAR2plus = (sum(num.inf.t>=2) / n() * base.inf) %>% round(., 3),
                           indvAR3plus = (sum(num.inf.t>=3) / n() * base.inf) %>% round(., 3),
                           # for severe outcomes
                           totHospRate = (sum(num.hosp.t) / n() * base.hosp) %>% round(., 1), 
                           indvHospRate1plus = (sum(num.hosp.t>=1) / n() * base.hosp) %>% round(., 1),
                           indvHospRate2plus = (sum(num.hosp.t>=2) / n() * base.hosp) %>% round(., 1),
                           indvHospRate3plus = (sum(num.hosp.t>=3) / n() * base.hosp) %>% round(., 1)
                 ) %>%
                 dplyr::mutate(mask.grp = 'all',
                        age.grp = 'all') %>%
                 ungroup(),
               rec1 %>% group_by(ens, masking) %>% # all age groups combined, by masking
                 summarise(totAR = (sum(num.inf.t) / n() * base.inf)  %>% round(., 3), 
                           indvAR1plus = (sum(num.inf.t>=1) / n() * base.inf) %>% round(., 3),
                           indvAR2plus = (sum(num.inf.t>=2) / n() * base.inf) %>% round(., 3),
                           indvAR3plus = (sum(num.inf.t>=3) / n() * base.inf) %>% round(., 3),
                           # for severe outcomes
                           totHospRate = (sum(num.hosp.t) / n() * base.hosp) %>% round(., 1), 
                           indvHospRate1plus = (sum(num.hosp.t>=1) / n() * base.hosp) %>% round(., 1),
                           indvHospRate2plus = (sum(num.hosp.t>=2) / n() * base.hosp) %>% round(., 1),
                           indvHospRate3plus = (sum(num.hosp.t>=3) / n() * base.hosp) %>% round(., 1)
                 ) %>%
                 dplyr::mutate(mask.grp = factor(masking, levels = c(0, 1), labels = c('no mask', 'use mask')),
                        age.grp = 'all', masking = NULL) %>%
                 ungroup(),
               rec1 %>% group_by(ens, age.grp) %>% # all mask groups combined, by age
                 summarise(totAR = (sum(num.inf.t) / n() * base.inf)  %>% round(., 3), 
                           indvAR1plus = (sum(num.inf.t>=1) / n() * base.inf) %>% round(., 3),
                           indvAR2plus = (sum(num.inf.t>=2) / n() * base.inf) %>% round(., 3),
                           indvAR3plus = (sum(num.inf.t>=3) / n() * base.inf) %>% round(., 3),
                           # for severe outcomes
                           totHospRate = (sum(num.hosp.t) / n() * base.hosp) %>% round(., 1), 
                           indvHospRate1plus = (sum(num.hosp.t>=1) / n() * base.hosp) %>% round(., 1),
                           indvHospRate2plus = (sum(num.hosp.t>=2) / n() * base.hosp) %>% round(., 1),
                           indvHospRate3plus = (sum(num.hosp.t>=3) / n() * base.hosp) %>% round(., 1)
                 ) %>%
                 dplyr::mutate(mask.grp = 'all', masking = NULL) %>%
                 ungroup(),
               rec1 %>% group_by(ens, age.grp, mask.grp) %>% # by age group, mask group, 
                 summarise(totAR = (sum(num.inf.t) / n() * base.inf)  %>% round(., 3),  
                           indvAR1plus = (sum(num.inf.t>=1) / n() * base.inf) %>% round(., 3),
                           indvAR2plus = (sum(num.inf.t>=2) / n() * base.inf) %>% round(., 3),
                           indvAR3plus = (sum(num.inf.t>=3) / n() * base.inf) %>% round(., 3),
                           # for severe outcomes
                           totHospRate = (sum(num.hosp.t) / n() * base.hosp) %>% round(., 1), 
                           indvHospRate1plus = (sum(num.hosp.t>=1) / n() * base.hosp) %>% round(., 1),
                           indvHospRate2plus = (sum(num.hosp.t>=2) / n() * base.hosp) %>% round(., 1),
                           indvHospRate3plus = (sum(num.hosp.t>=3) / n() * base.hosp) %>% round(., 1)
                 ) %>% ungroup()) %>% 
    # somehow did not work on the cluster...
    # dplyr::mutate(mask.grp = factor(mask.grp, levels = c('all', 'no mask', 'use mask', paste0('mask.', c('lowE', 'midE', 'highE')))),
    #               age.grp = factor(age.grp, levels = c('all', age.grps$label))) %>%
    dplyr::mutate(timeframe = 'this.year') %>%  
  data.table::setcolorder(., c('ens', 'timeframe', 'mask.grp', 'age.grp', 'totAR', 'indvAR1plus', 'indvAR2plus', 'indvAR3plus', 
                               'totHospRate', 'indvHospRate1plus', 'indvHospRate2plus', 'indvHospRate3plus')) # %>%
    # dplyr::arrange(., ens, age.grp, mask.grp)
  
  return(list(stat.cum = stat.cum, stat.this = stat.this))
}

fn_getStatYrlyByCohort = function(rec, rec.last, idx.t){
  # record detailed year 1 stats
  rec1 = rec %>% dplyr::select(-status, -tm.last.event, -tm.last.vax) %>% # , -num.nonHHctc, -num.HHctc
    left_join(., y = rec.last %>% dplyr::select(-age), 
              by = c('ens', 'ens.idx', 'idx', 'birth.year'), suffix = c('','.t0')) %>%
    dplyr::mutate(num.inf.t = case_when(!is.na(num.inf.t0) ~ num.inf - num.inf.t0,
                                 is.na(num.inf.t0) & (tm.last.inf %in% idx.t) ~ 1, # in case the individual age replaced
                                 T ~ 0),
           num.hosp.t = case_when(!is.na(num.hosp.t0) ~ num.hosp - num.hosp.t0,
                                  is.na(num.hosp.t0) & (tm.last.hosp %in% idx.t) ~ 1, # in case the individual age replaced
                                  T ~ 0)
    ) %>% # record whether a new infection/hosp occurs this year
    dplyr::mutate(mask.grp = case_when(masking == 1 & (p.redn.mask <mask.E.cut[1]) ~ 'mask.lowE',
                                masking == 1 & (p.redn.mask >=mask.E.cut[1] & p.redn.mask < mask.E.cut[2]) ~ 'mask.midE',
                                masking == 1 & (p.redn.mask >=mask.E.cut[2]) ~ 'mask.highE',
                                T ~ 'no mask')) %>%
    {eval(parse(text = paste('dplyr::mutate(., age.grp = case_when(', str_age_recode, ', T ~ NA_character_))', sep='')))}  %>%
    data.table()
  
  
  base.hosp = 1e5;  # hospitalization per 100K
  base.inf = 100; # infection, in %
  
  stat.cum = rbind(rec1 %>% 
                group_by(ens, birth.year, age) %>% # all mask groups combined
                summarise(totAR = (sum(num.inf) / n() * base.inf)  %>% round(., 3), 
                          indvAR1plus = (sum(num.inf>=1) / n() * base.inf) %>% round(., 3),
                          indvAR2plus = (sum(num.inf>=2) / n() * base.inf) %>% round(., 3),
                          indvAR3plus = (sum(num.inf>=3) / n() * base.inf) %>% round(., 3),
                          # for severe outcomes
                          totHospRate = (sum(num.hosp) / n() * base.hosp) %>% round(., 1), 
                          indvHospRate1plus = (sum(num.hosp>=1) / n() * base.hosp) %>% round(., 1),
                          indvHospRate2plus = (sum(num.hosp>=2) / n() * base.hosp) %>% round(., 1),
                          indvHospRate3plus = (sum(num.hosp>=3) / n() * base.hosp) %>% round(., 1)
                ) %>%
                dplyr::mutate(mask.grp = 'all') %>%
                ungroup(),
              rec1 %>% group_by(ens, birth.year, age, masking) %>% # by masking
                summarise(totAR = (sum(num.inf) / n() * base.inf)  %>% round(., 3), 
                          indvAR1plus = (sum(num.inf>=1) / n() * base.inf) %>% round(., 3),
                          indvAR2plus = (sum(num.inf>=2) / n() * base.inf) %>% round(., 3),
                          indvAR3plus = (sum(num.inf>=3) / n() * base.inf) %>% round(., 3),
                          # for severe outcomes
                          totHospRate = (sum(num.hosp) / n() * base.hosp) %>% round(., 1), 
                          indvHospRate1plus = (sum(num.hosp>=1) / n() * base.hosp) %>% round(., 1),
                          indvHospRate2plus = (sum(num.hosp>=2) / n() * base.hosp) %>% round(., 1),
                          indvHospRate3plus = (sum(num.hosp>=3) / n() * base.hosp) %>% round(., 1)
                ) %>%
                dplyr::mutate(mask.grp = factor(masking, levels = c(0, 1), labels = c('no mask', 'use mask')),
                       masking = NULL) %>%
                ungroup()
  )  %>% data.table() %>% 
    dplyr::mutate(timeframe = 'cum2date') %>%
    data.table::setcolorder(., c('ens', 'timeframe','birth.year', 'age', 'mask.grp', 'totAR', 'indvAR1plus', 'indvAR2plus', 'indvAR3plus', 
                     'totHospRate', 'indvHospRate1plus', 'indvHospRate2plus', 'indvHospRate3plus')) # %>%
    # dplyr::arrange(., ens, birth.year, age, mask.grp)
  
  stat.this = rbind(rec1 %>% 
                 group_by(ens, birth.year, age) %>% # all mask groups combined
                 summarise(totAR = (sum(num.inf.t) / n() * base.inf)  %>% round(., 3), 
                           indvAR1plus = (sum(num.inf.t>=1) / n() * base.inf) %>% round(., 3),
                           indvAR2plus = (sum(num.inf.t>=2) / n() * base.inf) %>% round(., 3),
                           indvAR3plus = (sum(num.inf.t>=3) / n() * base.inf) %>% round(., 3),
                           # for severe outcomes
                           totHospRate = (sum(num.hosp.t) / n() * base.hosp) %>% round(., 1), 
                           indvHospRate1plus = (sum(num.hosp.t>=1) / n() * base.hosp) %>% round(., 1),
                           indvHospRate2plus = (sum(num.hosp.t>=2) / n() * base.hosp) %>% round(., 1),
                           indvHospRate3plus = (sum(num.hosp.t>=3) / n() * base.hosp) %>% round(., 1)
                 ) %>%
                 dplyr::mutate(mask.grp = 'all') %>%
                 ungroup(),
               rec1 %>% group_by(ens, birth.year, age, masking) %>% # by masking
                 summarise(totAR = (sum(num.inf.t) / n() * base.inf)  %>% round(., 3), 
                           indvAR1plus = (sum(num.inf.t>=1) / n() * base.inf) %>% round(., 3),
                           indvAR2plus = (sum(num.inf.t>=2) / n() * base.inf) %>% round(., 3),
                           indvAR3plus = (sum(num.inf.t>=3) / n() * base.inf) %>% round(., 3),
                           # for severe outcomes
                           totHospRate = (sum(num.hosp.t) / n() * base.hosp) %>% round(., 1), 
                           indvHospRate1plus = (sum(num.hosp.t>=1) / n() * base.hosp) %>% round(., 1),
                           indvHospRate2plus = (sum(num.hosp.t>=2) / n() * base.hosp) %>% round(., 1),
                           indvHospRate3plus = (sum(num.hosp.t>=3) / n() * base.hosp) %>% round(., 1)
                 ) %>%
                 dplyr::mutate(mask.grp = factor(masking, levels = c(0, 1), labels = c('no mask', 'use mask')),
                        masking = NULL) %>%
                 ungroup()
  )  %>% data.table() %>% 
    dplyr::mutate(timeframe = 'this.year') %>%
    data.table::setcolorder(., c('ens', 'timeframe','birth.year', 'age', 'mask.grp', 'totAR', 'indvAR1plus', 'indvAR2plus', 'indvAR3plus', 
                     'totHospRate', 'indvHospRate1plus', 'indvHospRate2plus', 'indvHospRate3plus')) # %>%
    # dplyr::arrange(., ens, birth.year, age, mask.grp)
  
  return(list(stat.cum = stat.cum, stat.this = stat.this))
}


fn_getStat = function(rec){
  # NOTE: the tallies here are based on all those 'alive' (i.e., will miss those pass and get replace during population renewal)
  # for the total including the deceased, use the res.daily records
  rec1 = rec %>% dplyr::mutate(mask.grp = case_when(masking == 1 & (p.redn.mask <mask.E.cut[1]) ~ 'mask.lowE',
                                                   masking == 1 & (p.redn.mask >=mask.E.cut[1] & p.redn.mask < mask.E.cut[2]) ~ 'mask.midE',
                                                   masking == 1 & (p.redn.mask >=mask.E.cut[2]) ~ 'mask.highE',
                                                   T ~ 'no mask')) %>%
    {eval(parse(text = paste('dplyr::mutate(., age.grp = case_when(', str_age_recode, ', T ~ NA_character_))', sep='')))} 
  
  base.hosp = 1e5;  # hospitalization per 100K
  base.inf = 100; # infection, in %
  # aggregate across ens, since the res.yrly.stat already aggregate by year, by ens
  res.stat = rbind(rec1 %>% # combine all
                     # group_by(ens) %>%
                     summarise(totAR = (sum(num.inf) / n() * base.inf)  %>% round(., 3), 
                               indvAR1plus = (sum(num.inf>=1) / n() * base.inf) %>% round(., 3),
                               indvAR2plus = (sum(num.inf>=2) / n() * base.inf) %>% round(., 3),
                               indvAR3plus = (sum(num.inf>=3) / n() * base.inf) %>% round(., 3),
                               # for severe outcomes
                               totHospRate = (sum(num.hosp) / n() * base.hosp) %>% round(., 1), 
                               indvHospRate1plus = (sum(num.hosp>=1) / n() * base.hosp) %>% round(., 1),
                               indvHospRate2plus = (sum(num.hosp>=2) / n() * base.hosp) %>% round(., 1),
                               indvHospRate3plus = (sum(num.hosp>=3) / n() * base.hosp) %>% round(., 1)
                     ) %>%
                     dplyr::mutate(mask.grp = 'all',
                                   age.grp = 'all') %>%
                     ungroup(),
                   rec1 %>% group_by(masking) %>% # all age groups combined, by masking
                     summarise(totAR = (sum(num.inf) / n() * base.inf)  %>% round(., 3), 
                               indvAR1plus = (sum(num.inf>=1) / n() * base.inf) %>% round(., 3),
                               indvAR2plus = (sum(num.inf>=2) / n() * base.inf) %>% round(., 3),
                               indvAR3plus = (sum(num.inf>=3) / n() * base.inf) %>% round(., 3),
                               # for severe outcomes
                               totHospRate = (sum(num.hosp) / n() * base.hosp) %>% round(., 1), 
                               indvHospRate1plus = (sum(num.hosp>=1) / n() * base.hosp) %>% round(., 1),
                               indvHospRate2plus = (sum(num.hosp>=2) / n() * base.hosp) %>% round(., 1),
                               indvHospRate3plus = (sum(num.hosp>=3) / n() * base.hosp) %>% round(., 1)
                     ) %>%
                     dplyr::mutate(mask.grp = factor(masking, levels = c(0, 1), labels = c('no mask', 'use mask')),
                                   age.grp = 'all', masking = NULL) %>%
                     ungroup(),
                   rec1 %>% group_by(age.grp) %>% # all mask groups combined, by age
                     summarise(totAR = (sum(num.inf) / n() * base.inf)  %>% round(., 3), 
                               indvAR1plus = (sum(num.inf>=1) / n() * base.inf) %>% round(., 3),
                               indvAR2plus = (sum(num.inf>=2) / n() * base.inf) %>% round(., 3),
                               indvAR3plus = (sum(num.inf>=3) / n() * base.inf) %>% round(., 3),
                               # for severe outcomes
                               totHospRate = (sum(num.hosp) / n() * base.hosp) %>% round(., 1), 
                               indvHospRate1plus = (sum(num.hosp>=1) / n() * base.hosp) %>% round(., 1),
                               indvHospRate2plus = (sum(num.hosp>=2) / n() * base.hosp) %>% round(., 1),
                               indvHospRate3plus = (sum(num.hosp>=3) / n() * base.hosp) %>% round(., 1)
                     ) %>%
                     dplyr::mutate(mask.grp = 'all', masking = NULL) %>%
                     ungroup(),
                   rec1 %>% group_by(age.grp, mask.grp) %>% # by age group, mask group, 
                     summarise(totAR = (sum(num.inf) / n() * base.inf)  %>% round(., 3),  
                               indvAR1plus = (sum(num.inf>=1) / n() * base.inf) %>% round(., 3),
                               indvAR2plus = (sum(num.inf>=2) / n() * base.inf) %>% round(., 3),
                               indvAR3plus = (sum(num.inf>=3) / n() * base.inf) %>% round(., 3),
                               # for severe outcomes
                               totHospRate = (sum(num.hosp) / n() * base.hosp) %>% round(., 1), 
                               indvHospRate1plus = (sum(num.hosp>=1) / n() * base.hosp) %>% round(., 1),
                               indvHospRate2plus = (sum(num.hosp>=2) / n() * base.hosp) %>% round(., 1),
                               indvHospRate3plus = (sum(num.hosp>=3) / n() * base.hosp) %>% round(., 1)
                     ) %>% ungroup() 
                   ) %>% data.table() %>%
    # somehow did not work on the cluster...
    # dplyr::mutate(mask.grp = factor(mask.grp, levels = c('all', 'no mask', 'use mask', paste0('mask.', c('lowE', 'midE', 'highE')))),
    #               age.grp = factor(age.grp, levels = c('all', age.grps$label))) %>%
    data.table::setcolorder(., c('mask.grp', 'age.grp', 'totAR', 'indvAR1plus', 'indvAR2plus', 'indvAR3plus', 
                                 'totHospRate', 'indvHospRate1plus', 'indvHospRate2plus', 'indvHospRate3plus')) # %>%
  # dplyr::arrange(., age.grp, mask.grp)
  
  res.stat
}

fn_summary = function(x, vname){
  round.git = 3
  list(variable = vname, 
       avg = mean(x) %>% round(., round.git),
       med = median(x) %>% round(., round.git),
       mn = min(x) %>% round(., round.git),
       mx = max(x) %>% round(., round.git),
       sd = sd(x) %>% round(., round.git),
       q.025 = quantile(x, prob = .025) %>% round(., round.git),
       q.05 = quantile(x, prob = .05) %>% round(., round.git),
       q.1 = quantile(x, prob = .1) %>% round(., round.git),
       q.2 = quantile(x, prob = .2) %>% round(., round.git),
       q.3 = quantile(x, prob = .3) %>% round(., round.git),
       q.4 = quantile(x, prob = .4) %>% round(., round.git),
       q.5 = quantile(x, prob = .5) %>% round(., round.git),
       q.6 = quantile(x, prob = .6) %>% round(., round.git),
       q.7 = quantile(x, prob = .7) %>% round(., round.git),
       q.8 = quantile(x, prob = .8) %>% round(., round.git),
       q.9 = quantile(x, prob = .9) %>% round(., round.git),
       q.95 = quantile(x, prob = .95) %>% round(., round.git),
       q.975 = quantile(x, prob = .975) %>% round(., round.git)
       )
}
fn_getPopStat = function(rec){
  rec1 = rec %>% {eval(parse(text = paste('dplyr::mutate(., age.grp = case_when(', str_age_recode, ', T ~ NA_character_))', sep='')))} %>% data.table()
  # rec1 %>% group_by(ens, age.grp, masking) %>% # by masking
  #   summarise(num.nonHHctc.mn = num.nonHHctc %>% mean,num.nonHHctc.median = num.nonHHctc %>% median, 
  #             num.nonHHctc.lwr = num.nonHHctc %>% quantile(., prob = .05), num.nonHHctc.upr = num.nonHHctc %>% quantile(., prob = .95),
  #             num.HHctc.mn = num.HHctc %>% mean,num.HHctc.median = num.HHctc %>% median) %>% view
  
  
  rbind(rec1[, fn_summary(num.nonHHctc, 'num.nonHHctc'), by = list(ens, age.grp, masking)], 
        rec1[, fn_summary(num.nonHHctc, 'num.HHctc'), by = list(ens, age.grp, masking)])
}
