# Model code used in Yang & Shaman "Reconciling the efficacy and effectiveness of masking on epidemic outcomes"
# author: Wan Yang
# date: March 2023
# note: this code is for research use only; may not be entirely optimized nor neatly organized. 

# read daily tallies


out.daily <- lapply(files.t, function(x) {
  
  tmp = try(load(x), silent = T)
  
  if(any(class(tmp) == 'try-error')){
    tmp = try(res$res.daily, silent = T)
    if(any(class(tmp) == 'try-error')){
      print(x)
      d = NULL
    }
    
  } else {
    id <- gsub(dir_res,'', x) %>% 
      strsplit('/') %>% unlist %>% tail(.,1) %>% 
      strsplit('_') %>% unlist
    
    d = res$res.daily %>% # daily tallies
      reshape2::melt(., id.vars = c('measure', 'mask.grp', 'date')) 
    
    if(unique(d$variable) == 'ens1')
      d = d %>% mutate(variable = NULL)
    
    ens.t = id[grepl('ens',id)]
    if(length(ens.t) > 0){
      # it's run for each ens member and saved so
      ens.t = ens.t %>% strsplit('\\.') %>% unlist %>% head(.,1) %>% gsub('ens','',.) %>% as.numeric()
      # update the ens number
      d$ens = ens.t
    }
    
    for(parm.t in key.parms){
      eval(parse(text = paste0('d$', parm.t,'=res$parm.list$',parm.t)))
    }
  }
  
  d
}) %>%
  rbindlist() %>%
  (function(d) d[, list(mean = round(mean(value, na.rm = T), 3), 
                        median  = round(median(value, na.rm = T), 3), 
                        ci50lwr = round(quantile(value, prob = .25, na.rm = T), 3), 
                        ci50upr = round(quantile(value, prob = .75, na.rm = T), 3),
                        ci95lwr = round(quantile(value, prob = .025, na.rm = T), 3), 
                        ci95upr = round(quantile(value, prob = .975, na.rm = T), 3)), 
                 by = c(key.parms,
                        'measure','mask.grp','date')]) 
save(out.daily, file = paste0(dir_out, 'out.daily_',fname.t,'.RData'))
print('done reading out.daily')
