# Model code used in Yang & Shaman "Reconciling the efficacy and effectiveness of masking on epidemic outcomes"
# author: Wan Yang
# date: March 2023
# note: this code is for research use only; may not be entirely optimized nor neatly organized. 

# get main runs settings

# test/comparison:
diff.ctc.redn.bymask_vec = c(F, T)
diff.ctc.bymask_vec = c(F, T)

# sensitivity analysis
parm.type_vec = c('default', 'lwr', 'upr') # diff parm sets
cont_freq_vec = c('all', "monthly_plus", "weekly_plus", "daily"); # for the contact data, diff contact frequency
wt_nonphys_vec = c(.1, .2, .5, 1); # 0.2 = assume non-phyical contact has 20% of the chance resulting infection

frac.masking_vec = seq(0, 1, by = .1)  # fm
cut.mask.start_vec = c(0, .01/100, .02/100, .05/100, .1/100, .2/100, .5/100, 1/100, 2/100) # cm
cut.mask.start_vec.lab = c('always','0.01%','0.02%','0.05%','0.1%','0.2%','0.5%', '1%', '2%')

vsBaseCTCoMOBo = T # use the same baseline regardless of diff.ctc.bymask.t / diff.ctc.redn.bymask.t
inclVax.t = T # whether to include vaccination
UseContactDuration.t = T #  whether to use contact duration directly or use number of contacts
UseVxData.t = T # whether to use vax data for the vaccination module
UseAgeSpecSuscept.t = F # whether to assume diff susceptibility by age. set it to FALSE to reduce uncertainty
MASKATHOME.t = F # set this to TRUE only when testing a very specific scenario that ppl also mask at home
parm.type.t = 'default' # which parm set to use
cont_freq.t = 'all' # for the contact data, diff contact frequency
wt_nonphys.t = .2  # assume non-phyical contact has 20% of the chance resulting infection
prob.tx.per.ctc0.t = 0.5
