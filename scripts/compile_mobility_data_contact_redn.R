# Model code used in Yang & Shaman "Reconciling the efficacy and effectiveness of masking on epidemic outcomes"
# author: Wan Yang
# date: March 2023
# note: this code is for research use only; may not be entirely optimized nor neatly organized. 

# compile mobility data (using NYC as an example) 
# to represent reduction in transmission rate due to reduction in contact rate

# TO USE - FIRST SET THE WORKING DIRECTORY TO SOURCE FILE LOCATION
dir_data = '../data/'
dir_code = '../scripts/'

source(paste0(dir_code, 'loadPackages.R'))

# mobility 
mob.type.bus = c('retail_and_recreation_percent_change_from_baseline',
                 'grocery_and_pharmacy_percent_change_from_baseline',
                 'transit_stations_percent_change_from_baseline',
                 'workplaces_percent_change_from_baseline')  
mob.type.bus0 = c('retail_and_recreation_percent_change_from_baseline',
                  'transit_stations_percent_change_from_baseline',
                  'workplaces_percent_change_from_baseline')
# ,'parks_percent_change_from_baseline',
# 'residential_percent_change_from_baseline'
mob.type.full = c('retail_and_recreation_percent_change_from_baseline',
                  'grocery_and_pharmacy_percent_change_from_baseline',
                  'parks_percent_change_from_baseline',
                  'transit_stations_percent_change_from_baseline',
                  'workplaces_percent_change_from_baseline',
                  'residential_percent_change_from_baseline')

url.mob = 'https://www.gstatic.com/covid19/mobility/Global_Mobility_Report.csv'
d.mob = read_csv(url(url.mob)) %>% data.table() 

mob = d.mob[country_region %in% c("United States") & sub_region_1 == 'New York' &
              sub_region_2 %in% c("Bronx County", "New York County", "Queens County", "Kings County", "Richmond County")]
rm(d.mob)
mob = mob[, lapply(.SD, median, na.rm=T), by = c('date', "sub_region_1"), 
          .SD = mob.type.full]
mob$mob.bus = 1 + rowMeans(mob[,mob.type.bus,with=F])/ 100  # take the average for the 
mob$mob.full = 1 + rowMeans(mob[,mob.type.full,with=F])/ 100 
mob$mob.bus0 = 1 + rowMeans(mob[,mob.type.bus0,with=F])/ 100  # take the average for the 

mob.daily = mob %>% dplyr::select(date, mob.bus, mob.full, mob.bus0)
write.csv(mob.daily, paste0(dir_data, 'da_mobility_nyc_daily.csv'), row.names = F)

matplot(mob[,c('mob.bus0','mob.bus','mob.full')], type='l', lty =1, col = c('orange','blue','black'))
abline(h=0); legend('bottomleft',c('business','all'), lty=1, col = c('blue','black'), bty='n')

