library(scholar)
jn <- c("bioinformatics", "methods in ecology and evolution","molecular biosystems", 
       "molecular biology and evolution")
get_impactfactor(jn)

jn <- agrep("bioinformatics", scholar:::impactfactor$Journal, ignore.case=T, value=T, max.distance=0.05)
get_impactfactor(jn)



#get personal profile
library(tidyverse)
id <- get_publications("Xr9n9cUAAAAJ&hl")
ky <- get_profile(id)
ky$name 

y <- id %>% 
  select(year, author, journal, title) %>% 
  mutate(impactFactor = get_impactfactor(journal)$ImpactFactor) %>% 
  filter(grepl("^G Yu", author)) %>% group_by(year)

y %>% summarize(total_IF = sum(impactFactor, na.rm=T), mean_IF = mean(impactFactor, na.rm=T)) 



