library(dplyr)      # for data manipulation functions
library(tidyr)      # for data manipulation functions
library(data.table) # for function `fread`
library(broom)      # for function `tidy`

#overall test
shapirotest <- updated_207_quant %>% 
  gather(key = "variable_name", value = "value") %>% 
  group_by(variable_name)  %>% 
  do(tidy(shapiro.test(.$value))) %>% 
  ungroup()%>% 
  select(-method)
write.csv(shapirotest)

#test by cluster
shapirotestcluster <- 
  updated_207 %>% 
  gather(key = "variable_name", value = "value", -cluster) %>% 
  group_by(variable_name, cluster)  %>% 
  do(tidy(shapiro.test(.$value))) %>% 
  ungroup()%>% 
  select(-method)
write.csv(shapirotestcluster)
