library(maps)
library(ggplot2)
iowa_map <- map_data('county','iowa')
iowa_map$region <- iowa_map$subregion


mapfun <- function(group, countyname)
{
  groupdat <- data.frame(county = countyname, group = as.factor(group))
  gm <- ggplot()+geom_map(data = groupdat,aes(fill=group,map_id=countyname),map=iowa_map,color='grey')+expand_limits(x=iowa_map$long,y=iowa_map$lat)
  return(gm)
}

 