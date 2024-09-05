## 1. Beaver Creek Mark-Recapture data

## This script does all data import and data manipulation, and will be sourced
## in subsequent scripts.


library(tidyverse)  # for data manipulation
library(riverdist)  # for spatial stuff with the river network
library(dsftools)   # for some data manipulation stuff and ASL summaries


# a quick & dirty function to tabulate columns
column_tabler <- function(x) {
  x <- as.data.frame(x)
  for(j in 1:ncol(x)) {
    print(names(x)[j])
    if(!(class(x[,j]) %in% c("numeric","integer"))) {
      print(table(x[,j], useNA="always"))
    }
    print(table(is.na(x[,j])))
  }
}


### reading mark-recapture data and doing some basic reformatting

mark_float <- read_csv("Data/mark_float.csv") %>%
  select(which(colMeans(is.na(.)) < 1)) %>%
  filter(rowMeans(is.na(.)) < 1) %>%
  mutate(Length = as.numeric(Length)) %>%
  rename(Fish = `Fish #`) %>%
  rename(Tag = `Mark Tag#`) %>%
  rename(Comments = `...6`)
# head(mark_float)
# column_tabler(mark_float)

mark_hike <- read_csv("Data/mark_hike.csv") %>%
  select(which(colMeans(is.na(.)) < 1)) %>%
  filter(rowMeans(is.na(.)) < 1) %>%
  mutate(Length = as.numeric(Length)) %>%
  rename(Fish = `Fish #`) %>%
  rename(Tag = `Mark Tag#`) %>%
  rename(Comments = `...6`)
# head(mark_hike)
# column_tabler(mark_hike)


recap_float <- read_csv("Data/recap_float.csv") %>%
  select(which(colMeans(is.na(.)) < 1)) %>%
  filter(rowMeans(is.na(.)) < 1) %>%
  mutate(Length = as.numeric(Length)) %>%
  rename(Fish = `Fish #`) %>%
  rename(Tag = `Recap Tag#`) %>%
  rename(Comments = `...6`)
# head(recap_float)
# column_tabler(recap_float)


recap_hike <- read_csv("Data/recap_hike.csv") %>%
  select(which(colMeans(is.na(.)) < 1)) %>%
  filter(rowMeans(is.na(.)) < 1) %>%
  mutate(Length = as.numeric(Length)) %>%
  rename(Fish = `Fish #`) %>%
  rename(Tag = `Recap Tag#`) %>%
  rename(Comments = `...6`)
# head(recap_hike)
# column_tabler(recap_hike)




# making sure the column names match, then cbinding
names(mark_float)
names(mark_hike)
names(recap_float)
names(recap_hike)

# compiling all four dataframes into one master dataset bc_all
mark_all <- rbind(cbind(mark_float, sample="float"),
                  cbind(mark_hike, sample="hike"))
recap_all <- rbind(cbind(recap_float, sample="float"),
                   cbind(recap_hike, sample="hike"))
bc_all <- rbind(cbind(mark_all, event="mark"),
                cbind(recap_all, event="recap"))
summary(bc_all)



# Error checking
plot(bc_all$Length)   # two lengths under 250, did we want to censor?
head(sort(bc_all$Length))

tagtab <- table(bc_all$Tag, bc_all$event)
# tagtab
tagtab[tagtab[,1]>1 | tagtab[,2]>1,]
# tag 146 recorded twice (mark float)
# tag 772 recorded three times (recap hike)

bc_all %>% filter(Tag==146)
bc_all %>% filter(Tag==772)




# Data questions for Lisa:
# - repeat tags in both events (146 & 772)
#    + yes, keep only one observation.  Maybe average length & longlat
# - censor the small ones?
#    + yes, censor
# - confirm censor the ones that say "Recap?"
#    + yes, censor

# Additional data fixing per Lisa email:
bc_all <- bc_all %>%
  filter(is.na(Length) | Length >= 250) %>%
  filter(is.na(Tag) | Tag != "Recap?")

# dealing with fish 146
bc_all$Latitude[bc_all$Length==260 & bc_all$Tag==146 & bc_all$event=="mark"] <-
  mean(bc_all$Latitude[bc_all$Tag==146 & bc_all$event=="mark"], na.rm=TRUE)
bc_all$Longitude[bc_all$Length==260 & bc_all$Tag==146 & bc_all$event=="mark"] <-
  mean(bc_all$Longitude[bc_all$Tag==146 & bc_all$event=="mark"], na.rm=TRUE)
bc_all$Length[bc_all$Length==260 & bc_all$Tag==146 & bc_all$event=="mark"] <-
  mean(bc_all$Length[bc_all$Tag==146 & bc_all$event=="mark"], na.rm=TRUE)
bc_all <- bc_all[!(bc_all$Length==265 & bc_all$Tag==146 & bc_all$event=="mark"),]

# dealing with fish 772
bc_all$Latitude[bc_all$Length==362 & bc_all$Tag==772 & bc_all$event=="recap"] <-
  mean(bc_all$Latitude[bc_all$Tag==772 & bc_all$event=="recap"], na.rm=TRUE)
bc_all$Longitude[bc_all$Length==362 & bc_all$Tag==772 & bc_all$event=="recap"] <-
  mean(bc_all$Longitude[bc_all$Tag==772 & bc_all$event=="recap"], na.rm=TRUE)
bc_all$Length[bc_all$Length==362 & bc_all$Tag==772 & bc_all$event=="recap"] <-
  mean(bc_all$Length[bc_all$Tag==772 & bc_all$event=="recap"], na.rm=TRUE)
bc_all <- bc_all[!(bc_all$Length==370 & bc_all$Tag==772 & bc_all$event=="recap"),]


## quality control on coordinates (will use coords to assign sections)
with(bc_all, plot(x=Longitude, y=Latitude))  # ok some serious outliers
bc_all <- bc_all %>%
  mutate(Longitude = ifelse(Longitude < -150, Longitude+10, Longitude)) %>%
  mutate(Latitude = ifelse(Latitude==65.6534096, 65.34096, Latitude)) %>%
  mutate(Latitude = ifelse(Latitude==65.72964, 65.32964, Latitude)) %>%
  mutate(Longitude = ifelse(Longitude==-147.9696, -147.69616, Longitude)) %>%
  mutate(Longitude = ifelse(Longitude==-147.95062, -146.95062, Longitude)) %>%
  mutate(Longitude = ifelse(Longitude==-147.9976, -146.9976, Longitude)) %>%
  mutate(Longitude = ifelse(Longitude==-146.35953, -146.75953, Longitude)) %>%
  mutate(Longitude = ifelse(Longitude==-146.38034, -146.68034, Longitude))
with(bc_all, plot(x=Longitude, y=Latitude))  # much better!




# reading the .csv file with section starts/stops
site_boundaries <- read_csv("Data/site_boundaries.csv", skip=2)[1:40,] %>%
  rename(Latitude1 = Latitude...2) %>%
  rename(Latitude2 = Latitude...4) %>%
  rename(Longitude1 =Longitude...3)%>%
  rename(Longitude2 =Longitude...5)
with(site_boundaries, plot(Longitude1, Latitude1, pch="+", col=3,
     xlim=range(site_boundaries$Longitude1, site_boundaries$Longitude2, bc_all$Longitude),
     ylim=range(site_boundaries$Latitude1, site_boundaries$Latitude2, bc_all$Latitude)))
with(site_boundaries, points(Longitude2, Latitude2, pch="x", col=2))
with(bc_all, points(Longitude, Latitude, col=adjustcolor(1, alpha.f=.1)))


# loading the rivernetwork to make sure that points fall close enough to river
load(file="Data/beaver_cr_rivernetwork_op.Rdata")
plot(beaver_cr_op)

# reprojecting data locations to the same coord system as rivernetwork
AKalbers <- "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154
    +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80"
points_albers <- sf::sf_project(pts=bc_all[,3:2], to=AKalbers)
points(points_albers, pch=16, col=2)

# simplifying rivernetwork
bc_river <- trimriver(rivers=beaver_cr_op, trimto=c(5:8,13:16))
bc_river <- setmouth(seg=8, vert=211, rivers=bc_river)
plot(bc_river)
points(points_albers)

# snapping river coordinates
points_segvert <- xy2segvert(x=points_albers[,1],
                             y=points_albers[,2],
                             rivers=bc_river)
bc_all$seg <- points_segvert$seg
bc_all$vert <- points_segvert$vert

# adding upstream position to dataset (makes distance calculations easy!)
bc_all$upstream <- with(bc_all, mouthdist(seg=seg, vert=vert, rivers=bc_river))/1000


plot(bc_river)
points(sf::sf_project(pts=site_boundaries[,3:2], to=AKalbers))
lines(sf::sf_project(pts=site_boundaries[,3:2], to=AKalbers))


# assigning sites (sections of river network) to point data
# haha this is kludgy!!
bc_all$Site <- NA
for(i in 1:nrow(site_boundaries)) {
  bc_all$Site[bc_all$Longitude <= site_boundaries$Longitude1[i] &
                bc_all$Longitude >= site_boundaries$Longitude2[i]] <- i
}
bc_all$Site[bc_all$Longitude < min(site_boundaries$Longitude2)] <- 1+nrow(site_boundaries)
bc_all %>% mutate(Site=as.factor(Site)) %>%
  ggplot(aes(x=Longitude, y=Latitude, col=Site)) +
           geom_point()


# separating datasets for mark, recap, and both
bc_cap1 <- filter(bc_all, event=="mark")
bc_cap2 <- filter(bc_all, event=="recap")

# tags observed in both
recap_tags <- bc_cap1$Tag[!is.na(bc_cap1$Tag) & bc_cap1$Tag %in% bc_cap2$Tag]
bc_cap1_recaps <- filter(bc_cap1, Tag %in% recap_tags)
bc_cap2_recaps_justtags <- filter(bc_cap2, Tag %in% recap_tags)
bc_cap2_recaps <- filter(bc_cap2, Tag %in% c(recap_tags, "Recap"))


### MAKE SURE FINAL CALCS INCLUDE NON-NUMBERED FISH
### - Tag = "Recap"
### - Tag = "none"



## things that will need to be checked:
## * immigration emigration relative magnitude
##   - calculate movement by recaptured indivs (range of upstream!)
##   - see if this varies spatially, estimate near endpoints
##     + no wait, they can also immigrate/emigrate to/from other streams
##     + might be able to get around this with summer fidelity, see what is in OP
## * size selectivity (chi2 I believe) - KS
## * spatial selectivity (also chi2) - could even KS with upstream!!
