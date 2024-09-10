## 3. Mark-Recapture estimates and ASL summary

# sourcing data
# relevant data objects:
# - bc_all
# - bc_cap1 & bc_cap2
# - bc_cap1_recaps & bc_cap2_recaps
# - recaps
# - bc_river

source("R/2_BCMR_assumptions.R")
# source("R/1_BCMR_data.R")


## will make this its own file
# figure out where the hike/float break is
# make the fish that was captured in both one or the other
tapply(bc_cap1_recaps$upstream, bc_cap1_recaps$sample, summary)
tapply(bc_cap2_recaps$upstream, bc_cap2_recaps$sample, summary)
cumsum(rev(bc_river$lengths))/1000  # break at upstream dist of 42.81

subset(recaps, sample1 != sample2)
# Tag 515
# a little further into float, make it a float (or maybe try both)
recaps$sample1[recaps$Tag == 515] <- recaps$sample2[recaps$Tag == 515] <- "float"

# doing this with the object that will be used to tally!
bc_cap2_recaps$sample[bc_cap2_recaps$Tag == 515] <- "float"


# making double sure we're using the right counts
subset(bc_cap1, is.na(as.numeric(Tag)))
table(bc_cap2$Tag, is.na(as.numeric(bc_cap2$Tag)))
table(bc_cap2$Tag, is.na(as.numeric(bc_cap2$Tag))) %>% dim
table(bc_cap2_recaps$Tag, is.na(as.numeric(bc_cap2_recaps$Tag)))
table(bc_cap2_recaps$Tag, is.na(as.numeric(bc_cap2_recaps$Tag))) %>% dim
subset(bc_cap2, Tag == "Recap")


n1 <- table(bc_cap1$sample)    # this should work
n2 <- table(bc_cap2$sample)    # this should work
# m2 <- table(recaps$sample1)    # wrong, there are Recap $Tag's in bc_cap2
m2 <- table(bc_cap2_recaps$sample)  # this will work now




######### ----------------- abundance estimates ----------------- #########

library(recapr)
NBailey(n1=n1, n2=n2, m2=m2)   # estimated abundance
seBailey(n1=n1, n2=n2, m2=m2)  # SE of estimated abundance

seBailey(n1=n1, n2=n2, m2=m2)/NBailey(n1=n1, n2=n2, m2=m2)   # CV of estimated abundance

## both strata combined
Nstrat(n1=n1, n2=n2, m2=m2, estimator = "Bailey")   # estimated abundance
sestrat(n1=n1, n2=n2, m2=m2, estimator = "Bailey")  # SE of estimated abundance

sestrat(n1=n1, n2=n2, m2=m2, estimator = "Bailey") /
  Nstrat(n1=n1, n2=n2, m2=m2, estimator = "Bailey")   # CV of estimated abundance







######### ----------------- ASL estimates (length bins) ----------------- #########

lengthbreaks <- c(250, 270, 300, 350, 400, 450) # c(250, 270, 330, 400)

# just using first event sample
ASL_table(stratum = as.numeric(as.factor(bc_cap1$sample)),
          # length = bc_cap1$Length,
          age = cut(bc_cap1$Length, lengthbreaks, right=FALSE),
          Nhat = as.numeric(NBailey(n1=n1, n2=n2, m2=m2)),
          se_Nhat = as.numeric(seBailey(n1=n1, n2=n2, m2=m2)))

# just using second event sample
ASL_table(stratum = as.numeric(as.factor(bc_cap2$sample)),
          # length = bc_cap2$Length,
          age = cut(bc_cap2$Length, lengthbreaks, right=FALSE),
          Nhat = as.numeric(NBailey(n1=n1, n2=n2, m2=m2)),
          se_Nhat = as.numeric(seBailey(n1=n1, n2=n2, m2=m2)))

# pooling both events
# ASL_table(stratum = as.numeric(as.factor(c(bc_cap1$sample, bc_cap2$sample))),
#           # length = c(bc_cap1$Length, bc_cap2$Length),
#           age = cut(c(bc_cap1$Length, bc_cap2$Length), lengthbreaks, right=FALSE),
#           Nhat = as.numeric(NBailey(n1=n1, n2=n2, m2=m2)),
#           se_Nhat = as.numeric(seBailey(n1=n1, n2=n2, m2=m2)))

ASL_table(stratum = as.numeric(as.factor(bc_all$sample)),
          # length = bc_all$Length,
          age = cut(bc_all$Length, lengthbreaks, right=FALSE),
          Nhat = as.numeric(NBailey(n1=n1, n2=n2, m2=m2)),
          se_Nhat = as.numeric(seBailey(n1=n1, n2=n2, m2=m2)))


# pooling both events, but separating by capture stratum
with(subset(bc_all, sample == "float"),
     ASL_table(
       # length = Length,
       age = cut(Length, lengthbreaks, right=FALSE),
       Nhat = as.numeric(NBailey(n1=n1, n2=n2, m2=m2))[1],
       se_Nhat = as.numeric(seBailey(n1=n1, n2=n2, m2=m2))[1]))

with(subset(bc_all, sample == "hike"),
     ASL_table(
       # length = Length,
       age = cut(Length, lengthbreaks, right=FALSE),
       Nhat = as.numeric(NBailey(n1=n1, n2=n2, m2=m2))[2],
       se_Nhat = as.numeric(seBailey(n1=n1, n2=n2, m2=m2))[2]))


boilerplate <- capture.output(
  ASL_boilerplate(stratum = as.numeric(as.factor(bc_all$sample)),
                # length = bc_all$Length,
                age = cut(bc_all$Length, lengthbreaks, right=FALSE),
                Nhat = as.numeric(NBailey(n1=n1, n2=n2, m2=m2)),
                se_Nhat = as.numeric(seBailey(n1=n1, n2=n2, m2=m2))))
boilerplate_fix <- gsub(pattern="age", replacement="length", x=boilerplate)
cat(boilerplate_fix)
