## 3. Mark-Recapture estimates and ASL summary

# sourcing data
# relevant data objects:
# - bc_all
# - bc_cap1 & bc_cap2
# - bc_cap1_recaps & bc_cap2_recaps
# - recaps
# - bc_river

# source("R/2_BCMR_assumptions.R")
source("R/1_BCMR_data.R")


# ## will make this its own file
# # figure out where the hike/float break is
# # make the fish that was captured in both one or the other
# tapply(bc_cap1_recaps$upstream, bc_cap1_recaps$sample, summary)
# tapply(bc_cap2_recaps$upstream, bc_cap2_recaps$sample, summary)
# cumsum(rev(bc_river$lengths))/1000  # break at upstream dist of 42.81
#
# subset(recaps, sample1 != sample2)
# # Tag 515
# # a little further into float, make it a float (or maybe try both)
# recaps$sample1[recaps$Tag == 515] <- recaps$sample2[recaps$Tag == 515] <- "float"
#
# # doing this with the object that will be used to tally!
# bc_cap2_recaps$sample[bc_cap2_recaps$Tag == 515] <- "float"


# making double sure we're using the right counts
subset(bc_cap1, is.na(as.numeric(Tag)))
table(bc_cap2$Tag, is.na(as.numeric(bc_cap2$Tag)))
table(bc_cap2$Tag, is.na(as.numeric(bc_cap2$Tag))) %>% dim
table(bc_cap2_recaps$Tag, is.na(as.numeric(bc_cap2_recaps$Tag)))
table(bc_cap2_recaps$Tag, is.na(as.numeric(bc_cap2_recaps$Tag))) %>% dim
subset(bc_cap2, Tag == "Recap")


n1_Stratum1 <- table(bc_cap1$Stratum1)    # this should work
n2_Stratum1 <- table(bc_cap2$Stratum1)    # this should work
m2_Stratum1 <- table(bc_cap2_recaps$Stratum1)  # this will work now

n1_Stratum2 <- table(bc_cap1$Stratum2)    # this should work
n2_Stratum2 <- table(bc_cap2$Stratum2)    # this should work
m2_Stratum2 <- table(bc_cap2_recaps$Stratum2)  # this will work now




######### ----------------- abundance estimates ----------------- #########

library(recapr)

## estimates for Stratum1
nhat_Stratum1_num <- NBailey(n1=n1_Stratum1, n2=n2_Stratum1, m2=m2_Stratum1)   # estimated abundance
nhat_Stratum1 <- nhat_Stratum1_num %>%
  # unname %>%
  round %>%
  formatC(format="d", big.mark=",")

se_nhat_Stratum1_num <- seBailey(n1=n1_Stratum1, n2=n2_Stratum1, m2=m2_Stratum1)  # SE of estimated abundance
se_nhat_Stratum1 <- se_nhat_Stratum1_num  %>%
  # unname %>%
  round %>%
  formatC(format="d", big.mark=",")

seBailey(n1=n1_Stratum1, n2=n2_Stratum1, m2=m2_Stratum1) /
  NBailey(n1=n1_Stratum1, n2=n2_Stratum1, m2=m2_Stratum1)   # CV of estimated abundance


## estimates for Stratum2
nhat_Stratum2_num <- NBailey(n1=n1_Stratum2, n2=n2_Stratum2, m2=m2_Stratum2)   # estimated abundance
nhat_Stratum2 <- nhat_Stratum2_num %>%
  # unname %>%
  round %>%
  formatC(format="d", big.mark=",")

se_nhat_Stratum2_num <- seBailey(n1=n1_Stratum2, n2=n2_Stratum2, m2=m2_Stratum2)  # SE of estimated abundance
se_nhat_Stratum2 <- se_nhat_Stratum2_num  %>%
  # unname %>%
  round %>%
  formatC(format="d", big.mark=",")

seBailey(n1=n1_Stratum2, n2=n2_Stratum2, m2=m2_Stratum2) /
  NBailey(n1=n1_Stratum2, n2=n2_Stratum2, m2=m2_Stratum2)   # CV of estimated abundance


## estimates for both strata combined
nhat_combined_Stratum1_num <-
  Nstrat(n1=n1_Stratum1, n2=n2_Stratum1, m2=m2_Stratum1, estimator = "Bailey")   # estimated abundance
nhat_combined_Stratum1 <- nhat_combined_Stratum1_num  %>%
  # unname %>%
  round %>%
  formatC(format="d", big.mark=",")

se_nhat_combined_Stratum1_num <-
  sestrat(n1=n1_Stratum1, n2=n2_Stratum1, m2=m2_Stratum1, estimator = "Bailey")  # SE of estimated abundance
se_nhat_combined_Stratum1 <- se_nhat_combined_Stratum1_num  %>%
  # unname %>%
  round %>%
  formatC(format="d", big.mark=",")

sestrat(n1=n1_Stratum1, n2=n2_Stratum1, m2=m2_Stratum1, estimator = "Bailey") /
  Nstrat(n1=n1_Stratum1, n2=n2_Stratum1, m2=m2_Stratum1, estimator = "Bailey")   # CV of estimated abundance



Nstrat(n1=n1_Stratum2, n2=n2_Stratum2, m2=m2_Stratum2, estimator = "Bailey")   # estimated abundance
sestrat(n1=n1_Stratum2, n2=n2_Stratum2, m2=m2_Stratum2, estimator = "Bailey")  # SE of estimated abundance

sestrat(n1=n1_Stratum2, n2=n2_Stratum2, m2=m2_Stratum2, estimator = "Bailey") /
  Nstrat(n1=n1_Stratum2, n2=n2_Stratum2, m2=m2_Stratum2, estimator = "Bailey")   # CV of estimated abundance







######### ----------------- ASL estimates (length bins) ----------------- #########

### should recaps only count once???
bc_all_justonce <- rbind(bc_cap1,
                         subset(bc_cap2, Tag=="-"))

lengthbreaks <- c(250, 270, 300, 350, 450) # c(250, 270, 300, 350, 400, 450) # c(250, 270, 330, 400)

# Full study area - just using first event sample
ASL_firstevent <-
  ASL_table(stratum = as.numeric(as.factor(bc_cap1$Stratum1)),
          # length = bc_cap1$Length,
          age = cut(bc_cap1$Length, lengthbreaks, right=FALSE),
          Nhat = as.numeric(NBailey(n1=n1_Stratum1, n2=n2_Stratum1, m2=m2_Stratum1)),
          se_Nhat = as.numeric(seBailey(n1=n1_Stratum1, n2=n2_Stratum1, m2=m2_Stratum1)))

# # just using second event sample
# ASL_secondevent <-
#   ASL_table(stratum = as.numeric(as.factor(bc_cap2$Stratum1)),
#           # length = bc_cap2$Length,
#           age = cut(bc_cap2$Length, lengthbreaks, right=FALSE),
#           Nhat = as.numeric(NBailey(n1=n1_Stratum1, n2=n2_Stratum1, m2=m2_Stratum1)),
#           se_Nhat = as.numeric(seBailey(n1=n1_Stratum1, n2=n2_Stratum1, m2=m2_Stratum1)))
#
# # pooling both samples
# ASL_pooled <-
#   with(bc_all_justonce, # bc_all
#      ASL_table(stratum = as.numeric(as.factor(Stratum1)),
#                # length = bc_all$Length,
#                age = cut(Length, lengthbreaks, right=FALSE),
#                Nhat = as.numeric(NBailey(n1=n1_Stratum1, n2=n2_Stratum1, m2=m2_Stratum1)),
#                se_Nhat = as.numeric(seBailey(n1=n1_Stratum1, n2=n2_Stratum1, m2=m2_Stratum1))))



# # just using first event sample
# ASL_table(stratum = as.numeric(as.factor(bc_cap1$Stratum2)),
#           # length = bc_cap1$Length,
#           age = cut(bc_cap1$Length, lengthbreaks, right=FALSE),
#           Nhat = as.numeric(NBailey(n1=n1_Stratum2, n2=n2_Stratum2, m2=m2_Stratum2)),
#           se_Nhat = as.numeric(seBailey(n1=n1_Stratum2, n2=n2_Stratum2, m2=m2_Stratum2)))
#
# # just using second event sample
# ASL_table(stratum = as.numeric(as.factor(bc_cap2$Stratum2)),
#           # length = bc_cap2$Length,
#           age = cut(bc_cap2$Length, lengthbreaks, right=FALSE),
#           Nhat = as.numeric(NBailey(n1=n1_Stratum2, n2=n2_Stratum2, m2=m2_Stratum2)),
#           se_Nhat = as.numeric(seBailey(n1=n1_Stratum2, n2=n2_Stratum2, m2=m2_Stratum2)))
#
# # pooling both samples
# with(bc_all_justonce, # bc_all
#      ASL_table(stratum = as.numeric(as.factor(Stratum2)),
#                # length = bc_all$Length,
#                age = cut(Length, lengthbreaks, right=FALSE),
#                Nhat = as.numeric(NBailey(n1=n1_Stratum2, n2=n2_Stratum2, m2=m2_Stratum2)),
#                se_Nhat = as.numeric(seBailey(n1=n1_Stratum2, n2=n2_Stratum2, m2=m2_Stratum2))))




##  separating by capture stratum
ASL_Beaver <-
  with(subset(bc_all_justonce, Stratum1 == "Beaver"),    ###### using both events
     ASL_table(
       # length = Length,
       age = cut(Length, lengthbreaks, right=FALSE),
       Nhat = as.numeric(NBailey(n1=n1_Stratum1, n2=n2_Stratum1, m2=m2_Stratum1)["Beaver"]),
       se_Nhat = as.numeric(seBailey(n1=n1_Stratum1, n2=n2_Stratum1, m2=m2_Stratum1)["Beaver"])))

ASL_Nome <-
  with(subset(bc_cap1, Stratum1 == "Nome"),    ###### using just first event
     ASL_table(
       # length = Length,
       age = cut(Length, lengthbreaks, right=FALSE),
       Nhat = as.numeric(NBailey(n1=n1_Stratum1, n2=n2_Stratum1, m2=m2_Stratum1)["Nome"]),
       se_Nhat = as.numeric(seBailey(n1=n1_Stratum1, n2=n2_Stratum1, m2=m2_Stratum1)["Nome"])))



ASL_2000 <-
  with(subset(bc_all_justonce, Stratum2 == "2000 Study"),    ###### using both events
     ASL_table(
       # length = Length,
       age = cut(Length, lengthbreaks, right=FALSE),
       Nhat = as.numeric(NBailey(n1=n1_Stratum2, n2=n2_Stratum2, m2=m2_Stratum2)["2000 Study"]),
       se_Nhat = as.numeric(seBailey(n1=n1_Stratum2, n2=n2_Stratum2, m2=m2_Stratum2)["2000 Study"])))

# ASL_table(
#   # length = Length,
#   age = cut(bc_all$Length[bc_all$Stratum2 == "2000 Study"], lengthbreaks, right=FALSE),
#   Nhat = as.numeric(NBailey(n1=n1_Stratum2, n2=n2_Stratum2, m2=m2_Stratum2))[1],
#   se_Nhat = as.numeric(seBailey(n1=n1_Stratum2, n2=n2_Stratum2, m2=m2_Stratum2))[1])

# with(subset(bc_all_justonce, Stratum2 == "Beaver"),
#      ASL_table(
#        # length = Length,
#        age = cut(Length, lengthbreaks, right=FALSE),
#        Nhat = as.numeric(NBailey(n1=n1_Stratum2, n2=n2_Stratum2, m2=m2_Stratum2))[2],
#        se_Nhat = as.numeric(seBailey(n1=n1_Stratum2, n2=n2_Stratum2, m2=m2_Stratum2))[2]))
#
# with(subset(bc_all_justonce, Stratum2 == "Lower Nome"),
#      ASL_table(
#        # length = Length,
#        age = cut(Length, lengthbreaks, right=FALSE),
#        Nhat = as.numeric(NBailey(n1=n1_Stratum2, n2=n2_Stratum2, m2=m2_Stratum2))[3],
#        se_Nhat = as.numeric(seBailey(n1=n1_Stratum2, n2=n2_Stratum2, m2=m2_Stratum2))[3]))
#
# with(subset(bc_all_justonce, Stratum2 == "Upper Nome"),
#      ASL_table(
#        # length = Length,
#        age = cut(Length, lengthbreaks, right=FALSE),
#        Nhat = as.numeric(NBailey(n1=n1_Stratum2, n2=n2_Stratum2, m2=m2_Stratum2))[4],
#        se_Nhat = as.numeric(seBailey(n1=n1_Stratum2, n2=n2_Stratum2, m2=m2_Stratum2))[4]))




# boilerplate <- capture.output(
#   ASL_boilerplate(stratum = as.numeric(as.factor(bc_all$sample)),
#                 # length = bc_all$Length,
#                 age = cut(bc_all$Length, lengthbreaks, right=FALSE),
#                 Nhat = as.numeric(NBailey(n1=n1, n2=n2, m2=m2)),
#                 se_Nhat = as.numeric(seBailey(n1=n1, n2=n2, m2=m2))))
# boilerplate_fix <- gsub(pattern="age", replacement="length", x=boilerplate)
# cat(boilerplate_fix)
