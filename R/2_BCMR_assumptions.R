## 2. Checking Mark-Recapture Assumptions


# To summarize:
# - Size selectivity was detected when pooling both samples
#   . KS tests were satisfied by length stratification at 270mm FL
# - Spatial selectivity was also detected when pooling both samples (KS by upriver position)
# - **BUT** stratification by sample (float vs hike) alleviated selectivity
#   . both by length and by upriver position!
#   . no evidence of length selectivity in either event
#     (ks pvals = 0.74, 0.97, 0.08, and 0.13)
#     so I think we can pool both events for each sample for estimating proportions!
# - Not concerned about growth recruitment
# - We lose an estimated 1% of tagged fish out the bottom, so abundance est might be biased low by ~1%
# - One fish (tag 515) was marked in one stratum and recaptured in the other.
#   For the purpose of stratification, it was assigned to the "float" stratum, since
#   it was observed further in that direction.



# sourcing data
# relevant data objects:
# - bc_all
# - bc_cap1 & bc_cap2
# - bc_cap1_recaps & bc_cap2_recaps
source("R/1_BCMR_data.R")



## --------------- length selectivity ------------------ ##

## defining a function to illustrate differences in distribution, with ks test
ksplot <- function(x1, x2, legend=c("x1","x2"), main="", col=c(1,1), lty=c(1,1), xlab="") {
  d1 <- density(x1, na.rm=TRUE)
  d2 <- density(x2, na.rm=TRUE)
  plot(d1, main=main, col=col[1], lty=lty[1],
       xlim=range(d1$x, d2$x), ylim=range(d1$y, d2$y), xlab=xlab)
  lines(density(x2, na.rm=TRUE), col=col[2], lty=lty[2])
  legend("topright", lty=lty, col=col,
         legend=paste0(legend, " (n=",c(sum(!is.na(x1)), sum(!is.na(x2))),")"))

  ksks <- suppressWarnings(ks.test(x1, x2))
  plot(ecdf(x1), main=c(main,
                        paste0("D=", signif(ksks$statistic, digits=3), ", ",
                        "pval=", signif(ksks$p.value, digits=3), collapse=NULL)),
       col=col[1], lty=lty[1], xlab=xlab)
  plot(ecdf(x2), col=col[2], lty=lty[2], add=TRUE)
  legend("bottomright", lty=lty, col=col,
         legend=paste0(legend, " (n=",c(sum(!is.na(x1)), sum(!is.na(x2))),")"))
}

# # KS tests for size selectivity
# ks.test(bc_cap1$Length, bc_cap1_recaps$Length)  # cap vs recap for first event
# # D = 0.1758, p-value = 0.001329
# # evidence of size selectivity in second event
#
# # ksplot(x1=bc_cap1$Length, x2=bc_cap1_recaps$Length,
# #        main="Event 1", legend=c("All","Recaps"), col=c(1,2))

par(mfrow=c(2,2))
ksplot(x1=bc_cap1$Length, x2=bc_cap1_recaps$Length,
       main="Event 1", legend=c("All","Recaps"), col=c(1,2),
       xlab="Fork Length (mm)")
ksplot(x1=bc_cap2$Length, x2=bc_cap2_recaps$Length,
       main="Event 2", legend=c("All","Recaps"), col=c(1,4),
       xlab="Fork Length (mm)")



# actually very curious how cap-cap and recap-recap compare
par(mfrow=c(2,2))
ksplot(x1=bc_cap1$Length, x2=bc_cap2$Length,
       main="All Captures", legend=c("Event 1","Event 2"), lty=c(1,2),
       xlab="Fork Length (mm)")
ksplot(x1=bc_cap1_recaps$Length, x2=bc_cap2_recaps$Length,
       main="Recaptures", legend=c("Event 1","Event 2"), col=c(2,4),
       xlab="Fork Length (mm)")






# what happens at various breakpoints - all tests satisfied at 270
i <- 1
p1 <- p2 <- p3 <- p4 <- NA
possiblebreaks <- 260:300
for(breakpt in possiblebreaks) {
  p1[i] <- ks.test(bc_cap1$Length[bc_cap1$Length < breakpt],
          bc_cap1_recaps$Length[bc_cap1_recaps$Length < breakpt])$p.value
  p2[i] <- ks.test(bc_cap2$Length[bc_cap2$Length < breakpt],
          bc_cap2_recaps$Length[bc_cap2_recaps$Length < breakpt])$p.value

  p3[i] <- ks.test(bc_cap1$Length[bc_cap1$Length >= breakpt],
          bc_cap1_recaps$Length[bc_cap1_recaps$Length >= breakpt])$p.value
  p4[i] <- ks.test(bc_cap2$Length[bc_cap2$Length >= breakpt],
          bc_cap2_recaps$Length[bc_cap2_recaps$Length >= breakpt])$p.value
  i <- i+1
}
par(mfrow=c(1,1))
plot(possiblebreaks, p1, ylim=range(p1,p2,p3,p4), type='l', log="y", col=2,
     xlab="possible length breakpoints", ylab="p-value")
lines(possiblebreaks, p2, col=3)
lines(possiblebreaks, p3, col=4)
lines(possiblebreaks, p4, col=5)
abline(h=.05, lty=2)
text(x=min(possiblebreaks), y=0.05, labels="0.05", pos=3)
legend("topright",lty=1, col=2:5,
       legend=c("Event 1 All vs Recaps - small fish",
                "Event 2 All vs Recaps - small fish",
                "Event 1 All vs Recaps - large fish",
                "Event 2 All vs Recaps - large fish"))


## trying a break at 270

sum(bc_cap1_recaps$Length < 270)
sum(bc_cap2_recaps$Length < 270)



par(mfrow=c(2,2))
ksplot(x1=bc_cap1$Length %s_l% 270, x2=bc_cap1_recaps$Length %s_l% 270,
       main="Event 1, lengths < 270mm", legend=c("All","Recaps"), col=c(1,2),
       xlab="Fork Length (mm)")
ksplot(x1=bc_cap2$Length %s_l% 270, x2=bc_cap2_recaps$Length %s_l% 270,
       main="Event 2, lengths < 270mm", legend=c("All","Recaps"), col=c(1,4),
       xlab="Fork Length (mm)")

ksplot(x1=bc_cap1$Length %s_geq% 270, x2=bc_cap1_recaps$Length %s_geq% 270,
       main="Event 1, lengths >= 270mm", legend=c("All","Recaps"), col=c(1,2),
       xlab="Fork Length (mm)")
ksplot(x1=bc_cap2$Length %s_geq% 270, x2=bc_cap2_recaps$Length %s_geq% 270,
       main="Event 2, lengths >= 270mm", legend=c("All","Recaps"), col=c(1,4),
       xlab="Fork Length (mm)")

# par(mfrow=c(2,2))
# plot(density(bc_cap1$Length[bc_cap1$Length < 270], na.rm=T))
# lines(density(bc_cap1_recaps$Length[bc_cap1_recaps$Length < 270], na.rm=T), col=2)
# plot(ecdf(bc_cap1$Length[bc_cap1$Length < 270]))
# plot(ecdf(bc_cap1_recaps$Length[bc_cap1_recaps$Length < 270]), add=T, col=2)
#
# plot(density(bc_cap2$Length[bc_cap2$Length < 270], na.rm=T))
# lines(density(bc_cap2_recaps$Length[bc_cap2_recaps$Length < 270], na.rm=T), col=4)
# plot(ecdf(bc_cap2$Length[bc_cap2$Length < 270]))
# plot(ecdf(bc_cap2_recaps$Length[bc_cap2_recaps$Length < 270]), add=T, col=4)
#
# plot(density(bc_cap1$Length[bc_cap1$Length >= 270], na.rm=T))
# lines(density(bc_cap1_recaps$Length[bc_cap1_recaps$Length >= 270], na.rm=T), col=2)
# plot(ecdf(bc_cap1$Length[bc_cap1$Length >= 270]))
# plot(ecdf(bc_cap1_recaps$Length[bc_cap1_recaps$Length >= 270]), add=T, col=2)
#
# plot(density(bc_cap2$Length[bc_cap2$Length >= 270], na.rm=T))
# lines(density(bc_cap2_recaps$Length[bc_cap2_recaps$Length >= 270], na.rm=T), col=4)
# plot(ecdf(bc_cap2$Length[bc_cap2$Length >= 270]))
# plot(ecdf(bc_cap2_recaps$Length[bc_cap2_recaps$Length >= 270]), add=T, col=4)


## trying breaks at 270 and 300
## tests fail for 300-500 but not 270-500!  i hate ks tests

breaks <- c(250, 270, 500, 300)#
bc_cap1$Lengthbin <- cut(bc_cap1$Length, breaks, right=FALSE)
bc_cap2$Lengthbin <- cut(bc_cap2$Length, breaks, right=FALSE)
bc_cap1_recaps$Lengthbin <- cut(bc_cap1_recaps$Length, breaks, right=FALSE)
bc_cap2_recaps$Lengthbin <- cut(bc_cap2_recaps$Length, breaks, right=FALSE)
for(i in 1:(length(breaks)-1)) {
  x1 <- bc_cap1$Length[as.numeric(bc_cap1$Lengthbin) == i]
  x2 <- bc_cap1_recaps$Length[as.numeric(bc_cap1_recaps$Lengthbin) == i]
  x3 <- bc_cap2$Length[as.numeric(bc_cap2$Lengthbin) == i]
  x4 <- bc_cap2_recaps$Length[as.numeric(bc_cap2_recaps$Lengthbin) == i]

  # print(ks.test(x1, x2))  # cap vs recap for first event
  #
  # par(mfrow=c(2,2))
  # plot(density(x1, na.rm=T))
  # lines(density(x2, na.rm=T), col=2)
  # plot(ecdf(x1))
  # plot(ecdf(x2), add=T, col=2)
  #
  # print(ks.test(x3, x4))  # cap vs recap for first event
  #
  # plot(density(x3, na.rm=T))
  # lines(density(x4, na.rm=T), col=4)
  # plot(ecdf(x3))
  # plot(ecdf(x4), add=T, col=4)

  par(mfrow=c(2,2))
  ksplot(x1=x1, x2=x2,
         main=paste("Event 1,", levels(bc_cap1$Lengthbin)[i]),
         legend=c("All","Recaps"), col=c(1,2),
         xlab="Fork Length (mm)")
  ksplot(x1=x3, x2=x4,
         main=paste("Event 2,", levels(bc_cap1$Lengthbin)[i]),
         legend=c("All","Recaps"), col=c(1,4),
         xlab="Fork Length (mm)")
}


## -------------- spatial selectivity? ---------------- ##
ks.test(bc_cap1$upstream, bc_cap1_recaps$upstream)
# D = 0.22825, p-value = 8.778e-06
ksplot(bc_cap1$upstream, bc_cap1_recaps$upstream,
       main="Event 1", legend=c("All","Recaps"), col=c(1,2),
       xlab="Upstream Position (rkm)")

ks.test(bc_cap2$upstream, bc_cap2_recaps$upstream)
# D = 0.49899, p-value < 2.2e-16
ksplot(bc_cap2$upstream, bc_cap2_recaps$upstream,
       main="Event 2", legend=c("All","Recaps"), col=c(1,4),
       xlab="Upstream Position (rkm)")


## actually I don't think we need to worry about this one
### actually yes we do!!  This suggests much more fish in the lower river
### might even be float vs hike!   - IT IS, KEEP THIS STRATIFICATION!!!
# par(mfrow=c(2,2))
# ks.test(bc_cap1$upstream[bc_cap1$sample=="float"], bc_cap1_recaps$upstream[bc_cap1_recaps$sample=="float"])
# # D = 0.13742, p-value = 0.428
# ksplot(bc_cap1$upstream[bc_cap1$sample=="float"], bc_cap1_recaps$upstream[bc_cap1_recaps$sample=="float"],
#        main="Event 1 - float", legend=c("All","Recaps"), col=c(1,2),
#        xlab="Upstream Position (rkm)")
#
# ks.test(bc_cap2$upstream[bc_cap2$sample=="float"], bc_cap2_recaps$upstream[bc_cap2_recaps$sample=="float"])
# # D = 0.1145, p-value = 0.6445
# ksplot(bc_cap2$upstream[bc_cap2$sample=="float"], bc_cap2_recaps$upstream[bc_cap2_recaps$sample=="float"],
#        main="Event 2 - float", legend=c("All","Recaps"), col=c(1,4),
#        xlab="Upstream Position (rkm)")
#
# par(mfrow=c(2,2))
# ks.test(bc_cap1$upstream[bc_cap1$sample=="hike"], bc_cap1_recaps$upstream[bc_cap1_recaps$sample=="hike"])
# # D = 0.11804, p-value = 0.2432
# ksplot(bc_cap1$upstream[bc_cap1$sample=="hike"], bc_cap1_recaps$upstream[bc_cap1_recaps$sample=="hike"],
#        main="Event 1 - hike", legend=c("All","Recaps"), col=c(1,2),
#        xlab="Upstream Position (rkm)")
#
# ks.test(bc_cap2$upstream[bc_cap2$sample=="hike"], bc_cap2_recaps$upstream[bc_cap2_recaps$sample=="hike"])
# # D = 0.11701, p-value = 0.2823
# ksplot(bc_cap2$upstream[bc_cap2$sample=="hike"], bc_cap2_recaps$upstream[bc_cap2_recaps$sample=="hike"],
#        main="Event 2 - hike", legend=c("All","Recaps"), col=c(1,4),
#        xlab="Upstream Position (rkm)")


par(mfrow=c(2,2))
ks.test(bc_cap1$upstream[bc_cap1$Stratum1=="Beaver"], bc_cap1_recaps$upstream[bc_cap1_recaps$Stratum1=="Beaver"])
# D = 0.11782, p-value = 0.6284
ksplot(bc_cap1$upstream[bc_cap1$Stratum1=="Beaver"], bc_cap1_recaps$upstream[bc_cap1_recaps$Stratum1=="Beaver"],
       main="Event 1 - Beaver", legend=c("All","Recaps"), col=c(1,2),
       xlab="Upstream Position (rkm)")

ks.test(bc_cap2$upstream[bc_cap2$Stratum1=="Beaver"], bc_cap2_recaps$upstream[bc_cap2_recaps$Stratum1=="Beaver"])
# D = 0.12232, p-value = 0.5175
ksplot(bc_cap2$upstream[bc_cap2$Stratum1=="Beaver"], bc_cap2_recaps$upstream[bc_cap2_recaps$Stratum1=="Beaver"],
       main="Event 2 - Beaver", legend=c("All","Recaps"), col=c(1,4),
       xlab="Upstream Position (rkm)")

par(mfrow=c(2,2))
ks.test(bc_cap1$upstream[bc_cap1$Stratum1=="Nome"], bc_cap1_recaps$upstream[bc_cap1_recaps$Stratum1=="Nome"])
# D = 0.1425, p-value = 0.09111
ksplot(bc_cap1$upstream[bc_cap1$Stratum1=="Nome"], bc_cap1_recaps$upstream[bc_cap1_recaps$Stratum1=="Nome"],
       main="Event 1 - Nome", legend=c("All","Recaps"), col=c(1,2),
       xlab="Upstream Position (rkm)")

ks.test(bc_cap2$upstream[bc_cap2$Stratum1=="Nome"], bc_cap2_recaps$upstream[bc_cap2_recaps$Stratum1=="Nome"])
# D = 0.11333, p-value = 0.3215
ksplot(bc_cap2$upstream[bc_cap2$Stratum1=="Nome"], bc_cap2_recaps$upstream[bc_cap2_recaps$Stratum1=="Nome"],
       main="Event 2 - Nome", legend=c("All","Recaps"), col=c(1,4),
       xlab="Upstream Position (rkm)")

par(mfrow=c(2,2))
ks.test(bc_cap1$upstream[bc_cap1$Stratum2=="Lower Nome"], bc_cap1_recaps$upstream[bc_cap1_recaps$Stratum2=="Lower Nome"])
# D = 0.52273, p-value = 0.1898
ksplot(bc_cap1$upstream[bc_cap1$Stratum2=="Lower Nome"], bc_cap1_recaps$upstream[bc_cap1_recaps$Stratum2=="Lower Nome"],
       main="Event 1 - Lower Nome", legend=c("All","Recaps"), col=c(1,2),
       xlab="Upstream Position (rkm)")

ks.test(bc_cap2$upstream[bc_cap2$Stratum2=="Lower Nome"], bc_cap2_recaps$upstream[bc_cap2_recaps$Stratum2=="Lower Nome"])
# D = 0.43333, p-value = 0.3659
ksplot(bc_cap2$upstream[bc_cap2$Stratum2=="Lower Nome"], bc_cap2_recaps$upstream[bc_cap2_recaps$Stratum2=="Lower Nome"],
       main="Event 2 - Lower Nome", legend=c("All","Recaps"), col=c(1,4),
       xlab="Upstream Position (rkm)")

par(mfrow=c(2,2))
ks.test(bc_cap1$upstream[bc_cap1$Stratum2=="2000 Study"], bc_cap1_recaps$upstream[bc_cap1_recaps$Stratum2=="2000 Study"])
# D = 0.13215, p-value = 0.205
ksplot(bc_cap1$upstream[bc_cap1$Stratum2=="2000 Study"], bc_cap1_recaps$upstream[bc_cap1_recaps$Stratum2=="2000 Study"],
       main="Event 1 - 2000 Study", legend=c("All","Recaps"), col=c(1,2),
       xlab="Upstream Position (rkm)")

ks.test(bc_cap2$upstream[bc_cap2$Stratum2=="2000 Study"], bc_cap2_recaps$upstream[bc_cap2_recaps$Stratum2=="2000 Study"])
# D = 0.072168, p-value = 0.919
ksplot(bc_cap2$upstream[bc_cap2$Stratum2=="2000 Study"], bc_cap2_recaps$upstream[bc_cap2_recaps$Stratum2=="2000 Study"],
       main="Event 2 - 2000 Study", legend=c("All","Recaps"), col=c(1,4),
       xlab="Upstream Position (rkm)")

ks.test(bc_cap1$upstream[bc_cap1$Stratum2=="Upper Nome"], bc_cap1_recaps$upstream[bc_cap1_recaps$Stratum2=="Upper Nome"])
# D = 0.32609, p-value = 0.3363
ksplot(bc_cap1$upstream[bc_cap1$Stratum2=="Upper Nome"], bc_cap1_recaps$upstream[bc_cap1_recaps$Stratum2=="Upper Nome"],
       main="Event 2 - Upper Nome", legend=c("All","Recaps"), col=c(1,4),
       xlab="Upstream Position (rkm)")

ks.test(bc_cap2$upstream[bc_cap2$Stratum2=="Upper Nome"], bc_cap2_recaps$upstream[bc_cap2_recaps$Stratum2=="Upper Nome"])
# D = 0.27016, p-value = 0.5519
ksplot(bc_cap2$upstream[bc_cap2$Stratum2=="Upper Nome"], bc_cap2_recaps$upstream[bc_cap2_recaps$Stratum2=="Upper Nome"],
       main="Event 2 - Upper Nome", legend=c("All","Recaps"), col=c(1,4),
       xlab="Upstream Position (rkm)")


## confirming this with our standard chi^2 tests
## first using Stratum2 (effectively)
library(recapr)  # for automated consistency tests

with(subset(bc_all, event=="mark"),
     table(cut(Site, breaks=c(0,15,23,28,32,35,42))))

chisq_strat <- rep(NA, nrow(bc_all))
chisq_strat[bc_all$seg==7 | bc_all$seg==8] <- 1
chisq_strat[bc_all$seg==5 | bc_all$seg==6| bc_all$seg==4] <- 2
chisq_strat[bc_all$Stratum2=="Lower Nome"] <- 3
chisq_strat[bc_all$Stratum2=="2000 Study"] <- 4
chisq_strat[bc_all$Stratum2=="Upper Nome"] <- 5

table(chisq_strat)
n1 <- table(chisq_strat[bc_all$event=="mark"])
n2 <- table(chisq_strat[bc_all$event=="recap"])

recap_tags <- bc_cap1_recaps$Tag
m2strata1 <- chisq_strat[bc_all$Tag %in% recap_tags & bc_all$event=="mark"][order(bc_all$Tag[bc_all$Tag %in% recap_tags & bc_all$event=="mark"])]
m2strata2 <- chisq_strat[bc_all$Tag %in% recap_tags & bc_all$event=="recap"][order(bc_all$Tag[bc_all$Tag %in% recap_tags & bc_all$event=="recap"])]

consistencytest(n1=n1, n2=n2, m2strata1=m2strata1, m2strata2 = m2strata2)

with(bc_all, {
  n1 <- table(chisq_strat[event=="mark"])
  n2 <- table(chisq_strat[event=="recap"])

  m2strata1 <- chisq_strat[Tag %in% recap_tags & event=="mark"][order(Tag[Tag %in% recap_tags & event=="mark"])]
  m2strata2 <- chisq_strat[Tag %in% recap_tags & event=="recap"][order(Tag[Tag %in% recap_tags & event=="recap"])]

  consistencytest(n1=n1, n2=n2, m2strata1=m2strata1, m2strata2 = m2strata2)
})  # tests are not satisfied

bc_all$chisq_strat <- chisq_strat
with(subset(bc_all, Stratum1=="Beaver"), {
  n1 <- table(chisq_strat[event=="mark"])
  n2 <- table(chisq_strat[event=="recap"])

  m2strata1 <- chisq_strat[Tag %in% recap_tags & event=="mark"][order(Tag[Tag %in% recap_tags & event=="mark"])]
  m2strata2 <- chisq_strat[Tag %in% recap_tags & event=="recap"][order(Tag[Tag %in% recap_tags & event=="recap"])]

  consistencytest(n1=n1, n2=n2, m2strata1=m2strata1, m2strata2 = m2strata2)
}) # tests 2 and 3 are satisfied
with(subset(bc_all, Stratum1=="Nome"), {
  n1 <- table(chisq_strat[event=="mark"])
  n2 <- table(chisq_strat[event=="recap"])

  m2strata1 <- chisq_strat[Tag %in% recap_tags & event=="mark"][order(Tag[Tag %in% recap_tags & event=="mark"])] %>% as.factor %>% as.numeric
  m2strata2 <- chisq_strat[Tag %in% recap_tags & event=="recap"][order(Tag[Tag %in% recap_tags & event=="recap"])]  %>% as.factor %>% as.numeric

  consistencytest(n1=n1, n2=n2, m2strata1=m2strata1, m2strata2 = m2strata2)
}) # test 3 is satisfied


## trying the same tests (2 and 3 anyway) by hand, with different breakpoints
range(bc_all$upstream)
range(bc_all$upstream[bc_all$Stratum1=="Beaver"])
range(bc_all$upstream[bc_all$Stratum1=="Nome"])
range(bc_cap1_recaps$upstream[bc_cap1_recaps$Stratum1=="Beaver"])
range(bc_cap1_recaps$upstream[bc_cap1_recaps$Stratum1=="Nome"])
range(bc_cap2_recaps$upstream[bc_cap2_recaps$Stratum1=="Beaver"])
range(bc_cap2_recaps$upstream[bc_cap2_recaps$Stratum1=="Nome"])

par(mfrow=c(2,2))
beaverbreak <- seq(2, 40)
beaverpval <- minn <- NA*beaverbreak
for(i in seq_along(beaverbreak)) {
  row2 <- table(cut(bc_cap1_recaps$upstream[bc_cap1_recaps$Stratum1=="Beaver"], breaks=c(0, beaverbreak[i], 100)))
  row1 <- table(cut(bc_cap1$upstream[bc_cap1$Stratum1=="Beaver"], breaks=c(0, beaverbreak[i], 100))) - row2
  minn[i] <- min(rbind(row1,row2))
  beaverpval[i] <- suppressWarnings(chisq.test(rbind(row1,row2))$p.value)
}
plot(beaverbreak, beaverpval, type="b", log="y", main="Beaver Event 1")
abline(h=.05)
plot(minn, beaverpval, type="b", log="y", main="Beaver Event 1")
abline(h=.05)
abline(v=10)   # this case is completely fine

beaverpval <- minn <- NA*beaverbreak
for(i in seq_along(beaverbreak)) {
  row2 <- table(cut(bc_cap2_recaps$upstream[bc_cap2_recaps$Stratum1=="Beaver"], breaks=c(0, beaverbreak[i], 100)))
  row1 <- table(cut(bc_cap2$upstream[bc_cap2$Stratum1=="Beaver"], breaks=c(0, beaverbreak[i], 100))) - row2
  minn[i] <- min(rbind(row1,row2))
  beaverpval[i] <- suppressWarnings(chisq.test(rbind(row1,row2))$p.value)
}
plot(beaverbreak, beaverpval, type="b", log="y", main="Beaver Event 2")
abline(h=.05)
plot(minn, beaverpval, type="b", log="y", main="Beaver Event 2")
abline(h=.05)
abline(v=10)  # small pvals only come about with small sample sizes - probably fine

nomebreak <- seq(50, 76)
nomepval <- minn <- NA*nomebreak
for(i in seq_along(nomebreak)) {
  row2 <- table(cut(bc_cap1_recaps$upstream[bc_cap1_recaps$Stratum1=="Nome"], breaks=c(0, nomebreak[i], 100)))
  row1 <- table(cut(bc_cap1$upstream[bc_cap1$Stratum1=="Nome"], breaks=c(0, nomebreak[i], 100))) - row2
  minn[i] <- min(rbind(row1,row2))
  nomepval[i] <- suppressWarnings(chisq.test(rbind(row1,row2))$p.value)
}
plot(nomebreak, nomepval, type="b", log="y", main="Nome Event 1")
abline(h=.05)
plot(minn, nomepval, type="b", log="y", main="Nome Event 1")
abline(h=.05)
abline(v=10)  # this one is problematic, but it's the same stratum/event that was problematic with size comp

nomepval <- minn <- NA*nomebreak
for(i in seq_along(nomebreak)) {
  row2 <- table(cut(bc_cap2_recaps$upstream[bc_cap2_recaps$Stratum1=="Nome"], breaks=c(0, nomebreak[i], 100)))
  row1 <- table(cut(bc_cap2$upstream[bc_cap2$Stratum1=="Nome"], breaks=c(0, nomebreak[i], 100))) - row2
  minn[i] <- min(rbind(row1,row2))
  nomepval[i] <- suppressWarnings(chisq.test(rbind(row1,row2))$p.value)
}
plot(nomebreak, nomepval, type="b", log="y", main="Nome Event 2")
abline(h=.05)
plot(minn, nomepval, type="b", log="y", main="Nome Event 2")
abline(h=.05)
abline(v=10)  # small pvals only come about with small sample sizes - probably fine


## checking to see if there is still length selectivity detected after stratifying by sample
### THIS ACCOUNTS FOR MUCH OF THE LENGTH SELECTIVITY!!!
# par(mfrow=c(2,2))
# ks.test(bc_cap1$Length[bc_cap1$sample=="float"], bc_cap1_recaps$Length[bc_cap1_recaps$sample=="float"])
# # D = 0.10712, p-value = 0.7409
# ksplot(bc_cap1$Length[bc_cap1$sample=="float"], bc_cap1_recaps$Length[bc_cap1_recaps$sample=="float"],
#        main="Event 1 - float", legend=c("All","Recaps"), col=c(1,2),
#        xlab="Fork Length (mm)")
#
# ks.test(bc_cap2$Length[bc_cap2$sample=="float"], bc_cap2_recaps$Length[bc_cap2_recaps$sample=="float"])
# # D = 0.075859, p-value = 0.9699
# ksplot(bc_cap2$Length[bc_cap2$sample=="float"], bc_cap2_recaps$Length[bc_cap2_recaps$sample=="float"],
#        main="Event 2 - float", legend=c("All","Recaps"), col=c(1,4),
#        xlab="Fork Length (mm)")
#
# par(mfrow=c(2,2))
# ks.test(bc_cap1$Length[bc_cap1$sample=="hike"], bc_cap1_recaps$Length[bc_cap1_recaps$sample=="hike"])
# # D = 0.14723, p-value = 0.07561
# ksplot(bc_cap1$Length[bc_cap1$sample=="hike"], bc_cap1_recaps$Length[bc_cap1_recaps$sample=="hike"],
#        main="Event 1 - hike", legend=c("All","Recaps"), col=c(1,2),
#        xlab="Fork Length (mm)")
#
# ks.test(bc_cap2$Length[bc_cap2$sample=="hike"], bc_cap2_recaps$Length[bc_cap2_recaps$sample=="hike"])
# # D = 0.137, p-value = 0.1373
# ksplot(bc_cap2$Length[bc_cap2$sample=="hike"], bc_cap2_recaps$Length[bc_cap2_recaps$sample=="hike"],
#        main="Event 2 - hike", legend=c("All","Recaps"), col=c(1,4),
#        xlab="Fork Length (mm)")

par(mfrow=c(2,2))
ks.test(bc_cap1$Length[bc_cap1$Stratum1=="Beaver"], bc_cap1_recaps$Length[bc_cap1_recaps$Stratum1=="Beaver"])
# D = 0.10478, p-value = 0.7664
ksplot(bc_cap1$Length[bc_cap1$Stratum1=="Beaver"], bc_cap1_recaps$Length[bc_cap1_recaps$Stratum1=="Beaver"],
       main="Event 1 - Beaver", legend=c("All","Recaps"), col=c(1,2),
       xlab="Fork Length (mm)")

ks.test(bc_cap2$Length[bc_cap2$Stratum1=="Beaver"], bc_cap2_recaps$Length[bc_cap2_recaps$Stratum1=="Beaver"])
# D = 0.067667, p-value = 0.9869
ksplot(bc_cap2$Length[bc_cap2$Stratum1=="Beaver"], bc_cap2_recaps$Length[bc_cap2_recaps$Stratum1=="Beaver"],
       main="Event 2 - Beaver", legend=c("All","Recaps"), col=c(1,4),
       xlab="Fork Length (mm)")

par(mfrow=c(2,2))
ks.test(bc_cap1$Length[bc_cap1$Stratum1=="Nome"], bc_cap1_recaps$Length[bc_cap1_recaps$Stratum1=="Nome"])
# D = 0.16239, p-value = 0.03623
ksplot(bc_cap1$Length[bc_cap1$Stratum1=="Nome"], bc_cap1_recaps$Length[bc_cap1_recaps$Stratum1=="Nome"],
       main="Event 1 - Nome", legend=c("All","Recaps"), col=c(1,2),
       xlab="Fork Length (mm)")

ks.test(bc_cap2$Length[bc_cap2$Stratum1=="Nome"], bc_cap2_recaps$Length[bc_cap2_recaps$Stratum1=="Nome"])
# D = 0.14644, p-value = 0.09537
ksplot(bc_cap2$Length[bc_cap2$Stratum1=="Nome"], bc_cap2_recaps$Length[bc_cap2_recaps$Stratum1=="Nome"],
       main="Event 2 - Nome", legend=c("All","Recaps"), col=c(1,4),
       xlab="Fork Length (mm)")

par(mfrow=c(2,2))
ks.test(bc_cap1$Length[bc_cap1$Stratum2=="Lower Nome"], bc_cap1_recaps$Length[bc_cap1_recaps$Stratum2=="Lower Nome"])
# D = 0.31818, p-value = 0.755
ksplot(bc_cap1$Length[bc_cap1$Stratum2=="Lower Nome"], bc_cap1_recaps$Length[bc_cap1_recaps$Stratum2=="Lower Nome"],
       main="Event 1 - Lower Nome", legend=c("All","Recaps"), col=c(1,2),
       xlab="Fork Length (mm)")

ks.test(bc_cap2$Length[bc_cap2$Stratum2=="Lower Nome"], bc_cap2_recaps$Length[bc_cap2_recaps$Stratum2=="Lower Nome"])
# D = 0.3, p-value = 0.8104
ksplot(bc_cap2$Length[bc_cap2$Stratum2=="Lower Nome"], bc_cap2_recaps$Length[bc_cap2_recaps$Stratum2=="Lower Nome"],
       main="Event 2 - Lower Nome", legend=c("All","Recaps"), col=c(1,4),
       xlab="Fork Length (mm)")

par(mfrow=c(2,2))
ks.test(bc_cap1$Length[bc_cap1$Stratum2=="2000 Study"], bc_cap1_recaps$Length[bc_cap1_recaps$Stratum2=="2000 Study"])
# D = 0.15153, p-value = 0.1002
ksplot(bc_cap1$Length[bc_cap1$Stratum2=="2000 Study"], bc_cap1_recaps$Length[bc_cap1_recaps$Stratum2=="2000 Study"],
       main="Event 1 - 2000 Study", legend=c("All","Recaps"), col=c(1,2),
       xlab="Fork Length (mm)")

ks.test(bc_cap2$Length[bc_cap2$Stratum2=="2000 Study"], bc_cap2_recaps$Length[bc_cap2_recaps$Stratum2=="2000 Study"])
# D = 0.14634, p-value = 0.1611
ksplot(bc_cap2$Length[bc_cap2$Stratum2=="2000 Study"], bc_cap2_recaps$Length[bc_cap2_recaps$Stratum2=="2000 Study"],
       main="Event 2 - 2000 Study", legend=c("All","Recaps"), col=c(1,4),
       xlab="Fork Length (mm)")

par(mfrow=c(2,2))
ks.test(bc_cap1$Length[bc_cap1$Stratum2=="Upper Nome"], bc_cap1_recaps$Length[bc_cap1_recaps$Stratum2=="Upper Nome"])
# D = 0.17935, p-value = 0.9433
ksplot(bc_cap1$Length[bc_cap1$Stratum2=="Upper Nome"], bc_cap1_recaps$Length[bc_cap1_recaps$Stratum2=="Upper Nome"],
       main="Event 1 - Upper Nome", legend=c("All","Recaps"), col=c(1,2),
       xlab="Fork Length (mm)")

ks.test(bc_cap2$Length[bc_cap2$Stratum2=="Upper Nome"], bc_cap2_recaps$Length[bc_cap2_recaps$Stratum2=="Upper Nome"])
# D = 0.28629, p-value = 0.494
ksplot(bc_cap2$Length[bc_cap2$Stratum2=="Upper Nome"], bc_cap2_recaps$Length[bc_cap2_recaps$Stratum2=="Upper Nome"],
       main="Event 2 - Upper Nome", legend=c("All","Recaps"), col=c(1,4),
       xlab="Fork Length (mm)")




# ## curious to see how samples compare
# ksplot(bc_cap1$Length[bc_cap1$sample=="float"], bc_cap1$Length[bc_cap1$sample=="hike"],
#        main="Event 1 - all", legend=c("float","hike"), lty=c(1,2),
#        xlab="Fork Length (mm)")
# ksplot(bc_cap2$Length[bc_cap2$sample=="float"], bc_cap2$Length[bc_cap2$sample=="hike"],
#        main="Event 2 - all", legend=c("float","hike"), lty=c(1,2),
#        xlab="Fork Length (mm)")
# ksplot(bc_cap1_recaps$Length[bc_cap1_recaps$sample=="float"], bc_cap1_recaps$Length[bc_cap1_recaps$sample=="hike"],
#        main="Event 1 - recaps", legend=c("float","hike"), col=c(2,4),
#        xlab="Fork Length (mm)")
# ksplot(bc_cap2_recaps$Length[bc_cap2_recaps$sample=="float"], bc_cap2_recaps$Length[bc_cap2_recaps$sample=="hike"],
#        main="Event 2 - recaps", legend=c("float","hike"), col=c(2,4),
#        xlab="Fork Length (mm)")



## ---------------- growth recruitment? ----------------- ##
length1 <- bc_cap1_recaps$Length[order(bc_cap1_recaps$Tag)]
length2 <- bc_cap2_recaps_justtags$Length[order(bc_cap2_recaps_justtags$Tag)]

diffs <- length2 - length1

# checking that the order thing above actually worked - it does!
bc_cap1_recaps$Tag[order(bc_cap1_recaps$Tag)] ==
bc_cap2_recaps_justtags$Tag[order(bc_cap2_recaps_justtags$Tag)]

diffs <- diffs[abs(diffs) < 80]  # what happens when we remove that big outlier

plot(diffs)
plot(diffs,
     col=1+as.numeric(as.factor((paste(bc_cap1_recaps$event,
                                       bc_cap1_recaps$sample))[order(bc_cap1_recaps$Tag)])))
plot(diffs,
     col=1+as.numeric(as.factor((paste(bc_cap2_recaps_justtags$event,
                                       bc_cap2_recaps_justtags$sample))[order(bc_cap2_recaps_justtags$Tag)])))
median(diffs)
sd(diffs)
t.test(diffs)

# table combination of beaver/nome mark/recap, then plot diffs
recaps <- data.frame(Tag = bc_cap1_recaps$Tag[order(bc_cap1_recaps$Tag)],
                     length_cap1 = bc_cap1_recaps$Length[order(bc_cap1_recaps$Tag)],
                     length_cap2 = bc_cap2_recaps_justtags$Length[order(bc_cap2_recaps_justtags$Tag)],
                     sample_cap1 = bc_cap1_recaps$sample[order(bc_cap1_recaps$Tag)],
                     sample_cap2 = bc_cap2_recaps_justtags$sample[order(bc_cap2_recaps_justtags$Tag)],
                     upstream_cap1 = bc_cap1_recaps$upstream[order(bc_cap1_recaps$Tag)],
                     upstream_cap2 = bc_cap2_recaps_justtags$upstream[order(bc_cap2_recaps_justtags$Tag)],
                     Stratum1_cap1 = bc_cap1_recaps$Stratum1[order(bc_cap1_recaps$Tag)],
                     Stratum1_cap2 = bc_cap2_recaps_justtags$Stratum1[order(bc_cap2_recaps_justtags$Tag)],
                     Stratum2_cap1 = bc_cap1_recaps$Stratum2[order(bc_cap1_recaps$Tag)],
                     Stratum2_cap2 = bc_cap2_recaps_justtags$Stratum2[order(bc_cap2_recaps_justtags$Tag)])
recaps$diff <- recaps$length_cap2 - recaps$length_cap1
table(recaps$sample_cap1, recaps$sample_cap2)
par(mfrow=c(1,1))
recaps %>%
  ggplot(aes(y=diff, x=seq_along(diff),col=sample_cap1)) +
  geom_point()
recaps %>%
  ggplot(aes(y=diff, x=sample_cap1)) +
  geom_boxplot()



table(recaps$Stratum1_cap1, recaps$Stratum1_cap2)  # one problem fish
table(recaps$Stratum2_cap1, recaps$Stratum2_cap2)  # four problem fish

# endpoints of strata
riverlengths <- cumsum(rev(bc_river$lengths))/1000
Stratum1_bdy <- riverlengths[4]
Stratum2_bdy <- c(riverlengths[4], airstrip_upstream, riverlengths[6])

### looking at where the problem fish were with respect to boundaries
## one problem fish wrt Stratum1  ---------
subset(recaps, Stratum1_cap1=="Beaver" & Stratum1_cap2=="Nome")
recaps$upstream_cap1[recaps$Stratum1_cap1=="Beaver" & recaps$Stratum1_cap2=="Nome"] - Stratum1_bdy
recaps$upstream_cap2[recaps$Stratum1_cap1=="Beaver" & recaps$Stratum1_cap2=="Nome"] - Stratum1_bdy
# much further downstream than upstream: -0.46 vs 0.05 rkm
# call it Beaver (tag 515)

## four problem fish wrt Stratum2  ---------
subset(recaps, Stratum2_cap1=="2000 Study" & Stratum2_cap2=="Upper Nome")
recaps$upstream_cap1[recaps$Stratum2_cap1=="2000 Study" & recaps$Stratum2_cap2=="Upper Nome"] - Stratum2_bdy[3]
recaps$upstream_cap2[recaps$Stratum2_cap1=="2000 Study" & recaps$Stratum2_cap2=="Upper Nome"] - Stratum2_bdy[3]
# much further downstream than upstream: -10.36 vs 1.13
# call it 2000 Study (tag 1358)

subset(recaps, Stratum2_cap1=="Beaver" & Stratum2_cap2=="Lower Nome")
# much further downstream than upstream: -0.46 vs 0.05 rkm
# call it Beaver (tag 515)  - this is the same fish as earlier

subset(recaps, Stratum2_cap1=="Lower Nome" & Stratum2_cap2=="2000 Study")
recaps$upstream_cap1[recaps$Stratum2_cap1=="Lower Nome" & recaps$Stratum2_cap2=="2000 Study"] - Stratum2_bdy[2]
recaps$upstream_cap2[recaps$Stratum2_cap1=="Lower Nome" & recaps$Stratum2_cap2=="2000 Study"] - Stratum2_bdy[2]
# further upstream than downstream: 7.13 vs -1.6 rkm
# can call it 2000 Study (tag 1378)

subset(recaps, Stratum2_cap1=="Upper Nome" & Stratum2_cap2=="2000 Study")
recaps$upstream_cap1[recaps$Stratum2_cap1=="Upper Nome" & recaps$Stratum2_cap2=="2000 Study"] - Stratum2_bdy[3]
recaps$upstream_cap2[recaps$Stratum2_cap1=="Upper Nome" & recaps$Stratum2_cap2=="2000 Study"] - Stratum2_bdy[3]
# further upstream than downstream: 4.61 vs -0.59 rkm
# call it Upper Nome (tag 751)



# # should try a Darroch.  I don't want to use a Darroch, but should try it
# library(recapr)
# NDarroch(n1=unname(table(bc_cap1$Stratum2)),
#          n2=unname(table(bc_cap2$Stratum2)),
#          stratamat=as.matrix(table(recaps$Stratum2_cap1, recaps$Stratum2_cap2)))
# # does not seem to work and I don't know why


## ----------------- movement in/out? ------------------ ##
up1 <- bc_cap1_recaps$upstream[order(bc_cap1_recaps$Tag)]
up2 <- bc_cap2_recaps_justtags$upstream[order(bc_cap2_recaps_justtags$Tag)]
diffs <- up2-up1
plot(diffs)
sd(diffs)
mean(abs(diffs))
median(abs(diffs))
mean(abs(diffs)>.5)
par(mfrow=c(1,1))
plot(ecdf(abs(diffs)), xlim=c(0,5))

par(mfrow=c(1,1))
plot(x=up1, y=diffs, asp=1,
     main="Upriver Position - Mark Event (rkm)",
     xlab="", ylab="Movement (rkm)")
plot(x=up2, y=diffs, asp=1,
     main="Upriver Position - Recapture Event (rkm)",
     xlab="", ylab="Movement (rkm)")  # hardly any fish near endpoints, it's probably just fine

plot(x=up1, y=up2)

plot(NA,
     xlim=range(up1, up2),
     ylim=c(1, length(up1)))
segments(x0=up1[order((up1+up2)/2)],
         x1=up2[order((up1+up2)/2)],
         y0=seq_along(up1))
points(x=up1[order((up1+up2)/2)], y=seq_along(up1))
points(x=up2[order((up1+up2)/2)], y=seq_along(up1))

plot(NA,
     ylim=range(up1, up2),
     xlim=c(1, length(up1)))
segments(y0=up1[order((up1+up2)/2)],
         y1=up2[order((up1+up2)/2)],
         x0=seq_along(up1))
points(y=up1[order((up1+up2)/2)], x=seq_along(up1))
points(y=up2[order((up1+up2)/2)], x=seq_along(up1))

# try resampling distances and adding them to all capture locations to see how
# inferences are affected
# at the very least, I can see what happens to the tagged population subject to recapture:

nboot <- 10000
count_table <- matrix(nrow=nboot, ncol=4)
for(i in 1:nboot) {
  upstream_boot <- bc_cap1$upstream + sample(diffs, size=nrow(bc_cap1), replace=TRUE)
  # table(cut(upstream_boot,
  #           c(-100, 0, 42.81, 82.52, 200)))
  #           # levels=c("(-100,0]",   " (0,42.8]", "(42.8,82.5]",  "(82.5,200])")))
  count_table[i,1] <- sum(upstream_boot <= 0)
  count_table[i,2] <- length(upstream_boot %s_inside(]% c(0, 42.81))
  count_table[i,3] <- length(upstream_boot %s_inside(]% c(42.81, 82.5))
  count_table[i,4] <- sum(upstream_boot > 82.5)
}
prop_table <- count_table/nrow(bc_cap1)
colMeans(prop_table)
# we might lose 1% of tagged fish out the bottom, and the top is negligible









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
# m2 <- table(recaps$sample1)
# n1 <- table(bc_cap1$sample)
# n2 <- table(bc_cap2$sample)
#
# library(recapr)
# NBailey(n1=n1, n2=n2, m2=m2)
# seBailey(n1=n1, n2=n2, m2=m2)
#
# seBailey(n1=n1, n2=n2, m2=m2)/NBailey(n1=n1, n2=n2, m2=m2)
#
# # investigate the possibility of movement at ends & between strata?  simulate somehow
#
#
#
#
# # do the ASL_table thing
# library(dsftools)
#
# lengthbreaks <- c(250, 270, 330, 400)
#
# # just using first event sample
# ASL_table(stratum = as.numeric(as.factor(bc_cap1$sample)),
#           # length = bc_cap1$Length,
#           age = cut(bc_cap1$Length, lengthbreaks, right=FALSE),
#           Nhat = as.numeric(NBailey(n1=n1, n2=n2, m2=m2)),
#           se_Nhat = as.numeric(seBailey(n1=n1, n2=n2, m2=m2)))
#
# # just using second event sample
# ASL_table(stratum = as.numeric(as.factor(bc_cap2$sample)),
#           # length = bc_cap2$Length,
#           age = cut(bc_cap2$Length, lengthbreaks, right=FALSE),
#           Nhat = as.numeric(NBailey(n1=n1, n2=n2, m2=m2)),
#           se_Nhat = as.numeric(seBailey(n1=n1, n2=n2, m2=m2)))
#
# # pooling both events
# # ASL_table(stratum = as.numeric(as.factor(c(bc_cap1$sample, bc_cap2$sample))),
# #           # length = c(bc_cap1$Length, bc_cap2$Length),
# #           age = cut(c(bc_cap1$Length, bc_cap2$Length), lengthbreaks, right=FALSE),
# #           Nhat = as.numeric(NBailey(n1=n1, n2=n2, m2=m2)),
# #           se_Nhat = as.numeric(seBailey(n1=n1, n2=n2, m2=m2)))
#
# ASL_table(stratum = as.numeric(as.factor(bc_all$sample)),
#           # length = bc_all$Length,
#           age = cut(bc_all$Length, lengthbreaks, right=FALSE),
#           Nhat = as.numeric(NBailey(n1=n1, n2=n2, m2=m2)),
#           se_Nhat = as.numeric(seBailey(n1=n1, n2=n2, m2=m2)))
#
#
# # pooling both events, but separating by capture stratum
# with(subset(bc_all, sample == "float"),
#      ASL_table(
#        # length = Length,
#        age = cut(Length, lengthbreaks, right=FALSE),
#        Nhat = as.numeric(NBailey(n1=n1, n2=n2, m2=m2))[1],
#        se_Nhat = as.numeric(seBailey(n1=n1, n2=n2, m2=m2))[1]))
#
# with(subset(bc_all, sample == "hike"),
#      ASL_table(
#        # length = Length,
#        age = cut(Length, lengthbreaks, right=FALSE),
#        Nhat = as.numeric(NBailey(n1=n1, n2=n2, m2=m2))[2],
#        se_Nhat = as.numeric(seBailey(n1=n1, n2=n2, m2=m2))[2]))
