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

# plot(density(bc_cap1$Length, na.rm=T),
#      main="Event 1: All vs Recaptures", xlab="")
# legend("topright", lty=1, col=c(1,2), legend=c("All","Recaps"))
# lines(density(bc_cap1_recaps$Length, na.rm=T), col=2)
# plot(ecdf(bc_cap1$Length),
#      main="Event 1: All vs Recaptures", xlab="")
# legend("bottomright", lty=1, col=c(1,2), legend=c("All","Recaps"))
# plot(ecdf(bc_cap1_recaps$Length), add=T, col=2)


# ks.test(bc_cap2$Length, bc_cap2_recaps$Length)  # cap vs recap for second event
# # D = 0.16438, p-value = 0.002293
# # evidence of size selectivity in first event
#
# plot(density(bc_cap2$Length, na.rm=T),
#      main="Event 2: All vs Recaptures", xlab="")
# legend("topright", lty=1, col=c(1,4), legend=c("All","Recaps"))
# lines(density(bc_cap2_recaps$Length, na.rm=T), col=4)
# plot(ecdf(bc_cap2$Length),
#      main="Event 2: All vs Recaptures", xlab="")
# legend("bottomright", lty=1, col=c(1,4), legend=c("All","Recaps"))
# plot(ecdf(bc_cap2_recaps$Length), add=T, col=4)


# actually very curious how cap-cap and recap-recap compare
par(mfrow=c(2,2))
ksplot(x1=bc_cap1$Length, x2=bc_cap2$Length,
       main="All Captures", legend=c("Event 1","Event 2"), lty=c(1,2),
       xlab="Fork Length (mm)")
ksplot(x1=bc_cap1_recaps$Length, x2=bc_cap2_recaps$Length,
       main="Recaptures", legend=c("Event 1","Event 2"), col=c(2,4),
       xlab="Fork Length (mm)")

# ks.test(bc_cap1$Length, bc_cap2$Length)  # cap for first event vs cap for second
# D = 0.057789, p-value = 0.03195

# plot(density(bc_cap1$Length, na.rm=T),
#      main="Event 1 All vs Event 2 All", xlab="")
# legend("topright", lty=c(1,2), legend=c("Event 1", "Event 2"))
# lines(density(bc_cap2$Length, na.rm=T), lty=2)
# plot(ecdf(bc_cap1$Length),
#      main="Event 1 All vs Event 2 All", xlab="")
# legend("bottomright", lty=c(1,2), legend=c("Event 1", "Event 2"))
# plot(ecdf(bc_cap2$Length), add=T, lty=2)
#
#
# ks.test(bc_cap1_recaps$Length, bc_cap2_recaps$Length)  # recap for first event vs recap for second
# # D = 0.051852, p-value = 0.9934
#
# plot(density(bc_cap1_recaps$Length, na.rm=T),
#      main="Event 1 Recaps vs Event 2 Recaps", xlab="", col=2)
# legend("topright", col=c(2,4), legend=c("Event 1", "Event 2"))
# lines(density(bc_cap2_recaps$Length, na.rm=T), col=4)
# plot(ecdf(bc_cap1_recaps$Length), col=2,
#      main="Event 1 Recaps vs Event 2 Recaps", xlab="")
# legend("bottomright", col=c(2,4), lty=1, legend=c("Event 1", "Event 2"))
# plot(ecdf(bc_cap2_recaps$Length), add=T, col=4)





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
par(mfrow=c(2,2))
ks.test(bc_cap1$upstream[bc_cap1$sample=="float"], bc_cap1_recaps$upstream[bc_cap1_recaps$sample=="float"])
# D = 0.13742, p-value = 0.428
ksplot(bc_cap1$upstream[bc_cap1$sample=="float"], bc_cap1_recaps$upstream[bc_cap1_recaps$sample=="float"],
       main="Event 1 - float", legend=c("All","Recaps"), col=c(1,2),
       xlab="Upstream Position (rkm)")

ks.test(bc_cap2$upstream[bc_cap2$sample=="float"], bc_cap2_recaps$upstream[bc_cap2_recaps$sample=="float"])
# D = 0.1145, p-value = 0.6445
ksplot(bc_cap2$upstream[bc_cap2$sample=="float"], bc_cap2_recaps$upstream[bc_cap2_recaps$sample=="float"],
       main="Event 2 - float", legend=c("All","Recaps"), col=c(1,4),
       xlab="Upstream Position (rkm)")

par(mfrow=c(2,2))
ks.test(bc_cap1$upstream[bc_cap1$sample=="hike"], bc_cap1_recaps$upstream[bc_cap1_recaps$sample=="hike"])
# D = 0.11804, p-value = 0.2432
ksplot(bc_cap1$upstream[bc_cap1$sample=="hike"], bc_cap1_recaps$upstream[bc_cap1_recaps$sample=="hike"],
       main="Event 1 - hike", legend=c("All","Recaps"), col=c(1,2),
       xlab="Upstream Position (rkm)")

ks.test(bc_cap2$upstream[bc_cap2$sample=="hike"], bc_cap2_recaps$upstream[bc_cap2_recaps$sample=="hike"])
# D = 0.11701, p-value = 0.2823
ksplot(bc_cap2$upstream[bc_cap2$sample=="hike"], bc_cap2_recaps$upstream[bc_cap2_recaps$sample=="hike"],
       main="Event 2 - hike", legend=c("All","Recaps"), col=c(1,4),
       xlab="Upstream Position (rkm)")




## checking to see if there is still length selectivity detected after stratifying by sample
### THIS ACCOUNTS FOR MUCH OF THE LENGTH SELECTIVITY!!!
par(mfrow=c(2,2))
ks.test(bc_cap1$Length[bc_cap1$sample=="float"], bc_cap1_recaps$Length[bc_cap1_recaps$sample=="float"])
# D = 0.10712, p-value = 0.7409
ksplot(bc_cap1$Length[bc_cap1$sample=="float"], bc_cap1_recaps$Length[bc_cap1_recaps$sample=="float"],
       main="Event 1 - float", legend=c("All","Recaps"), col=c(1,2),
       xlab="Fork Length (mm)")

ks.test(bc_cap2$Length[bc_cap2$sample=="float"], bc_cap2_recaps$Length[bc_cap2_recaps$sample=="float"])
# D = 0.075859, p-value = 0.9699
ksplot(bc_cap2$Length[bc_cap2$sample=="float"], bc_cap2_recaps$Length[bc_cap2_recaps$sample=="float"],
       main="Event 2 - float", legend=c("All","Recaps"), col=c(1,4),
       xlab="Fork Length (mm)")

par(mfrow=c(2,2))
ks.test(bc_cap1$Length[bc_cap1$sample=="hike"], bc_cap1_recaps$Length[bc_cap1_recaps$sample=="hike"])
# D = 0.14723, p-value = 0.07561
ksplot(bc_cap1$Length[bc_cap1$sample=="hike"], bc_cap1_recaps$Length[bc_cap1_recaps$sample=="hike"],
       main="Event 1 - hike", legend=c("All","Recaps"), col=c(1,2),
       xlab="Fork Length (mm)")

ks.test(bc_cap2$Length[bc_cap2$sample=="hike"], bc_cap2_recaps$Length[bc_cap2_recaps$sample=="hike"])
# D = 0.137, p-value = 0.1373
ksplot(bc_cap2$Length[bc_cap2$sample=="hike"], bc_cap2_recaps$Length[bc_cap2_recaps$sample=="hike"],
       main="Event 2 - hike", legend=c("All","Recaps"), col=c(1,4),
       xlab="Fork Length (mm)")


## curious to see how samples compare
ksplot(bc_cap1$Length[bc_cap1$sample=="float"], bc_cap1$Length[bc_cap1$sample=="hike"],
       main="Event 1 - all", legend=c("float","hike"), lty=c(1,2),
       xlab="Fork Length (mm)")
ksplot(bc_cap2$Length[bc_cap2$sample=="float"], bc_cap2$Length[bc_cap2$sample=="hike"],
       main="Event 2 - all", legend=c("float","hike"), lty=c(1,2),
       xlab="Fork Length (mm)")
ksplot(bc_cap1_recaps$Length[bc_cap1_recaps$sample=="float"], bc_cap1_recaps$Length[bc_cap1_recaps$sample=="hike"],
       main="Event 1 - recaps", legend=c("float","hike"), col=c(2,4),
       xlab="Fork Length (mm)")
ksplot(bc_cap2_recaps$Length[bc_cap2_recaps$sample=="float"], bc_cap2_recaps$Length[bc_cap2_recaps$sample=="hike"],
       main="Event 2 - recaps", legend=c("float","hike"), col=c(2,4),
       xlab="Fork Length (mm)")



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

# table combination of float/hike mark/recap, then plot diffs
recaps <- data.frame(Tag = bc_cap1_recaps$Tag[order(bc_cap1_recaps$Tag)],
                     length1 = bc_cap1_recaps$Length[order(bc_cap1_recaps$Tag)],
                     length2 = bc_cap2_recaps_justtags$Length[order(bc_cap2_recaps_justtags$Tag)],
                     sample1 = bc_cap1_recaps$sample[order(bc_cap1_recaps$Tag)],
                     sample2 = bc_cap2_recaps_justtags$sample[order(bc_cap2_recaps_justtags$Tag)],
                     upstream1 = bc_cap1_recaps$upstream[order(bc_cap1_recaps$Tag)],
                     upstream2 = bc_cap2_recaps_justtags$upstream[order(bc_cap2_recaps_justtags$Tag)])
recaps$diff <- recaps$length2 - recaps$length1
table(recaps$sample1, recaps$sample2)
par(mfrow=c(1,1))
recaps %>%
  ggplot(aes(y=diff, x=seq_along(diff),col=sample1)) +
  geom_point()
recaps %>%
  ggplot(aes(y=diff, x=sample1)) +
  geom_boxplot()



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
