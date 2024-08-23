## 2. Checking Mark-Recapture Assumptions

# sourcing data
# relevant data objects:
# - bc_all
# - bc_cap1 & bc_cap2
# - bc_cap1_recaps & bc_cap2_recaps
source("R/1_BCMR_data.R")



## --------------- length selectivity ------------------ ##

# KS tests for size selectivity
ks.test(bc_cap1$Length, bc_cap1_recaps$Length)  # cap vs recap for first event
# D = 0.1758, p-value = 0.001329
# evidence of size selectivity in second event

par(mfrow=c(2,2))
plot(density(bc_cap1$Length, na.rm=T))
lines(density(bc_cap1_recaps$Length, na.rm=T), col=2)
plot(ecdf(bc_cap1$Length))
plot(ecdf(bc_cap1_recaps$Length), add=T, col=2)


ks.test(bc_cap2$Length, bc_cap2_recaps$Length)  # cap vs recap for second event
# D = 0.16438, p-value = 0.002293
# evidence of size selectivity in first event

plot(density(bc_cap2$Length, na.rm=T))
lines(density(bc_cap2_recaps$Length, na.rm=T), col=4)
plot(ecdf(bc_cap2$Length))
plot(ecdf(bc_cap2_recaps$Length), add=T, col=4)


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
plot(possiblebreaks, p1, ylim=range(p1,p2,p3,p4), type='l', log="y", col=2)
lines(possiblebreaks, p2, col=3)
lines(possiblebreaks, p3, col=4)
lines(possiblebreaks, p4, col=5)
abline(h=.05)


## trying a break at 270

sum(bc_cap1_recaps$Length < 270)
sum(bc_cap2_recaps$Length < 270)

par(mfrow=c(2,2))
plot(density(bc_cap1$Length[bc_cap1$Length < 270], na.rm=T))
lines(density(bc_cap1_recaps$Length[bc_cap1_recaps$Length < 270], na.rm=T), col=2)
plot(ecdf(bc_cap1$Length[bc_cap1$Length < 270]))
plot(ecdf(bc_cap1_recaps$Length[bc_cap1_recaps$Length < 270]), add=T, col=2)

plot(density(bc_cap2$Length[bc_cap2$Length < 270], na.rm=T))
lines(density(bc_cap2_recaps$Length[bc_cap2_recaps$Length < 270], na.rm=T), col=4)
plot(ecdf(bc_cap2$Length[bc_cap2$Length < 270]))
plot(ecdf(bc_cap2_recaps$Length[bc_cap2_recaps$Length < 270]), add=T, col=4)

plot(density(bc_cap1$Length[bc_cap1$Length >= 270], na.rm=T))
lines(density(bc_cap1_recaps$Length[bc_cap1_recaps$Length >= 270], na.rm=T), col=2)
plot(ecdf(bc_cap1$Length[bc_cap1$Length >= 270]))
plot(ecdf(bc_cap1_recaps$Length[bc_cap1_recaps$Length >= 270]), add=T, col=2)

plot(density(bc_cap2$Length[bc_cap2$Length >= 270], na.rm=T))
lines(density(bc_cap2_recaps$Length[bc_cap2_recaps$Length >= 270], na.rm=T), col=4)
plot(ecdf(bc_cap2$Length[bc_cap2$Length >= 270]))
plot(ecdf(bc_cap2_recaps$Length[bc_cap2_recaps$Length >= 270]), add=T, col=4)


## trying breaks at 270 and 300
## tests fail for 300-500 but not 270-500!  i hate ks tests

breaks <- c(250, 270, 500)# , 300
bc_cap1$Lengthbin <- cut(bc_cap1$Length, breaks, right=FALSE)
bc_cap2$Lengthbin <- cut(bc_cap2$Length, breaks, right=FALSE)
bc_cap1_recaps$Lengthbin <- cut(bc_cap1_recaps$Length, breaks, right=FALSE)
bc_cap2_recaps$Lengthbin <- cut(bc_cap2_recaps$Length, breaks, right=FALSE)
for(i in 1:(length(breaks)-1)) {
  x1 <- bc_cap1$Length[as.numeric(bc_cap1$Lengthbin) == i]
  x2 <- bc_cap1_recaps$Length[as.numeric(bc_cap1_recaps$Lengthbin) == i]
  x3 <- bc_cap2$Length[as.numeric(bc_cap2$Lengthbin) == i]
  x4 <- bc_cap2_recaps$Length[as.numeric(bc_cap2_recaps$Lengthbin) == i]

  print(ks.test(x1, x2))  # cap vs recap for first event

  par(mfrow=c(2,2))
  plot(density(x1, na.rm=T))
  lines(density(x2, na.rm=T), col=2)
  plot(ecdf(x1))
  plot(ecdf(x2), add=T, col=2)

  print(ks.test(x3, x4))  # cap vs recap for first event

  plot(density(x3, na.rm=T))
  lines(density(x4, na.rm=T), col=4)
  plot(ecdf(x3))
  plot(ecdf(x4), add=T, col=4)
}


## -------------- spatial selectivity? ---------------- ##
ks.test(bc_cap1$upstream, bc_cap1_recaps$upstream)
# D = 0.22825, p-value = 8.778e-06

ks.test(bc_cap2$upstream, bc_cap2_recaps$upstream)
# D = 0.49899, p-value < 2.2e-16

## actually I don't think we need to worry about this one




## ---------------- growth recruitment? ----------------- ##
length1 <- bc_cap1_recaps$Length[order(bc_cap1_recaps$Tag)]
length2 <- bc_cap2_recaps$Length[order(bc_cap2_recaps$Tag)]

# checking that the order thing above actually worked - it does!
bc_cap1_recaps$Tag[order(bc_cap1_recaps$Tag)] ==
bc_cap2_recaps$Tag[order(bc_cap2_recaps$Tag)]

plot(length2 - length1)
plot(length2 - length1,
     col=1+as.numeric(as.factor((paste(bc_cap1_recaps$event, bc_cap1_recaps$sample))[order(bc_cap1_recaps$Tag)])))
plot(length2 - length1,
     col=1+as.numeric(as.factor((paste(bc_cap2_recaps$event, bc_cap2_recaps$sample))[order(bc_cap2_recaps$Tag)])))
median(length2-length1)
sd(length2-length1)
t.test(length2-length1)




## ----------------- movement in/out? ------------------ ##
up1 <- bc_cap1_recaps$upstream[order(bc_cap1_recaps$Tag)]
up2 <- bc_cap2_recaps$upstream[order(bc_cap2_recaps$Tag)]
plot(up2-up1)
sd(up2-up1)
mean(abs(up2-up1))
median(abs(up2-up1))
