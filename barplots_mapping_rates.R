library(dplyr)
library(ggplot2)

# reading in the mapping stats
dat = read.table("mapping_stats_all_isolates.txt")
str(dat)
dat = dat[c(2,3)]
head(dat)
colnames(dat) <- c("isolate","mapping to MG8")
dat$`mapping to MG8` = as.double(sub("%","",dat$`mapping to MG8`))

# to plot them in a certain order, we must transforme into a factor and fix the factor levels
dat$isolate <- factor(dat$isolate, levels = dat$isolate[order(dat$`mapping to MG8`,decreasing = T)])

dat %>% ggplot(aes(y=`mapping to MG8`, x = isolate)) + geom_bar(stat = "identity") + labs(y="percentage of reads", title = "Overall read mapping rates to the MG8 reference")
