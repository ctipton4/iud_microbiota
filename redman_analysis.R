
#####
#colors
#####
cb <-  scale_color_manual(values=c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))
col = scale_fill_manual(values = c(brewer_pal(8, "Dark2"), brewer_pal(8, "Accent")))
col1 = scale_colour_manual(values = c(brewer_pal(8, "Dark2"), brewer_pal(8, "Accent")))

#hexcodes
col2 = scale_fill_manual(values = c("#FFFF00", "#996633", "#33CC33"))
col3 = scale_colour_manual(values = c("#FFFF00", "#996633", "#33CC33"))
#by name
colors()
col4 = scale_fill_manual(values = c("thistle", "violet", "wheat"))
col5 = scale_colour_manual(values = c("thistle", "violet", "wheat"))

library(gridExtra)
cbbFill2 <- scale_fill_manual(values=c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00","#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "thistle4", "yellowgreen", "violetred1"))
cbbColour2 <- scale_colour_manual(values=c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00"))
gray2 = scale_colour_manual(values = c("black", "gray50"))
gray2_fill = scale_fill_manual(values = c("black", "gray50"))

col30 = scale_fill_manual(values = c(rev(c(brewer_pal(palette = "Dark2")(8))), c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FF0000", "#B3DE69", "#FCCDE5", "#00008B", "#BC80BD", "#CCEBC5", "#FFED6F","#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#333333", "#66C2A5", "#3288BD", "#5E4FA2", "#CCCCCC")))
col30c = scale_colour_manual(values = c(rev(c(brewer_pal(palette = "Dark2")(8))), c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FF0000", "#B3DE69", "#FCCDE5", "#00008B", "#BC80BD", "#CCEBC5", "#FFED6F","#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#333333", "#66C2A5", "#3288BD", "#5E4FA2", "#CCCCCC")))


sample_data(d)$Group = ordered(sample_data(d)$Group, levels = c("Normal", "Abnormal"))
sample_data(d)$Patient = ordered(sample_data(d)$Patient, levels = sort(as.numeric(as.character(unique(met$Patient)))))

### this bit of code is useful for subsetting, so not really needed here
# d = d 
# d_sub = subset_samples(d, sample_data(d)$Complete == "yes")

met = data.frame(sample_data(d))


list = row.names(met)

otu = subset(otu, row.names(otu) %in% list)
gg = subset(gen_p, row.names(gen_p) %in% list)
sp = subset(spp_p, row.names(spp_p) %in% list)
# ff = subset(fam_p, row.names(fam_p) %in% list)
# pp = subset(phy_p, row.names(phy_p) %in% list)
unif_uw = unif_uw[list,list]
unif_w = unif_w[list,list]

reads = rowSums(otu)
# write.csv(file = "reads.csv", reads)
summary(reads)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 29963   39149   45677   43142   47409   52506

### gtsummary for demographic table
library(gtsummary)
library(kableExtra)
library(flextable)

met$bmi30 = 0
met$bmi30[met$BMI > 30] = 1
saveRDS(met, file = "dat.rda")

tab1 = met %>% select(Group, age, bmi30, Gravida, Parity, Brand, DurationYears)
tbl_summary(tab1, by = Group, label =
              list(age ~ "Age", bmi30 ~ "BMI >30",
                   DurationYears ~ "Duration (Years)"),
            type = list(Gravida ~ "continuous",
                             Parity ~ "continuous",
                        DurationYears ~ "continuous")) %>%
  bold_labels() %>%
  add_p() %>% 
  add_overall() %>%
  modify_caption(caption = "Table 1. Patient Demographics and Reproductive History") %>%
  modify_header(all_stat_cols() ~ "**{level}**\n*N = {n}*")

tab1 = met %>% select(Group, Bleeding, PelvicPain, 
                      DeviceComplication, ExpiredIUD, HxVaginitis, ToConcieve, Other, Complaints)
tbl_summary(tab1, by = Group, label =
              list(Bleeding ~ "Abnormal Bleeding", ExpiredIUD ~ "Expired IUD",
                   Complaints ~ "Total Reasons/Complaints",
                   PelvicPain ~ "Pelvic Pain",
                   HxVaginitis ~ "History of Vaginitis",
                   ToConcieve ~ "To Concieve"),
            type = list(Complaints ~ "continuous")) %>%
  bold_labels() %>%
  add_overall() %>%
  modify_caption(caption = "Table 2. Reported Reasons and Complaints Related to IUD Removal") %>%
  modify_header(all_stat_cols() ~ "**{level}**\n*N = {n}*")

### a few packages I use frequently but may not call below
library(ggplot2)
library(reshape2)
library(vegan)

R = estimate_richness(d) # calculates a suite of alpha div measures. 
dtemp = cbind(R, sample_data(d))
dtemp$Hill = exp(dtemp$Shannon) # for evenness based testing, I tend to prefer Hill1 Diversity which is calculated from Shannon

write.csv(file = "alphadiv.csv", dtemp) # you can use this output to make figures in other programs if you like
#####
#Rarefaction
#####
# r = rarerichness(otu, "sample", from = 1000, to = 20000, by = 1000)
# r = join(r, met, by = "sample")
# write.csv(file="rarefaction.csv", r)
# p = ggplot(r, aes(variable, value, group = sample)) + geom_point()  + theme_bw(base_size = 10) + 
#   geom_line()+ scale_x_discrete("Number of Sequences") + scale_y_continuous("OTU Richness") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + cb
# png('rare.png', res = 600, width = 6, height = 3, units = 'in')
# print(p)
# dev.off()

#####
# demographic table
#####
library(gtsummary)
## gtsummary table for demographics
tab1 = met %>% select(Group, age, BMI, Brand, DurationYears, Parity)
tbl_summary(tab1, by = Group,
                   type = c(DurationYears ~ "continuous",
                            Parity ~ "continuous")) %>%
  bold_labels() %>%
  add_p() %>% 
  add_overall() %>%
  bold_p() %>%
  modify_caption(caption = "Table 1. Patient Demographics per each Group") %>%
  modify_header(all_stat_cols() ~ "**{level}**\n*N = {n}*")

tab1 = met %>% select(Group, Bleeding, PelvicPain, DeviceComplication, ExpiredIUD, HxVaginitis, ToConcieve,
                      Other)
tbl_summary(tab1, by = Group) %>%
  bold_labels() %>%
  add_p() %>% 
  bold_p() %>%
  modify_caption(caption = "Table 2. Cited Reasons and Complaints Leading to IUD Removal") %>%
  modify_header(all_stat_cols() ~ "**{level}**\n*N = {n}*")


#####
# bacterial load analysis
#####
bl = dtemp
bl = bl[!is.na(bl$Load),]

m = aov(log(Load) ~ Group, bl)
summary(m)
t.test(log10(Load) ~ Group, bl, alternative = "less", var.equal = T)

library(scales)
library(ggpubr)
p4 = ggplot(bl, aes(Group, Load)) + geom_boxplot() + 
  geom_jitter(height = 0, width = 0.25, size = 2, alpha = 0.6) +
  theme_bw() + labs(x = "IUD Removal Cause", y = "Bacterial Load (16S copies)") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
    theme(axis.title = element_text(size = 9, face = "bold"), axis.text = element_text(size = 9)) +
  stat_compare_means(method = 't.test', method.args = list(alternative = "greater", var.equal = T),
                     vjust = 2)

# png('load_group.png', res = 600, width = 3, height = 3, units = 'in')
# print(p)
# dev.off()

## sum genus table to Lacto and non-lacto, then transform by load
lact = gg
lact$NL = rowSums(lact) - lact$g__Lactobacillus
lact = lact[,c("g__Lactobacillus", "NL")]
names(lact) = c("Lactobacillus", "NL")

# do groups segregate by lacto/non-lacto relative abundance?
ll = lact[which(row.names(lact) != "7.357wf.806r"),] # for absolute abundance below
lact = cbind(met, lact)

wilcox.test(Lactobacillus ~ Group, lact)
wilcox.test(NL ~ Group, lact)

p = ggplot(lact, aes(Group, Lactobacillus)) + geom_boxplot() + 
  geom_jitter(height = 0, width = 0.25, size = 2, alpha = 0.6) +
  theme_bw() + labs(x = "IUD Removal Cause", y = "Lactobacillus Relative Abundance") +
  theme(axis.title = element_text(size = 9, face = "bold"), axis.text = element_text(size = 9)) +
  stat_compare_means(method = 't.test', method.args = list(alternative = "greater", var.equal = T),
                     vjust = 3)
p5 = p
# png('lacto_ra_group.png', res = 600, width = 3, height = 3, units = 'in')
# print(p)
# dev.off()

p = ggplot(lact, aes(Group, NL)) + geom_boxplot() + 
  geom_jitter(height = 0, width = 0.25, size = 2, alpha = 0.6) +
  theme_bw() + labs(x = "IUD Removal Cause", y = "Non-Lactobacillus Relative Abundance") +
  theme(axis.title = element_text(size = 9, face = "bold"), axis.text = element_text(size = 8))
# png('nonlacto_ra_group.png', res = 600, width = 3, height = 3, units = 'in')
# print(p)
# dev.off()

# ## ratio ra plot?
# p = ggplot(lact, aes(Group, (NL/Lactobacillus))) + geom_boxplot() + 
#   geom_jitter(height = 0, width = 0.25, size = 2, alpha = 0.6) +
#   theme_bw() + labs(x = "IUD Removal Cause", y = "Lactobacillus Relative Abundance") +
#   theme(axis.title = element_text(size = 9, face = "bold"), axis.text = element_text(size = 8))
# p5 = p
# png('lacto_ra_group.png', res = 600, width = 3, height = 3, units = 'in')
# print(p)
# dev.off()

# do groups segregate by lacto/non-lacto absolute abundance?

ll = ll*bl$Load
ll = cbind(bl, ll)

t.test(log10(Lactobacillus+1) ~ Group, ll, alternative = "less", var.equal = T)
t.test(log10(NL) ~ Group, ll, alternative = "less", var.equal = T)

p = ggplot(ll, aes(Group, Lactobacillus+1)) + geom_boxplot(outlier.shape = NA) + 
  geom_jitter(height = 0, width = 0.25, size = 2, alpha = 0.6) +
  theme_bw() + labs(x = "IUD Removal Cause", y = "Estimated Lactobacillus Abundance") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme(axis.title = element_text(size = 9, face = "bold"), axis.text = element_text(size = 9)) +
  stat_compare_means(method = 't.test', method.args = list(alternative = "greater", var.equal = T),
                     vjust = 1.5)
p6 = p
# png('lacto_group.png', res = 600, width = 3, height = 3, units = 'in')
# print(p)
# dev.off()

p = ggplot(ll, aes(Group, NL)) + geom_boxplot() + 
  geom_jitter(height = 0, width = 0.25, size = 2, alpha = 0.6) +
  theme_bw() + labs(x = "IUD Removal Cause", y = "Estimated Non-Lactobacillus Abundance") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme(axis.title = element_text(size = 9, face = "bold"), axis.text = element_text(size = 9)) +
  stat_compare_means(method = 't.test', method.args = list(alternative = "greater", var.equal = T),
                     vjust = 2)
p7 = p
# png('nonlacto_group.png', res = 600, width = 3, height = 3, units = 'in')
# print(p)
# dev.off()

## merge figure 3 together
library(patchwork)
tiff('Fig3.tiff', res = 600, width = 6, height = 6, units = 'in', compression = 'lzw')
print((p4 + p5)/(p6 + p7) + plot_annotation(tag_levels = "a"))
dev.off()


####
#Alpha diversity
#####

### formal analysis

## quick exploratory model
library(car)
m = aov(Observed ~ Group + age + BMI + Duration + Mirena + Complaints, dtemp)
m = aov(Observed ~ age, dtemp)
Anova(m)

m = aov(Hill ~ Group + age + BMI + Duration + Mirena + Complaints, dtemp)
m = aov(Hill ~ BMI, dtemp)
Anova(m)

# briefly check assumptions of normality and variance structure
hist(dtemp$Observed)
hist(dtemp$Hill)
# sample size small, but may violate normality assumptions
plot(dtemp$age, dtemp$Observed)
plot(dtemp$Duration, dtemp$Observed)
plot(dtemp$Complaints, dtemp$Observed)
plot(log(dtemp$Load), dtemp$Observed)
ggplot(dtemp, aes(Group, Observed)) + geom_boxplot()
# ggplot(dtemp, aes(Ethnicity, Observed)) + geom_boxplot()
# ggplot(dtemp, aes(Group, Observed)) + geom_boxplot() + facet_grid(~Ethnicity)
## literature suggests ethnicity may drive large pH and thus microbiome differences
# hispanics overrepresented compared to others. drop white and black (n =2)?

plot(log2(dtemp$Load), dtemp$Shannon)

plot(dtemp$age, dtemp$Hill)
plot(dtemp$Duration, dtemp$Hill)
plot(dtemp$Complaints, dtemp$Hill)

## variance may scale with group and/or age. 

summary(aov(age ~ Group, dtemp))
ggplot(dtemp, aes(Group, age)) + geom_boxplot()
summary(aov(Duration ~ Group, dtemp))
ggplot(dtemp, aes(Group, Duration)) + geom_boxplot() + geom_jitter(height = 0)

library(lmtest)
# m = lm(Observed ~ Group, dtemp)
# Anova(m)
# bptest(m)
# dtemp$otu_res = m$residuals
# 
# m = lm(Hill ~ Group, dtemp)
# Anova(m)
# bptest(m)
# dtemp$hill_res = m$residuals

## age seems to be more important factor, and confounded with group. Weight age in regression model

m = lm(Observed ~ age, dtemp)
Anova(m)
bptest(m) # significant result indicates that variance in Observed scales with Age
dtemp$otu_res = m$residuals

# m = lm(Hill ~ age, dtemp)
# Anova(m)
# bptest(m)
# dtemp$hill_res = m$residuals

## calculate weight for regression model
w = 1/abs(dtemp$otu_res)
m2 = lm(Observed ~ age, dtemp, weights = w)
summary(m2)

## make new dtemp for dummy coding
dtemp2 = dtemp
dtemp2$Group = gsub("Abnormal", 1, dtemp2$Group)
dtemp2$Group = gsub("Normal", 0, dtemp2$Group)
dtemp2$Group = as.numeric(dtemp2$Group)

## repeat backward selection with new weighted model
# Group + age + BMI + Duration + Mirena + Complaints

m2 = lm(Observed ~ age + Group, dtemp2, weights = w)
summary(m2)
## less informative, use age alone
m2 = lm(Observed ~ age, dtemp, weights = w)
summary(m2)

### plotting


p = ggplot(dtemp, aes(Group, Observed)) + geom_jitter(size = 2, height = 0, width = 0.3) +
  theme_bw() + geom_boxplot(alpha = 0.4, outlier.shape = NA) +
  theme(axis.text.x = element_text(size = 9, face = "bold"), axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 9, face = "bold")) +
  labs(x = "", y = "OTUs Observed")
p1 = p
# png('otu_boxplot_group.png', res = 600, width = 6, height = 8, units = 'in')
# print(p)
# dev.off()

p = ggplot(dtemp, aes(age, abs(otu_res))) + geom_point(size = 2) + 
  geom_smooth(method = "lm", color = "black", linetype = "dashed", size = 0.75) + 
  theme_bw() + labs(x = "Patient Age (years)", y = "Residual OTU Variance") + 
  theme(axis.title = element_text(size = 9, face = "bold"), axis.text = element_text(size = 8)) +
  geom_hline(yintercept = 0, linetype = "dashed")
png('residual_plot.png', res = 600, width = 3, height = 3, units = 'in')
print(p)
dev.off()

p = ggplot(dtemp, aes(age, Observed)) + geom_point(size = 2) + 
  geom_smooth(method = "lm", color = "black", linetype = "dashed", size = 0.75, se = F) + 
  theme_bw() + labs(x = "Patient Age (years)", y = "OTUs Observed") + 
  theme(axis.title = element_text(size = 9, face = "bold"), axis.text = element_text(size = 8))
p3 = p

## merge plots together
# will need to have p1-3 made, use Ctrl+f feature to quickly find each plot and make adjustments 
library(patchwork)
png('Fig2.png', res = 600, width = 8, height = 3.5, units = 'in')
print(p1 + p2 + p3 + plot_annotation(tag_levels = "a"))
dev.off()


## Hill
p = ggplot(dtemp, aes(Group, Hill)) + geom_jitter(aes(color = Group), size = 4, height = 0, width = 0.3) +
  theme_bw((base_size = 20)) + geom_boxplot(aes(fill = Group), alpha = 0.4, outlier.shape = NA) +
  cbbColour2 + cbbFill2  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16), legend.position = "none")
p = p + xlab("") + ylab("Hill Diversity")
png('hill_boxplot_group.png', res = 600, width = 7, height = 8, units = 'in')
print(p)
dev.off()


#####
#barplots
#####
# pd = pd[,-c(9:24)]
pd = met
pd$Other = NULL

p = ra_bar(gg, pd, 19, "set6", "Genus", c("Group"), "Group")
p = p + theme_bw(base_size = 16) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 20)) + xlab("") +
  # facet_grid(Body_Site~Time, scales = "free")
  png('bar_genus_group.png', res = 600, width = 10, height = 10, units = 'in')
print(p)
dev.off()

p = ra_bar(gg, pd, 19, "set6", "Genus", c("Group", "Patient"), "Patient")
p = p + theme_bw(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16), legend.key.size = unit(0.75, "line")) + 
  xlab("") +
  facet_wrap(~Group, scales = "free", nrow = 1)
png('bar_genus_Patient.png', res = 600, width = 11, height = 6, units = 'in')
print(p)
dev.off()

## pub version
# cleanup genus level name for plot; too long as is.
gg2 = gg
names(gg)[grep("Rhodospirillales", names(gg))]
names(gg2) = gsub("o__Rhodospirillales; f__Unclassified; g__Unclassified", "o__Rhodospirillales", names(gg2))
names(gg)[grep("Clostridiales", names(gg))]
names(gg2) = gsub("o__Clostridiales; f__Unclassified; g__Unclassified", "o__Clostridiales", names(gg2))
names(gg)[grep("p__Unclassified", names(gg))]
names(gg2) = gsub("p__Unclassified; c__Unclassified; o__Unclassified; f__Unclassified; g__Unclassified", "p__Unclassified", names(gg2))
# use gg2 instead of gg for genus level barplot
p = ra_bar(gg2, pd, 19, "set6", "Genus", c("Group", "Patient"), "Patient")
p = p + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10), legend.key.size = unit(0.75, "line"),
        strip.text = element_text(size = 10, face = "bold")) + 
  xlab("") +
  facet_wrap(~Group, scales = "free", nrow = 1)
png('Fig1.png', res = 600, width = 6, height = 3.75, units = 'in')
print(p)
dev.off()

##

p = ra_bar(sp, pd, 19, "set6", "Species", c("Group"), "Group")
p = p + theme_bw(base_size = 16) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 20)) + xlab("") +
  # facet_grid(Body_Site~Time, scales = "free")
  png('bar_species_group.png', res = 600, width = 10, height = 10, units = 'in')
print(p)
dev.off()

p = ra_bar(sp, pd, 19, "set6", "Species", c("Group", "Patient"), "Patient")
p = p + theme_bw(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16), legend.key.size = unit(0.75, "line")) + 
  xlab("") +
  facet_wrap(~Group, scales = "free", nrow = 1)
png('bar_species_Patient.png', res = 600, width = 11, height = 6, units = 'in')
print(p)
dev.off()

## pub version
# cleanup genus level name for plot; too long as is.
sp2 = sp
names(sp)[grep("Rhodospirillales", names(sp))]
names(sp2) = gsub("o__Rhodospirillales; f__Unclassified; g__Unclassified; s__Unclassified", "o__Rhodospirillales", names(sp2))
# use sp2 instead of sp for genus level barplot
p = ra_bar(sp2, pd, 19, "set6", "Species", c("Group", "Patient"), "Patient")
p = p + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10), legend.key.size = unit(0.75, "line"),
        strip.text = element_text(size = 10, face = "bold")) + 
  xlab("") +
  facet_wrap(~Group, scales = "free", nrow = 1)
png('FigS3.png', res = 600, width = 6, height = 3.75, units = 'in')
print(p)
dev.off()

#####
#Beta diversity
#####
# Group + age + BMI + Duration + Mirena + Complaints

#WEIGHTED UNIFRAC
wuf = as.dist(unif_w)
a = adonis(as.dist(wuf) ~ met$age)
a
# no variables significant here.

pc = capscale(wuf~1, comm = gg)
cap = data.frame(pc$CA$u)
ar = cap[,1:2]
cap = merge(met, cap, by = "row.names")
rownames(cap) = cap[,1]
cap = cap[,2:ncol(cap)]
s = summary(eigenvals(pc))

p = ggplot(cap, aes(MDS1, MDS2)) + geom_point(aes(color = Group), size = 2, alpha = 0.7) + theme_bw() +
  gray2 + gray2_fill + stat_ellipse(aes(linetype = Group, color = Group), size = 0.75) + 
  labs(x = paste("Axis 1 (", round(s[2,1]*100, digits =2), "%)", sep = ''),
       y = paste("Axis 2 (", round(s[2,2]*100, digits =2), "%)", sep = ''),
       color = "", linetype = "") +
  theme(legend.position = "bottom", legend.text = element_text(size = 9, face = "bold"),
        axis.title = element_text(size = 9, face = "bold"), axis.text = element_text(size = 8))
p2 = p
# png(file= "pcoa_w_group.png", width = 10, height = 8, units = 'in', res = 600)
# print(p)
# dev.off()

# ### species biplot
# library(ggrepel)
# 
# c = t(cor(ar, gg )) #correlation of axes with variables of interest
# c = as.data.frame(na.omit(c))
# # dat = data.frame(dat, gender = p$data$gender, group = p$data$group) #add factors to axes data.frame (e.g. get from metadata)
# dists =   sqrt(c[,1]^2 + c[,2]^2) #see which genera show maximal combined correlation with axes
# keep =  order(dists, decreasing = T)[1:5] #keep the top 5
# postn = rep(0.5, length(keep)) #this is used to offset the labels- if x cor is negative, the hjust = 1, else hjust = 0 
# postn[c[keep,1]<0] = 1
# postn[c[keep,1]>0] = 0
# #xend and yend are multiplies by a fraction that is the scaling factor you'd like for plotting
# d.arrows = data.frame(x = rep(0,length(keep)),y = rep(0,length(keep)), xend = c[keep,1]*.5, yend = c[keep,2]*.5, labels = rownames(c)[keep], p = postn)
# d.arrows$labels = gsub("\\;.*", "", d.arrows$labels)
# 
# p = ggplot(cap, aes(MDS1, MDS2)) + geom_point(aes(color = Group), size = 5, alpha = 0.7) + theme_bw(base_size = 20) + 
#   cbbFill2 + cbbColour2 
# # p = p + labs(x = paste("Axis 1 (", round(s$importance[2,1]*100, digits =2), "%)", sep = ''),
# #              y = paste("Axis 2 (", round(s$importance[2,2]*100, digits =2), "%)", sep = ''))
# p = p + labs(x = paste("Axis 1 (", round(s[2,1]*100, digits =2), "%)", sep = ''),
#              y = paste("Axis 2 (", round(s[2,2]*100, digits =2), "%)", sep = ''))
# 
# find_hull <- function(df) df[chull(df$MDS1, df$MDS2), ]
# hulls <- ddply(p$data, c("Group"), find_hull)
# p = p + geom_polygon(data = hulls, aes(fill = Group), alpha = 0.15, linetype = 0)
# p = p + annotate("segment", x = d.arrows$x, xend = 0.4*d.arrows$xend, y = d.arrows$y, yend = 0.4*d.arrows$yend, size = 1.25, arrow = arrow(length=unit(0.2, "cm")))
# # p = p + geom_label(data = d.arrows, aes(label = str_wrap(labels, 10), x = (0 + 0.45*d.arrows$xend), y = (0 + 0.45*d.arrows$yend)))
# p = p + geom_label_repel(data = d.arrows, aes(x = 0.4*xend, y = 0.4*yend, label = str_wrap(labels)) , point.padding = NA)
# 
# png(file= "pcoa_w_group_ggarrows.png", width = 10, height = 8, units = 'in', res = 600)
# print(p)
# dev.off()


### UNWEIGHTED UNIFRAC
uwuf = as.dist(unif_uw)
a = adonis(uwuf ~ met$Group)
a


pc = capscale(uwuf~1, comm = gg)
cap = data.frame(pc$CA$u)
ar = cap[,1:2]
cap = merge(met, cap, by = "row.names")
rownames(cap) = cap[,1]
cap = cap[,2:ncol(cap)]
s = summary(eigenvals(pc))


p = ggplot(cap, aes(MDS1, MDS2)) + geom_point(aes(color = Group), size = 5, alpha = 0.7) + theme_bw(base_size = 20) +
  cbbFill2 + cbbColour2
p = p + labs(x = paste("Axis 1 (", round(s[2,1]*100, digits =2), "%)", sep = ''),
             y = paste("Axis 2 (", round(s[2,2]*100, digits =2), "%)", sep = ''))
# p = p + labs(x = paste("Axis 1 (", round(s$importance[2,1]*100, digits =2), "%)", sep = ''),
#              y = paste("Axis 2 (", round(s$importance[2,2]*100, digits =2), "%)", sep = ''))

find_hull <- function(df) df[chull(df$MDS1, df$MDS2), ]
hulls <- ddply(p$data, c("Group"), find_hull)
p = p + geom_polygon(data = hulls, aes(fill = Group), alpha = 0.15, linetype = 0)

png(file= "pcoa_uw_group.png", width = 10, height = 8, units = 'in', res = 600)
print(p)
dev.off()

# #####
# # ANCOM
# #####
# # d_s <- tax_glom(d, taxrank = "Species")
# # pd = sample_data(d_s)
# # #dp = transform_sample_counts(dp, function(x) x/sum(x))
# # gc = data.frame(t(otu_table(d_s)))
# # 
# # n = paste(tax_table(d_s)[,1]," ; ", 
# #           tax_table(d_s)[,2]," ; ",
# #           tax_table(d_s)[,3]," ; ",
# #           tax_table(d_s)[,4]," ; ",
# #           tax_table(d_s)[,5]," ; ",
# #           tax_table(d_s)[,6], sep = '')
# # n = clean_blast2_names(n, tax_level="species")
# # 
# # names(gc) = n
# # gc = subset(gc, row.names(gc) %in% list)
# 
# 
# write.csv(res, "ancom_results_injection.csv")
# print(xtable(res, caption="Taxa Detected by ANCOM for Injection difference"), 
#       file = 'ancom_res_injection.tex', caption.placement = "top")
# 
# a = pgc[,c(abs((length(names(amet))-1)-dim(pgc)[2]):dim(pgc)[2], which(names(pgc) %in% res$otu_ids))]
# m = melt(a, id.vars = names(amet))
# # m$variable = gsub("_Unknown", "", m$variable)
# m$variable = gsub("_g__Unclassified", "", m$variable)
# p = ggplot(m, aes(x = Injection, y = value)) + geom_boxplot(aes(fill = Injection), alpha = 0.4, outlier.shape = NA) +
#   geom_jitter(aes(color = Injection), size = 4, height = 0, width = 0.25, alpha = 0.6) + theme_bw((base_size = 13)) + 
#   theme(legend.position = "none") + 
#   cbbColour2 + cbbFill2 + ggtitle("ANCOM Detections for Differentially Abundant Taxa\nInjection") + 
#   theme(plot.title = element_text(hjust = 0.5))
# p = p + facet_wrap(~variable, scales = "free_y", labeller = label_wrap_gen(width = 10, multi_line = T), nrow = 6) + 
#   labs(x = NULL, y = "Relative Abundance")
# png('ancom_boxplot_injection.png', res = 600, width = 13, height = 12, units = 'in')
# print(p)
# dev.off()


#####
# Does within Lactobacillus species composition vary between groups?
#####
## create new object with Lactobacillus species
sp2 = sp
list = names(sp)[grep("Lactobacillus", names(sp))]
sp2 = sp2[,which(names(sp2) %in% list)]
slac = merge(sp2, met, by = "row.names")
row.names(slac) = slac$Row.names
slac$Row.names = NULL
## ancom to test for differences in species between groups
gc = sp2
keep3 = data.frame()
for (i in unique(names(gc))){
  if(length(which(gc[,i] > 0)) > 0.1*dim(gc)[1]){ keep1 = "yes"} else {
    keep1 = "no"
  }
  occurrences = length(which(gc[,i] > 0))
  keep2 = data.frame(i, keep1, occurrences)
  names(keep2) = c("otu_ids", "keep", "Observed")
  keep3 = rbind(keep3, keep2)
}
# res = merge(res, keep3, by = "otu_ids")
keep = keep3[which(keep3$keep == "yes"),]
gc = gc[,which(names(gc) %in% keep$otu_ids)]
names(gc) = gsub("; ", "_", names(gc))
names(gc) = gsub(" ", "_", names(gc))
keep$otu_ids = gsub("; ", "_", keep$otu_ids)
keep$otu_ids = gsub(" ", "_", keep$otu_ids)

# pgc = data.frame(merge(gc, amet ,by="row.names",all.y=T))

res = ancom.W(otu_data = gc, var_data = met, adjusted = F,repeated = F,
        main.var = "Group",adj.formula = NULL,
        repeat.var = NULL,long = F,rand.formula = NULL ,
        multcorr = 3,sig = 0.05)

res = merge(res, keep, by = "otu_ids")
res = res[,-4]

### no differences in relative abundance detected by ancom using moderate correction
### L. reuteri detected to be different using no correction and makes sense with plots. 

## plot to assess differences in relative abundance, absolute abundance and detection (0,1)
library(reshape2)
m = melt(slac, measure.vars = list)
names(m)[16] = "Species"
## absolute
m2 = m[!is.na(m$Load),]
m2$abs = (m2$Load*m2$value)
m2$abs2 = log10(m2$abs+1)
# ## detection
# sp3 = sp2
# sp3[(sp3 != 0)] = 1
# colSums(sp3)
# sp3 = merge(met, sp3, by = "row.names")
# m3 = melt(sp3, measure.vars = list)

## name cleanup for plots
m$Species = gsub("s__", "", m$Species)
m2$Species = gsub("s__", "", m2$Species)
m$Species = gsub("Lactobacillus", "L.", m$Species)
m2$Species = gsub("Lactobacillus", "L.", m2$Species)
## plotting lactobacillus species
p1 = ggplot(m, aes(Group, value)) + geom_boxplot(outlier.shape = NA, width = 0.5) + geom_jitter(width = 0.2, alpha = 0.7) +
  facet_wrap(~Species, scales = "free_y") + theme_bw() + labs(x = "", y = "Relative Abundance") +
  theme(axis.title = element_text(size = 9, face = "bold"), axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 9, face = "bold"), strip.text = element_text(size = 9, face = "bold.italic"))

p2 = ggplot(m2, aes(Group, abs2)) + geom_boxplot(outlier.shape = NA, width = 0.5) + geom_jitter(width = 0.2, alpha = 0.7) +
  facet_wrap(~Species, scales = "free_y") + theme_bw() + labs(x = "", y = "Estimated Abundance (log10 scale)") +
  theme(axis.title = element_text(size = 9, face = "bold"), axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 9, face = "bold"), strip.text = element_text(size = 9, face = "bold.italic"))
png('lacto_species.png', res = 600, width = 5.5, height = 6, units = 'in')
print(p1/p2 + plot_annotation(tag_levels = "a"))
dev.off()

#####
# heatmap - qualitative analysis on IUD types and symptoms
#####

who = names(sort(colSums(sp), decreasing = TRUE))[1:20]
f = sp[,names(sp) %in% who]

abx_cols = c("dodgerblue4", "dodgerblue", "#9ecae1", "#feeda0","#fec44f","#fdae6b", "#fe9929",
                          "darkorange2", "darkorange4"  ,"#9e9ac8", "lightcoral", "#B78560", "darkolivegreen",
                          "grey20",
                          "darkviolet","#74c476", "red2")
                          
library(ComplexHeatmap)
library(pheatmap)
library(viridis)
# ## pcr2 missing one name
# row.names(f)[!row.names(f) %in% row.names(pcr2)]
# pcr2
## names must be in consistent order going in
f = f[order(row.names(f)),]
names(f) = gsub("; .__Unclass.*", "", names(f))
# pcr2 = pcr2[order(row.names(pcr2)),]
tmp = met[,c(11:17)]
# arg3_short = arg3_short[order(row.names(arg3_short)),]

## wuf clustering?
hc = hclust(wuf)
plot(hc)
coldend = as.dendrogram(hc)

## species
ht1 = Heatmap(t(f), column_dend_height = unit(0.6, 'cm'),
              row_dend_width = unit(0.2, 'cm'),
              name = "Relative Abundance",
              column_names_rot = 45, column_names_gp = gpar(fontsize = 0),
              cluster_columns = coldend,
              row_names_gp = gpar(fontsize = 9),
              col = viridis(n = 4, option = "mako", direction = -1),
              ## annotations
              heatmap_legend_param = list(direction = 'horizontal',
                                          legend_width = unit(1.3, "in")),
              top_annotation = HeatmapAnnotation(Brand = met$Brand,
                                                 IUD_Removal = met$Group,
                                                 simple_anno_size = unit(0.3, 'cm'),
                                                 col = list(
                                                   Brand = c("Kyleena" = "lightcoral", "Lyletta" = "#9ecae1",
                                                             "Mirena" = "dodgerblue", "Paraguard" = "darkorange2", 
                                                             "Skyla" = "darkolivegreen"),
                                                   IUD_Removal = c("Normal" = "gray60", "Abnormal" = "red")),
                                                 annotation_name_gp = gpar(fontsize = 9)),
              row_title = "Top 20 Bacteria")

ht2 = Heatmap(t(tmp), 
        column_names_rot = 45, column_names_gp = gpar(fontsize = 0),
        cluster_rows = F,
        row_names_gp = gpar(fontsize = 9),
        col = viridis(n = 2, option = "mako", direction = -1),
        row_title = "Reasons for IUD Removal",
        name = "Reasons")

#####
# ht1 = Heatmap(t(f), column_split = met$Group,
#               column_dend_height = unit(0.4, 'cm'),
#               row_dend_width = unit(0.2, 'cm'),
#               name = "Relative Abundance",
#               column_names_rot = 45, column_names_gp = gpar(fontsize = 0),
#               row_names_gp = gpar(fontsize = 8),
#               col = viridis(n = 4, option = "mako", direction = -1),
#               ## annotations
#               heatmap_legend_param = list(direction = 'horizontal',
#                                           legend_width = unit(1.3, "in")),
#               top_annotation = HeatmapAnnotation(Infection = cmet3$PostOp_Infection,
#                                                  NGS = cmet3$NGS,
#                                                  Gender = cmet3$Gender,
#                                                  col = list(
#                                                    Infection = c("Yes" = "red"),
#                                                    NGS = c("Negative" = "skyblue", "Positive" = "dodgerblue4"),
#                                                    Gender = c("Female" = "lightcoral", "Male" = "darkcyan", 
#                                                               "Other" = "gray80")),
#                                                  simple_anno_size = unit(0.3, 'cm'),
#                                                  annotation_name_gp = gpar(fontsize = 9)),
#               row_title = "Species Detected by NGS")
# 
# ht4 = Heatmap(t(log10(pcr2+1)), column_split = cmet3$arm,
#               name = "qPCR Log Abundance",
#               row_dend_width = unit(0.2, 'cm'),
#               column_names_rot = 45, column_names_gp = gpar(fontsize = 0),
#               row_names_gp = gpar(fontsize = 8),
#               col = viridis(n = 4, option = "magma", direction = 1),
#               heatmap_legend_param = list(direction = 'horizontal',
#                                           legend_width = unit(1.3, 'in')),
#               row_title = "qPCR")
# # ## antibiotic determinations
# 
# abx_colors = structure(abx_cols[1:17], names = unique(abx_list[1:17]))
# ht2 = Heatmap(t(abx), column_split = cmet3$arm, 
#               column_names_rot = 45, column_names_gp = gpar(fontsize = 0),
#               row_names_gp = gpar(fontsize = 8),
#               col = abx_colors,
#               heatmap_legend_param = list(ncol = 3, title_position = "leftcenter-rot",
#                                           at = abx_list,
#                                           labels = abx_list[1:17], 
#                                           legend_gp = gpar(fill = 1:17)),
#               row_title = "",
#               name = "Antibiotic",
#               ## annotations
#               top_annotation = HeatmapAnnotation(ABXMatch = cmet3$ABXMatch,
#                                                  col = list(
#                                                    ABXMatch = c("Yes" = "red", "No" = "gray20")
#                                                  ),
#                                                  simple_anno_size = unit(0.3, 'cm'),
#                                                  annotation_name_gp = gpar(fontsize = 9))
# )
# ## qpcr resistance profile
# ht3 = Heatmap(t(arg3_short), column_split = cmet3$arm, 
#               column_names_rot = 45, column_names_gp = gpar(fontsize = 0),
#               row_names_gp = gpar(fontsize = 8),
#               col = viridis(n = 2, option = "mako", direction = -1),
#               row_title = "ARG",
#               name = "ARG Detection")
#####

ht_list = ht1 %v% ht2

tiff('heatmap_qual_analysis.tiff', res = 900, width = 6.75, height = 6, units = 'in', compression = 'lzw')
print(draw(ht_list, heatmap_legend_side = 'bottom'))
dev.off()








