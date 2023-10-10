######### Test comparing GE of top GST gene vs settlement #########

library(tidyverse)
library(readxl)
library(bayesnec)
library(ggpubr)
library(multcomp) # Dunnett's

####----Load and mutate data----####
data <- read_excel("GEvsSettlement_BayesnecInputData.xlsx", sheet="GE") %>%
  data.frame() 


datGE <- data %>% 
  mutate(x = as.numeric(as.character(x)), 
         x = ifelse(x==0, 0.1, x),     
         #log_x = log(x),
         sqrt_x = sqrt(x), 
         Tr = as.factor(Tr),
         GST_prop = (median(data[1:6,4])/GST_7641),
         pirin_prop = (median(data[1:6,5])/pirin_38837),
         epid_prop = (median(data[1:6,7])/epid_7393),
         prost_prop = (median(data[1:6,11])/prost_14113),
         acylCoA_prop = (median(data[1:6,13])/acylCoA_3029))

data2 <- read_excel("GEvsSettlement_BayesnecInputData.xlsx", sheet="settlement") %>%
  data.frame() 


datSettle <- data2 %>% 
  mutate(TPAH = as.numeric(as.character(TPAH)), 
         #TPAH = ifelse(TPAH==0, 0.1, TPAH),     
         Tr = as.numeric(as.character(Tr)),
         log_TPAH = log(TPAH),
         sqrt_TPAH = sqrt(TPAH),
         sqrt_Tr = sqrt(Tr),
         No.settled = as.integer(No.settled),
         Total = as.integer(Total))


# best genes to try deriving thresholds for are Glutathione-S-Transferase (GST), pirin, epimerase dehydratase (epid), prostaglandin reductase 1-like (prost), and acylCoA

pGE <- ggplot(datGE, aes(x=as.factor(x), y=GST_7641)) + geom_boxplot(outlier.colour="red", outlier.shape=8,
                                                                                 outlier.size=2, fill="#0000ff")+ theme_classic() +   
  labs(x="", y = "") +
  ggtitle("Glutathione S Transferase") +
  theme(plot.title = element_text(size = 10))
pGE

pGE2 <- ggplot(datGE, aes(x=as.factor(x), y=GST_prop)) + geom_boxplot(outlier.colour="red", outlier.shape=8,
                                                                      outlier.size=2, fill="#0000ff")+ theme_classic() +   
  labs(x="Treatment (% WAF)", y = "Proportion expression rel. to ctrl (inv.)")
pGE2

pGE3 <- ggplot(datGE, aes(x=as.factor(x), y=pirin_38837)) + geom_boxplot(outlier.colour="red", outlier.shape=8,
                                                                       outlier.size=2, fill="#990066")+ theme_classic() +   
  labs(x="", y = "") +
  ggtitle("Pirin") +
  theme(plot.title = element_text(size = 10))
  
pGE3

pGE4 <- ggplot(datGE, aes(x=as.factor(x), y=pirin_prop)) + geom_boxplot(outlier.colour="red", outlier.shape=8,
                                                                       outlier.size=2, fill="#990066")+ theme_classic() +   
  labs(x="Treatment (% WAF)", y = "Proportion expression rel. to ctrl (inv.)")
pGE4

pGE5 <- ggplot(datGE, aes(x=as.factor(x), y=friz_16133)) + geom_boxplot(outlier.colour="red", outlier.shape=8,
                                                                          outlier.size=2, fill="#86b1da")+ theme_classic() +   
  labs(x="Treatment (% WAF)", y = "Expression")
pGE5

pGE6 <- ggplot(datGE, aes(x=as.factor(x), y=friz_prop)) + geom_boxplot(outlier.colour="red", outlier.shape=8,
                                                                         outlier.size=2, fill="#86b1da")+ theme_classic() +   
  labs(x="Treatment (% WAF)", y = "Proportion expression rel. to ctrl (inv.)")
pGE6

pGE7 <- ggplot(datGE, aes(x=as.factor(x), y=epid_7393)) + geom_boxplot(outlier.colour="red", outlier.shape=8,
                                                                         outlier.size=2, fill="#009900")+ theme_classic() +   
  labs(x="", y = "") +
  ggtitle("Epimerase dehydratase") +
  theme(plot.title = element_text(size = 10))
pGE7

pGE8 <- ggplot(datGE, aes(x=as.factor(x), y=epi_prop)) + geom_boxplot(outlier.colour="red", outlier.shape=8,
                                                                        outlier.size=2, fill="#009900")+ theme_classic() +   
  labs(x="Treatment (% WAF)", y = "Proportion expression rel. to ctrl (inv.)")
pGE8

pGE9 <-  ggplot(datGE, aes(x=as.factor(x), y=prost_14113)) + geom_boxplot(outlier.colour="red", outlier.shape=8,
                                                                        outlier.size=2, fill="#6666cc")+ theme_classic() +   
  labs(x="", y = "") +
  ggtitle("Prostaglandin reductase 1-like") +
  theme(plot.title = element_text(size = 10))

pGE10 <- ggplot(datGE, aes(x=as.factor(x), y=prost_prop)) + geom_boxplot(outlier.colour="red", outlier.shape=8,
                                                                            outlier.size=2, fill="#6666cc")+ theme_classic() +   
  labs(x="Treatment (% WAF)", y = "Proportion expression rel. to ctrl (inv.)")

pGE11 <- ggplot(datGE, aes(x=as.factor(x), y=acylCoA_3029)) + geom_boxplot(outlier.colour="red", outlier.shape=8,
                                                                       outlier.size=2, fill="#ff9999")+ theme_classic() +   
  labs(x="", y = "") +
  ggtitle("Acyl-CoA synthetase bubblegum family member") +
  theme(plot.title = element_text(size = 10))

pGE12 <- ggplot(datGE, aes(x=as.factor(x), y=acylCoA_prop)) + geom_boxplot(outlier.colour="red", outlier.shape=8,
                                                                        outlier.size=2, fill="#ff9999")+ theme_classic() +   
  labs(x="Treatment (% WAF)", y = "Proportion expression rel. to ctrl (inv.)")


pSettle <- ggplot(datSettle2, aes(x=as.factor(Tr), y=p.settled)) + geom_boxplot(outlier.colour="red", outlier.shape=8,
                                                                      outlier.size=2, fill="#a6a6a6")+ theme_classic() +   
  labs(x="", y = "Percent Settlement")
pSettle

###---- Combined boxplots ----####

pGSTsettle <- ggarrange(pGE ,pGE2, pSettle, ncol = 1, nrow = 3,
                  align = "v")
pGSTsettle

pGEcomb <- ggarrange(pGE, pGE3, pGE7, pGE9, pGE11, pSettle, ncol = 2, nrow = 3,
                 align = "v")
ggpubr::annotate_figure(
  pGEcomb, left = text_grob(expression("Expression"), rot = 90),
  bottom = text_grob(expression("Treatment (% WAF)")))

pGEcomb

###---- Linear models and Dunnett's ----####

iso7641.lm <- lm(GST_prop ~ Tr, data = datGE) #using proportion from control data
plot(iso7641.lm, ask = F, which = 1:2) #Q-Q plot
plot(rstandard(iso7641.lm)) #residuals
anova(iso7641.lm) #F val = 117.01, p 2.2e-16
summary(glht(iso7641.lm, linfct = mcp(Tr="Dunnett"))) 


iso7641.lm2 <- lm(GST_7641 ~ Tr, data = datGE) #using raw vst expression data
plot(iso7641.lm2, ask = F, which = 1:2) #Q-Q plot
plot(rstandard(iso7641.lm2)) #residuals
anova(iso7641.lm2) #F val = 163.02, p 2.2e-16
summary(glht(iso7641.lm2, linfct = mcp(Tr="Dunnett"))) 



datSettle2 <- data2 %>% 
  mutate(Tr2=c(rep("A_0_WAF", 12), rep("B_2_WAF", 8), rep("C_5_WAF", 8), rep("D_10_WAF", 4), rep("E_20_WAF", 6), rep("F_50_WAF", 6), rep("G_70_WAF", 6), rep("H_100_WAF", 6)),
         Tr2=as.factor(Tr2))

settle.lm <- lm(p.settled ~ Tr2, data = datSettle2) 
plot(settle.lm, ask = F, which = 1:2) #Q-Q plot
plot(rstandard(settle.lm)) #residuals
anova(settle.lm) #F val = 12.295, p 6.957e-09
summary(glht(settle.lm, linfct = mcp(Tr2="Dunnett"))) 


# remove dilutions not used in GE assays

datSettle2 <- datSettle %>% 
  filter(Tr!=5, Tr!=70) %>% 
  mutate(Tr2=c(rep("A_0_WAF", 12), rep("B_2_WAF", 8), rep("C_10_WAF", 4), rep("D_20_WAF", 6), rep("E_50_WAF", 6), rep("F_100_WAF", 6)),
         Tr2=as.factor(Tr2))

settle.lm2 <- lm(p.settled ~ Tr2, data = datSettle2) 
plot(settle.lm2, ask = F, which = 1:2) #Q-Q plot
plot(rstandard(settle.lm2)) #residuals
anova(settle.lm2) #F val = 11.937, p 7.488e-07
summary(glht(settle.lm2, linfct = mcp(Tr2="Dunnett"))) 


pSettle2 <- ggplot(datSettle2, aes(x=as.factor(Tr), y=p.settled)) + geom_boxplot(outlier.colour="red", outlier.shape=8,
                                                                               outlier.size=2, fill="#a6a6a6")+ theme_classic() +   
  labs(x="Treatment (% WAF)", y = "Percent Settlement")

pSettle2

pGSTsettle2 <- ggarrange(pGE ,pGE2, pSettle, pSettle2, ncol=2, nrow = 2, 
                         align = "v")
pGSTsettle2

###---- Top DE GST gene isogroup7641 ----####

# variance stabilized data from DESeq2 as input expression - proportion from control 

set.seed(333)

mod1 <- bnec(GST_prop ~ crf(sqrt_x, model = "all"), data = datGE)


save(mod1, file = "GE.GST_bayesmodfits_all_sqrtx.RData")

fail.mods <- rhat(mod1)$failed #none

check_chains(mod1)

summary(mod1) #necsigm, followed by nechorme4pwr and nechorme4

#EC10_1 <- ecx(mod1, ecx_val = 10)
#EC50_1 <- ecx(mod1, ecx_val = 50)


plot(mod1, add_nec = TRUE, lxform=function(x){x^2}, ylim = c(0.6, 1.1), main = "GST - isogroup7641_c0_g1",
     ylab = expression("Proportion of gene expression"), xlab = expression("Concentration (% WAF)"))


nsecGST <- nsec(mod1)

###---- pirin_38837 ----####

mod2 <- bnec(pirin_prop ~ crf(sqrt_x, model = "all"), data = datGE)

fail.mods <- rhat(mod2)$failed #nec4param

mod2.good <- amend(mod2, drop = fail.mods)

save(mod2.good, file="pirin_bayesmodfits_all.good_sqrtPercentWAF.RData")

check_chains(mod2.good)

summary(mod2.good) #nechorme4pwr highest, followed by nechorme4, ecxhormebc5, ecxwb2 


par(mfrow = c(1,1))
plot(mod2.good, add_nec = TRUE, lxform=function(x){x^2}, ylim = c(0.4, 1.3), main = "pirin_isogroup38837",
     ylab = expression("Prop. expression rel. to ctrl (inv.)"), xlab = expression("Concentration (% WAF)"))


###---- epid_7393 ----####

mod3 <- bnec(epid_prop ~ crf(sqrt_x, model = "all"), data = datGE)

fail.mods <- rhat(mod2)$failed #none

save(mod3, file="epid_bayesmodfits_all_sqrtPercentWAF.RData")

check_chains(mod3)

summary(mod3) #ecxsigm highest, followed by nec4param, nechorme4pwr

par(mfrow = c(1,1))
plot(mod3, add_nec = TRUE, lxform=function(x){x^2}, ylim = c(0.4, 1.3), main = "epid_isogroup7393",
     ylab = expression("Prop. expression rel. to ctrl (inv.)"), xlab = expression("Concentration (% WAF)"))


###---- prost_14113 ----####

mod4 <- bnec(prost_prop ~ crf(sqrt_x, model = "all"), data = datGE)

fail.mods <- rhat(mod4)$failed #none

save(mod4, file="prost_bayesmodfits_all_sqrtPercentWAF.RData")

plot(mod4, add_nec = TRUE, lxform=function(x){x^2}, ylim = c(0.4, 1.3), main = "prost_isogroup14113",
     ylab = expression("Prop. expression rel. to ctrl (inv.)"), xlab = expression("Concentration (% WAF)"))

###---- acylCoA_3029 ----####

mod5 <- bnec(acylCoA_prop ~ crf(sqrt_x, model = "all"), data = datGE)

fail.mods <- rhat(mod5)$failed #none

save(mod5, file="acylCoA_bayesmodfits_all_sqrtPercentWAF.RData")

plot(mod5, add_nec = TRUE, lxform=function(x){x^2}, ylim = c(0.8, 1.05), main = "acylCoA_3029",
     ylab = expression("Prop. expression rel. to ctrl (inv.)"), xlab = expression("Concentration (% WAF)"))

####---- Settlement ----####	 
# remove dilution treatments I didn't have in GE assay

datSettle2 <- datSettle %>% 
  filter(Tr!=5, Tr!=70)


mod6 <- bnec(No.settled | trials(Total) ~ crf((sqrt_Tr), model = "all"), data = datSettle2)

save(mod6, file="settlementNo5_70_bayesmodfits_all_sqrtPercentWAF.RData")

fail.mods <- rhat(mod6)$failed #none

check_chains(mod6)

summary(mod6) #ecxexp highest weight followed by ecxll3/necsigm


#EC10_5 <- ecx(mod6, ecx_val = 10)

plot(mod6, add_nec = TRUE, lxform=function(x){x^2}, ylim = c(0, 1), main = "No 5 or 70%",
     ylab = expression("Proportion of settlement"), xlab = expression("Concentration (% WAF)"))


autoplot(mod6) +
  labs(x="Percent WAF", y="Proportion of settlement") +
  scale_x_continuous(breaks = sqrt(c(0, 6.6, 25, 56.3, 100)),
                     labels = c("0", "6.6", "25", "56.3", "100"))

    scale_x_continuous(trans = "log", labels = scales::number_format(accuracy = 0.1)) + #0.1 = 1 decimal place 0.01 = 2, 1=none

###---- Load Mods ----####

load("settlementNo5_70_bayesmodfits_all_sqrtPercentWAF.RData")
load("GE.GST_bayesmodfits_all_sqrtx.RData")
load("pirin_bayesmodfits_all.good_sqrtPercentWAF.RData")
load("epid_bayesmodfits_all_sqrtPercentWAF.RData")
load("acylCoA_bayesmodfits_all_sqrtPercentWAF.RData")
load("prost_bayesmodfits_all_sqrtPercentWAF.RData")


###---- Posteriors ----####

GE_Settle_List <- list(GST = mod1, pirin = mod2.good, epid = mod3, prost = mod4, acylCoA = mod5, settlement =mod6)

post_comp <- compare_posterior(GE_Settle_List, comparison = "nec") 

all_plots <- lapply(GE_Settle_List, function(x) {
  autoplot(x) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_blank()) +
    ggtitle("")
})
figure <- ggpubr::ggarrange(
  plotlist = all_plots, nrow = 3, ncol = 2, labels = names(all_plots),
  font.label = list(color = "black", size = 12, face = "plain"), align = "hv"
)
ggpubr::annotate_figure(
  figure, left = text_grob(expression("Prop. expression rel. to ctrl (inv.)"), rot = 90),
  bottom = text_grob(expression("Concentration (% WAF)"))
)

ggplot(data = post_comp$posterior_data, mapping = aes(x = value)) + 
  geom_density(mapping = aes(group = model, colour = model, fill = model),
               alpha = 0.3) +
  theme_classic()

post_comp$prob_diff

#comparison      prob
#           GST-pirin 0.7321830
#            GST-epid 0.8052013
#           GST-prost 0.7494374
#         GST-acylCoA 0.8543386
#      GST-settlement 0.9009752
#          pirin-epid 0.6164041
#         pirin-prost 0.5277569
#       pirin-acylCoA 0.7475619
#    pirin-settlement 0.8172043
#         epid-prost 0.4026007
#       epid-acylCoA 0.6644161
#    epid-settlement 0.7598150
#      prost-acylCoA 0.7475619
#   prost-settlement 0.8099525
# acylCoA-settlement 0.6291573

#nec_prost <- pull_out(mod10, model = "nec")
#ecx_prost <- pull_out(mod10, model= "ecx")

#post_comp_prost <- compare_posterior(list("all" = mod3.2, "ecx" = ecx_prost, "nec" = nec_prost), comparison = "nsec")

#ggplot(data = post_comp_prost$posterior_data, mapping = aes(x = value)) + 
#  geom_density(mapping = aes(group = model, colour = model, fill = model),
#               alpha = 0.3) +
#  theme_classic()

#ggplot(data = post_comp_prost$diff_data, mapping = aes(x = diff)) + 
#  geom_density(mapping = aes(group = comparison, colour = comparison, fill = comparison),
#               alpha = 0.3) +
#  theme_classic()