##==== Effect of sainfoin (Onobrychis viciifolia) ====
##==== on cyathostomin fitness, community structure, and management ====
##==== with ivermectin in horses (data analysis) ===

# == Cyathostomin fitness before and after IVM treatment and 
# == IVM pharmacokinetic analysis (Fig. 2, 3 and 5) ==

## == Packages ==
require(dplyr)
require(ggplot2)
require(RColorBrewer)
theme_set(theme_bw())
require(nlme)
require(geepack)
require(rstatix)
require(ggpubr)
require(PKNCA)
require(reshape2)

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)}


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL){
  require(grid)
# Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
# If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
# Make the panel
# ncol: Number of columns of plots
# nrow: Number of rows needed, calculated from # of cols
  layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
  print(plots[[1]])
  } else {
# Set up the page
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
# Make each plot, in the correct location
  for (i in 1:numPlots) {
# Get the i,j matrix positions of the regions that contain this subplot
    matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
    print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))}}}


## === Work file ===
setwd("D:/R - (pc)")

###============================================================================
###=== Satisfaction of energy (UFC) and protein (MADC) requirements analysis ===
###============================================================================
nrj=read.csv(file='Besoin couvert.CSV', header=TRUE, sep=';', dec=',')
ufc<-nrj %>% select(-MADC_cover) #Keep only the data concerning UFC cover percentage
madc<-nrj %>% select(-UFC_cover) #Keep only the data concerning MADC cover percentage
### === Graphs ===
Graph_ufc<-ggplot(ufc, aes(x = Date,y= UFC_cover, colour = Group))+
  geom_boxplot(alpha=1, outlier.shape = NA, size=0.5, width=0.5)+
  geom_point(aes(x = Date,y= UFC_cover, colour = Group), size = 1, shape = 1,position = position_jitterdodge(0))+
  labs(y='UFC requirement covered (%)', x="")+
  scale_colour_manual(values = c("#737373","#005a32"))+
  theme(legend.title = element_blank(),legend.position="top",
        legend.text = element_text(size=15),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))+
  theme(axis.title.y = element_text(size=10.5),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=10, face = "bold"))

Graph_madc<-ggplot(madc, aes(x = Date,y= MADC_cover, colour = Group))+
  geom_boxplot(alpha=1, outlier.shape = NA, size=0.5, width=0.5)+
  geom_point(aes(x = Date,y= MADC_cover, colour = Group), size = 1, shape = 1,position = position_jitterdodge(0))+
  labs(y='MADC requirement covered (%)', x="Date")+
  scale_colour_manual(values = c("#737373","#005a32"))+
  theme(legend.title = element_blank(),legend.position="top",
        legend.text = element_text(size=15),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))+
  theme(axis.title.y = element_text(size=10),
        axis.text.y = element_text(size=10.5, face = "bold"),
        axis.title.x = element_text(size=14),
        axis.text.x = element_text(size=10, face = "bold",angle = 45, hjust = 1))

### === Analysis of UFC requirements between the control and sainfoin group ==
mod <- lme(UFC_cover~ Group, random = ~1|Name, data = ufc) 
summary(mod)

# Linear mixed-effects model fit by REML
# Data: ufc 
#     AIC      BIC   logLik
# 1487.37 1503.699 -739.685
# 
# Random effects:
# Formula: ~1 | Name
#           (Intercept) Residual
# StdDev:   0.6866302  1.24045
# 
# Fixed effects: UFC_cover ~ Group 
#                    Value Std.Error  DF  t-value p-value
# (Intercept)   103.37182 0.2326806 420 444.2649       0
# GroupSainfoin   5.42182 0.3290601  18  16.4767       0
# Correlation: 
#              (Intr)
# GroupSainfoin -0.707
# 
# Standardized Within-Group Residuals:
#   Min           Q1          Med           Q3          Max 
# -10.96092519  -0.05039286   0.02544692   0.16273606   2.22318722 
# 
# Number of Observations: 440
# Number of Groups: 20 

### === Analysis of MADC requirements between the control and sainfoin group ==
mod <- lme(MADC_cover~ Group,random = ~1|Name, data = madc) 
summary(mod)

# Linear mixed-effects model fit by REML
# Data: madc 
#      AIC      BIC    logLik
# 2498.377 2514.706 -1245.189
# 
# Random effects:
#   Formula: ~1 | Name
#          (Intercept) Residual
# StdDev:     2.80067 3.896049
# 
# Fixed effects: MADC_cover ~ Group 
#                      Value Std.Error  DF   t-value p-value
#  (Intercept)   252.13455 0.9237812 420 272.93751  0.0000
# GroupSainfoin  -2.06636 1.3064239  18  -1.58169  0.1311
# Correlation: 
#   (Intr)
# GroupSainfoin -0.707
# 
# Standardized Within-Group Residuals:
#   Min          Q1         Med          Q3         Max 
# -9.95492990 -0.04134008  0.03137638  0.13591603  2.60811459 
# 
# Number of Observations: 440
# Number of Groups: 20 


###=======================================================================
###=== FEC percentage between two groups before IVM treatment analysis ===
###=======================================================================
xpsainfoin=read.csv(file='sainfoin_data.csv', header=TRUE, sep=';', dec=',')
xpsainfoin=xpsainfoin[xpsainfoin$Day!='d-7',]
xpsainfoin$id=factor(xpsainfoin$ID) %>% as.numeric(xpsainfoin$id)

### === Graphs ===
#Boxplot
ggplot(xpsainfoin, aes(x = Day ,y= EPG, fill=Group))+
  geom_boxplot(alpha=1, outlier.shape = NA, size=0.9, width=0.8)+
  geom_point(aes(x = Day, y = EPG, fill = Group), size = 1.5, shape = 1,position = position_jitterdodge(0))+
  labs(y='Fecal egg count (eggs per gram)', x="Days")+
  theme(axis.line = element_line(size = 1, linetype = "solid"),
        legend.title = element_blank(),legend.position="bottom",
        legend.text = element_text(size=14),
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=14, face = "bold"),
        axis.title.x = element_text(size=18), 
        axis.text.y = element_text(size=14, face = "bold"))+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))+
  theme(legend.text = element_text(size=20))+
  scale_y_continuous(breaks=c(500,1000, 1500, 2000,2500))+
  scale_fill_manual(values = c("#bdbdbd","#016450" ))+
  scale_colour_manual(values = c("#252525"))+ 
  scale_x_discrete(labels=c("0","7", "14", '21'))

### === Statistical analysis of the effects of the diet on the FEC ====
gee.fit.EPG <- geeglm(EPG ~ Group * Day,id = id, data = xpsainfoin, family = poisson,
                  corstr = "ar1", scale.fix = TRUE, std.err = "san.se")
summary(gee.fit.EPG)
# Coefficients:
#                      Estimate Std.err    Wald Pr(>|W|)    
# (Intercept)            7.0733  0.1360 2704.32   <2e-16 ***
# GroupSainfoin         -0.0258  0.2042    0.02    0.900    
# Dayd07                -0.2820  0.1758    2.57    0.109    
# Dayd14                -0.2877  0.1955    2.16    0.141    
# Dayd21                -0.4532  0.1891    5.75    0.017 *  
# GroupSainfoin:Dayd07   0.1015  0.3212    0.10    0.752    
# GroupSainfoin:Dayd14  -0.0505  0.3151    0.03    0.873    
# GroupSainfoin:Dayd21   0.0650  0.3033    0.05    0.830    
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Correlation structure = ar1 
# Scale is fixed.
# Link = identity 
# Estimated Correlation Parameters:
# Estimate Std.err
# alpha        0       0
# Number of clusters:   80  Maximum cluster size: 1

###====================================================================================
###=== Eggs development percentage between two groups before IVM treatment analysis ===
###====================================================================================
### === Graphs ===
#Boxplot
ggplot(xpsainfoin, aes(x = Day ,y= Dev, fill=Group))+
  geom_boxplot(alpha=1, outlier.shape = NA, size=0.9, width=0.8)+
  geom_point(aes(x = Day, y = Dev, fill = Group), size = 1.5, shape = 1,position = position_jitterdodge(0))+
  labs(y='Development percentage (%)', x="Days")+
  theme(axis.line = element_line(size = 1, linetype = "solid"),
        legend.title = element_blank(),legend.position="bottom",
        legend.text = element_text(size=14),
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=14, face = "bold"),
        axis.title.x = element_text(size=18), 
        axis.text.y = element_text(size=14, face = "bold"))+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))+
  theme(legend.text = element_text(size=20))+
  scale_fill_manual(values = c("#bdbdbd","#016450"))+
  scale_colour_manual(values = c("#252525"))+ 
  scale_x_discrete(labels=c("0","7", "14", '21'))

### === Statistical analysis of the effects of the diet on the Cyathotomins development ====
D1421=xpsainfoin[xpsainfoin$Day=='d14' | xpsainfoin$Day=='d21',]
geeglm(formula = Dev ~ Group * Day, family = poisson, data = D1421, 
         id = id, corstr = "ar1", scale.fix = TRUE, std.err = "san.se")
# Call:
#   geeglm(formula = Dev ~ Group * Day, family = poisson, data = D1421, 
#          id = id, corstr = "ar1", scale.fix = TRUE, std.err = "san.se")
# 
# Coefficients:
#                       Estimate  Std.err     Wald Pr(>|W|)    
# (Intercept)           2.85209  0.08555 1111.506  < 2e-16 ***
# GroupSainfoin         0.39704  0.17275    5.282   0.0215 *  
# Dayd21                0.85458  0.18551   21.221 4.09e-06 ***
# GroupSainfoin:Dayd21 -0.53438  0.26113    4.188   0.0407 *  
#   ---
gee.fit.Dev <- geeglm(Dev ~ Group * Day,id = id, data = xpsainfoin, family = poisson,
                  corstr = "ar1", scale.fix = TRUE, std.err = "san.se")
summary(gee.fit.Dev)
# Coefficients:
#                      Estimate Std.err   Wald Pr(>|W|)    
# (Intercept)            3.5052  0.1441 591.37  < 2e-16 ***
# GroupSainfoin         -0.0638  0.1706   0.14    0.708    
# Dayd07                -0.4036  0.1963   4.23    0.040 *  
# Dayd14                -0.6531  0.1676  15.18  9.8e-05 ***
# Dayd21                 0.2014  0.2188   0.85    0.357    
# GroupSainfoin:Dayd07   0.3000  0.2753   1.19    0.276    
# GroupSainfoin:Dayd14   0.4609  0.2428   3.60    0.058 .  
# GroupSainfoin:Dayd21  -0.0735  0.2597   0.08    0.777    
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Correlation structure = ar1 
# Scale is fixed.
# Link = identity 
# Estimated Correlation Parameters:
# Estimate Std.err
# alpha        0       0
# Number of clusters:   80  Maximum cluster size: 1 

###========================================
###=== FEC after IVM treatment analysis ===
###========================================
PostIvr=read.csv(file='Suivipostiver.csv', header=TRUE, dec=',', sep=';')
PostIvr$id=factor(PostIvr$ID) %>% as.numeric(PostIvr$id)
### === Graphs ===
#Boxplot
PlotFECpost<-ggplot(PostIvr, aes(x = day ,y= EPG, fill=Group))+
  geom_boxplot(alpha=1, outlier.shape = NA, size=0.9, width=0.8)+
  geom_point(aes(x = day, y = EPG, fill = Group), size = 1.5, shape = 1,position = position_jitterdodge(0))+
  labs(y='FEC after IVM treatment \n (eggs per gram)', x="Days")+
  theme(axis.line = element_line(size = 1, linetype = "solid"),
        legend.title = element_blank(),legend.position="bottom",
        legend.text = element_text(size=14),
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=14, face = "bold"),
        axis.title.x = element_text(size=18), 
        axis.text.y = element_text(size=14, face = "bold"))+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))+
  theme(legend.text = element_text(size=20))+
  scale_fill_manual(values = c("#bdbdbd","#016450"))+
  scale_colour_manual(values = c("#252525"))+ 
  scale_x_discrete(labels=c("36", "50", "63", '71', '78'))

### === Statistical analysis of the effects of the diet on the FEC after IVM treatment ====
mod.IVM <- lme(EPG~ day*Group, random = ~1|ID, data = PostIvr) 
summary(mod.IVM)
#Linear mixed-effects model fit by REML
# Data: PostIvr 
#      AIC      BIC    logLik
# 768.9599 795.9418 -372.4799
# 
# Random effects:
#   Formula: ~1 | ID
#          (Intercept) Residual
# StdDev:    18.29815 39.69751
# 
# Fixed effects: EPG ~ day * Group 
# Value Std.Error DF    t-value p-value
# (Intercept)            0.00  15.45443 56  0.0000000  1.0000
# dayd+29                0.00  19.84876 56  0.0000000  1.0000
# dayd+42               12.50  19.84876 56  0.6297624  0.5314
# dayd+50               18.75  19.84876 56  0.9446436  0.3489
# dayd+57               31.25  19.84876 56  1.5744060  0.1210
# GroupSainfoin          0.00  21.85586 14  0.0000000  1.0000
# dayd+29:GroupSainfoin  0.00  28.07038 56  0.0000000  1.0000
# dayd+42:GroupSainfoin -6.25  28.07038 56 -0.2226546  0.8246
# dayd+50:GroupSainfoin 50.00  28.07038 56  1.7812370  0.0803
# dayd+57:GroupSainfoin 56.25  28.07038 56  2.0038917  0.0499

Sainfoin=PostIvr[PostIvr$Group=='Sainfoin',]
mod.IVM.Sainfoin <- lme(EPG ~ day, random = ~1|ID, data = Sainfoin)
summary(mod.IVM.Sainfoin)
# Linear mixed-effects model fit by REML
# Data: Sainfoin 
#      AIC      BIC    logLik
# 400.2181 411.1055 -193.1091
# 
# Random effects:
#   Formula: ~1 | ID
#          (Intercept) Residual
# StdDev:    15.38204 49.95534
# 
# Fixed effects: EPG ~ day 
#              Value Std.Error DF  t-value p-value
# (Intercept)  0.00  18.48020 28 0.000000  1.0000
# dayd+29      0.00  24.97767 28 0.000000  1.0000
# dayd+42      6.25  24.97767 28 0.250224  0.8042
# dayd+50     68.75  24.97767 28 2.752459  0.0103
# dayd+57     87.50  24.97767 28 3.503129  0.0016

Control=PostIvr[PostIvr$Group=='Control',]
mod.IVM.Control <- lme(EPG ~ day,random = ~1|ID, data = Control)
summary(mod.IVM.Control)
# Linear mixed-effects model fit by REML
# Data: Control 
#      AIC      BIC    logLik
# 360.9611 371.8485 -173.4805
# 
# Random effects:
#   Formula: ~1 | ID
#          (Intercept) Residual
# StdDev:    20.80952 25.61738
# 
# Fixed effects: EPG ~ day 
#             Value Std.Error DF   t-value p-value
# (Intercept)  0.00  11.66879 28 0.0000000  1.0000
# dayd+29      0.00  12.80869 28 0.0000000  1.0000
# dayd+42     12.50  12.80869 28 0.9759001  0.3375
# dayd+50     18.75  12.80869 28 1.4638502  0.1544
# dayd+57     31.25  12.80869 28 2.4397504  0.0213

###==========================================
###=== IVM concentration and AUC analysis ===
###==========================================
PKivr=read.csv(file='PK Ivr.csv',header=TRUE, dec=',', sep=';')

### === Graphs ===
Pk<- data_summary(PKivr, varname="Dosage",groupnames=c("Assay", "Time_hour"))

PlotPK<-ggplot(Pk, aes(x=Time_hour, y=Dosage, colour=Assay))+ 
  geom_line(linetype = "solid", size = 1)+
  geom_point()+
  geom_pointrange(aes(ymin=Dosage, ymax=Dosage+sd), width=1.5)+
  labs(y='Ivermectin concentration \n (ng/mL)', x="Times post-treatment (hour)")+
  scale_colour_manual(values=c("#bdbdbd","#016450"))+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", linetype = "solid"))+
  theme(legend.title=element_blank(),legend.position="bottom",
        legend.text=element_text(size=16),
        axis.title.x = element_text(size = 18),
        axis.text.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 14,face = "bold"))+
  scale_x_continuous(breaks=c(0,24,48,72, 96))


PlotPk<- multiplot(PlotFECpost,PlotPK)

### === Statistical analysis of the effects of the diet on the IVM concentratration ====
#Data distribution 
bxp <- ggboxplot(
  PKivr, x = "Assay", y = "Dosage", 
  ylab = "Dosage", xlab = "Groups", add = "jitter"
)
bxp

#Identify outliers
PKivr %>%
  group_by(Assay) %>%
  identify_outliers(Dosage)
# Assay   Horses       Time_hour Dosage is.outlier is.extreme
# 1 Control Iegadelavega      24   10.2 TRUE       FALSE  

#Check for normality
PKivr %>%
  group_by(Assay) %>%
  shapiro_test(Dosage)
# Assay    variable statistic        p
# 1 Control  Dosage       0.932 0.000953
# 2 Sainfoin Dosage       0.941 0.00245 
ggqqplot(PKivr, x = "Dosage", facet.by = "Assay")

#Mann-Whitney U test
test=PKivr[PKivr$Time_hour!="0",]

stat.test <- PKivr %>%
  group_by(Time_hour) %>%
  wilcox_test(Dosage ~ Assay)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test
#   Time_hour .y.    group1  group2      n1    n2 statistic        p   p.adj p.adj.signif
# *     <int> <chr>  <chr>   <chr>    <int> <int>     <dbl>    <dbl>   <dbl> <chr>       
# 1         1 Dosage Control Sainfoin    10    10      50   1        1       ns          
# 2         2 Dosage Control Sainfoin    10    10      59.5 0.496    1       ns          
# 3        24 Dosage Control Sainfoin    10    10      95   0.000206 0.00124 **          
# 4        48 Dosage Control Sainfoin    10    10      91   0.00105  0.0063  **          
# 5        72 Dosage Control Sainfoin    10    10      84   0.0113   0.0678  ns          
# 6        96 Dosage Control Sainfoin    10    10      85   0.00908  0.0545  ns

### === Areas under the concentration-time curves (AUC) determined ====
PKivrO<-PKivr[order(PKivr$Horses),]
df = PKivrO %>% data.table::dcast(Time_hour ~ Horses, value.var = 'Dosage')
df$Time_hour=NULL
dfl<-as.list(df)

AUC<-list()
mytime<-c(0,1,2,24,48,72,96)

DF=data.frame(Horses=PKivr[c(1:20),1],Group=PKivr[c(1:20),2],AUC="x")
DF<-DF[order(DF$Horses),]

d=0
for(i in dfl){
  K=pk.calc.auc(i, mytime, interval=c(0, Inf), method='lin up/log down')
  DF[d+1,3]<-K
  d=d+1
}

DF$AUC=as.numeric(DF$AUC)
DF
#            Horses    Group              AUC
# 2     Icilifairy Sainfoin 203.495717032329
# 8     Idomissyou Sainfoin 248.878097781186
# 17  Iegadelavega  Control 514.456255418553
# 4           Iena Sainfoin 155.044309044932
# 14       Ilander  Control 409.386865216929
# 3        Ilivala Sainfoin 168.933858936632
# 11       Ilsouri  Control 266.897043333099
# 5         Imbala Sainfoin 363.505442132615
# 9     Imeethappy Sainfoin 239.661120882686
# 1            Ine Sainfoin  164.66016308388
# 10       Inechai Sainfoin 374.304750876822
# 18        Ipflap  Control 462.688404111873
# 15         Iplat  Control 224.125938928953
# 20   Iquissyoubi  Control 352.167128096689
# 6     Iritabella Sainfoin  269.76836508442
# 13       Iscashe  Control 544.237669600048
# 16 Ispiritusantu  Control  494.59805685902
# 12        Itarki  Control 392.226158036979
# 7         Itropi Sainfoin 216.079225731677
# 19     Ivaliente  Control 536.894860201937

### === Statistical analysis of the AUC between the two diet ====
#Data distribution 
bxp <- ggboxplot(
  DF, x = "Group", y = "AUC", 
  ylab = "AUC", xlab = "Groups", add = "jitter"
)
bxp

#Identify outliers
DF %>%
  group_by(Group) %>%
  identify_outliers(AUC)

#Check for normality
DF %>%
  group_by(Group) %>%
  shapiro_test(AUC)
#  Group    variable statistic     p
# 1 Control  AUC          0.920 0.355
# 2 Sainfoin AUC          0.889 0.166
ggqqplot(DF, x = "AUC", facet.by = "Group")

#Equality of variances
DF %>% levene_test(AUC ~ Group) #No equality

#Welch t_test (unpaired t-test)
stat.test <- DF %>% 
  t_test(AUC ~ Group) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") 
stat.test
# .y.   group1  group2      n1    n2 statistic    df        p    p.adj p.adj.signif
# * <chr> <chr>   <chr>    <int> <int>     <dbl> <dbl>    <dbl>    <dbl> <chr>       
# 1 AUC   Control Sainfoin    10    10      4.17  16.0 0.000727 0.000727 ***   
#=============================================
      




