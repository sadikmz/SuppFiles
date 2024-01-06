# <h1 style="text-align: center;"> Evaluation of 20 enset (*Ensete ventricosum*) landraces for response to *Xanthomonas vasicola* pv. *musacearum* infection</h1>


# 1. General information 

This supplementary information mainly describes data analysis part used in this manuscript. The script used has been modified from [Schandry (2017)](https://www.frontiersin.org/articles/10.3389/fpls.2017.00623/full).

A common garden experiment was conducted at a single site on 20 enset (*Ensete ventricosum*) landraces of 14 months old seedlings to assess their resistance / tolerance against Xanthomonas wilt disease caused by *Xanthomonas vasicola pv. musacearum* (Xvm). Disease infection data was recorded as count of symptomatic and total leaves per plant for each landrace.  Data recorded at 15 days interval starting from 21st days post inoculation (dpi) to 155th day for 10 time points. Percentage of symptomatic leaves calculated and converted into disease index of 0-4 range numeric score multiplying each record by 4 and named as ‘disease index’ or DI.

# 2. Installing packages required for data managment and analysis

If the required list of packages any are not available in the machine, then the  missing package(s) will be installed from CRAN.

```{r echo=TRUE, warning=FALSE}
# List of package required
packages = c("tidyr", "broom", "tibble", "dplyr", "ggplot2", "lme4", "lmerTest", "multcomp", "lsmeans", "MESS", 'agricolae', 'magrittr', "rms","coxme", "survival", "stringr","grid","tinytex")

# Install packages missing packages.
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
```

#  3. Data entry and specifying contrasts varaibles

## 3.1 Importing and formatting 

The input data was named as enset_xvm.csv and need to be placed in the working directory so the the script will be able to read it.

```{r import_data, echo=TRUE,warning=FALSE}

# get path to the current working direcotry
work_dir <- getwd()

path_to_enset_xvm <- file.path(work_dir, basename("enset_xvm.csv"))
  
enset_xvm<- read.csv(path_to_enset_xvm)
head(enset_xvm) # view the first six row of imported data
```

In the dataframe variable dss, leaf and di refers to count of diseased leaves, count of total leaves, and disease index, respectively. The numbers trailing dss, leaf, and di variables indicate the day on which data was recorded. Disease index is calculated as 4 times percentage of diseased leaves (4*dss/leaf). For details see Material and Method section of the manuscript.

Inspect data structure
```{r data_structure, echo=FALSE,warning=FALSE}
str(enset_xvm)
```

The wide formatted enset_xvm was converted to long format for ease of analysis. In order to keep the original data as it was we have grouped identical variable together for all landraces and then regrouped into single dataframe to convert it to long formatted data.

```{r reformat_data, echo=TRUE,warning=FALSE}
#extract  number of leaves showing disease symptoms at each evaluation day
leaf_dss=enset_xvm[,c(1:6,9,12,15,18,21,24,27,30,33)]

#extract  total leaves at each evaluation day
leaf_num <- enset_xvm[,c(1:5,7,10,13,16,19,22,25,28,31,34)]

#extract  disease index values at each evaluation day
di_wk <- enset_xvm[,c(1:5,8,11,14,17,20,23,26,29,32,35)]

library(tidyr)
#convert the extracted data into long format
leaf_dss_long <- gather(leaf_dss,key = "dpi",
                     value = "dss_leaf", 
                     dss21:dss155,na.rm = F)# converting data of disease leaf (this is 3 to 13) into long format
leaf_num_long <- gather(leaf_num,key = "leaf_number",
                     value = "leaf_num", 
                     leaf21:leaf155,na.rm = F)#converting total leaf data (this is 3 to 13) into long format
di_wk_long <- gather(di_wk,key = "di",
                  value = "disease_index", 
                  di21:di155,na.rm = F)# converting data disease index (this is 3 to 13) into long format

#Combine long formatted subset data 
enset_xvm_long <- cbind(leaf_dss_long,leaf_num_long$leaf_num,di_wk_long$disease_index)

enset_xvm_long$subject <- c(1:nrow(enset_xvm_long))# add subject variable. This variable helps to assign a unique, numeric identifier to each individual.
enset_xvm_long$count_dslf <- c(rep(1,6000))# add variable count_dslf: a repeats of 1 from 1 to the length of data set (6000). This variable will be used to count the number of observations in landrace at each evaluation period that show disease. It might be used to calculate incubation period of landraces for the disease under study.     
colnames(enset_xvm_long) <- c("landrace", "plant", "ptno", "block", "ind", "dpi", "diseased", "leaf", "di", "subject", "count_dslf")# rename column headings
enset_xvm_long$id <- interaction(enset_xvm_long$landrace,enset_xvm_long$ptno,enset_xvm_long$block,enset_xvm_long$ind) # create id column of landrace, plant number and block.

enset_xvm_long <- na.omit(enset_xvm_long) #Generate a table that does not contain missing observations
enset_xvm_long$dpi <- na.omit(as.numeric (unlist( strsplit( as.character( enset_xvm_long$dpi ), "dss" ))))  #discard the "X" and save the number as a number (instead of a factor).

# Factorize variables landrace, plant and block 
enset_xvm_long$block <- as.factor(enset_xvm_long$block)
enset_xvm_long$plant <- as.factor(enset_xvm_long$plant)
enset_xvm_long$landrace <- as.factor(enset_xvm_long$landrace)

# # view structure of long formatted data
# str(enset_xvm_long)
```

Checking overall summary for missing values, or out of range misrecorded values. Missing values coded as NA and will be observed print in summary function if there are any. Out of range values can be checked by looking at the min and max values of the summary table and make sure if they're within the expected range.
```{r data summary,echo=TRUE,warning=FALSE}
summary(enset_xvm_long)
```

## 3.2 Defining the contrasts
Assigning and setting proper contrast is in analysis of variance. See for detail [here](https://cran.r-project.org/web/packages/codingMatrices/vignettes/codingMatrices.pdf). "Treatment" contrasts specify that the first alphabetical level will be used as a reference for all others (see landrace below), while a "sum" constrast means that the reference value is the mean across all levels of that variable (the grand mean).
```{r contrasts,echo=TRUE,warning=FALSE}
####Specify what should be "appropriate" contrasts
contrasts(enset_xvm_long$landrace) <- "contr.treatment" #First alphabetical landrace will be the baseline!
contrasts(enset_xvm_long$block) <- "contr.sum" # blockes will be averaged to generate the baseline!
```

## 3.3 Add a censoring variable
A new variable called "Useful" is added the data and is a binary yes/no variable. The purpose of this variable is marking observations, that are re-observations of a subject that has previously reached disease index 4. As disease index of 4 means 100% wilting, the plant either dead or considered as dead when it reached this degree of wilting. Since death is permanent, continuing to observe this plant is unlikely to provide new information. To more accurately capture the linear phase of disease development, all observations where no symptoms are visible, except the one on the day before disease onset, will also be assigned to Useful = No.

```{r censor_dead_plants, warning=FALSE, message=FALSE, error=FALSE, echo=T}
###Interesting R-related observation, the below does not work when subsetting is done with filter(),
###because filter does not retain rownames!
for (i in 1:max(enset_xvm_long$subject)) { ###Go by subject
  dummy1 <- enset_xvm_long[enset_xvm_long$subject==i,] ##Create a first dummy object, that is a subset of the full data containing the current subject
  if(min(dummy1$di) == 0){
  #remove those observations that are before disease onset, except the one directly before disease onset.
  dummy1 <- dummy1[dummy1$dpi %in% (max(dummy1$dpi[dummy1$di==0]):max(dummy1$dpi)),]
  }
  if (max(dummy1$di) == 4) { ###If this subject dies at some point
    dummy2 <- dummy1[dummy1$di==4,] ###Make a new dummy object, that only contains those recordings where di=4
    NEW <- enset_xvm_long[enset_xvm_long$subject==i & (enset_xvm_long$dpi %in% min(dummy1$dpi):min(dummy2$dpi)),] ###Generate data subset "NEW", which contains those observations for a subject, that are between (including) the last recording where di=0 and the first recording where disease index is 4.
  } else { ###If this dubject does not die
    NEW <- dummy1 ###New is the same as dummy1
  }
  enset_xvm_long$Useful[rownames(enset_xvm_long) %in% rownames(NEW)] <- c("Yes") ###All of those row(names) that are part of the "NEW" object are useful. Therefore these receive status "Yes" in column "Useful"
}
enset_xvm_long$Useful[which(is.na(enset_xvm_long$Useful))] <- c("No") ###Those that are not yes (and therefore are NA) become No
rm(dummy1,dummy2)
```

Inspecting the final data before analysis
```{r inspect_table, echo=TRUE,warning=FALSE,message=FALSE}
library(broom)
library(tibble)
as_tibble(enset_xvm_long)
# tidy(enset_xvm_long) ###This can be used to assess the descriptive statistics of the data.
```

# 4. Analysis of the area under the disease progression curve (AUDPS)

Here the data is summarized to visualize area under the curve using disease index (y) and time (x).

Summarizing disease index recordings and generating disease area plot
```{r audps_summary, warning=FALSE, message=FALSE, error=FALSE,echo=TRUE}
library(dplyr)
di_summary <- enset_xvm_long %>% 
  group_by(plant, landrace, block, dpi) %>% 
  summarise(mean=mean(di),sd=sd(di),se=sd(di)/sqrt(length(di))) ##calculate within block mean, sd, and se for each plant/landrace combination.
##Averages within each replicate. This summary table is mainly helpful for plotting, not used for analysis...
colnames(di_summary) <- c("plant", "landrace", "block", "dpi", "mean", "sd", "se") ### rename column names

library(ggplot2)
ggplot(di_summary,aes(x=dpi,y=mean)) + ###Use data di_summary
  geom_area(aes(color=block,fill=block))+
  # geom_line() +#show line by each block
  #aes(x=dpi,y=mean,color=landrace) + ###Color by landrace, specify x and y.
  # geom_area(aes(fill=block),position="identity",alpha=0.15) + ###Area plot, colored by landrace
  facet_wrap(~landrace) + ###One plot per landrace
  labs(x = "Days post infection", y = "Average disease index") + #Labels
  ggtitle("Disease areas per landrace, for blocks represented by three filled color")+ #Title
  theme_bw()+
  theme(strip.text.x = element_text(size = 14,face="bold"),
        axis.text.x = element_text(size = 10,face = "bold"),
        axis.title.x = element_text(size = 12,face = "bold"),
        axis.title.y = element_text(size = 12,face = "bold"),
        axis.text.y =element_text(size = 12,face ="bold" ),
        legend.position = "none")
```  

From this plot, one can see that in the example dataset the areas differ quite drastically among different landraces. To calculated the actual audps for each individual in the dataset a new data frame is created. AUDCP is calculated using audps function in agricolae package.

```{r calculate_AUDPS, message=FALSE, echo=TRUE,warning=FALSE}
library("MESS")
library('agricolae') # Used to calculated audps
library('magrittr')  # Adds %>% operator for chaining commands

# Calculate AUDPS
audps_data <- enset_xvm_long %>%
  group_by(landrace, ptno, block) %>% # do any following calculations on every combination of treatments
  summarise(audps = round(audps(di, dpi),0))
```

## 4.1 Model
Linear model and linear mixed effects model were fitted to identify and select the better fit model to area under disease progress data to be used for final analysis.

```{r AUDPS_analysis,message=FALSE,warning=FALSE,echo=TRUE}
library("lme4")
library("lmerTest")
auc_lm <- lm(audps~landrace+block,data=audps_data) #Linear model to identify landrace and block effect.
auc_lm1<-lm(audps~landrace,data=audps_data) ### linear model to identify landrace effects.
auc_lmer <- lmer(audps ~ landrace + (1|block),data = audps_data  ) #Linear mixed effects model. AUC modeled as a function of landrace and random effects of block.  
AIC(auc_lm,auc_lm1,auc_lmer) #Lower AIC, better fit, linear mixed model is slightly better. #Hence, mixed model with block as random effect will be used for further analysis.
```

Summary of mixed linear model fit to AUDPS
```{r summary auc_lmer,echo=TRUE,warning=FALSE,message=FALSE}
summary(auc_lmer)  ### A model summary, containing information on the model.
```

Model summaries contain information on the estimated mean (Estimate) and standard errors, together with a t - and corresponding p-value in the output of the summary function. Note, that the above only contains information on differences between different levels and the reference, which is called (Intercept) in this case it's landrace Abatemerza. The reference is contrast dependent.


```{r test_audps_random_effect,echo=TRUE,warning=FALSE,message=FALSE}
ranova(auc_lmer) 
#The result displays block has a significant random effect explaining the importance of including block variable in the model.
```

## 4.2 Post-hoc of audps

As observed in the model summary R by default does contrast with the alphabetically first landrace. This can be analyzed using a generalized linear hypothesis test (glht), while adjusting for multiple comparisons using Tukey's method. Below is shown glht analysis using Tukey's method for area under disease progress stairs (AUDPS).

```{r audps post-hocs,message=FALSE, error=FALSE,warning=FALSE,echo=TRUE,warning=FALSE,message=FALSE}
library("multcomp")
library("lsmeans")
##carry out mean separation and adjust for compact letter desply (cld)###
cld_audps_lmer <- cld(glht(auc_lmer, linfct = mcp(landrace = "Tukey")), test = adjusted("holm"))###Save letters
cld_audps_cbind <- cbind(levels(enset_xvm_long$landrace),cld_audps_lmer$mcletters$Letters) ###bind letters to columns of respective landraces
cld_audps_df <- data.frame(cld_audps_cbind)#format extracted cld with landraces from mean separation to data.frame
colnames(cld_audps_df) <- c("landrace","aud_let")###Name columns
rownames(cld_audps_df) <- NULL #remove repeated rownames
cld_audps_df$landrace <- as.factor(cld_audps_df$landrace)

#Extrating estimated mean and confidence interval of landraces using lsmeanLT
audps_lsmeansLT <- lsmeansLT(auc_lmer,test.effs="landrace")#lsmeansLT calculates Least Squares Means and Confidence Intervals for the factors of a fixed part of mixed effects model of lmer object
#audps_lsmeansLT=audps_lsmeansLT$lsmeans.table# extract table containing means with confidence interval
rownames(audps_lsmeansLT) <- NULL #remove row names containing duplicates of landrace lists
audps_lsmeansLT<-as.data.frame(audps_lsmeansLT)[,-1]

library (dplyr)

cld_audps <- cbind(audps_lsmeansLT,cld_audps_df[,2]) #combine extracted mean data from lmer using lsmeansLT and dataframe of letters from cld
rownames(cld_audps) <- NULL #remove repeated rownames
colnames(cld_audps)=c("landrace", "audps", "Stand_Error","df","t-value", "lower", "upper", "p_value", "audp_let")# rename column of combined orinal data. Here landrace is duplicated and can be remove but left as landrace
cld_audps <- as.data.frame(cld_audps) #Coerce data containing means with confidence interval and letter display to dataframe
cld_audps <- cld_audps[order(cld_audps$audps),]# reorder by mean value to cheat ggplot
#cld_audps_lmerer$landrace<-factor(cld_audps_lmerer$landrace,levels=cld_audps_lmerer$landrace)# put factors in level

audps_data__cld <- left_join(audps_data,cld_audps_df, by ="landrace",copy=T) ###Add letter information to raw data - audps
head(audps_data__cld)
```

This information can, for example, be integrated into a boxplot of the individual disease areas. For example, using a compact letter display, in combination with audps, combines statistical and visual information. First, however, the grouping letters need to be calculated. Then, these letters are added to the boxplot generated earlier.

## 4.3 Ploting audps

```{r ploting_audps,echo=TRUE,warning=FALSE,message=FALSE}
####Generate plot of meanCI of AUDCP with significance letters and raw data as jittered points
cld_audps <- cld_audps[order(cld_audps$audps),]# reorder by mean value to cheat ggplot
cld_audps$landrace <- factor(cld_audps$landrace,levels = cld_audps$landrace)#put factors in level

ggplot(aes(x=landrace, y=audps, color=landrace),data=audps_data__cld) + ###Plot the audps_data
  scale_shape_manual(values=c(9, 19, 22,26))+
  geom_crossbar(data = cld_audps, aes(x = landrace, y = audps, ymin = lower, ymax = upper,fill=landrace), alpha=0.3) +
  geom_jitter(aes(shape=block),data = audps_data__cld) + ###with jitter overplotted, symbol shape defined by block
  geom_text(aes(x=landrace, y=0.039, label=aud_let),color="black", data=audps_data__cld) + ###Get the letters from auc_cld
  #and write those to position y=-3
  labs(y="audps") + #Y-Axis label
  ggtitle("Mean and raw audps values per landrace and grouping letters") +#Title
  theme_bw()+
  theme(axis.title.x = element_text(face="bold", size=12), axis.text.x  =  element_text(angle=90, vjust=0.5, size=12))+
  theme(legend.position="none")
```

# 5. Analysis of disease development
Linear mixed model [(Bates, 2010)](https://www.jstatsoft.org/article/view/v067i01/) used to handle the repeated measure data collected in this experiment for testing landrace specific influence on the area under the disease progression curve. Here, individual data of disease index which formatted earlier will be used. Diease index (di) is modeled on the fixed effects "landrace" and an interaction of landrace and time. Block is included as a random effect, meaning it is not of direct interest, but assumed to capture introduced variation, specifically by affecting the intercept.


```{r di_set_constrast,echo=TRUE,warning=FALSE,message=FALSE}
library("lmerTest")
###Define contrasts for lmer
contrasts(enset_xvm_long$landrace) <- "contr.treatment"###First alphabetical landrace will be the baseline!

contrasts(enset_xvm_long$block) <- "contr.sum" ###blockes will be averaged
###Drop things that are not "Useful"
enset_xvm_long_useful <- filter(enset_xvm_long, Useful=="Yes")
# str(enset_xvm_long_useful)
```

## 5.1 Build linear mixed effect model(s) 
```{r model_di,message=FALSE,error=FALSE,warning=FALSE, echo=TRUE}
## only landrace as fixed effect (a simple model),
## block and subject are assumed to have random effects.
## landrace and time are assumed to interact, meaning that for each day, the landrace specific slope is estimated.
## See page 6 (Table2) of the lme4 vignette on how to construct error terms
## This model is supposed to capture disease development, so it will work with the data previously deemed "useful"
disease_lmer <- lmer(di~landrace+(1|dpi/landrace)+(1|block/plant),data=enset_xvm_long_useful,REML = F)# Here the final optimum model is that best fits the data is used and model simplication steps not shown here.
```

View summary of ananlysis of variance
```{r anova_di,echo=TRUE,message=FALSE,warning=FALSE}
#analysis of variance
anova(disease_lmer)# Analysis of variance shows highly significant difference among landrace for disease index
```
The model can be investigated using summary functions.
```{r summary_di_lmer,echo=TRUE,warning=FALSE,message=FALSE}
###Check model summary
summary(disease_lmer)
```
## 5.2 Testing random effects 
It's important to check statistically importance random effects in the model so that decision to include or remove them will be made before proceeding to any analysis.   

```{r test_random_effect_di, echo=TRUE,warning=FALSE,message=FALSE}
rand(disease_lmer)# the random effects have shown highly significance suggesting that it wise to include these random effects in the model
```

## 5.3 Post hoc analysis of disease index
Here function lsmeanLT from lmer package to extract means with their confidence interval. Then multiple comparison of landrace for disease index can be carried using general linear hypotheses (glht)and multiple comparisons for parametric models from multicomp package. The result adjusted to letter display for difference among landrace using compact letter display (cld) function. ggplot2 used to plot the result combining with original disease index data.

```{r di_post_hoc,error=FALSE,warning=FALSE,message=FALSE, echo=TRUE}
di_lsmeansLT=lsmeansLT(disease_lmer,test.effs="landrace")#lsmeansLT calculates Least Squares Means and Confidence Intervals for the factors of a fixed part of mixed effects model of lmer object

di_lsmeansLT<-as.data.frame(di_lsmeansLT)[,-1]
rownames(di_lsmeansLT) <- NULL #remove row names containing duplicates of landrace lists

library (dplyr)
di_lsmeansLT <- di_lsmeansLT %>% mutate_if(is.numeric, funs(round(., 2)))#round data.frame to two digist using dplyr

##carry out mean separation and adjust for compact letter desply (cld)###
cld_di_lmer <- cld(glht(disease_lmer, linfct = mcp(landrace = "Tukey")), test = adjusted("holm"))###Save letters
cld_di_lmer <- cbind(levels(enset_xvm_long$landrace),cld_di_lmer$mcletters$Letters) ###bind letters to columns of respective landraces
cld_di_lmer <- data.frame(cld_di_lmer)#formatextracted cld with landraces from mean separation to data.frame
colnames(cld_di_lmer) <- c("landrace","di_let")###Name columns
rownames(cld_di_lmer) <- NULL #remove repeated rownames
di_lmer_cld <- cbind(di_lsmeansLT,cld_di_lmer[,2]) #combine extracted mean data from lmer using lsmeansLT and dataframe of letters from cld
rownames(di_lmer_cld) <- NULL #remove repeated rownames
colnames(di_lmer_cld) <- c("landrace", "mean", "Stand_Error", "DF", "t-value", "lower", "upper", "p-value", "di_let")# rename column of combined orinal data.
di_lmer_cld <- as.data.frame(di_lmer_cld) #Coerce data containing means with confidence interval and letter display to dataframe
di_lmer_cld <- di_lmer_cld[order(di_lmer_cld$mean),]# reorder by mean value to cheat ggplot
di_lmer_cld$landrace<-factor(di_lmer_cld$landrace,levels=di_lmer_cld$landrace)# put factors in level
di_long_cld <- left_join(enset_xvm_long,cld_di_lmer,by="landrace",copy=T) ###Add letter information to raw data - disease index
```

## 5.4 Ploting disease index lmer
```{r di_post_hoc_plot,echo=TRUE,warning=FALSE,message=FALSE}
di_long_cld_copy <- di_long_cld
di_long_cld_copy <- di_long_cld_copy[order(di_long_cld_copy$di),]# reorder by mean value to cheat ggplot
di_lmer_cld=di_lmer_cld[order(di_lmer_cld$mean),]# reorder by mean value to cheat ggplot
di_lmer_cld$landrace<-factor(di_lmer_cld$landrace, levels = di_lmer_cld$landrace)# put factors in level
di_long_cld$block <- as.factor(di_long_cld$block)
####Generate plot of mean+confidence interval of disease index with significance letters and raw data as jittered points#
ggplot(aes(x=landrace, y=di, color=landrace),data=di_long_cld) + ###Plot the origial disease index data
  geom_crossbar(data = di_lmer_cld, aes(x = landrace, y = mean, ymin = lower, ymax = upper,fill=landrace), alpha=0.3) +
    geom_jitter(aes(shape=block)) + ###with jitter overplotted, symbol shape defined by block
  scale_shape_manual(values=c(9, 19, 22,23))+
  geom_text(aes(x=landrace, y=-0.25, label=di_let,size=10),color="black", data=di_long_cld) + ###Get the letters from di_long_cld
  #and write those to position y=-3
  labs(y="disease index") + #Y-Axis label
  ggtitle("disease index raw values and mean from the LMM\ nwith 95% CI per landrace and grouping letters") +#Title
 theme_bw()+
  theme(axis.title.x = element_text(face="bold", size=14), axis.text.x  =  element_text(angle=90, vjust=0.5, size=12,face = "bold"))+
  theme(axis.title.y= element_text(size = 12,face = "bold"))+
  theme(legend.position="none")
```

# 6. Analysis of apparent infection rate (AIR)

Apparent infection rate which is the speed at which an epidemic develops [(Meena et al., 2011)](https://www.tandfonline.com/doi/abs/10.1080/03235400903345281), was calculated as the slope of disease index development.

Disease index data though best fits to mixed model we could not able to extract rate of disease development value for each landrace using mixed model approach. Instead we considered linear development of disease and used disease index as function of time and intercept for analysis (i.e y=ax+b, where y is disease index, x is time and b is intercept). In this way a function that goes through each landrace is constructed to extract the slop of disease index from each landrace.

```{r generate_air_data,echo=TRUE,warning=FALSE,message=FALSE}
slop_data=di_wk# take copy of the original wide formated data for analysis 

slop_api<-numeric(600)
time<-c(21,35,50,65,80,95,110,125,140,155)# disease evaluation days after inoculation of the pathogen.
for(i in 1:600) # go row by row from the first to the last row (60th)
{
  y<-as.numeric(slop_data[i,6:15])#drop column 1-3 from the data disease index (di)
  l<-lm(y~time)  # construct linear model for each row data (each landrace in blocks so that the coeffiecent of linear model will extract for slop analysis)
  slop_api[i]<- l$coefficients[2]# extract slop of the model 
}
di_wk$slop_api=round(slop_api,3)# add extract slop to original disease index data
slop=di_wk[,c(1:6,16)]# take columns containing landrace, plant and slop 
# head(slop)#view the first six rows
```

```{r air_contrast,echo=TRUE,warning=FALSE,message=FALSE}
slop$block=as.factor(slop$block)#factorize block
slop$landrace=as.factor(slop$landrace)#factorize landrace
enset_xvm$id=interaction(enset_xvm$landrace,enset_xvm$ptno,enset_xvm$block)# create combination id of landrace plant number and block which it belongs
contrasts(slop$landrace) <- "contr.treatment"###First alphabetical landrace will be the baseline!
contrasts(slop$block) <- "contr.sum" ###blockes will be averaged to generate the baseline!

```

## 6.1 Model comparison for apparent infection rate
```{r air_test_models, echo=TRUE,warning=FALSE,message=FALSE}
#analyse data use
slop_lm=lm(slop_api~ block+landrace, data=slop)
slop_lmer=lmer(slop_api~ landrace +(1|block), data=slop)
AIC(slop_lm,slop_lmer)# The lower AIC the better the model is and when comparing negative AIC model with, model with high negative value is better as indicated before in anakysis of area under disease progress. So linear model is the choice for analysis.
summary(slop_lm)# view model summary
```

## 6.2 Post hoc analysis of slop
```{r air_post_hoc, error=FALSE,warning=FALSE,message=FALSE, echo=TRUE}
library(lsmeans)

slop_lsmeans=lsmeans(slop_lm,"landrace") #calculate estimated mean, standard error, DF, and confidence intervals
slop_lsmeans=summary(slop_lsmeans)
slop_lsmeans_df=data.frame(slop_lsmeans)# converting it data frame
str(slop_lsmeans_df)# view data structure
colnames(slop_lsmeans_df)=c("landrace", "Estimate", "Stand_Error", "DF", "lower","upper")#rename column headings
slop_lsmeans_df=slop_lsmeans_df %>% mutate_if(is.numeric, funs(round(., 3)))#round data.frame to two digits using dplyr

###Extract estimated mean and confidence interval of landraces
di_lsmeansLT=di_lsmeansLT %>% mutate_if(is.numeric, funs(round(., 2)))#round data.frame to two digist using dplyr
slop_lsmeans=lsmeans(slop_lm,"landrace") #calculate estimated mean, standard error, DF, and confidence intervals
slop_lsmeans=summary(slop_lsmeans)
slop_lsmeans_df=data.frame(slop_lsmeans)# converting it data frame
str(slop_lsmeans_df)# view data structure
colnames(slop_lsmeans_df)=c("landrace", "Estimate", "Stand_Error", "DF", "lower","upper")#rename column headings
slop_lsmeans_df=slop_lsmeans_df %>% mutate_if(is.numeric, funs(round(., 3)))#round data.frame to two digist using dplyr

##carry out mean separation and adjust for compact letter desply (cld)###
cld_slop_lm=cld(glht(slop_lm, linfct = mcp(landrace = "Tukey")), test = adjusted("holm"))###Save letters
cld_slop_lm <- cbind(levels(enset_xvm_long_useful$landrace),cld_slop_lm$mcletters$Letters) ###bind letters to columns of respective landraces
cld_slop_lm=data.frame(cld_slop_lm)#formatextracted cld with landraces from mean separation to data.frame
colnames(cld_slop_lm)=c("landrace","Letter")###Name columns
rownames(cld_slop_lm) <- NULL #remove repeated rownames
slop_lm_cld <- cbind(slop_lsmeans_df,cld_slop_lm$Letter) #combine extracted mean data from lmer using lsmeansLT and dataframe of letters from cld
rownames(slop_lm_cld)=NULL
colnames(slop_lm_cld)=c("landrace", "mean", "std_err", "DF", "lower", "upper", "Letters")# rename column of combined orinal data. slop_lmer_cld <- as.data.frame(slop_lmer_cld) #Coerce data containing means with confidence interval and letter display to dataframe
slop_lm_cld=slop_lm_cld[order(slop_lm_cld$mean),]# reorder by mean value to cheat ggplot
slop_lm_cld$landrace<-factor(slop_lm_cld$landrace,levels=slop_lm_cld$landrace)# put factors in level
slop_cld <- left_join(slop,cld_slop_lm,by="landrace",copy=T) ###Add letter information to raw data - disease index
slop_cld=slop_cld %>% mutate_if(is.numeric, funs(round(., 3)))#round data.frame to two digist using dplyr
# head(slop_lm_cld)
```

## 6.3 Plotting anlaysis for apparent infection rate

```{r air_plotting, echo=TRUE,warning=FALSE, message=FALSE}
#Generate plot of mean+confidence interval of slop with significance letters and raw data as jittered points##
ggplot(aes(x=landrace, y=slop_api, color=landrace),data=slop_cld) + ###Plot the origial disease index data
  geom_crossbar(data = slop_lm_cld, aes(x = landrace, y = mean, ymin = lower, ymax = upper,fill=landrace), alpha=0.3) +
  geom_jitter(aes(shape=block)) + ###with jitter overplotted, symbol shape defined by block
  scale_shape_manual(values=c(9, 19, 22,23))+
  geom_text(aes(x=landrace, y=0.04, label=Letters,size=10),color="black", data=slop_lm_cld)  +###Get the letters from cld_di_lmer
  #and write those to position y=-3
  labs(y="disease index") + #Y-Axis label
  ggtitle("Mean and raw slop of disease index per landrace and grouping letters") + #Title  
  theme_bw()+
  theme(axis.title.x = element_text(face="bold", size=10), axis.text.x  =  element_text(angle=90, vjust=0.5, size=12,face = "bold"))+
  theme(axis.title.y= element_text(size = 12,face = "bold"))+
  theme(legend.position="none")
```


# 7. Survival Analysis

Survival analysis, also known as time-to-event analysis, builds on a different dataset, which can be generated from raw disease indices. Survival analysis provides a way to analyze the time-to-event recordings within one population. In the case of the data used here, what is of interest is an event that can be referred to as *death*. Correctly defining *death* in this context is crucial for the outcome of the analysis. *Death*, here, is defined as a subject reaching a certain disease index, from which it cannot recover. In this script the disease index that defines the threshold to *death* is called "cutoff".  We used a cutoff value of 2.23 (i.e 55.75% of the leaves showing disease symptoms due to Xvm) which is the high percent infection recorded for the least susceptible enset landrace tested here 'Mazia'.

In the following plot all of the observations that are above the dotted line in the plot below are *dead*, at the day they cross that line. Those that never cross the line are alive until the end of observations, and are "right censored", meaning that their event was not observed during the time this subject was observed.
```{r plot_data_cutoff, echo=TRUE,warning=FALSE,message=FALSE}
cutoff <- c(2.23)
ggplot(data=enset_xvm_long) +  
  geom_jitter(aes(x=dpi, y=di,color=landrace,shape=block)) +
  geom_segment(aes(x=0, xend=max(dpi)+0.5,y=cutoff, yend=cutoff), linetype="dotted") +
                  labs(x = "Days post infection", 
                       y = "disease indices", 
                       title="Scatterplot of disease indices\n cutoff plotted as dotted line") +
                  coord_cartesian(xlim=c(19,66))
```
## 7.1 Generate data for survival analysis
A "survival table" can be generated using the following code. This code works on the long table, generated in the beginning, and the cutoff variable defined above.

```{r generate_survival_table, message = FALSE,error=FALSE,warning=FALSE,echo=TRUE}
###Generate survival table
surv_from_di <- data.frame(Subject=enset_xvm_long$subject,
landrace=enset_xvm_long$landrace,
plant=enset_xvm_long$plant,
block=enset_xvm_long$block)
###Fill survival table based on the enset_xvm_long table. This generates warnings. These can be ignored and come from the min()
for (i in 1:max(enset_xvm_long$subject)) { #Go by subject
  dummy <- enset_xvm_long[enset_xvm_long$subject==i,] #generate dummy for the subject
  if (is.infinite(min(dummy$dpi[which(dummy$di >= cutoff)]))) { #If none of the di is greater than the cutoff (this is where warnings are generated, min on an empty object returns infinite and a warning!)
    surv_from_di$End[i] <- max(dummy$dpi) #Generate a  observation, censoring at the maximum dpi recorded
    surv_from_di$Death[i] <- 0 #Still alive, because it did not pass the cutoff
  } else { #If more than zero di are greater than the cutoff
    surv_from_di$End[i] <- min(dummy$dpi[which(dummy$di >= cutoff)]) #Use the lowest dpi where condition is met
    surv_from_di$Death[i] <- 1 #record as dead
  }
}
rm(dummy)
```

## 7.2 Kaplan-Meier estimates of survival

Kaplan-Meier estimates of survival are the basic tool of survival analysis. These can be estimated using the survfit function from the "survival" package.

```{r survfit,echo=TRUE,warning=FALSE,message=FALSE}
library(survival)
library(stringr)
library(tidyr)
surv_di_fit <- survfit(Surv(End, Death) ~landrace +strata(block), data=surv_from_di)
```
The survminer package provides the ggsurvplot() function. This works nicely on datasets with few treatments. However, for the data presented here, it is easier to initially generate a data frame that contains the fit and plot with "normal" ggplot2

```{r generate_KM_dataframe, message=FALSE, warning=FALSE, error=FALSE, echo=TRUE}
###Strata dummy generation, modified from kmiddleton / rexamples
strata_dummy <-NULL
for(i in 1:length(surv_di_fit$strata)){
      # add vector for one strata according to number of rows of strata
      strata_dummy <- c(strata_dummy, rep(names(surv_di_fit$strata)[i], surv_di_fit$strata[i]))
}
###Data frame generation inspired by a post by Hadley Wickham to the ggplot2 Googlegroup
surv_di_fit.df <- data.frame(
  time = surv_di_fit$time,
  n.risk = surv_di_fit$n.risk,
  n.event = surv_di_fit$n.event,
  surv = surv_di_fit$surv,
  strata = strata_dummy,
  upper = surv_di_fit$upper,
  lower = surv_di_fit$lower
)
zeros <- data.frame(time = 0, surv = 1, strata = names((surv_di_fit$strata)),
  upper = 1, lower = 1)

surv_di_fit.df <- plyr::rbind.fill(zeros, surv_di_fit.df) ###I dont want to load plyr because i guess it will interfere with dplyr...
rm(strata_dummy)
rm(zeros)

##The block and landrace columns are regenerated from the strata field, there are probably more elegant ways to do this

surv_di_fit.df$block <- as.factor( str_split_fixed(
  matrix( nrow=length(surv_di_fit.df$strata),ncol=2, unlist(strsplit(as.character(surv_di_fit.df$strata),", ")), byrow=T )[,2],"=",2)[,2])
surv_di_fit.df$landrace <- as.factor( str_split_fixed(
    matrix( nrow=length(surv_di_fit.df$strata),ncol=2, unlist(strsplit(as.character(surv_di_fit.df$strata),", ")), byrow=T )[,1],"=",2)[,2])
###End of data frame generation
```

```{r plot_KM, echo=TRUE,warning=FALSE,message=FALSE}
###Start plotting
ggplot(surv_di_fit.df, aes(time, surv, colour = landrace)) +
  geom_step(aes(y = surv*100,linetype=block)) +
  facet_wrap(~landrace) +
  ggtitle("Survival estimates for all blockes")+
  theme(legend.position="none")
```

## 7.3 Cox-Proportional hazards and hazard ratios

Different approaches to survival analysis, are based on analysing the hazards. The hazard is the probability of experiencing an event at a given timepoint. Many hazard based analysis assume that hazards are proportional between treatments, meaning that they differ by a fixed factor. Hazards were strongly influenced by Cox and a basic model is the cox proportional hazards model.

```{r coxph and PH assumption,warning=FALSE,message=FALSE}
###Cox-Proportional hazards####
#Build model
srv_coxph <- coxph(Surv(End, Death) ~landrace + strata(block), data=surv_from_di)
summary(srv_coxph)
```

```{r check proportionality of hazard,warning=FALSE,message=FALSE}
###Check porportionality of hazards
cox.zph(srv_coxph)
```
From the above output, one can see that the p-levels several landraces such as landrace Arkia, Bededet, Dirbo, Lemat, Unjame and Wachiso show statistical significance (p<0.05) including the global statistics. This suggest a clear violation of proportionality of hazard and  one cannot assume the proportional hazards for further analysis. Alternative analysis would be uisng Cox-mixed model effects similar to linear mixed-effects model, but with a different type of response variable.  

## 7.4 Hazard analysis using cox mixed model

```{r Cox-mixed model,error=FALSE,warning=FALSE, message=FALSE}
library("coxme")
cxm=coxme(Surv(End,Death)~landrace+(1|block),data = surv_from_di)
```

```{r cox-mixed-model-anova,warning=FALSE,message=FALSE}
anova(cxm)
```

```{r cox-mixed-model-summary, warning=FALSE, message=FALSE}
summary(cxm)
```
```{r cox-mixed model post-hoc analysis,error=FALSE,warning=FALSE}
cxme_ph=cld(glht(cxm,linfct=mcp(landrace="Tukey")))
cxme_ph.df=data.frame(cbind(levels(surv_from_di$landrace),cxme_ph$mcletters$Letters))
colnames(cxme_ph.df)=c("Landrace","Letters")
rownames(cxme_ph.df)=NULL
cxme_ph.df
write.csv(cxme_ph.df,file = "cxme_ph.df.csv")
```
## 7.5 Analysis of survival curves and fits
Comparing the Kaplan-Meier survival estimates can be done in different ways. Here non-parametric log-rank testing and parametric survival regression will be used.

### Logrank testing
The below produces all pairwise comparisons of the Kaplan Meier estimate of survival using a logrank test.
```{r pairwise_logranks, message=FALSE,error=FALSE,warning=FALSE, echo=TRUE}
###Make a table of pairwise chisq pvalues, for the logrank test.
#Based on a post to the R Mailing list by T. Therneau
pw_logrank_test_type <- 0 ###0 for logrank, 1 for peto and peto
pw_logrank <- matrix(0., nlevels(surv_from_di$landrace),nlevels(surv_from_di$landrace))
for (i in 1:nlevels(surv_from_di$landrace)) {
  for (j in (1:nlevels(surv_from_di$landrace))[-i]) {
    datasubset <- droplevels(subset( surv_from_di,
                                     surv_from_di$landrace %in% (unique(surv_from_di$landrace))[c(i,j)]))
    temp <- survdiff(Surv(End, Death)~landrace+strata(block), data=datasubset, rho=pw_logrank_test_type)

    pw_logrank[i,j] <- pchisq(temp$chisq, df=1, lower=F) ##df will always be 1 because this is pairwise
  }
}
colnames(pw_logrank) <- levels(surv_from_di$landrace)
rownames(pw_logrank) <- levels(surv_from_di$landrace)
#Make dummy adjustment table
pw_logrank_adjBon <- pw_logrank
#Fill adjusted pvalue table.
for (i in 1:ncol(pw_logrank)) {
  pw_logrank_adjBon[,i] <- cbind(p.adjust(pw_logrank[,i], method="bonferroni"))
}
stargazer::stargazer(pw_logrank_adjBon,type="text",title="Pairwise Chisq p-values (Bonferroni adjusted)")
```
Hazards are found to be non-proportional, as is observed here, it might be a good idea to perform survival regression analysis, or pairwise log-rank testing (see below) instead of hazard ratio tests.

### Regressions
Generally a survival regression does not assume proportionality of hazards. A survival regression is fit to a distribution, defined by dist="".

```{r survival_regression,error=FALSE,warning=FALSE,echo=TRUE,message=FALSE}
library(rms)
####Survival Regression###
###This is done using functions from rms.
###psm is a survival::survreg wrapper. but the output is more handle-able for some other functions.
ddist <- datadist(surv_from_di)
options(datadist="ddist")
psurv_gaus <- psm(Surv(End, Death) ~landrace, data=surv_from_di, dist="gaussian")
psurv_logistic <- psm(Surv(End, Death) ~landrace, data=surv_from_di, dist="logistic")
psurv_lnorm <- psm(Surv(End, Death) ~landrace, data=surv_from_di, dist="lognormal")
psurv_wei <- psm(Surv(End, Death) ~landrace, data=surv_from_di, dist="weibull")
###Same with survreg()
s_reg_gaus <- survreg(Surv(End, Death) ~landrace, data=surv_from_di, dist="gaussian")
s_reg_logistic <- survreg(Surv(End, Death) ~landrace, data=surv_from_di, dist="logistic")
s_reg_lnorm <- survreg(Surv(End, Death) ~landrace, data=surv_from_di, dist="lognormal")
s_reg_wei <- survreg(Surv(End, Death) ~landrace, data=surv_from_di, dist="weibull")  
aic.scores.psurv <- rbind(
  extractAIC(psurv_wei),
  extractAIC(psurv_gaus),
  extractAIC(psurv_logistic),
  extractAIC(psurv_lnorm))
###Make useable AIC table
rownames(aic.scores.psurv) <- c("Weibull", "Gaussian", "Logist", "Lognorm")
colnames(aic.scores.psurv) <- c("df", "AIC")
stargazer::stargazer(aic.scores.psurv,type="text",title="AIC Scores")
###Call table
```
From the table above, the model with the lowest AIC score can be chosen. For this analysis, Gaussian model will be the choice for analysis. Then, one can inspect this model for significant differences.

### Plotting of parametric survival regression
It is possible, but not really easy, to plot the curves generated using parametric survival regression. These curves are the result of fitting the data to a distribution in the earlier section.
Doing this in a manner that is compatible with ggplot2 requires a nice dataframe. Below is code to generate plots of the KM estimates per block and the generated fit. This is performed for the four distributions above, and can be adapted to different distributions if necessary.

```{r generate_curves_survreg,message=FALSE,error=FALSE,warning=FALSE, echo=TRUE}

###Step 1, extract the coefficients. These are relative to landrace1 because landrace is treatment contrasted.

for (i in 1:nlevels(surv_di_fit.df$landrace)) { #For loop through landraces
  if(i==1) { #landrace1 is relative to itself, so no change
  coef_wei <- list()
  coef_logistic <- list()
  coef_gaus <- list()
  coef_lnorm <- list()
  coef_wei[i] <- coef(s_reg_wei)[i]
  coef_logistic[i] <- coef(s_reg_logistic)[i]
  coef_gaus[i] <- coef(s_reg_gaus)[i]
  coef_lnorm[i] <- coef(s_reg_lnorm)[i]
  } else { ###Other landraces are relative to 1
  coef_wei[i] <- coef(s_reg_wei)[1] + coef(s_reg_wei)[i]
  coef_logistic[i] <- coef(s_reg_logistic)[1] + coef(s_reg_logistic)[i]
  coef_gaus[i] <- coef(s_reg_gaus)[1] + coef(s_reg_gaus)[i]
  coef_lnorm[i] <- coef(s_reg_lnorm)[1] + coef(s_reg_lnorm)[i]
  }
}
##Step 2
####Store the coefficients and the scale in a new data frame, of parameters
### Keep in mind that survreg.distributions$weibull is different from rweibull, hence the difference in names.
sregparams <- data.frame(
  landrace = rep(levels(surv_from_di$landrace),4 ), #Fill with landraces
  scale.wei = exp(unlist(coef_wei)), #weibull fit scale parameters
  scale.logistic = rep(s_reg_logistic$scale, nlevels(surv_from_di$landrace)), #fill with logis scales
  scale.gaus = rep(s_reg_gaus$scale, nlevels(surv_from_di$landrace)), #fill with gaus scales
  scale.lnorm = rep(s_reg_lnorm$scale, nlevels(surv_from_di$landrace)), #fill with lnorm scale
  shape.wei =  rep(1/s_reg_wei$scale, nlevels(surv_from_di$landrace)), #shape for weibull
  shape.logistic = unlist(coef_logistic), #shape for logistic
  shape.gaus =  unlist(coef_gaus), #shape for gaus
  shape.lnorm =   unlist(coef_lnorm) #shape for lnorm
  )
##Step 3
###Calculate the "daily" value of each curve
for (i in 1:nlevels(surv_di_fit.df$landrace)){
  if(i==1) {
    wei <- list()
    logis <- list()
    gaus <- list()
    lnorm <- list()
  }
  x <- levels(surv_di_fit.df$landrace)[i]
  n <- c(1:max(surv_from_di$End))
  data <- filter(sregparams, landrace==x)
  time <- n
  wei <- cbind(wei, pweibull(
    q=n,
    scale=data$scale.wei,
    shape=data$shape.wei,
    lower.tail=FALSE))
  logis <- cbind(logis,plogis(
    q=n,
    scale=data$scale.logistic,
    location=data$shape.logistic,
    lower.tail=FALSE  ))
  gaus <- cbind(gaus,pnorm(
   q=n,
   sd=data$scale.gaus,
   mean=data$shape.gaus,
   lower.tail = F))
  lnorm <- cbind(lnorm, plnorm(
    q=n,
    sd=data$scale.lnorm,
    mean=data$shape.lnorm,
    lower.tail=F))
}

##Step 4
###Put all the curves into a data.frame that contains information on "time" and also "landrace", for compatibility with other data.frames
sreg_curves <- data.frame(
  wei.sreg = cbind(unlist(wei)),
  logis.sreg = cbind(unlist(logis)),
  gaus.sreg = cbind(unlist(gaus)),
  lnorm.sreg = cbind(unlist(lnorm)),
  landrace = rep(unlist(levels(surv_di_fit.df$landrace)),each=max(surv_from_di$End)),
  time = rep(c(1:max(surv_from_di$End)), nlevels(surv_di_fit.df$landrace))
)
##Step 5
###Turn that data.frame into a long data.frame (not used here but could be handy.)
sreg_long <- sreg_curves %>%
  group_by(landrace) %>%
  gather(key="distribution", value = "value", logis.sreg, wei.sreg, gaus.sreg, lnorm.sreg )
sreg_long$distribution <- as.factor(sreg_long$distribution)
##Levels: gaus.sreg lnorm.sreg logis.sreg wei.sreg
levels(sreg_long$distribution) <- c("Gaussian","Lognormal","Loglogistic","Weibull")
```

Now, these can be plotted and inspected visually..

```{r plot of KM+Gaussian, echo=TRUE,message=FALSE,warning=FALSE}
###Plot of KM+Gaussian  
  ggplot(surv_di_fit.df, aes(time, surv, colour = landrace)) +
    geom_step(aes(linetype=block)) +
    #geom_smooth(aes(linetype=block)) +
    geom_line(data=sreg_curves,aes(y=gaus.sreg),color="black") +
    facet_wrap(~landrace) +
    ggtitle("Kaplan-Meier estimates and fit to\nGaussian distribution")+
    theme(legend.position = "none")
```

# 8. Manuscript figures

The following scripts was used to generate figures in the manuscript.  
## Fig 5: Plot of actual AUDPS per landrace and, mean and CI of AUDPS 
```{r Fig5, echo=TRUE, error=FALSE, message=FALSE, warning=FALSE,fig.fullwidth = TRUE,eval=FALSE}

###For the purpose of visualization calculate summaries across all replicates
di_summary.noblock <- enset_xvm_long %>% group_by(plant, landrace, dpi) %>% summarise(mean(di),sd(di),sd(di)/sqrt(length(di)))
colnames(di_summary.noblock) <- c("plant", "landrace", "dpi", "mean", "sd", "se")

###Generate mean and standard deviation of disease index using the pipe function ( %>% ).
di_summary <- enset_xvm_long %>% group_by(plant, landrace, block, dpi) %>% summarise(mean(di),sd(di),sd(di)/sqrt(length(di)))
colnames(di_summary) <- c("plant", "landrace", "block", "dpi", "mean", "sd", "se") ###Assign correct columnnames

####Fig7a: Plot the actual areas per landrace
Fig5a <-
  ggplot(filter(di_summary),aes(x=dpi,y=mean)) + ###Use data di_summary
  geom_area(aes(color=block,fill=block))+
  scale_fill_manual(values = c("gray1","#999999","gray50"))+
#geom_line(aes(linetype=block)) +#show line by each block
#aes(x=dpi,y=mean,color=landrace) + ###Color by landrace, specify x and y.
#geom_area(aes(fill=landrace),position="identity",alpha=0.15) + ###Area plot, colored by landrace
facet_wrap(~landrace) + ###One plot per landrace
  labs(x = "Days after infection", y = "Mean disease index") + #Labels
  ggtitle("5A) Area under disease progression of landraces across blocks")+
  theme_bw()+
  theme(strip.text.x = element_text(size = 12,face = "bold"),
        axis.text.x =  element_text(size = 11,face = "bold"),
        axis.title.x = element_text(size = 12,face = "bold"),
        axis.title.y = element_text(size = 12,face = "bold"),
        axis.text.y =  element_text(size = 11,face = "bold"),
        legend.position = "none")

####Fig5b: Plot mean and confidence interval AUDPS per landrace###

color_fill_gray <- c("#999999","#999999","#999999","#999999","#999999",                               "#999999","#999999","#999999","#999999","#999999", "#999999","#999999","#999999","#999999","#999999",
         "#999999","#999999","#999999","#999999","#999999") # fill grey color for all landraces

####Generate plot of meanCI of AUDPS with significance letters and raw data as jittered points
cld_audps=cld_audps[order(cld_audps$audps),]# reorder by mean value to cheat ggplot
cld_audps$landrace =factor(cld_audps$landrace,levels = cld_audps$landrace)#put factors in level

Fig5b<-
  ggplot(aes(x=landrace, y=audps, fill=landrace),data=audps_data__cld) + ###Plot the audps_data
  scale_shape_manual(values=c(9, 19, 22))+
  geom_crossbar(data = cld_audps, aes(x = landrace, y = audps, ymin = lower, ymax = upper,fill=landrace), alpha=0.3,fill="gray") +
  scale_fill_manual(values = color_fill_gray)+
  geom_jitter(aes(shape=block),data = audps_data__cld) + ###with jitter overplotted, symbol shape defined by block
  geom_text(aes(x=landrace, y=-20, label=aud_let),color="black", size=3.5, data=audps_data__cld) + ###Get the letters from auc_cld
  #and write those to position y=-3
  labs(y="AUDPS") + #Y-Axis label
  ggtitle("5B) Mean and raw audps values per landrace and grouping letters") +#Title
  theme_bw()+
  theme(axis.title.y = element_text(face="bold", size=12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 11,face = "bold"),
        axis.text.x  = element_text(angle=90, vjust=0.5, size=12,face = "bold"),
        legend.position="none")

library(grid)
pushViewport(viewport(layout= grid.layout(nrow=1, ncol=2)))
print(Fig5a, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(Fig5b, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
### Fig5 done ###
# This figure might not appear in below and can be view running in R console on non-Rmakrdown Rscript file. 

```


##Fig 6 - Survival Regression
```{r Fig6, echo=TRUE, error=FALSE, message=FALSE, warning=FALSE,fig.fullwidth=TRUE,eval=FALSE}

# Script from Schandry (2017) used with modfication
survfit_to_df <- function(x) {
  strata_dummy <-NULL
  for(i in 1:length(x$strata)){
    # add vector for one strata according to number of rows of strata
    strata_dummy <- c(strata_dummy, rep(names(x$strata)[i], x$strata[i]))
  }
  #make x.df from x..
  x.df <- data.frame(
    time = x$time,
    n.risk = x$n.risk,
    n.event = x$n.event,
    surv = x$surv,
    strata = strata_dummy,
    upper = x$upper,
    lower = x$lower
  )
  zeros <- data.frame(time = 0, surv = 1, strata = names((x$strata)),
                      upper = 1, lower = 1)
  x.df <- plyr::rbind.fill(zeros, x.df)
  rm(strata_dummy)
  rm(zeros)
  return(x.df)
}

###Generate intermediate survival tables##
surv_table1 <- data.frame(subject=enset_xvm_long$subject,
                          landrace=enset_xvm_long$landrace,
                          plant=enset_xvm_long$plant,
                          block=enset_xvm_long$block)
###Fill survival table based on the enset_xvm_long table. This generates warnings. These can be ignored and come from the min()
cutoff <- c(2.23)
for (i in 1:max(enset_xvm_long$subject)) { #Go by subject
  dummy <- enset_xvm_long[enset_xvm_long$subject==i,] #generate dummy for the subject
  if (is.infinite(min(dummy$dpi[which(dummy$di >= cutoff)]))) { #If none of the di is greater than the cutoff (this is where warnings are generated, min on an empty object returns infinite!)
    surv_table1$End[i] <- max(dummy$dpi) #Generate a  observation, censoring at the maximum dpi recorded
    surv_table1$Death[i] <- 0 #Still alive, because it did not pass the cutoff
  } else { #If more than zero di are greater than the cutoff
    surv_table1$End[i] <- min(dummy$dpi[which(dummy$di >= cutoff)]) #Use the lowest dpi where condition is met
    surv_table1$Death[i] <- 1 #record as dead
  }
}
###Make survfits with the four survival tables..
surv_di_fit1 <- with(surv_table1, survfit(Surv(End, Death) ~landrace +strata(block), data=surv_from_di))

###Make survfits into surv_df
surv_di_fit1.df <- survfit_to_df(surv_di_fit1)

###However the strata field still needs to be split manually....
surv_di_fit1.df$block <- as.factor( str_split_fixed(
  matrix( nrow=length(surv_di_fit1.df$strata),ncol=2, unlist(strsplit(as.character(surv_di_fit1.df$strata),", ")), byrow=T )[,2],"=",2)[,2])
surv_di_fit1.df$landrace <- as.factor( str_split_fixed(
  matrix( nrow=length(surv_di_fit1.df$strata),ncol=2, unlist(strsplit(as.character(surv_di_fit1.df$strata),", ")), byrow=T )[,1],"=",2)[,2])
###End of intermediate table generation##

###Landraces that showed non-significant coxme hazard ration were removed.

surv_landraces0=c("Ado", "Alagena", "Arkia", "Bededet", "Bota Arkia",
                  "Gezewet",  "Godere  ", "Hae'la", "Kuro", "Lemat", "Mazia","Yesha")


kmdata <- filter(surv_di_fit1.df, landrace %in% surv_landraces0)
curvdata <- filter(sreg_curves, landrace %in% surv_landraces0)
#Labels

###Plot of KM+Gaussian  
Fig6a <- ggplot(kmdata, aes(time, surv, colour = landrace)) +
  geom_step(aes(linetype=block)) +
  geom_line(data=curvdata,aes(y=gaus.sreg),color="black") +
  facet_wrap(~landrace) + theme(legend.position="none")+
  ggtitle("Fig 6A)")+
  theme_bw()+
  labs(x = "Days post infection", y = "Fraction of plants alive") +
  theme(legend.position="none")+
  theme(axis.title.x = element_text(size = 12,face = "bold"),
        axis.title.y = element_text(size = 12,face = "bold"),
        axis.text.x = element_text(size = 8,face = "bold"),
        axis.text.y = element_text(size = 8,face = "bold"),
        strip.text = element_text(size = 12))


cb_palette <- c("yellow", "darkorchid", "red", "darkgrey", "#CDCD00",
                "#CDCD00", "peachpuff", "#228B22","dodgerblue4","#54FF9F",
                "paleturquoise2","indianred")

Fig6b <- ggplot(data=filter(surv_di_fit1.df), aes(x=time, y=surv,color=landrace))+
  #geom_step(aes(linetype=block),alpha=1)  +
  geom_line(data=curvdata,
            aes(y=lnorm.sreg,color=landrace),  size=1,  alpha=1) +
  scale_color_manual(values=cb_palette)+
  geom_point(data=curvdata,aes(y=lnorm.sreg,shape=landrace))+
  scale_shape_manual(values = c(0:5,8,9,11,25,15,16))+
  labs(x="Days post infection",
       y="Estimated survival",
       title=("Fig 6B)")) +
  coord_cartesian(xlim=c(54,155)) +
  scale_x_continuous(breaks=seq(50, 155, 15))+
  #labs(color="Legend")
  theme_bw()+
  theme(legend.position=c(0.01,1),
        legend.justification=c(0,1.6),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        axis.title.x = element_text(size = 12,face = "bold"),
        axis.title.y = element_text(size = 12,face = "bold"),
        axis.text.x = element_text(size = 8,face = "bold"),
        axis.text.y = element_text(size = 8,face = "bold"))


pushViewport(viewport(layout= grid.layout(nrow=1, ncol=2)))
print(Fig6a, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(Fig6b, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
# This figure might not appear in below and can be view running in R console on non-Rmakrdown Rscript file. 
####clear device
```

# 9. Session Info
```{r session info, echo=TRUE,warning=FALSE}
sessionInfo()
dev.off()
```
