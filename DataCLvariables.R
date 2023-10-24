
## Select medical variables 

#We are mostly interested to see the effect of medical covariates on the transition probabilities. 
#Therefore, we will be selecting which variables have enough measurements to be accounted 
#for in the model. We will be taking the measurement that was taken at the first moment of diagnosis.

# take out the baseline columns
names.df <- names(adnimerge)
names.df <- names.df[substr(names.df, nchar(names.df)-2, nchar(names.df)) != ".bl"]
takeout <- c("COLPROT", "ORIGPROT", "PTID", "SITE", "VISCODE",
             "AGE", "PTGENDER", "PTEDUCAT", "PTETHCAT", "PTRACCAT", "PTMARRY",
             "Month", "M", "DX", "FSVERSION", "IMAGEUID")
names.df <- setdiff(names.df, takeout)
print(names.df)
adni_medical <- adnimerge[,names.df]

rm(names.df, takeout)

# look at the missingness (%) per column
missing <- as.data.frame(sapply(adni_medical, 
                                function(x) sum(is.na(x))/length(x)))

# add medical information to relevant follow-up moments
long4 <- long4 %>%
  mutate(RID = as.numeric(RID),
         EXAMDATE = as.Date(EXAMDATE))

adni_medical <- adni_medical %>%
  mutate(RID = as.numeric(RID),
         EXAMDATE = as.Date(EXAMDATE))

joined <- left_join(long4, adni_medical, by = c("RID", "EXAMDATE"))

# look at the missingness (%) per column
missing2 <- as.data.frame(sapply(joined[, 6:ncol(joined)], 
                                 function(x) sum(is.na(x))/length(x)))

# only values at first moment of diagnosis
adni_medical_first <- joined %>%
  arrange(RID, EXAMDATE) %>%
  mutate(interim = paste(RID, DX)) %>%
  filter(!duplicated(interim))

# look at the missingness (%) per column
missing3 <- as.data.frame(sapply(adni_medical_first[, 6:(ncol(adni_medical_first)-1)], 
                                 function(x) sum(is.na(x))/length(x)))  
names(missing3) <- "FreqMissing"
missing3$Var <- rownames(missing3)

rm(joined, missing, missing2)


### Consider only columns with enough values 

#We will proceed only with the medical variables that have no more than $1/3$ of data missing.

VarNames_EnoughInfo <- missing3$Var[missing3$FreqMissing < 1/3]
print(VarNames_EnoughInfo)
adni_medical_first2 <- adni_medical_first %>%
  select(RID, EXAMDATE, DX, VISCODE, Month.bl, all_of(VarNames_EnoughInfo))

rm(adni_medical_first, VarNames_EnoughInfo, missing3)


# Data preparation 

#Since a lot of the measured medical variables are of similar nature, 
#they will surely show collinearity. Therefore, we will be removing variables 
#that are strongly correlated with (an)other variable(s). For this, we'll use a 
#cut-off value of $0.65$ for the correlation coefficient. 


# calculate correlations
cormatrix <- cor(adni_medical_first2[, 6:ncol(adni_medical_first2)], use = "pairwise.complete.obs")

# generate correlation plot 
pdf("corplot.pdf")
corrplot(cormatrix, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, tl.cex = 0.6)
dev.off()

# extract pairs with > 0.65 correlation coefficient 
k <- 1 
pairsmatrix <- matrix(NA, nrow(cormatrix)*ncol(cormatrix), ncol = 2)

for(i in 1:nrow(cormatrix)){
  for(j in 1:ncol(cormatrix)){
    if(abs(cormatrix[i,j]) > 0.65){
      pairsmatrix[k,1] <- rownames(cormatrix)[i]
      pairsmatrix[k,2] <- colnames(cormatrix)[j]
      k <- k + 1 
    }
  }
}

pairsmatrix <- as.data.frame(pairsmatrix) %>%
  filter(!is.na(V1)) %>%
  filter(V1 != V2)

# keep only unique pairs 
for(i in 1:nrow(pairsmatrix)){
  if(pairsmatrix[i,1] > pairsmatrix[i,2]) pairsmatrix[i,] = pairsmatrix[i,c(2,1)]
}

pairsmatrix <- pairsmatrix %>%
  distinct(V1, V2)

# look at missingness in these variables 
v <- unique(c(pairsmatrix$V1, pairsmatrix$V2))
sort(apply(adni_medical_first2[, v], 2, function(x) sum(is.na(x))))

# remove variables that have high correlations (> 0.65)
med_variables <- colnames(adni_medical_first2[, 6:ncol(adni_medical_first2)])
med_variables <- setdiff(med_variables, c("ADAS11", "ADAS13", "ADASQ4", "CDRSB",
                                          "mPACCdigit", "mPACCtrailsB", "LDELTOTAL",
                                          "RAVLT.learning", "RAVLT.perc.forgetting",
                                          "WholeBrain", "Entorhinal", "MidTemp"))

# clear environment
rm(i,j, k, v, pairsmatrix, cormatrix)


# Final dataframe 

#This dataframe contains only the variables we are interested in to possibly 
#use for analysis. 


# select non-medical final variables 
nonmed_variables <- c("RID", "DX", "Month.bl", 
                      "PTGENDER", "PTETHCAT", "PTRACCAT", 
                      "PTMARRY", "AGE", "PTEDUCAT")

# create final data frame 
adnimerge <- adnimerge %>%
  mutate(RID = as.numeric(RID),
         DX = as.character(DX))

adni_medical_first3 <- adni_medical_first2[,c("RID", "DX", 
                                              "Month.bl", med_variables)]

adni_final <- left_join(long4, adnimerge[,nonmed_variables])
adni_final <- left_join(adni_final, adni_medical_first3)
adni_final$VISCODE <- adni_final$EXAMDATE <- NULL

# dataframe for table 1 
baseline <- left_join(adni_final[adni_final$Month.bl == 0,], adnimerge) %>%
  select(colnames(adni_final), ADAS11, ADAS13, ADASQ4, CDRSB, mPACCdigit, 
         mPACCtrailsB, LDELTOTAL, RAVLT.learning, RAVLT.perc.forgetting,
         WholeBrain, Entorhinal, MidTemp)



# Clean environment

#Save only the objects that are necessary for the files that follow hereafter.

rm(adni_medical, adni_medical_first2, adni_medical_first3, adnimerge,
   deaths, long3, long4, trans3, trans4)



rm(Dxconfr, countpatients, countpatients1, countpatients2, countpatients4, countpatients5, lon3couter, longcout)
