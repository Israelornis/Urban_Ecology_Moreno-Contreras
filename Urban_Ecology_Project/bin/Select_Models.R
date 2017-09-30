#####################################################
# Script to select models employing MAMs and VIFs   #
#####################################################

### Script elaborated by Israel Moreno-Contreras (Ornithology, UNAM)

### Call libraries
require(MASS)
require(clusterGeneration)
library(car) # Evaluate multicollinearity in models

### Read csv file of the clean variables
clean_data <- read.csv("..\\data\\input_data\\clean_variables.csv", header = TRUE, row.names = 1)

### Minimum Adequate Models ###
# We excluded a priori Annual Precipitation variable due it is highly correlated with Elevation

### Model for Richness
ols_Richness <- lm(Richness ~ tem + veg + ele + hfp + abs, data = clean_data)# AIC=74.27041
step(ols_Richness, direction = "backward") 
ols_Richness <- lm(Richness ~ tem + ele + abs, data = clean_data)# AIC=70.98916

### Model for the ses.PD
ols_ses.PD <- lm(ses.PD ~ tem + veg + ele + hfp + abs, data = clean_data) # AIC=21.34279
step(ols_ses.PD, direction = "backward")
ols_ses.PD <- lm(ses.PD ~ ele + abs, data = clean_data) # AIC=18.07445

### Model for the NRI
ols_NRI <- lm(NRI ~ tem + veg + ele + hfp + abs, data = clean_data) # AIC=15.84582
step(ols_NRI, direction = "backward")
ols_NRI <- lm(NRI ~ tem + ele + hfp + abs, data = clean_data) # 15.46387

### Model for the NTI
ols_NTI <- lm(NTI ~ tem + veg + ele + hfp + abs, data = clean_data) # AIC=1.592076
step(ols_NTI, direction = "backward")
ols_NTI <- lm(NTI ~ tem + veg + ele + hfp + abs, data = clean_data) # AIC=1.592076

### Model for the ses. FD
ols_ses.FD <- lm(ses.FD ~ tem + veg + ele + hfp + abs, data = clean_data) # AIC=12.99652
step(ols_ses.FD, direction = "backward")
ols_ses.FD <- lm(ses.FD ~ veg + ele + hfp + abs, data = clean_data) # AIC=12.53493

### Model for the FNRI
ols_FNRI <- lm(FNRI ~ tem + veg + ele + hfp + abs, data = clean_data) # AIC=20.00647
step(ols_FNRI, direction = "backward")
ols_FNRI <- lm(FNRI ~ ele + abs, data = clean_data) # AIC=15.58166

### Model for the FNTI
ols_FNTI <- lm(FNTI ~ tem + veg + ele + hfp + abs, data = clean_data) # AIC=7.794455
step(ols_FNTI, direction = "backward")
ols_FNTI <- lm(FNTI ~ tem + ele + hfp + abs, data = clean_data) # AIC=6.469958

### Obtain a vector for each explanatory variable
ele <- clean_data$ele
pre <- clean_data$pre
tem <- clean_data$tem
veg <- clean_data$veg
hfp <- clean_data$hfp
abs <- clean_data$abs

### Obtain a vector for each response variable
Richness <- clean_data$Richness
ses.PD <- clean_data$ses.PD
NRI <- clean_data$NRI
NTI <- clean_data$NTI ### This model is overparametrized
ses.FD <- clean_data$ses.FD
FNRI <- clean_data$FNRI
FNTI <- clean_data$FNTI

### Obtain a data frame containing explanatory variables included inn the VIF analysis
explanatory_variables <- as.data.frame(cbind(tem, veg, ele, hfp, abs))

### Vif function elaborated by https://www.r-bloggers.com/collinearity-and-stepwise-vif-selection/
vif_func<-function(in_frame,thresh=10,trace=T,...){
  
  require(fmsb)
  
  if(class(in_frame) != 'data.frame') in_frame<-data.frame(in_frame)
  
  #get initial vif value for all comparisons of variables
  vif_init<-NULL
  var_names <- names(in_frame)
  for(val in var_names){
    regressors <- var_names[-which(var_names == val)]
    form <- paste(regressors, collapse = '+')
    form_in <- formula(paste(val, '~', form))
    vif_init<-rbind(vif_init, c(val, VIF(lm(form_in, data = in_frame, ...))))
  }
  vif_max<-max(as.numeric(vif_init[,2]), na.rm = TRUE)
  
  if(vif_max < thresh){
    if(trace==T){ #print output of each iteration
      prmatrix(vif_init,collab=c('var','vif'),rowlab=rep('',nrow(vif_init)),quote=F)
      cat('\n')
      cat(paste('All variables have VIF < ', thresh,', max VIF ',round(vif_max,2), sep=''),'\n\n')
    }
    return(var_names)
  }
  else{
    
    in_dat<-in_frame
    
    #backwards selection of explanatory variables, stops when all VIF values are below 'thresh'
    while(vif_max >= thresh){
      
      vif_vals<-NULL
      var_names <- names(in_dat)
      
      for(val in var_names){
        regressors <- var_names[-which(var_names == val)]
        form <- paste(regressors, collapse = '+')
        form_in <- formula(paste(val, '~', form))
        vif_add<-VIF(lm(form_in, data = in_dat, ...))
        vif_vals<-rbind(vif_vals,c(val,vif_add))
      }
      max_row<-which(vif_vals[,2] == max(as.numeric(vif_vals[,2]), na.rm = TRUE))[1]
      
      vif_max<-as.numeric(vif_vals[max_row,2])
      
      if(vif_max<thresh) break
      
      if(trace==T){ #print output of each iteration
        prmatrix(vif_vals,collab=c('var','vif'),rowlab=rep('',nrow(vif_vals)),quote=F)
        cat('\n')
        cat('removed: ',vif_vals[max_row,1],vif_max,'\n\n')
        flush.console()
      }
      
      in_dat<-in_dat[,!names(in_dat) %in% vif_vals[max_row,1]]
      
    }
    
    return(names(in_dat))
    
  }
  
}

### Explore VIF values in each minimum adequate model ###

### Richness
ols_Richness <- lm(Richness ~ tem + ele + abs, data = clean_data) # AIC=70.98916
Richness_explanatory_variables <- as.data.frame(cbind(tem, ele, abs))

keep.dat <- vif_func(in_frame = Richness_explanatory_variables,thresh = 7.5, trace = F)
form.in<-paste('Richness ~',paste(keep.dat,collapse='+'))
mod2<-lm(form.in,data=clean_data)
summary(mod2) # AIC=71.9505 excluding temperature lm(Richness ~ ele + abs, data = clean_data)

# MAM Richness
ols_Richness <- lm(Richness ~ tem + ele + abs, data = clean_data)


### ses.PD
ols_ses.PD <- lm(ses.PD ~ ele + abs, data = clean_data) # AIC=18.07445
ses.PD_explanatory_variables <- as.data.frame(cbind(ele, abs))

keep.dat <- vif_func(in_frame = ses.PD_explanatory_variables,thresh = 7.5, trace = F)
form.in<-paste('ses.PD ~',paste(keep.dat,collapse='+'))
mod2<-lm(form.in,data=clean_data)
summary(mod2) # AIC=18.07445

# MAM ses.PD
ols_ses.PD <- lm(ses.PD ~ ele + abs, data = clean_data) 


### NRI
ols_NRI <- lm(NRI ~ tem + ele + hfp + abs, data = clean_data) # 15.46387
NRI_explanatory_variables <- as.data.frame(cbind(tem, ele, hfp, abs))

keep.dat <- vif_func(in_frame = NRI_explanatory_variables,thresh = 7.5, trace = F)
form.in<-paste('NRI ~',paste(keep.dat,collapse='+'))
mod2<-lm(form.in,data=clean_data)
summary(mod2) # AIC=30.11055 excluding temperature lm(NRI ~ ele + hfp + abs, data = clean_data)

### MAM NRI
ols_NRI <- lm(NRI ~ tem + ele + hfp + abs, data = clean_data)


### NTI
ols_NTI <- lm(NTI ~ tem + veg + ele + hfp + abs, data = clean_data) # AIC=1.592076
NTI_explanatory_variables <- as.data.frame(cbind(tem, veg, ele, hfp, abs))

keep.dat <- vif_func(in_frame = NTI_explanatory_variables,thresh = 7.5, trace = F)
form.in<-paste('NTI ~',paste(keep.dat,collapse='+'))
mod2<-lm(form.in,data=clean_data)
summary(mod2) # AIC=24.02259

# MAM NTI
ols_NTI <- lm(NTI ~ veg + ele + hfp + abs, data = clean_data)


### ses.FD
ols_ses.FD <- lm(ses.FD ~ veg + ele + hfp + abs, data = clean_data) # AIC=12.53493
ses.FD_explanatory_variables <- as.data.frame(cbind(veg, ele, hfp, abs))

keep.dat <- vif_func(in_frame = ses.FD_explanatory_variables,thresh = 7.5, trace = F)
form.in<-paste('ses.FD ~',paste(keep.dat,collapse='+'))
mod2<-lm(form.in,data=clean_data)
summary(mod2) # AIC=12.53493

### MAM ses.FD
ols_ses.FD <- lm(ses.FD ~ veg + ele + hfp + abs, data = clean_data)


### FNRI
ols_FNRI <- lm(FNRI ~ ele + abs, data = clean_data) # AIC=15.58166
FNRI_explanatory_variables <- as.data.frame(cbind(ele, abs))

keep.dat <- vif_func(in_frame = FNRI_explanatory_variables,thresh = 7.5, trace = F)
form.in<-paste('FNRI ~',paste(keep.dat,collapse='+'))
mod2<-lm(form.in,data=clean_data)
summary(mod2) # AIC=15.58166

# MAM FNRI
ols_FNRI <- lm(FNRI ~ ele + abs, data = clean_data)


### FNTI
ols_FNTI <- lm(FNTI ~ tem + ele + hfp + abs, data = clean_data) # AIC=6.469958
FNTI_explanatory_variables <- as.data.frame(cbind(tem, ele, hfp,abs))


keep.dat <- vif_func(in_frame = FNTI_explanatory_variables,thresh = 7.5, trace = F)
form.in<-paste('FNTI ~',paste(keep.dat,collapse='+'))
mod2<-lm(form.in,data=clean_data)
summary(mod2) # AIC=6.545912 excluding temperature lm(FNTI ~ ele + hfp + abs, data = clean_data)

# MAM FNTI
ols_FNTI <- lm(FNTI ~ tem + ele + hfp + abs, data = clean_data)
