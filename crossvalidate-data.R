

# function validate, called by validatemain.R
# running the "source-call"

# pass to val: b=number of bootstraps; dat=data table
val<-function(b,dat,more,maxgrid=10000){

  # generate a collector for plotting variables with 3 columns (predictor, sens, spec) and rows = (predictors*b)

  #vplot<-matrix(nrow=b*9,ncol=3) #9 = ncol(dat)-1
#colnames(vplot)<-c("predictor","Sensitivity","Specificity")

#lis<-vector("list", 2*(ncol(dat)-3)) #-1 statt -3
#names(lis)<-c("quality.sens", "quality.spec", rep(colnames(dat[2:(ncol(dat)-3)]),2)) #-1 statt -3
#lis[[1]]<-c(rep("Sensitivity",b))
#lis[[2]]<-c(rep("Specificity",b))

#  mata<-matrix(,nrow=b,ncol=(ncol(dat)-3))
#colnames(mata)<-c("quality",colnames(dat[2:(ncol(dat)-3)]))
  
#pairwise analysis for all prediction data in table=second column onwards
  for(m in 2:(ncol(dat)-2)){ #ncol(dat)
    
# select column
dati<-dat[,c(1,m)]

# size of training dataset:
dimi<-dim(dati)
n1<-round(0.66*dimi[1])
# for bootstrap cross validation:
# training data is a random sample of 2/3 of the data
# test data is the complementary dataset

# set sensitive grid size for thresholds 
thresh<-c(seq(0,100,0.05),seq(100,maxgrid,1))
# thresh<-dat[,1]

lt<-length(thresh)

# to store sensitivity and specificity at increasing values 
# for the threshold, which is needed to apply Maria`s optimality 
# criteria on those values afterwards.

sens<-matrix(nrow=b,ncol=lt)
spec<-matrix(nrow=b,ncol=lt)

# to store sensitivity and specificity at the optimal threshold, at the 0.5%, 2%, 50nM, 500nM and 5000nM threshold
# for each bootstrap sample
optithresh<-rep(0,b)
sens.optithresh<-rep(0,b)
spec.optithresh<-rep(0,b)


# to store validation values for sensitivity and specificity at different thresholds
sens.val<-rep(0,b)
spec.val<- rep(0,b) 
acc.val<-rep(0,b)
sens.05p.val<-rep(0,b)
spec.05p.val<-rep(0,b)
acc.05p.val<-rep(0,b)
sens.2p.val<-rep(0,b)
spec.2p.val<-rep(0,b)
acc.2p.val<-rep(0,b)
sens.50.val<-rep(0,b)
spec.50.val<-rep(0,b)
acc.50.val<-rep(0,b)
sens.500.val<-rep(0,b)
spec.500.val<-rep(0,b)
acc.500.val<-rep(0,b)
sens.5000.val<-rep(0,b)
spec.5000.val<-rep(0,b) 
acc.5000.val<-rep(0,b)
sens.meanopti.val<-rep(0,b)
spec.meanopti.val<-rep(0,b)
acc.meanopti.val<-rep(0,b)
# loop for bootstrap samples
for(j in 1:b){
  # get bootstrap sample used as training sample in the following:
  svec<-sample(1:dimi[1],n1,replace=FALSE)
  sami<-dati[svec,]
  # get complementary test dataset:
  tami<-dati[-svec,]
  
  k<-1
  # loop for increasing thresholds
  # aim: detect optimal threshold
  for(i in thresh){
    TP<-which(sami[,2]<=i & sami[,1]=="binder")
    TN<-which(sami[,2]>i & sami[,1]=="nb")
    sens[j,k]<-length(TP)/length(which(sami[,1]=="binder"))
    spec[j,k]<-length(TN)/length(which(sami[,1]=="nb"))
    k<-k+1
  }
  
  # Maria`s optimality criteria for trade off between sensitivity and specificity
  
  #  "Die Kriterien die ich für die Grenzwert-Empfehlung verwende sind folgende:
  
  #    1)  Spezifität mindestens 2/3 = FPR maximal 1/3
  #    2)  Sensitivität mindestens 2x(1-Spezifität)
  #    3)  Innerhalb 1) und 2) die höchstmögliche Sensitivität"
  
  # apply for each bootstrap sample
  cond1<-which(spec[j,]>=2/3 & sens[j,]>=2*(1-spec[j,]))
  # which.max returns index of the first maximum
  cond2<-which.max(sens[j,cond1])
  optithresh[j]<-thresh[cond1[cond2]]
  # attained sensitivity and specificity at optimal threshold: 
  sens.optithresh[j]<-sens[j,cond1[cond2]]
  spec.optithresh[j]<-spec[j,cond1[cond2]]
  
  # now validation of the optimal threshold based on the test dataset:
  
  TP.val<-which(tami[,2]<=optithresh[j] & tami[,1]=="binder")
  TN.val<-which(tami[,2]>optithresh[j] & tami[,1]=="nb")

  sens.val[j]<-length(TP.val)/length(which(tami[,1]=="binder"))
  spec.val[j]<-length(TN.val)/length(which(tami[,1]=="nb"))
  acc.val[j]<-(length(TP.val)+length(TN.val))/(length(which(tami[,1]=="binder"))+length(which(tami[,1]=="nb")))
  
  # now validation of the common threshold based on the test dataset:    
  TP.val<-which(tami[,2]<=0.5 & tami[,1]=="binder")
  TN.val<-which(tami[,2]>0.5 & tami[,1]=="nb")
  
  sens.05p.val[j]<-length(TP.val)/length(which(tami[,1]=="binder"))
  spec.05p.val[j]<-length(TN.val)/length(which(tami[,1]=="nb"))
  acc.05p.val[j]<-(length(TP.val)+length(TN.val))/(length(which(tami[,1]=="binder"))+length(which(tami[,1]=="nb")))

  TP.val<-which(tami[,2]<=2 & tami[,1]=="binder")
  TN.val<-which(tami[,2]>2 & tami[,1]=="nb") 
  
  sens.2p.val[j]<-length(TP.val)/length(which(tami[,1]=="binder"))
  spec.2p.val[j]<-length(TN.val)/length(which(tami[,1]=="nb"))
  acc.2p.val[j]<-(length(TP.val)+length(TN.val))/(length(which(tami[,1]=="binder"))+length(which(tami[,1]=="nb")))
  
  TP.val<-which(tami[,2]<=50 & tami[,1]=="binder")
  TN.val<-which(tami[,2]>50 & tami[,1]=="nb")
  
  sens.50.val[j]<-length(TP.val)/length(which(tami[,1]=="binder"))
  spec.50.val[j]<-length(TN.val)/length(which(tami[,1]=="nb"))
  acc.50.val[j]<-(length(TP.val)+length(TN.val))/(length(which(tami[,1]=="binder"))+length(which(tami[,1]=="nb")))
  
  TP.val<-which(tami[,2]<=500 & tami[,1]=="binder")
  TN.val<-which(tami[,2]>500 & tami[,1]=="nb")
  
  sens.500.val[j]<-length(TP.val)/length(which(tami[,1]=="binder"))
  spec.500.val[j]<-length(TN.val)/length(which(tami[,1]=="nb"))
  acc.500.val[j]<-(length(TP.val)+length(TN.val))/(length(which(tami[,1]=="binder"))+length(which(tami[,1]=="nb")))
  
  TP.val<-which(tami[,2]<=5000 & tami[,1]=="binder")
  TN.val<-which(tami[,2]>5000 & tami[,1]=="nb")
  
  sens.5000.val[j]<-length(TP.val)/length(which(tami[,1]=="binder"))
  spec.5000.val[j]<-length(TN.val)/length(which(tami[,1]=="nb"))
  acc.5000.val[j]<-(length(TP.val)+length(TN.val))/(length(which(tami[,1]=="binder"))+length(which(tami[,1]=="nb")))
  
  TP.val<-which(tami[,2]<=more[m-1] & tami[,1]=="binder")
  TN.val<-which(tami[,2]>more[m-1] & tami[,1]=="nb")
  
  sens.meanopti.val[j]<-length(TP.val)/length(which(tami[,1]=="binder"))
  spec.meanopti.val[j]<-length(TN.val)/length(which(tami[,1]=="nb"))
  acc.meanopti.val[j]<-(length(TP.val)+length(TN.val))/(length(which(tami[,1]=="binder"))+length(which(tami[,1]=="nb")))
  
  }


# plot for sensitivity and specificity values depending on optimal threshold for IC50.
# optimal thresholds with respect to Marias cirteria based on b bootstrap samples.
# That means each triple of the following plot shows
# the optimal IC50 Threshold together with attained sensitivity and specificity
# at this threshold.

sens.optithresh.mean<-tapply(sens.optithresh,optithresh,mean)
spec.optithresh.mean<-tapply(spec.optithresh,optithresh,mean)
optithresh.mean<-as.numeric(names(sens.optithresh.mean))

xs<-"optimaler IC50 Bootstrap-Wert"
ys<-"Sensitivität (grün) und Spezifität (blau)"

# create a directory in the current working directory that has the name of the currently analyzed column of dat
dir.create(paste(colnames(dati[2])))

# the working directory is set to the newly created directory
setwd(paste(colnames(dati[2])))

# the following command starts to put output into a new svg-file
svg("Bootstrap-Plot.svg")
plot(optithresh.mean,sens.optithresh.mean,type="l",col="green",xlab=xs,ylab=ys, main=colnames(dati[2]))
lines(optithresh.mean,spec.optithresh.mean,type="l",col="blue")
abline(h=2/3,v=500)
dev.off() # ends the stream to svg file


# boxplots of validated sensitivity and specificity 
# (depending on computed bootstrap threshold)

# the following command starts to put output into a new svg-file
svg("expectedSensSpec.svg")
# generate boxplot, as boxplot(vector of x-values, [vector of y-values], parameters)
boxplot(c(sens.val,spec.val)~c(rep("Sensitivity",length(sens.val)),rep("Specificity",length(spec.val))),par(cex.lab=1.5, cex.axis=1.5), ylim=c(0,1), boxwex=0.5, main=colnames(dati[2]))
dev.off() # ends the stream to svg file


#print(m)
#lis[[m+1]]<-sens.val
#lis[[m+(ncol(dat)-3)]]<-spec.val


# bootstrap estimate for the 95% confidence interval of sensitivity and specificity
CI.sens.val<-quantile(sens.val,probs=c(0.025,0.5,0.975))
CI.spec.val<-quantile(spec.val,probs=c(0.025,0.5,0.975))

optithresh<-median(optithresh)

# store output variables of function val

out.val<-list(optithresh,CI.sens.val,CI.spec.val)
names(out.val)<-c("median of optimal threshold","bootstrap quantiles for sensitivity",
                  "bootstrap quantiles for specificity")

# go back to parent directory
setwd("..")



# start writing following output into a text-file that can be added up by following calls
sink(file="output.txt", append=TRUE)
cat("=========================\n")
cat(colnames(dati[2])) 
cat("\n")
cat("=========================\n")
cat("outcome\n")
print(out.val) # return(out.val)
cat("bootstrap Sensitivity Optimal threshold\n")
print(sens.val)
cat("bootstrap Specificity Optimal threshold\n")
print(spec.val)
cat("bootstrap Accuracy Optimal threshold\n")
print(acc.val)
cat("bootstrap Sensitivity 0.5% threshold\n")
print(sens.05p.val)
cat("bootstrap Specificity 0.5% threshold\n")
print(spec.05p.val)
cat("bootstrap Accuracy 0.5% threshold\n")
print(acc.05p.val)
cat("bootstrap Sensitivity 2% threshold\n")
print(sens.2p.val)
cat("bootstrap Specificity 2% threshold\n")
print(spec.2p.val)
cat("bootstrap Accuracy 2% threshold\n")
print(acc.2p.val)
cat("bootstrap Sensitivity 50nM threshold\n")
print(sens.50.val)
cat("bootstrap Specificity 50nM threshold\n")
print(spec.50.val)
cat("bootstrap Accuracy 50nM threshold\n")
print(acc.50.val)
cat("bootstrap Sensitivity 500nM threshold\n")
print(sens.500.val)
cat("bootstrap Specificity 500nM threshold\n")
print(spec.500.val)
cat("bootstrap Accuracy 500nM threshold\n")
print(acc.500.val)
cat("bootstrap Sensitivity 5000nM threshold\n")
print(sens.5000.val)
cat("bootstrap Specificity 5000nM threshold\n")
print(spec.5000.val)
cat("bootstrap Accuracy 5000nM threshold\n")
print(acc.5000.val)
cat(paste("bootstrap Sensitivity ",more[m-1],"nM threshold\n",sep=""))
print(sens.meanopti.val)
cat(paste("bootstrap Specificity ",more[m-1],"nM threshold\n",sep=""))
print(spec.meanopti.val)
cat(paste("bootstrap Accuracy ",more[m-1],"nM threshold\n",sep=""))
print(acc.meanopti.val)
sink() # end writing
print(out.val)
} #end of for loop

  
# the following command starts to put all plots into a new svg-file
#svg(paste("overviewSensSpec.svg",sep=""))  

#boxplot(vplot[,2]~vplot[,1])    

#dev.off() # ends the stream to svg file  
} # end of function val



