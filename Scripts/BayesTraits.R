library(reshape2)
library(coda)

filelist = read.csv('Data/BayesTrait_ListFiles.txt', sep='\t')
names(filelist) = 'FileName'

filelist$diff = NA
filelist$LHdep = NA
filelist$LHind = NA
filelist$AICdep = NA
filelist$AICind = NA

for (i in seq(1, length(filelist$FileName))){
  foc = filelist[i,1]
  dep = scan(paste0('Data/', foc, '_dependent_MC_Log.txt'), what='numeric', sep='\n')
  startH = grep("Iteration", dep)[3]
  dep = dep[-c(1:startH)]
  dep2 = colsplit(dep, "\t", c("Iteration", "Lh", "treeno", "q12","q13","q21","q24","q31","q34","q42","q43", "Rest"))[,-c(3,12)]
  dep2 = dep2[-dim(dep2)[1],]
  
  dep2= tail(dep2, n=50000)
  thin= seq(1,dim(dep2)[1],2)
  dep2 = dep2[thin,]
  
  ESS = effectiveSize(dep2[,-1])
  if (all(ESS > 200)) {
    dep_mean= mean(dep2$Lh)
    filelist$LHdep[i]=dep_mean
    AIC_dep = (2*8)-2*dep_mean
    filelist$AICdep[i]=AIC_dep
    
  } else {
    print(i); print( ESS )
    AIC_dep = NA
    warning("LH of the dependend model has not converged") 
  }
  
  inddep = scan(paste0('Data/', foc, '_independent_MC_Log.txt'), what='numeric', sep='\n')
  startH = grep("Iteration", inddep)[3]
  ind = inddep[-c(1:startH)]
  ind2 = colsplit(ind, "\t", c("Iteration", "Lh", "treeno", "a1", "b1", "a2", "b2", "Rest"))[,c(-3,-8)]
  ind2 = ind2[-dim(ind2)[1],]
  
  ind2= tail(ind2, n=50000)
  thin= seq(1,dim(ind2)[1],2)
  ind2= ind2[thin,]
  
  ESS = effectiveSize(ind2[,-1])
  if (all(ESS > 200)) {
    ind_mean= mean(ind2$Lh)
    filelist$LHind[i]= ind_mean
    AIC_ind = (2*4)-2*ind_mean
    filelist$AICind[i]=AIC_ind
    
  } else {
    print(i); print( ESS )
    AIC_ind = NA
    warning("LH of the independend model has not converged") 
  }
  
  filelist$diff[i] =   AIC_ind - AIC_dep
  
  print(i)
}

candsMC = filelist[ filelist$diff > 0 , ]
