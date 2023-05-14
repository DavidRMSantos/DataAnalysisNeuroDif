missimput=function(datamat,repvec,minrep,nneib){
  repgroups=unique(repvec)
  repvec=c(0,repvec)
  validcount=matrix(data=0,nrow=nrow(datamat),ncol=length(repgroups))
  refsets=list()
  nreps=vector()
  minval=vector()
  for (i in 1:length(repgroups)){
    nreps[i]=sum(repvec==repgroups[i])
    validcount[,i]=rowSums(!is.na(datamat[,repvec==repgroups[i]]))
    refsets[[i]]=datamat[validcount[,i]==nreps[i],repvec==repgroups[i]]
    minval[i]=min(datamat[,repvec==repgroups[i]],na.rm=T)
  }
  maxvalid=apply(validcount,1,max)
  datamat=datamat[maxvalid>=minrep,]
  validcount=validcount[maxvalid>=minrep,]
  for (i in 1:nrow(validcount)){
    for (j in 1:ncol(validcount)){
      if (validcount[i,j]>=(nreps[j]/2)){
        datamat[i,repvec==repgroups[j]]=findknn(as.numeric(datamat[i,repvec==repgroups[j]]),as.matrix(refsets[[j]]),3)
        #knn
      } else {
        values=datamat[i,repvec==repgroups[j]]
        values[is.na(values)]=minval[j]-1
        #min-const
        datamat[i,repvec==repgroups[j]]=values
      }
    }
  }
  list(data=datamat,detected=validcount)
  
}

findknn=function(seed,lookuptab,nneib){
  seedmat=matrix(seed,nrow=nrow(lookuptab),ncol=ncol(lookuptab),byrow = T)
  dist2seed=rowSums((seedmat-lookuptab)^2,na.rm=T)
  lookupord=lookuptab[order(dist2seed),]
  knnavg=colMeans(lookupord[1:nneib,])
  seed[is.na(seed)]=knnavg[is.na(seed)]
  seed
}
geral=read.csv("geral_samples_id.csv",header=T,stringsAsFactors = F)
geral[geral==0]=NA
loggeral=geral
loggeral[,2:21]=log2(geral[,2:21])



loggeralm=missimput(loggeral,c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4),2,3)

loggeralv=loggeralm[[2]]
loggeralm=loggeralm[[1]]
totdetected=rowSums(loggeralv)
table(totdetected)

summary(loggeralm[totdetected==20,2:21])
samplemed=apply(loggeralm[totdetected==20,2:21],2,median)

colSums(loggeralm[,2:21])
summary(loggeralm[,2:21])

t3pca <- prcomp(t(loggeralm[,2:21]),center=T,scale.=T)

library(factoextra)
fviz_eig(t3pca)

fviz_pca_ind(t3pca, repel = TRUE)

t3cor=cor(loggeralm[,2:21],use="pairwise")
col=colorRampPalette(c("blue","white","red"))(20)

library("gplots")

heatmap.2(t3cor, col=col, trace="none", key = TRUE, keysize = 1.5, margins = c(2,10))

write.csv(loggeralm,file="loggeralm.csv",row.names=F)

unloggeralm=loggeralm
unloggeralm[,2:21]=2^loggeralm[,2:21]

write.csv(unloggeralm,file="unloggeralm.csv",row.names=F)
