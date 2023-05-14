#Linear model, significance and fold change
#RA
limma_RA_t0=read.csv("limmT0_RA.csv",header=T,stringsAsFactors = F)
limma_RA_dmso=read.csv("limmRA_DMSO.csv",header=T,stringsAsFactors = F)
sig_RA_t0=limma_RA_t0[limma_RA_t0$adj.P.Val<=0.05,]
sig_RA_dmso=limma_RA_dmso[limma_RA_t0$adj.P.Val<=0.05,]
library(dplyr)
sig_RA <-sig_RA_t0[is.element(sig_RA_t0$X, sig_RA_dmso$X),]
sig_RA=sig_RA[abs(sig_RA$logFC)>0.3,]

#NS
limma_NS_t0=read.csv("limmT0_SF.csv",header=T,stringsAsFactors = F)
sig_NS=limma_NS_t0[limma_NS_t0$adj.P.Val<=0.05,]
sig_NS=sig_NS[abs(sig_NS$logFC)>0.3,]

#DMSO
limma_DMSO_t0=read.csv("limmT0_DMSO.csv",header=T,stringsAsFactors = F)
sig_DMSO=limma_DMSO_t0[limma_DMSO_t0$adj.P.Val<=0.05,]
sig_DMSO=sig_DMSO[abs(sig_DMSO$logFC)>0.3,]


#KEGGIDs
peak2kegg=read.csv("peak2kegg.csv",header=T,stringsAsFactors = F)

kegg_RA=peak2kegg[is.element(peak2kegg$peakid,sig_RA$X),]
kegg_NS=peak2kegg[is.element(peak2kegg$peakid,sig_NS$X),]
kegg_DMSO=peak2kegg[is.element(peak2kegg$peakid,sig_DMSO$X),]

#Para ver o número de m/z identificados sem as condições do T1
#Metabolitos significamente alterados
Numbsig=full_join(sig_RA,sig_NS)
Numbsig=full_join(Numbsig,sig_DMSO)
Numbsig=unique(Numbsig[1])
#m/z com um respetivo KEGG id associado
AllData=read.csv("unloggeralm.csv",header=T,stringsAsFactors = F)
Numbid=peak2kegg[is.element(peak2kegg$peakid,AllData$peakID),]
sum(Numbid$allkegg=="")

#KEGGIDs, FC and adj.p.val (Só tem os closests KEGG ids)
#Para analise complementar quanto ao aumento ou diminuição de determinado metabolito
#RA
metabolite_RA=inner_join(sig_RA,peak2kegg, by=c("X"="peakid"))
metabolite_RA=select(metabolite_RA,X,allkegg,allkeggnames,adj.P.Val,logFC)
metabolite_RA=metabolite_RA[!metabolite_RA$allkegg=="",]
write.csv(metabolite_RA,file="metabolite_RA.csv",row.names = F)
#NS
metabolite_NS=inner_join(sig_NS,peak2kegg, by=c("X"="peakid"))
metabolite_NS=select(metabolite_NS,X,allkegg,allkeggnames,adj.P.Val,logFC)
metabolite_NS=metabolite_NS[!metabolite_NS$allkegg=="",]
write.csv(metabolite_NS,file="metabolite_NS.csv",row.names = F)
#DMSO
metabolite_DMSO=inner_join(sig_DMSO,peak2kegg, by=c("X"="peakid"))
metabolite_DMSO=select(metabolite_DMSO,X,allkegg,allkeggnames,adj.P.Val,logFC)
metabolite_DMSO=metabolite_DMSO[!metabolite_DMSO$allkegg=="",]
write.csv(metabolite_DMSO,file="metabolite_DMSO.csv",row.names = F)
#Comuns
metabolite_common0=as.data.frame(intersect(metabolite_RA$allkeggnames,metabolite_NS$allkeggnames))
metabolite_common=as.data.frame(intersect(metabolite_common0$`intersect(metabolite_RA$allkeggnames, metabolite_NS$allkeggnames)`,metabolite_DMSO$allkeggnames))
write.csv(metabolite_common,file="metabolite_common",row.names = F)

#Tem todos os KEGGIDs
#Este aumentou, perguntar
kegg_RAid=data.frame(id=unlist(strsplit(kegg_RA$allkegg,split="|",fixed=T)))
kegg_RAid=unique(kegg_RAid)

kegg_NSid=data.frame(id=unlist(strsplit(kegg_NS$allkegg,split="|",fixed=T)))
kegg_NSid=unique(kegg_NSid)

kegg_DMSOid=data.frame(id=unlist(strsplit(kegg_DMSO$allkegg,split="|",fixed=T)))
kegg_DMSOid=unique(kegg_DMSOid)

kuniverse=data.frame(id=unlist(strsplit(peak2kegg$allkegg,split="|",fixed=T)))

#Pathway analysis do Metaboanalist
write.csv(kegg_RAid,file="kegg_RAid.csv",row.names = F)
write.csv(kegg_NSid,file="kegg_NSid.csv",row.names = F)
write.csv(kegg_DMSOid,file="kegg_DMSOid.csv",row.names = F)

write.csv(kuniverse,file="kegg_univ.csv",row.names = F)

#rwr analysis

library(igraph)

load("ghsa_main.RData")
load("hsa_met_tab.RData")
load("hsa_path_sets.RData")
load("hsa_path_freq.RData")


RAnodes=1*(is.element(hsa_met_tab$id,kegg_RAid$id))
NSnodes=1*(is.element(hsa_met_tab$id,kegg_NSid$id))
DMSOnodes=1*(is.element(hsa_met_tab$id,kegg_DMSOid$id))

RArwr=page_rank(ghsa_main,directed=F,damping=0.5,personalized=RAnodes)

plot(degree(ghsa_main),RArwr$vector)

NSrwr=page_rank(ghsa_main,directed=F,damping=0.5,personalized=NSnodes)

DMSOrwr=page_rank(ghsa_main,directed=F,damping=0.5,personalized=DMSOnodes)

pathrwr=function(rwrvec,pathsets){
  pathscores=vector()
  for (i in 1:length(pathsets)){
    pathscores[i]=sum(rwrvec[pathsets[[i]]])
  }
  pathscores
}

pathrwrp=function(g,seedvec,pathsets,nrep,d){
  rwrvec=page_rank(g,directed=F,damping=d,personalized=seedvec)$vector
  pathscores=pathrwr(rwrvec,pathsets)
  pvec=rep(0,length(pathscores))
  cumscore=rep(0,length(pathscores))
  for (i in 1:nrep){
    rseed=sample(seedvec)
    randrwrvec=page_rank(g,directed=F,damping=d,personalized=rseed)$vector
    randpathscores=pathrwr(randrwrvec,pathsets)
    pvec=pvec+1*(randpathscores>pathscores)
    cumscore=cumscore+randpathscores
  }
  pvec=pvec/nrep
  fcvec=pathscores/(cumscore/nrep)
  list(p=pvec,fc=fcvec)
}

pRA=pathrwrp(ghsa_main,RAnodes,hsa_path_sets,5000,0.5)

hsa_path_freq$pRA=pRA$p
hsa_path_freq$fcRA=pRA$fc

pNS=pathrwrp(ghsa_main,NSnodes,hsa_path_sets,5000,0.5)

hsa_path_freq$pNS=pNS$p
hsa_path_freq$fcNS=pNS$fc

pDMSO=pathrwrp(ghsa_main,DMSOnodes,hsa_path_sets,5000,0.5)


hsa_path_freq$pDMSO=pDMSO$p
hsa_path_freq$fcDMSO=pDMSO$fc

#Pathways
alldiff = hsa_path_freq[hsa_path_freq$pRA<0.05 & hsa_path_freq$pNS<0.05 & hsa_path_freq$pDMSO<0.05,]
RAdiff = hsa_path_freq[hsa_path_freq$pRA<0.05 & hsa_path_freq$pNS>0.05 & hsa_path_freq$pDMSO>0.05,]
NSdiff = hsa_path_freq[hsa_path_freq$pRA>0.05 & hsa_path_freq$pNS<0.05 & hsa_path_freq$pDMSO>0.05,]
DMSOdiff = hsa_path_freq[hsa_path_freq$pRA>0.05 & hsa_path_freq$pNS>0.05 & hsa_path_freq$pDMSO<0.05,]

RA_NSdiff = hsa_path_freq[hsa_path_freq$pRA<0.05 & hsa_path_freq$pNS<0.05 & hsa_path_freq$pDMSO>0.05,]
RA_DMSOdiff = hsa_path_freq[hsa_path_freq$pRA<0.05 & hsa_path_freq$pNS>0.05 & hsa_path_freq$pDMSO<0.05,]
NS_DMSOdiff = hsa_path_freq[hsa_path_freq$pRA>0.05 & hsa_path_freq$pNS<0.05 & hsa_path_freq$pDMSO<0.05,]

write.csv(alldiff,file="alldiff.csv", row.names = F)
write.csv(RAdiff,file="RAdiff.csv", row.names = F)
write.csv(NSdiff,file="NSdiff.csv", row.names = F)
write.csv(DMSOdiff,file="DMSOdiff.csv", row.names = F)
write.csv(RA_NSdiff,file="RA_NSdiff.csv", row.names = F)
write.csv(RA_DMSOdiff,file="RA_DMSOdiff.csv", row.names = F)
write.csv(NS_DMSOdiff,file="NS_DMSOdiff.csv", row.names = F)
#KEGG Color map
#RA
keggRAmap1=anti_join(kegg_RAid,kegg_NSid)
keggRAmap=anti_join(keggRAmap1,kegg_DMSOid)
keggRAmap=cbind(keggRAmap, "red")
#NS
keggNSmap1=anti_join(kegg_NSid,kegg_RAid)
keggNSmap=anti_join(keggNSmap1,kegg_DMSOid)
keggNSmap=cbind(keggNSmap, "blue")
#DMSO
keggDMSOmap1=anti_join(kegg_DMSOid,kegg_NSid)
keggDMSOmap=anti_join(keggDMSOmap1,kegg_RAid)
keggDMSOmap=cbind(keggDMSOmap, "yellow")
#Comuns em todos
keggALLmap1=inner_join(kegg_RAid,kegg_NSid)
keggALLmap=inner_join(keggALLmap1,kegg_DMSOid)
keggALLmap=cbind(keggALLmap, "black")
#RA e NS
keggRA_NSmap1=inner_join(kegg_RAid,kegg_NSid)
keggRA_NSmap=anti_join(keggRA_NSmap1,kegg_DMSOid)
keggRA_NSmap=cbind(keggRA_NSmap, "purple")
#RA e DMSO
keggRA_DMSOmap1=inner_join(kegg_RAid,kegg_DMSOid)
keggRA_DMSOmap=anti_join(keggRA_DMSOmap1,kegg_NSid)
keggRA_DMSOmap=cbind(keggRA_DMSOmap, "orange")
#NS e DMSO
keggNS_DMSOmap1=inner_join(kegg_NSid,kegg_DMSOid)
keggNS_DMSOmap=anti_join(keggNS_DMSOmap1,kegg_RAid)
keggNS_DMSOmap=cbind(keggNS_DMSOmap, "green")

write.csv(keggRAmap,file="keggRAmap.csv",row.names = F)
write.csv(keggNSmap,file="keggNSmap.csv",row.names = F)
write.csv(keggDMSOmap,file="keggDMSOmap.csv",row.names = F)
write.csv(keggALLmap,file="keggALLmap.csv",row.names = F)
write.csv(keggRA_NSmap,file="keggRA_NSmap.csv",row.names = F)
write.csv(keggRA_DMSOmap,file="keggRA_DMSOmap.csv",row.names = F)
write.csv(keggNS_DMSOmap,file="keggNS_DMSOmap.csv",row.names = F)
