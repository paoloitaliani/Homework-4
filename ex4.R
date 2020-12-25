library(mclust)
library(cluster)
library(mclust)
library(smacof)
library(factoextra)
library(dbscan)
library(knitr)
setwd("~/Desktop/Università/MODERN STATISTICS AND BIG DATA ANALYTICS ")
d=read.table("wdbc.txt",header = F,sep=",")
wdbcc <- d[,3:12]
wdbcdiag <- as.integer(d[,2])
wdbcc.2=wdbcc^2
wdbcc.2s=scale(wdbcc.2)


###Exercise 1
 
ARI<- matrix(0, nrow = 1, ncol = 9)
colnames(ARI)=c("Kmeans","Average Linkage","Complete Linkage","Single Linkage","Ward's Method", "pam","VVV","EVI","VEI")
rownames(ARI)=c("ARI")

kmeans=kmeans(wdbcc.2s,2,nstart=100)
ARI[1,1]=adjustedRandIndex(kmeans$cluster,wdbcdiag)



dist.m=dist(wdbcc.2s,method = "euclidean")
clust.avg <- hclust(as.dist(dist.m),method="average")
clust.avg=cutree(clust.avg,2)
ARI[1,2]=adjustedRandIndex(clust.avg,wdbcdiag)

clust.c <- hclust(as.dist(dist.m),method="complete")
clust.c=cutree(clust.c,2)
ARI[1,3]= adjustedRandIndex(clust.c,wdbcdiag)


clust.s <- hclust(as.dist(dist.m),method="single")
clust.s=cutree(clust.s,2)
ARI[1,4]= adjustedRandIndex(clust.s,wdbcdiag)

clust.w <- hclust(as.dist(dist.m),method="ward.D2")
clust.w=cutree(clust.w,2)
ARI[1,5]= adjustedRandIndex(clust.w,wdbcdiag)

pam=pam(as.dist(dist.m),2)
ARI[1,6] = adjustedRandIndex(pam$cluster,wdbcdiag)




clust.m=Mclust(wdbcc,modelName="VVV",G=2)
ARI[1,7]= adjustedRandIndex(clust.m$classification,wdbcdiag)


clust.m=Mclust(wdbcc,modelName="EVI",G=2)
ARI[1,8]=adjustedRandIndex(clust.m$classification,wdbcdiag)

clust.m=Mclust(wdbcc,modelName="VEI",G=2)
ARI[1,9]=adjustedRandIndex(clust.m$classification,wdbcdiag)

pairs(wdbcc,col=wdbcdiag)
pairs(wdbcc,col=clust.m$classification)

clust.m$classification[clust.m$classification==1]=3
clust.m$classification[clust.m$classification==2]=1
clust.m$classification[clust.m$classification==3]=2

pairs(wdbcc,col=wdbcdiag)
pairs(wdbcc,col=clust.m$classification)





#Exercise 2

ASW<- matrix(0, nrow = 2, ncol = 5)
colnames(ASW)=c("Single Linkage","Complete Linkage","Average Linkage","Ward's Method", "pam")
rownames(ASW)=c("ASW","Number of clusters")

#Single
clust.s <- hclust(as.dist(dist.m),method="single")
tasw <- NA
tclusk <- list()
tsil <- list()
for (k in 2:30){
  tclusk[[k]] <- cutree(clust.s,k)
  tsil[[k]] <- silhouette(tclusk[[k]],dist=as.dist(dist.m))
  tasw[k] <- summary(silhouette(tclusk[[k]],dist=as.dist(dist.m)))$avg.width
}
ASW[2,1]=match(max(tasw[-1]),tasw[-1])+1
ASW[1,1]=max(tasw[-1])
ASW[2,]=as.integer(ASW[2,])
#Complete
clust.c <- hclust(as.dist(dist.m),method="complete")
tasw <- NA
tclusk <- list()
tsil <- list()
for (k in 2:30){
  tclusk[[k]] <- cutree(clust.c,k)
  tsil[[k]] <- silhouette(tclusk[[k]],dist=as.dist(dist.m))
  tasw[k] <- summary(silhouette(tclusk[[k]],dist=as.dist(dist.m)))$avg.width
}
ASW[2,2]=match(max(tasw[-1]),tasw[-1])+1
ASW[1,2]=max(tasw[-1])


#Average
clust.avg <- hclust(as.dist(dist.m),method="average")
tasw <- NA
tclusk <- list()
tsil <- list()
for (k in 2:30){
  tclusk[[k]] <- cutree(clust.avg,k)
  tsil[[k]] <- silhouette(tclusk[[k]],dist=as.dist(dist.m))
  tasw[k] <- summary(silhouette(tclusk[[k]],dist=as.dist(dist.m)))$avg.width
}
ASW[2,3]=match(max(tasw[-1]),tasw[-1])+1
ASW[1,3]=max(tasw[-1])

#Ward
clust.w <- hclust(as.dist(dist.m),method="ward.D2")
tasw <- NA
tclusk <- list()
tsil <- list()
for (k in 2:30){
  tclusk[[k]] <- cutree(clust.w,k)
  tsil[[k]] <- silhouette(tclusk[[k]],dist=as.dist(dist.m))
  tasw[k] <- summary(silhouette(tclusk[[k]],dist=as.dist(dist.m)))$avg.width
}
ASW[2,4]=match(max(tasw[-1]),tasw[-1])+1
ASW[1,4]=max(tasw[-1])

#Pam

pasw <- NA
pclusk <- list()
psil <- list()
for (k in 2:30){
  pclusk[[k]] <-pam(as.dist(dist.m),k)
  # Computation of silhouettes:
  psil[[k]] <- silhouette(pclusk[[k]],dist=as.dist(dist.m))
  # ASW needs to be extracted:
  pasw[k] <- summary(psil[[k]])$avg.width
}

ASW[2,5]=match(max(pasw[-1]),pasw[-1])+1
ASW[1,5]=max(pasw[-1])




clust=Mclust(wdbcc)




fviz_mclust_bic(clust, model.names = NULL, shape = 1,
                color = "model", palette = NULL, legend = NULL,
                main = "Model selection", xlab = "Number of components",
                ylab = "BIC")


#####EXERCISE 4
d=read.table("bundestag.txt",header = T,sep="")
data=d[,1:5]

kNNdistplot(data, k = 5)
abline(h=.06, col = "red", lty=2)

res <- dbscan(data, eps = 0.06, minPts = 5)
pairs(data,cex=0.3,col=res$cluster)
pairs(data,cex=0.3,col=d$ewb)

adjustedRandIndex(res$cluster,d$ewb)
