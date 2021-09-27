library(GENIE3)

obj <- load(hmrds)
regulators <- load(hmTF)

# calculate type by type
celltype <- subset(obj,idents= 'celltype')

celltype.2 <- rownames(celltype)[which(rowSums(as.matrix(celltype@assays$RNA@counts)>0)/ncol(celltype@assays$RNA@counts)>0.1)]

TF <- list()
TF <- list(AC.2,BC.2,HC.2,Rod.2,Cone.2,Muller.2,RGC.2)
retina <- list()
retina <- c(AC,BC,HC,Rod,Cone,Muller,Micro,RGC,AST)

name <- c("AC","BC","HC","Rod","Cone","Muller","RGC")
for (i in 1:length(retina)) {
  data <- retina[[i]]@assays$RNA@counts
  regulators2 <- TF[[i]][which(TF[[i]]%in%regulators)]
  regulators2 <- unique(regulators2)
  target <- TF[[i]]
  target=unique(target)
  weightMat <- matrix()
  weightMat <- GENIE3(as.matrix(data),regulators=regulators2,target=target,nCores=4,nTrees = 1000)
  
  linklist <- getLinkList(weightMat) #,reportMax = 400
  linklist$species <- rep("species",nrow(linklist))# 
  linklist$type <- rep(name[i],nrow(linklist)) 
  linklist[is.na(linklist)] <- 0
  retina[[i]] <- linklist
}
b <- data.frame()
for (i in 1:length(retina)) {
  a <- as.data.frame(retina[[i]])
  b <- rbind(b,a)
}
