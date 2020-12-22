source("D:/论文复现/Alignet/Alignet/Analysis.R")

mus_blast=read.matrix("D:/论文复现/Alignet/Alignet/MoreData/matrices/musBlast.tab",mode="col3")
cel_blast=read.matrix("D:/论文复现/Alignet/Alignet/MoreData/matrices/celBlast.tab",mode="col3")
mus_cel_blast=read.matrix("http://bioinfo.uib.es/~recerca/AligNet/Data/mus-cel-blast.tab",mode="col3")
mus_net=read.network("http://bioinfo.uib.es/~recerca/AligNet/Data/mus-net.tab",mode="edges",sep="\t")
cel_net=read.network("http://bioinfo.uib.es/~recerca/AligNet/Data/cel-net.tab",mode="edges",sep="\t")

#Some auxiliar functions
make.clusters <- function(net, blast.net) {
  ## ----matrices------------------------------------------------------------
  blast = matrix(0, nrow = vcount(net), ncol = vcount(net))

  ## ----adapt---------------------------------------------------------------
  dimnames(blast) = list(V(net)$name,V(net)$name)
  protsr = intersect(V(net)$name, rownames(blast.net))
  protsc = intersect(V(net)$name, colnames(blast.net))
  for(p in protsr){
    blast[p,protsc] = blast.net[p,protsc]
  }
  rm(protsr)
  rm(protsc)
  rm(blast.net)

  ## ----clusters------------------------------------------------------------
  blast = blast/max(blast)
  sigma = (as.matrix(blast) + compute.matrix(net))/2
  rm(blast)
  q3 = fivenum(sigma)[3]
  clust = cluster.network(sigma,q3,20)
  rm(sigma)
  clusters = extract.clusters(net,clust)
  rm(clust)
  names(clusters)=V(net)$name

  return(clusters)
}


# Clusters and global alignement
## ----clusters networks------------------------------------------------------------
clusters1 = make.clusters(mus_net, mus_blast)
rm(mus_blast)

clusters2 = make.clusters(cel_net, cel_blast)
rm(cel_blast)

## ----localalign----------------------------------------------------------

blast = matrix(0, nrow = vcount(mus_net), ncol = vcount(cel_net))
dimnames(blast) = list(V(mus_net)$name,V(cel_net)$name)
prots1 = intersect(V(mus_net)$name, rownames(mus_cel_blast))
prots2 = intersect(V(cel_net)$name, colnames(mus_cel_blast))
for(p1 in prots1){
  blast[p1,prots2] = mus_cel_blast[p1,prots2]
}
rm(prots1)
rm(prots2)
rm(mus_cel_blast)

blast = blast/max(blast)

multithread.enable = Sys.getenv('ALIGNET_LOCAL_MULTITHREADED') == '1'
numCores =
  if (multithread.enable) {
    library(parallel)
    detectCores()-1
  } else {
    1
  }
localAligns = align.local.all(clusters1,clusters2,blast,0,cores=numCores,1-blast)

rm(clusters1)
rm(clusters2)

## ----globalalign---------------------------------------------------------
global = align.global(localAligns,blast)

# Save alignement and EC and FC scores
#save(global, file='global.Rdata')
write.table(global[[2]], file='net1-net2-alignment.tab', col.names=FALSE, sep='\t')
data(go)

## ----scores--------------------------------------------------------------
ec = EC.score(global[[2]], mus_net, cel_net)
fc = FC.score(global[[2]], go)

ec
fc