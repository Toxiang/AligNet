source("D:/论文复现/Alignet/Alignet/Analysis.R")      

#网络1
edges1 = matrix(c(
  "85962.HP0109", "85962.HP0136",
"85962.HP0109", "85962.HP0137",
"85962.HP0136", "85962.HP0247",
"85962.HP0136", "85962.HP0303",
"85962.HP0137", "85962.HP0247",
"85962.HP0137", "85962.HP0853",
"85962.HP0247", "85962.HP1316"
), ncol = 2, byrow = TRUE)
hpy <- read.network(edges1, mode = "edges")

#网络2
edges2 = matrix(c(
"DBP2_YEAST", "RL2A_YEAST",
"HAS1_YEAST", "MAK5_YEAST",
"NOP10_YEAST", "DBP2_YEAST",
"NOP10_YEAST", "HAS1_YEAST",
"NOP10_YEAST", "MAK5_YEAST",
"NOP10_YEAST", "RL2A_YEAST",
"TSA1_YEAST", "HSP7F_YEAST",
"TSA1_YEAST", "TSA2_YEAST"
), ncol = 2, byrow = TRUE)
sce <- read.network(edges2, mode = "edges")

#论文中例子
#网络一 dm
e <- c(

    "dm11644", "dm10450",
    "dm247", "dm6389",
    "dm247", "dm11454",
    "dm6389", "dm11070",
    "dm6389", "dm11454",
    "dm11070", "dm2171",
    "dm2171", "dm11644",
    "dm11454", "dm11644",
    "dm10450", "dm8158"
)
edges_dm = matrix(e, ncol = 2, byrow = TRUE)
net_dm <- read.network(edges_dm, mode = "edges")
ver <- get.vertex.attribute(net_dm)$name
sim_dm <- dme.dme[rownames=ver,colnames=ver]
                    
#网络2 --hs
e2 <- c(
    "hs59", "hs3857",
    "hs3857", "hs553",
    "hs3857", "hs6992",
    "hs3857", "hs12566",
    "hs3857", "hs5433",
    "hs553", "hs6992",
    "hs553", "hs12566",
    "hs553", "hs5433",
    "hs6992", "hs12566",
    "hs6992", "hs5433",
    "hs6992", "hs5638",
    "hs12566", "hs5433",
    "hs12566", "hs5638",
    "hs399", "hs3857",
    "hs5433", "hs5638",
    "hs5638", "hs3857",
    "hs12206", "hs12566"
    )
edges_hs <- matrix(e2, ncol = 2, byrow = TRUE)
net_hs <- read.network(edges_hs, mode = "edges")
ver2 <- get.vertex.attribute(net_hs)$name
sim_hs <- hsa_hsa[rownames = ver2, colnames = ver2]

lambda_dm <- as.numeric(quantile(sim_dm)[3])
lambda_hs <- as.numeric(quantile(sim_hs)[3])
lambda_dm <- round(lambda_dm, 2)
lambda_hs <- round(lambda_hs,2)                          

#生成簇
clust_dm <- cluster.network(sigma = sim_dm, lambda = lambda_dm, k = 8)
clusters_dm <- extract.clusters(Net = net_dm, clust_dm)
for (i in 1:length(clusters_dm)) {
    print(names(clusters_dm)[i])
    plot(clusters_dm[[i]])
}

clust_hs <- cluster.network(sigma = sim_hs, lambda = lambda_hs, k = 8)
clusters_hs <- extract.clusters(net_hs, clust_hs)
for (i in 1:length(clusters_hs)) {
    print(names(clusters_hs)[i])
    plot(clusters_hs[[i]])
}

sim_all <- matrix(0, nrow = vcount(net_dm), ncol = vcount(net_hs))
dimnames(sim_all) <- list(V(net_dm)$name, V(net_hs)$name)
sim_all[1:7, 1:8] = Sim
sim_all[, 9] = sim_all[, 2]
sim_all[8,] = sim_all[3,]





#一对簇比对测试

localAligns <- align.local.all(clusters_dm, clusters_hs, sim_all, 0)

#align.local.test(clusters_dm[[1]], clusters_hs[[7]], "dm11644", "hs5433", mat = sim_all)


