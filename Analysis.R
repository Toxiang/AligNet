#'读取矩阵
#'
#'读取几种格式的矩阵
#'@param 矩阵或data.frame的fileName目录（对于模式col3）
#'@param mode table, col3, RData
#'@param sym如果矩阵必须对称，则为true。在mode table和RData中，如果矩阵非对称，则M [i，j]的值为（M [i，j] + M [j，i]）/ 2。
#'在3col模式下，M [i，j]的值将是文件中最后出现的值
#'@param def the default value in mode 3col
#'@return matrix
library(igraph)
library(Matrix)
library(data.table)
library(plot3D)
library(clue)
library(plyr)
library(parallel)
library(lpSolveAPI)



read.matrix <- function(fileName, mode = "table", sym = "TRUE", def = 0) {
    switch(mode,
           table = return(read.matrix.table(fileName, sym)),
           col3 = return(read.matrix.col3(fileName, def = def)),
           RData = return(read.matrix.RData(fileName, sym)),
           stop("Enter a valid mode"))
}

#对于table形式的文件使用as.matrix读取
read.matrix.table <- function(fileName, sym = "FALSE") {
    mat <- as.matrix(read.table(fileName))
    if (sym) {
        return((mat + t(mat)) / 2)
    }
    else
        return(mat)
    }
#对于RData形式的文件使用load读取
read.matrix.RData <- function(fileName, sym = "FALSE") {
    load(fileName)
    mat <- get(ls()[ls() != "fileName"]) #使用 ls()函数列出所有变量。ls()[ls()!="Sim1.RData"] [out]:"Sim1"
    if (sym) {
        return((mat + t(mat)) / 2)
    }
    return(mat)
}

read.matrix.col3 <- function(fileName, def = 0) {
    if (class(fileName) == "character") {
        MMM <-  read.table(fileName)
    }
    else {
        MMM <- fileName
    }
    colnames(MMM) <- c("v1","v2","v3")
    mat2 = with(MMM, sparseMatrix(
                i = as.numeric(v1), j = as.numeric(v2), x = v3,
                dimnames = list(levels(v1),levels(v2))
    ))
    return(mat2)
}


#'读网络
#'
#'从文件中读取网络
#'@param filename网络目录或data.frame
#'@param mode tab :DIP或BIOGRID中类型为tab或mitab的文件；
#'@param edges: 具有属性的边列表（不必要）
#'@param db 数据库以选择模式选项卡的蛋白质（uniprot，DIP，refseq，BIOGRID，gene）
#'@param cols 必须使用模式边缘选择哪些列，您可以使用“ all”选择所有列
#'@param sep 边缘模式下列之间的分隔符
#'@param int.type 布尔值，确定是否必须在选项卡模式下选择交互类型。

read.network <- function(fileName, mode = "tab", db = "uniprot", cols = "all"
                         , sep = " ", int.type = "TRUE") {
    switch(mode,
           tab = simplify(return(
           read.network.tab(fileName, db, int.type = int.type)
           )),
           edges = simplify(return(
           read.network.edges(fileName, cols, sep = sep)
           )),
           stop("error input"))
}

#'network tab
read.network.tab <- function(fileName, db = "uniprot", int.type = "TRUE") {
    if (class(fileName) == "character") {
        tab <- read.table(
            fileName,sep = "\t",row.names = NULL,stringsAsFactors = FALSE
       )
    }
    else {
        tab <- fileName
    }
    tab2 <- tab[, c[1, 2]]
    get.init.db <- function(i) {
        switch(db,
            refseq = {
                x=splitrefseq(tab2[i,1])
            },
            uniprot = {
            x = splituniprot(tab2[i, 1])
            },
            DIP = {
            x = splitDIP(tab2[i, 1])
            },
            BIOGRID = {
            x = splitBIOGRID(tab2[i, 1])
            },
            gene = {
            x = splitgene(tab2[i, 1])
            },
            stop("Enter valid database")
       )
        switch(db,
             refseq = {
             y = splitrefseq(tab2[i, 2])
             },
             uniprot = {
             y = splituniprot(tab2[i, 2])
             },
             DIP = {
             y = splitDIP(tab2[i, 2])
             },
             BIOGRID = {
             y = splitBIOGRID(tab2[i, 2])
             },
             gene = {
             y = splitgene(tab2[i, 2])
             },
             stop("Enter valid database")
             )
        return(list(x,y))
    }
    tab3 <- list()
    type <- c()
    for (i in l:dim(tab2)[1]) {
        x <- get.init.db(i)
        if (!is.na(x[[1]]) && !is.na(x[[2]])) {
            tab3[[i]] <- x[[1]]
            tab3[[i + dim(tab2)[1]]] <- x[[2]]
            type[i] = tab[i,12]
        }
        else {
            tab3[[i]] <- NA
            tab3[[i + dim(tab2)[1]]] <- NA
        }
    }
    tab3 <- matrix(unlist(tab3)[!is.na(tab3)], ncol = 2, byrow = FALSE)
    Net <- graph.data.frame(tab3, directed = FALSE)
    if (int.type) {
        type <- type[!is.na(type)]
        Net <- set.edge.attribute(Net,name = "Interaction type",value = type)
    }
    return(Net)
}

#' split refseq
#' @keywords internal
splitrefseq <- function(x) {
    x1 <- strsplit(x, split = "refseq:", fixed = TRUE)[[1]]
    if (length(x1) == 2) {
        return(strsplit(x1[2], split = "|", fixed = TRUE)[[1]][1])
    }
    return(NA)
}

#' split uniprot
#' @keywords internal
splituniprot <- function(x) {
    x1 <- strsplit(x, split = "uniprotkb:", fixed = TRUE)[[1]]
    if (length(x1) == 2) {
        return(x1[2])
    }
    return(NA)
}

#' split DIP
#' @keywords internal
splitDIP <- function(x) {
    x1 <- strsplit(x, split = "DIP-", fixed = TRUE)[[1]]
    if (length(x1) == 2) {
        return(paste("DIP-", strsplit(x1[2], split = "|", fixed = TRUE)[[1]][1], sep = ""))

    }
    return(NA)
}

#' split BIOGRID
#' @keywords internal
splitBIOGRID <- function(x) {
    x1 <- strsplit(x, split = "BIOGRID:", fixed = TRUE)[[1]]
    if (length(x1) == 2) {
        return(x1[2])
    }
    return(NA)
}

#' split gene
#' @keywords internal
splitgene <- function(x) {
    x1 <- strsplit(x, split = "locuslink:", fixed = TRUE)[[1]]
    if (length(x1) == 2) {
        return(strsplit(x1[2], split = "|", fixed = TRUE)[[1]][1])

    }
    return(NA)
}


read.network.edges <- function(fileName, cols, sep) {
    if (class(fileName) == "character") {
        tab <- read.table(
            fileName, sep = sep, header = TRUE,
            stringsAsFactors=FALSE
        )
    }
    else {
        tab <- fileName
    }
    if (cols[1] == "all") {
        return(graph.data.frame(tab,directed=FALSE))
    }
    return(graph.data.frame(tab[, cols], directed = FALSE))
}



random.k <- function(net,mn,mx) {
    options(digits = 1)
    l <- length(choose_k)
    x <- runif(1, 1, l-1) + 1
    eee <- get.edgelist(net)
    if(eee[choose_k[x]] != eee[choose_k[x - 1]]) {
        x <- x-1
    }
    return(x)
}

split.net <- function(net,k, size) {
    i <- 0
    #k <- 1
    eee <- get.edgelist(net)
    l <- length(eee)/2
    sub_edges <- c()
    before_k <- list()
    while (i < size) {
        i <- i + 1
        sub_edges <- append(sub_edges, eee[k,])
        x <- eee[k, 2]
        for (j in k:l) {
            if (eee[j, 1] == x) {
                choose_k <<- append(choose_k,k)
                k <- j
                break
            }
        }
    }
    print(choose_k)
    return(sub_edges)
}

sub.all <- function(num, net) {
    sub <- c()
    subs <- c()
    k <- 1
    for (i in 1:num) {
        sub <- split.net(net, k, 6)
        subs <- append(subs, sub)
        k<-random.k(net,1,10000)
    }
    return(subs)
}
#sub <- matrix(c(sub.all(8, dme)), ncol = 2, byrow = TRUE)
#net <- read.network(sub, mode = "edges")
#net <- minimum.spanning.tree(net)
#plot(net)

#' 计算矩阵
#' 计算几个矩阵
#' @param net1是一个igraph实例化对象
#' @param net2 an igraph object for modes \code{BLAST} and \code{FC} if Net2
#' is not NULL, 计算Net1和Net2中节点之间的矩阵。
#' @param 输入您要计算的矩阵(BLAST, DSD, Distance,FC, Degree)
#' @param mode for type \code{BLAST} : pident or bitscore.
#' For \code{DSD} and \code{Distance} distance,  similarity or similarity by component组件之间的相似性或相似性。
#' For \code{FC} 本体列表
#' @param byComp for \code{DSD} and \code{Distance} if the similarity or
#' distance must be normalized by components 相似性或距离必须通过组件归一化
#' @param database for \code{BLAST} and \code{FC} the database which
#' proteins belongs.
#' @param database2 for \code{BLAST} and \code{FC} the database which
#' proteins belongs.
#' @param normalized if true矩阵将被归一化
#' @return The matrix
compute.matrix <- function(net1, net2 = NULL, type = "Distance", mode = "Similarity",
                           database = NULL, database2 = NULL, byComp = TRUE, normalized = TRUE) {
    switch(
       type,
       BLAST = return(
        compute.matrix.Blast(net1, net2, mode, database, database2, normalized)
        ),
        DSD = return(compute.matrix.DSD(net1, mode, byComp, normalized)),
        Distance = return(compute.matrix.Distance(net1, mode, byComp, normalized)),
        FC = return(compute.matrix.FC(net1, net2, mode)),
        Degree = return(compute.matrix.Degree(net1, net2)),
        stop("error input")
        
      )
}

#DSD
#'@keyword interal
compute.matrix.DSD <- function(net, mode = "Similarity", byComp = TRUE, normalized = TRUE) {
    path <- paste(system.file(package = AligNet), "DSDmain.py", sep = "/")
    n <- length(V(net))
    prots <- V(net)$name
    DSD <- matrix(Inf, nrow = n, ncol = n)
    dimnames(DSD) <- list(prots, prots)
    cc <- decompose.graph(net)
    for (net in cc) {
        tmp <- tempfile()
        tmp2 <- tempfile()
        write.table(
            get.edgelist(net), quote = FALSE, file = tmp, row.names = FALSE, col.names = FALSE
       )
        command <- paste("python", path, "-m 1 -o", tmp2, tmp)
        response <- system(command, intern = T)
        table <- as.matrix(read.table(paste(tmp2, "DSD1", sep = ".")))
        diam <- max(table) + 1
        if (byComp) {
            if (mode == "Similarity") {
                DSD[rownames(table), colnames(table)] <- (max(table) + 1 - table) / (max(table) + 1)
            }
            else {
                DSD[rownames(table), colnames(table)] <- table / max(table)
            }
        }
        else {
            DSD[rownames(table), colnames(table)] <- table
        }
    }
    mmm <- max(DSD[DSD < Inf])
    if (!byComp) {
        if (mode == "Similarity") {
            DSD <- (mmm+1-DSD)/(mmm+1)
        }
        else {
            if (normalized) {
                DSD <- DSD / mmm
            }
        }
        if (mode == "Similarity") {
            DSD[DSD ==Inf] <- 0
        }
        DSD[DSD == -Inf] <- 0
        return (DSD)
    }
}

#'Distance
#'@keywords internal
compute.matrix.Distance <- function(net, mode = "Similarity",
                                    byComp = TRUE, normalized = TRUE) {
    n <- length(V(net))
    prots <- V(net)$name #网络名
    if (!byComp) {  #归一化
        dist <- shortest.paths(net)
        mmm <- max(dist[dist < Inf])
        if (mode == "Similarity") {
            dist <- (mmm + 1 - dist) / (mmm + 1)
            dist[dist == -nf] <- 0
            return(dist)
        }
        if (normalized) {
            dist <- dist / mmm
        }
        return(dist)
    }
    if (mode == "Similarity") {
        dist <- matrix(0, nrow = n, ncol = n)
    }
    else {
        dist <- matrix(Inf, nrow = n, ncol = n)
    }
    dimnames(dist) <- list(prots, prots) #为dist行列命名
    cc <- decompose.graph(net)#igraph对象
    for (nnn in cc) {#
        dist2 <- shortest.paths(nnn)#计算出每个节点到其他节点的最短距离
        mmm2 <- max(dist2[dist2 < Inf])#任意两个节点之间的最长距离
        if (mode == "Similarity") {
            dist2 <- (mmm2 + 1 - dist2) / (mmm2 + 1)#距离计算D(G)+1-dg(u,v)/D(g)+1
            dist2[dist2 == -Inf] <- 0
            dist2[dist2 == Inf] <- 0
            dist[V(nnn)$name, V(nnn)$name] <- dist2
        }
        if (normalized) {
            dist2 <- dist2 / max(dist2[dist2 < Inf]) #归一化
        }
        dist[V(nnn)$name, V(nnn)$name] <- dist2
    }
    return(dist)
}

#'FC
#'@keywords interal
compute.matrix.FC <- function(net1, net2 = NULL, gos) {
    onenet <- FALSE
    if (is.null(net2)) {
        net2 <- net1
        onenet <- TRUE
    }
    prots1 <- V(net1)$name
    prots2 <- V(net2)$name
    if (onenet) {
        FSim <- diag(1, nrow = length(prots1), ncol = length(prots2))
    }
    else {
        FSim <- matrix(0,nrow = length(prots1),ncol = length(prots2))
    }
    diameter(FSim) <- list(prots1, prots2)
    if (onenet) {
        for (i in l:(length(prots1) - 1)) {
            for (j in (i + 1):length(prots1)) {
                fc <- lenth(intersect(gos[[prots1[i]]],
                                       , gos[[prots1[j]]])) / length(union(gos[[prots1[i]]],
                                                                            gos[[prots1[j]]]))
                FSim[i, j] <- fc
                FSim[j, i] <- fc
            }
        }
    }
    else {
        for (i in l:length(prots1)) {
            for (j in l:length(prots2)) {
                fc <- length(intersect(gos[[prots1[i]]],
                                        gos[[prots2[j]]])) / length(union(gos[[prots1[i]]],
                                                                          gos[[prots2[j]]]))
                FSim[i,j] <- fc
            }
        }
    }
    return(FSim)
}

#'Degree
#'@keywords internal
compute.matrix.Degree <- function(net1, net2 = NULL) {
    #如果net2没有输入，则将net1赋值给2
    if (is.null(net2)) {
        net2 <- net1
    }
    deg <- degree(net1)
    deg2 <- degree(net2)
    return(matrix(unlist(lapply(deg, function(i)
        lapply(deg2, function(j)
            abs(i - j)))), nrow = length(deg), byrow = TRUE))

}

#'IC
#'@keywords internal
computeIC <- function(net1, net2) {
    degree1 <- degree(net1)
    degree2 <- degree(net2)
    neigh1 <- neighborhood(net1, order = 1)
    neigh2 <- neighborhood(net2, order = 1)
    lens1 <- unlist(lapply(neigh1, length))
    lens2 <- unlist(lapply(neigh2, length))
    sums1 <- unlist(lapply(neigh1, function(i)
        sum(1 / lens1[i])))
    sums2 <- unlist(lapply(neigh2, function(i)
        sum(1 / lens2[i])))
    mat <- matrix(0, nrow = vcount(net1), ncol = vcount(net2))
    maxneig <- max(degree1, degree2)
    for (i in 1:vcount(net1)) {
        for (j in 1:vcount(net2)) {
            mat[i, j] <- min(sums1[i], sums2[i]) / maxneig
        }
    }
    dimnames(mat) <- list(V(net1)$name, V(net2)$name)
    return(mat)
}

#'No funciona
#'@keywords internal
compute.matrix.Blast <- function(Net1, Net2, mode, database, database2, normalized) {
    tmp <- tempfile()
    buildBlast(Net1, Net2, mode, database, database2, tmp)
    BlastM <- read.matrix.table(tmp)
    if (normalized) {
        return(BlastM / max(BlastM))
    }
    return(BlastM)
}

#'No funciona
#'@keywords internal
buildBlast <- function(Net1, Net2, mode, database, database2, tmp) {
    path <-
    paste(system.file(package = "AligNet"), "Blast.jar", sep = "/")
    tmp1 <- tempfile()
    write(V(Net1)$name, file = tmp1)
    tmp2 <- tempfile()
    if (is.null(Net2)) {
        tmp2 <- tmp1
    }
    else {
        write(V(Net2)$name, file = tmp2)
    }
    if (is.null(database2)) {
        database <- database2
    }
    command <- paste(
    "java", path, "-prot1", tmp1, "-prot2", tmp2, "-db1",
    database, "-db2", database2, "-m", mode, "-outfmt mat -o", tmp
  )
    response <- system(command, intern = T)

}



#'Cluster Network
#'
#'根据相似矩阵sigma计算聚类矩阵，其中所有节点都有更多
#'相似性比λ大 ，并且簇的大小小于k。
#'
#'@param sigma 相似矩阵
#'@param lambda 相似度阈值 第三个四分位数
#'@param k 阈值大小
#'@return 聚类矩阵，如果节点i属于节点j的集群，则M [i，j]为TRUE 否则为FALSE
cluster.network <- function(sigma, lambda = 0, k = dim(sigma)[1]) {
    if (is.data.frame(sigma)) {
        sigma <- read.matrix.col3(sigma, def = 0)

    }
    n <- dim(sigma)[1]#lambda为相似度阈值
    clustmatrix <- matrix(c(as.matrix(sigma) > lambda), nrow = n, byrow = TRUE)#true/false矩阵
    sums <- apply(clustmatrix, 2, sum)#将cm矩阵按照列相加
    ind <- which(sums > k)
    if (length(ind) > 0) { #如果有大于k个节点，则按照s大小进行排序
        lapply(ind, function(i)
            clustmatrix[sort(as.matrix(sigma)[i,], index.return = TRUE)$ix[1:(n - k)], i] <<- FALSE)
        }
    dimnames(clustmatrix) <- dimnames(sigma)
    clustmatrix <- as(clustmatrix, "sparseMatrix")
    return(clustmatrix)
}

#'Extract clusters
#'
#'Compute the subnetworks of \code{Net} from a cluster matrix
#'@param Net an igraph object
#'@param ClustMat a cluster matrix (output of \code{cluster.matrix})
#'@return Clusters in igraph format

extract.clusters <- function(Net, ClustMat) {
    prots <- rownames(ClustMat)
    clusts <- lapply(1:dim(ClustMat)[2], function(i)
        induced.subgraph(Net, prots[ClustMat[, i] == 1]))
    names(clusts) = prots
    return(clusts)

}
#'Display Clusters
#'
#'Given a cluster matrix, see \code{cluster.network} and a network
#'用以下颜色显示位置（i，j）的聚类矩阵：
#'
#'- 如果蛋白质和蛋白质i不属于蛋白质j的簇并且蛋白质在网络中不相互作用，则为黄色
#'
#'-  如果蛋白质和蛋白质i属于蛋白质j的簇，则为黑色，但是蛋白质在网络中不相互作用
#'
#'- 如果蛋白质和蛋白质i不属于蛋白质j的簇但蛋白质在网络中相互作用，则为红色
#'
#'-  如果蛋白质和蛋白质i属于蛋白质j的簇并且蛋白质在网络中相互作用，则为绿色
#'
#'@param clust a matrix which is the output1 of \code{cluster.network}
#'@param cols a list of 4 colors if you want to change the default colors
#'@param zoom an integer to define the size of plot or NA, to plot all clusters
#'@param type 0 = "yellow", 1 = "black", 2="red", 3 = "green"
#'@param col color to use with zoom
#'@param Net an igraph object
#'@param ... Additional plotting parameters

display.clusters <- function(clust, Net, zoom = NA, type = 1
                             , cols = c("yellow", "black", "red", "green"),
                             col = cols[type + 1,], ...) {
    clust <- as.matrix(clust)
    if (is.na(zoom)) {
        mmm <- clust + 2 * as.matrix(get.adjacency(Net))
        dimnames(mmm) <- list(1:dim(mmm)[1], 1:dim(mmm)[2])
        # par(oma=c(1,1,1,1))
        image2D(
            mmm,
            col = cols,
        axex = F,
        xlab = NA,
        ylab = NA,
        colkey = FALSE,
        ...
        )
        axis(
        1, at = seq(0, 1, by = 1 / (dim(clust)[2] - 1)),
        labels = colnames(clust), las = 2
        , pos = -1 / (2 * (dim(clust)[1] - 1))

    )
        axis(
        2, at = seq(0, 1, by = 1 / (dim(clust)[1] - 1)),
        labels = rownames(clust), las = 2, pos = -1 / (2 * (dim(clust)[2] - 1))
        )
        legend(
        x = 0.2, y = 1.4, legend = 0:3, fill = cols, horiz = TRUE, bty = "n", xpd = TRUE)
    }
    else {
        mat <- clust + 2 * as.matrix(get.adjacency(Net))
        k <- floor(length(V(Net)) / zoom)
        conts <- unlist(lapply(1:zoom, function(i)
            lapply(1:zoom, function(j)
                conts.type(mat, (i - 1) * k + 1, i * k, (j - 1) * k + 1, j * k, type))))
            mat2 <- matrix(conts, nrow = zoom, byrow = TRUE)
        print("Min")
        print(min(mat2))
        print("Max")
        print(max(mat2))
        print(max(mat2))
        image2D(
        mat2, col = ramp.col(col = c("white", col), n = k * k),
        axes = F,
        xlab = NA,
        ylab = NA,
        ...
      )
    }
}
#' Cont type
#' @keywords internal
cont.type <- function(mat, init1, fin1, init2, fin2, type) {
    return(length(which(mat[init1:fin1, init2:fin2] == type)))
}


##########################################################
###########################################################
############################################################

#'Hungraian 匈牙利算法
#'@keywords interal
HungarianFinal <- function(mat, maxim = TRUE) {
    nodesnet1 <- rownames(mat)
    nodesnet2 <- colnames(mat)
    if (dim(mat)[1] > dim(mat)[2]) {
        alin <- solve_LSAP(t(mat), maximum = maxim)
        alin <- cbind(nodesnet1[alin],nodesnet2[seq_along(alin)])
    }
    else {
        alin <- solve_LSAP(mat, maximum = maxim)
        alin <- cbind(nodesnet1[seq_along(alin)],nodesnet2[alin])
    }
    return(alin)
}

#'EC分数
#'
#' Given two networks, \code{net1}=(\eqn{V_1,E_1}) , \code{net2}=(\eqn{V_2,E_2})
#' , and an alignment \eqn{g:V_1 \rightarrow V_2}, then the edge correctnees is
#' defined as: \deqn{\frac{|\{(u,v)\in E_1 \: : \;
#' (g(u),g(v))\in E_2\}|}{|E_1|}}
#'@param alin  alignment, with the format of \code{alin.local} or
#'\code{alin.global}
#'@param net1  network
#'@param net2  network
#'@return EC score
EC.score <- function(alin, net1, net2) {
    alini <- function(x, y) {
        matrix(c(as.character(alin[x]),
        as.character(alin[y])),nrow = 1)
    }
    E1 <- get.edgelist(net1)
    if (dim(E1)[1] == 0) {
        return(0)
    }
    eds <- t(mapply(alini, E1[, 1], E1[, 2]))
    nas <- unique(c(which(is.na(eds[, 2])), which(is.na(eds[, 1]))))
    if (length(nas) > 0) {
        eds <- eds[-nas,]
    }
    if (is.null(dim(eds))) {
        if (length(eds) == 2) {
            eds <- cbind(eds[1],eds[2])
        }
    }
    if (dim(eds)[1] == 0) {
        return(0)
    }
    else {
        if (dim(eds)[2] == 1) {
            return(0)
        }
    }
    Gnet3 <- graph.edgelist(eds, directed = FALSE)
    ggg4 <- graph.intersection(Gnet3, net2, byname = TRUE)
    Gnet1 <- induced.subgraph(net1, vids = names(alin))
    Gnet2 <- induced.subgraph(net2, vids = alin)
    if (min(length(E(Gnet1)), length(E(Gnet2))) == 0) {
        return(0)
    }
    return(length(E(ggg4))/min(length(E(Gnet1)),length(E(Gnet2))))
}

#'Functional Coherence Score
#'
#'Given an alignment \eqn{g:V_1 \rightarrow V_2}, and an ontologies \eqn{gos},
#'then the functional coherence score is defined as:
#'\deqn{\frac{1}{|V_1|} \sum_{v \in V_1} \frac{|gos(v) \cap gos(g(v))|}
#'{|gos(v)\cup gos(g(v))|}}
#'@param alin alignment, with the format of \code{alin.local} or
#'\code{alin.global}
#'@param gos a list of ontologies
#'@return FC score
FC.score <- function(alin, gos) {
    aa1 <- names(alin)
    aa2 <- alin

    fcs <- unlist(lapply(1:length(aa1), function(i)
        length(intersect(gos[[aa1[i]]], gos[[aa2[i]]])) /
            length(union(gos[[aa1[i]]], gos[[aa2[i]]]))))
    return(mean(fcs,na.rm=TRUE))
}

#'Local alignment
#'
#'Compute the local alignment between the networks \code{net1} and
#'\code{net2} with centers \code{p1} and \code{p2}
#'@param net1 network
#'@param net2 network
#'@param p1 center of net1
#'@param p2 center of net2
#'@param compute.ec compute ecscore (TRUE/FALSE)
#'@param mat 非相似矩阵
#'@return alignment


#输入两个聚类矩阵 将他们进行匹配
align.local <- function(net1, net2,p1,p2, compute.ec = FALSE, mat = NULL) {
    #计算两个网络的度之差
    mat1 <- compute.matrix.Degree(net1, net2)
    #mat作用未知
    if (is.null(mat)) {
        mat <- mat1
    }
    else {
        mat <- mat1+1-mat[V(net1)$name,V(net2)$name]
    }
    dimnames(mat) <- list(V(net1)$name, V(net2)$name)
    #neigh1 <- neighborhood(graph = net1, order = 1)
    #找到簇的邻居节点
    neigh1 <- neighborhood(graph = net1, order = 1)
    neigh1 <-
      lapply(1:length(V(net1)),function(i)
        V(net1)$name[neigh1[[i]]])
    names(neigh1) <- V(net1)$name
    neigh2 <- neighborhood(graph = net2,order = 1)
    neigh2 <-
      lapply(1:length(V(net2)),function(i)
        V(net2)$name[neigh2[[i]]])
    names(neigh2) <- V(net2)$name
    #################
    assign <- p2
    names(assign) <- c(p1)
    #85962. HP0109
    #"DBP2_YEAST"

    #暂定 分配节点列表
    completes <- c()
    incomplets <- c(p1)
    assignats <- c(p2)
    i <- 1
    while (length(assign) < length(neigh1)) {
        #print(paste("第", i, "次循环"))
        i <- i + 1
        
      q1 <- incomplets[1]
      q2 <- assign[q1]

      #print("q1:")
        #print(q1)
        #print("q2:")
        #print(q2)
      #setdiff 找neigh1[[q1]]中不同于names(assign)的元素
      n1 <- setdiff(neigh1[[q1]],names(assign))
      n2 <- setdiff(neigh2[[q2]], assign)
        #print("n1:")
        #print(n1)
        #print("n2")
        #print(n2)
      assign2 <- c()

      if (length(n2) > 1 && length(n1) > 1) {
        mat2 <- mat[n1,n2]
        dimnames(mat2) <- list(n1,n2)

      }
      else{
        mat2 <- c()
        if (length(n1) * length(n2) > 0) {
          if (length(n1) == 1 && length(n2) == 1) {
            assign2 <- n2
            names(assign2) <- n1
            assign <- c(assign,assign2)

          }
          else{
            if (length(n1) == 1) {
              assign2 <- names(which.min(mat[n1, n2]))
                #print("assign2:")
                #print(assign2)

              names(assign2) <- n1
                #print(paste("name:",names(assign2)))
              assign <- c(assign,assign2)
            }
            if (length(n2) == 1) {
              assign2 <- n2
              names(assign2) <- names(which.min(mat[n1,n2]))
              assign <- c(assign,assign2)
            }
          }
        }
      }

      if (!is.null(dim(mat2))) {
        if (dim(mat2)[1] > dim(mat2)[2]) {
          hung <- HungarianFinal(as.matrix(t(mat2)), maxim = FALSE)
            #print("hung:")
            #print(hung)
          assign2 <- hung[,1]

          names(assign2) <- hung[,2]
          inds <- sort(unlist(lapply(1:dim(hung)[1],
                                     function(i)
                                       t(mat2)[hung[i,1],hung[i,2]])),
                       decreasing = TRUE,index.return = TRUE)$ix
          assign2 <- assign2[inds]
            #print("assign2:")
            #print(assign2)
          assign <- c(assign,assign2)
          
        }
        else{
          hung <- HungarianFinal(as.matrix(mat2), maxim = FALSE)
            #print("hung:")
            #print(hung)
          assign2 <- hung[,2]
          names(assign2) <- hung[,1]
          inds <- sort(unlist(lapply(1:dim(hung)[1],
                                     function(i)
                                       (mat2)[hung[i,1],hung[i,2]])),
                       decreasing = TRUE,index.return = TRUE)$ix
          assign2 <- assign2[inds]
            #print("assign2:")
            #print(assign2)
          assign <- c(assign,assign2)
        }
      }
        #print("此时Assign")
        #print(assign)

      incomplets <- incomplets[ - 1] #保留除第一个元素外的其他元素
      incomplets <- c(incomplets,names(assign2))
        #print("incompletes:")
        #print(incomplets)

      if (length(incomplets) == 0) {
        break
      }
    }
    #print("Assign")
    #print(assign)
    if (length(assign) == 0) {
      assign <- p2
      names(assign) <- p1
    }
    if (compute.ec) {
      return(list(align = assign,ec = EC.score(assign,net1,net2)))
    }
    return(list(align = assign))
}
#'alin_aux
#'alin_aux
#'@keywords internal
alin_aux <- function(p1, mat, ll, clust1, clust2, mat2 = NULL) {
    pp2 <- intersect(names(which(mat[p1,] > ll)), names(clust2))
    protsc1 <- V(clust1[[p1]])$name
    if (length(pp2) == 0) {
        return(list())
    }
    return(lapply(pp2,
                function(p2)
                    align.local(clust1[[p1]], clust2[[p2]], p1, p2, mat = mat2)))

}
#'All local alignments
#'Compute a list of local alignments between \code{clust1} and
#'\code{clust2}. The function compute all the local alignments of all pairs
#'of clusters whose centers have a similarity greather than \code{ll}
#'@param mat a similarity matrix
#'@param threshold a threshold
#'@param cores number of cores
#'@param dismat a disimilarity matrix to use in the local aligments
#'@return a list of alignments

align.local.all <-
    function(clust1, clust2, mat, threshold, cores = 1, dismat = NULL) {
        prots1 <- intersect(names(clust1), rownames(mat))
        aligns <- mclapply(prots1,
                       function(i)
                           alin_aux(i, mat, threshold, clust1, clust2, dismat),
                       mc.cores = cores)
        return(aligns)
    }

#'Size score
#'
#'Compute the size score of the list of local alignments. Given an alignment
#'\eqn{g:V_1 \rightarrow V_2}, the size score of \eqn{g} is defined as
#'\eqn{|V_1|}
#'@param localAligns a list of local alignments
#'@return sizeScore a table with 3 columns, the first and second ones
#'represents the local alignment and the third the score of alignment

size.score.all <- function(localAligns) {
    #sim=Sim
    als <- localAligns
    #计算每个单独的align中有几个对应的匹配
    align.size <- unlist(lapply(als, function(i)
        length(i)))

    align.name1 <- unlist(lapply(als, function(i)
        names(i)[1]))
    align.name2 <- unlist(lapply(als, function(i)
        i[1]))

    return(cbind(align.name1, align.name2, align.size))
}

#'Similarity Score
#'
#'Compute the similarity score of the list of local alignments. Given an
#'alignment \eqn{g:V_1\rightarrow V_2}, and a similarity matrix \eqn{sim}, the
#'similarity score of the alignment \eqn{g} is defined as:
#'\deqn{\frac{1}{|V_1|} \sum_{v \in V_1} sim[v,g(v)]}
#'@param localAligns a list of local alignments
#'@param sim a similarity matrix
#'@return simScore a table with 3 columns, the first and second ones
#'represents the local alignment and the third the score of alignment
sim.score.all <- function(localAligns, sim) {
    als <- localAligns
    align.sim <- unlist(lapply(als, function(i)
        sim.score(i, sim)))
    align.name1 <- unlist(lapply(als, function(i)
        names(i)[1]))
    align.name2 <- unlist(lapply(als, function(i)
        i[1]))

    return(cbind(align.name1, align.name2, align.sim))
}

#'simScore.aux
#'@keywords internal
sim.score <- function(align, sim) {
    #diag()提取对角矩阵
    sco <- sum(diag(sim[names(align), align]))
    n <- length(align)
    return(sco / n)
}

#'Edges score
#'Compute the EC.score of the list of local alignments
#'@param localAligns, a list of local alignments
#'@param net1 an igraph object
#'@param net2 an igraph object
#'@return edgesScore  a table with 3 columns, the first and second ones
#'represents the local alignment and the third the score of alignment

EC.score.all <- function(localAligns, net1, net2) {
    als <- localAligns
    align.ec <- unlist(lapply(als, function(i)
        EC.score(i, net1, net2)))
    align.name <- unlist(lapply(als, function(i)
        names(i)[1]))
    align.name2 <- unlist(lapply(als, function(i)
        i[1]))
    return(cbind(align.name,align.name2,align.ec))
}

#'Select alignments
#'
#'Select a list of aligments that recover all the proteins, based on the
#'list of scores
#'@param localAligns a list of local alignments
#'@param scores a table of scores
#'@return selectAligns a list of local alignments
select.aligns <- function(localAligns, scores) {
    als <- localAligns
    als2 <- list()
    protsin <- 0
    prots <- c()
    n <- length(unique(scores[, 1]))
    while (protsin < n) {
        #找到score分数最大的簇对
        i <- which.max(scores[,2])
        al1 <- als[[i]]
        prots <- append(prots, names(al1))
        prots <- unique(prots)
        protsin <- length(prots)
        for (p1 in names(al1)) {
            if (p1 %in% scores[, 1]) {
                scores[which(scores[, 1] == p1), 2] <- - 1
            }
        }
        als2 <- append(als2,list(al1))
    }
    return(als2)
}
#'is.element2
#'@keywords internal
is.element2 <- function(i, j) {
    return(is.element(j, i))
}

#'Compute score
#'@keywords internal
compute.score <- function(als2, Sim) {
    blasts <- sim.score.all(als2, Sim)
    tamanys <- size.score.all(als2)
    score <- tamanys
    score[, 3] <- as.numeric(score[, 3]) / max(as.numeric(score[, 3]))
    if (max(as.numeric(blasts[, 3])) > 0) {
        score[, 3] <-
        as.numeric(score[, 3]) + as.numeric(blasts[, 3]) /
        max(as.numeric(blasts[, 3]))
    }
    colnames(score) <- c("V1", "V2", "V3")
    rownames(score) <- NULL
    sc3 <- floor(100 * as.numeric(score[, 3]))
    score <- as.data.frame(score)
    score$V3 <- sc3
    return(score)
}

#'Global Alignment of two Protein Interaction Network
#'
#'Return a global alignment, from the list of local alignments and the table
#'of scores. The function first calculate a list of alignments with
#'\code{selectAligns}, then found a solution of the hypergraph mathching
#'problem. And finally extend this alignment to a global alignment
#'@param localAligns a list of local alignments
#'@param Sim a similarity matrix
#'@param AllSteps a boolean to determine whether return all intermediate
#'global alignments
#'@param getGlobal a boolean to set if the return is a global alignment or a better aligment but not strictly global
#'@return Global a global alignment and, if AllSteps is TRUE all intermediate alignments
align.global <- function(localAligns,Sim, AllSteps = TRUE, getGlobal = FALSE ) {
  global <- c()
  als <-
    unlist(unlist(localAligns,recursive = FALSE),recursive = FALSE)
  scores <- compute.score(als, Sim)
  Mat <-
    with(scores, sparseMatrix(
      i = as.numeric(V1), j = as.numeric(V2),
      x = V3, dimnames = list(levels(V1), levels(V2))
    ))
  globals <- list()
  Mat <- as.matrix(Mat)
  while (max(Mat) > 0) {
    hg <- HungarianFinal(Mat)
    #得到使用匈牙利算法得到的匹配
    als2 <- get.aligns(als, hg, Mat)
    score <- compute.score(als2, Sim)
    als2 <- select.aligns(als2, score[, c(1, 3)])
    hy <- hypergraph.solve(als2)
    global <- c(global, hy)
    globals <- append(globals, list(global))
    als <- aligns.update(als, global)
    if (length(als) > 0) {
      als <- remove.global(als, global)
      scores <- compute.score(als, Sim)
    }
    else {
      scores <- matrix(1:3, nrow = 1, ncol = 3)[ - 1,]
    }
    if (dim(scores)[1] > 0) {
      Mat <- with(scores, sparseMatrix(
        i = as.numeric(V1),
        j = as.numeric(V2),
        x = V3,
        dimnames = list(levels(V1), levels(V2))
      ))
      Mat <- as.matrix(Mat)
    }
    else {
      Mat <- matrix(0, nrow = 1, ncol = 1)
    }
  }
  if (getGlobal) {
    global2 <- align.end(localAligns, global)
  } else {
    global2 = global
  }
  if (AllSteps){
  return(list(globals, global2))
  }
  else{
    return(global2)
  }
}

#'Update aligns
#'@keywords internal
aligns.update <- function(als,global) {
  if (length(global) == 0) {
    return(als)
  }
  prots1 <- names(global)
  prots2 <- global
  als2 <- list()
  for (al in als) {
    if (!al[1] %in% prots2) {
      if (!names(al)[1] %in% prots1) {
        als2 <- append(als2,list(al))
      }
    }
  }
  return(als2)
}

#'get aligns
#'@keywords internal
get.aligns <- function(als,hg,Mat) {
  als2 <- list()
  prots1 <- hg[,1]
  prots2 <- hg[,2]
    for (al in als) {
      #判断al[1]是否在prots2中
    if (al[1] %in% prots2) {
      i <- which(al[1] == prots2)
      if (names(al)[1] == prots1[i]) {
        if (Mat[names(al[1]),al[1]] > 0) {
          als2 <- append(als2,list(al))
        }
      }
    }
  }
  return(als2)
}

#'remove global
#'@keywords internal
remove.global <- function(als2,global) {
  for (i in 1:length(als2)) {
    als2[[i]] <- als2[[i]][which(!als2[[i]] %in% global)]
    als2[[i]] <-
      als2[[i]][which(!names(als2[[i]]) %in% names(global))]
  }
  return(als2)
}

#'Solve hypergraph
#'@keywords internal
hypergraph.solve <- function(als) {
  E2 <- als
  E1 <- lapply(E2,names)
  scores <- unlist(lapply(E1, length))
  #统计节点出现的次数
  cprots1 <- count(unlist(E1))
  prots1 <- cprots1[,1]
  cprots2 <- count(unlist(E2))
  prots2 <- cprots2[,1]
  vars <- length(E1)
  constr1 <- length(prots1)
  constr2 <- length(prots2)
  constr <- constr1 + constr2
  lprec <- make.lp(constr,vars)
  #检查E1向量中是否含有i E1共有三列向量
  constr11 <-
    lapply(prots1, function(i)
      as.numeric(unlist(lapply(E1,is.element2,i))))
  constr22 <-
    lapply(prots2, function(i)
      as.numeric(unlist(lapply(E2,is.element2,i))))

  ##Set constraints
  fcon1 <- function(i) {
    s <- sum(constr11[[i]])
    if (s > 0) {
      set.row(lprec,i,xt = rep(1,s),indices = which(constr11[[i]] == 1))
    }
  }

  aa <- lapply(1:constr1, fcon1)
  fcon2 <- function(i) {
    s <- sum(constr22[[i]])
    if (s > 0) {
      set.row(lprec,i + constr1,xt = rep(1,s),
              indices = which(constr22[[i]] == 1))
    }
  }

  bb <- lapply(1:constr2, fcon2)
  #添加Minimize
  set.objfn(lprec, unlist(scores))
    #<== instead free
  set.constr.type(lprec,c(rep("<=",constr)))
  set.rhs(lprec,c(rep(1,constr)))
  set.bounds(
    lprec,lower = rep(0,vars),upper = rep(1,vars),columns = 1:vars
  )

  set.type(lprec,1:vars,"binary")
  break.value <- min(length(prots1),length(prots2))
  lp.control(
    lprec,sense = "max",verbose = "neutral",break.at.first = TRUE
  )
  solve(lprec)
  sols <- which(get.variables(lprec) > 0)

  getprots <- function(i,E2) {
    matrix(c(names(E2[[i]]),E2[[i]]),ncol = 2)
  }
  mmm <- getprots(sols[[1]],E2)
  nsol <- length(sols)
  if (nsol > 1) {
    for (i in 2:nsol) {
      mmm <- rbind(mmm,getprots(sols[[i]],E2))
    }
  }

  global <- mmm[,2]
  names(global) <- mmm[,1]
  return(global)
}


#'Update Matrix
#'@keywords internal
matrix.update <- function(mat,global) {
  mat2 <- mat
  rows.delete <- which(rownames(mat2) %in% names(global))
  cols.delete <- which(colnames(mat2) %in% global)
  mat2 <- mat2[ - rows.delete,]
  if (is.null(dim(mat2))) {
    if (length(mat2) > 0) {
      mat2 <- matrix(mat2,nrow = 1)
      rownames(mat2) <- setdiff(rownames(mat),names(global))
      colnames(mat2) <- colnames(mat)
      mat3 <- mat2[, - cols.delete]
      mat3 <- matrix(mat3,nrow = 1)
      colnames(mat3) <- setdiff(colnames(mat2),global)
      rownames(mat3) <- rownames(mat2)
      return(mat3)

    }
    else{
      return(matrix(0, nrow = 1,ncol = 1))
    }
  }
  mat3 <- mat2[, - cols.delete]
  if (is.null(dim(mat3))) {
    if (length(mat3) > 0) {
      mat3 <- matrix(mat3,ncol = 1)
      colnames(mat3) <- setdiff(colnames(mat2),global)
      rownames(mat3) <- rownames(mat2)
    }
    else{
      mat3 <- matrix(0, nrow = 1,ncol = 1)
    }
  }
  return(mat3)
}

#'End alignment
#'@keywords internal
align.end <- function(localAligns,global) {
  als <-
    unlist(unlist(localAligns,recursive = FALSE),recursive = FALSE)
  mmm <- cbind(names(global),global)
  E2 <- lapply(seq(1,length(als),2), function(i)
    als[i][[1]])
  E1 <- lapply(E2,names)

  prots1 <- unique(unlist(E1))

  rest2 <- unlist(E2)
  rest1 <- names(rest2)
  restm <- count(matrix(c(rest1,rest2),byrow = FALSE,ncol = 2))
  restm <- data.frame(restm,row.names = 1:dim(restm)[1])
  colnames(restm) <- c("V1","V2","V3")
  mat2 <- with(restm, sparseMatrix(
    i = as.numeric(V1),
    j = as.numeric(V2),
    x = V3,
    dimnames = list(levels(V1), levels(V2))
  ))
  mat2 = as.matrix(mat2)
  mat2[mat2 > 0] <- 1

  mat2[intersect(mmm[,1],rownames(mat2)),] <- -1

  mat2[,intersect(mmm[,2],colnames(mat2))] <- -1
  mat2 <- as.matrix(mat2) + 1
  hg <- HungarianFinal(as.matrix(mat2))


  for (i in 1:dim(hg)[1]) {
    if (mat2[hg[i,1],hg[i,2]] > 0) {
      mmm <- rbind(mmm,hg[i,])
    }
  }

  global2 <- mmm[,2]
  names(global2) <- mmm[,1]
  return(global2)
}

#'Plot the alignment
#'@param net1 an igraph object
#'@param net2 an igrpah object
#'@param global an alignment
#'@param k1 the width of the new edges of the alignment
#'@param k2 the width of the old edges
#'@param edge.curved Specifies whether to draw curved
#'edges, or not. This can be a logical or a numeric vector or scalar.
#'@param ... further arguments to be passed to igraph plotting
align.plot <-
  function(net1,net2,global,k1=1, k2=1, edge.curved = 0.5, ...) {
    coms1 <- fastgreedy.community(net1)
    coms2 <- fastgreedy.community(net2)
    newedges <- cbind(names(global),global)
    net3 <- graph.data.frame(newedges,directed = FALSE)
    net4 <- graph.union(net1,net2)
    net5 <- graph.union(net3,net4)
    num.eds <- ecount(net5)
    eds <- get.edgelist(net5)
    E(net5)$weight <- k1
    E(net5)$edge.curved <- 0
    eds1 <- get.edgelist(net4)
    ids <- get.edge.ids(net5,t(eds1),FALSE)

    net5 <- set.edge.attribute(net5,name = "weight",index = ids,k2)
    net5 <- set.edge.attribute(net5,name = "edge.curved",
                               index = ids,edge.curved)

    inds1 <- sort(coms1$membership,index.return = TRUE)

    lay1 <- 1:vcount(net1)
    lay2 <- rep(0,vcount(net2))
    lay12 <- lay1[inds1$ix]
    names(lay12) <- V(net1)$name
    lay <- cbind(10,rep(0,vcount(net5)))
    for (p1 in V(net1)$name) {
      i1 <- which(V(net5)$name == p1)
      p2 <- global[p1]
      i2 <- which(V(net5)$name == p2)
      lay[i1,] <- c(0,lay12[p1])
      lay[i2,] <- c(10,lay12[p1])
    }
    print("plot")
    print(net5)
    plot(
      net5,layout = layout.norm(lay),rescale = TRUE,
      edge.curved = E(net5)$edge.curved,edge.width = E(net5)$weight,...
    )
  }

#'Search alignments
#' Search alignments in als from p1 to p2
#' @param als a list of local alignments
#' @param p1 a protein
#' @param p2 a protein
#' @return list of alignments that includes p1 and p2
search.aligns <- function(als,p1,p2) {
  clustp1 <- c()
  clustp2 <- c()
  for (al in als) {
    if (p1 == names(al)[1]) {
      if (al[p1] == p2) {
        clustp1 <- names(al)
        clustp2 <- als
      }
    }
  }
  als2 <- list()
  for (al in als) {
    al.aux <- al[intersect(names(al),clustp1)]
    if (length(al.aux) > 0) {
      als2 <- append(als2,list(al.aux))
    }
  }
  return(als2)

}

#'Local alignment plot
#'plot a local alignment and all the local alignments that
#'intersect with him
#'@param localAligns a list of local alignments
#'@param global an alignment
#'@param p1 the center of cluster1
#'@param p2 the center of cluster2
#'@param net1 the first network
#'@param net2 the second network
#'@param ... further arguments to be passed to igraph plotting
align.local.plot <- function(localAligns,global,p1,p2,net1,net2,...) {
  als <- unlist(unlist(localAligns,recursive = FALSE),recursive = FALSE)
  als <- search.aligns(als,p1,p2)
  alp1p2 <- NULL
  for (al in als) {
    if (p1 == names(al)[1]) {
      if (al[p1] == p2) {
        alp1p2 <- al
      }
    }
  }
  if (is.null(alp1p2)) {
    print("Local alignment not found")
    return(NULL)
  }
  alini <- function(x,y) {
    matrix(c(as.character(alp1p2[x])
             ,as.character(alp1p2[y])),nrow = 1)
  }
  E1 <- get.edgelist(net1)
  if (dim(E1)[1] == 0) {
    return(0)
  }
  eds <- t(mapply(alini,E1[,1],E1[,2]))
  nas <- unique(c(which(is.na(eds[,2])),which(is.na(eds[,1]))))
  if (length(nas) > 0) {
    eds <- eds[ - nas,]
    eds1 <- E1[ - nas,]
  }
  if (is.null(dim(eds))) {
    if (length(eds) == 2) {
      eds <- cbind(eds[1],eds[2])
      eds1 <- cbind(eds1[1],eds1[2])
    }
  }
  eds <- rbind(eds,eds1)
  Gnet3 <- graph.edgelist(eds,directed = FALSE)

  cols <- rainbow(length(als))
  new.edges <- cbind(names(alp1p2),alp1p2)
  net3 <- graph.data.frame(new.edges,directed = FALSE)

  nets <- list(net3)
  for (al in als) {
    new.edges <- cbind(names(al),al)
    net3 <- graph.data.frame(new.edges,directed = FALSE)
    nets <- append(nets,list(net3))
  }
  net12 <- induced.subgraph(net1,vids = unique(unlist(lapply(als,names))))
  net22 <- induced.subgraph(net2,vids = unique(unlist(als)))
  net4 <- graph.union(net12,net22)
  net5 <- net4
  for (i in 1:length(nets)) {
    net5 <- graph.union(net5,nets[[i]])
  }
  num.eds <- ecount(net5)
  eds <- get.edgelist(net5)
  newedges <- cbind(names(global),global)
  net3 <- graph.data.frame(newedges,directed = FALSE)
  net3 <- induced.subgraph(net3, vids = intersect(V(net3)$name,V(net5)$name))
  E(net5)$color <- "black"
  E(net5)$lty <- 3
  E(net5)$width <- 0.5
  E(net5)$edge.curved <- 1
  ids <- c()
  eds1 <- get.edgelist(nets[[1]])
  ids2 <- get.edge.ids(net5,t(eds1),FALSE)
  ids2 <- setdiff(ids2,ids)
  net5 <- set.edge.attribute(net5,name = "color",index = ids2,cols[1])
  net5 <- set.edge.attribute(net5,name = "edge.curved",index = ids2,0)
  net5 <- set.edge.attribute(net5,name = "lty",index = ids2,2)

  ids <- c(ids,ids2)

  for (i in 2:length(nets)) {
    eds1 <- get.edgelist(nets[[i]])
    ids2 <- get.edge.ids(net5,t(eds1),FALSE)
    ids2 <- setdiff(ids2,ids)
    net5 <- set.edge.attribute(net5,name = "color",index = ids2,cols[i])
    net5 <- set.edge.attribute(net5,name = "edge.curved",index = ids2,0)
    ids <- c(ids,ids2)
  }
  eds1 <- get.edgelist(net3)
  ids <- get.edge.ids(net5,t(eds1),FALSE)

  net5 <- set.edge.attribute(net5,name = "lty",index = ids,1)
  net5 <- set.edge.attribute(net5,name = "edge.curved",index = ids,0)
  net5 <- set.edge.attribute(net5,name = "width",index = ids,1)

  eds1 <- get.edgelist(Gnet3)
  ids <- get.edge.ids(net5,t(eds1),FALSE)
  net5 <- set.edge.attribute(net5,name = "lty",index = ids,1)
  net5 <- set.edge.attribute(net5,name = "width",index = ids,1)
  net5 <- set.edge.attribute(net5,name = "color",index = ids,cols[1])

  eds1 <- get.edgelist(net22)
  ids <- get.edge.ids(net5,t(eds1),FALSE)
  net5 <- set.edge.attribute(net5,name = "edge.curved",index = ids, - 1)
  net5 <- induced.subgraph(net5,vids = union(alp1p2,names(alp1p2)))
  lay <- cbind(10,rep(0,vcount(net5)))
  cont1 <- 1
  cont2 <- 1
  for (i in 1:vcount(net5)) {
    p <- V(net5)$name[i]
    if (is.element(p,V(net12)$name)) {
      lay[i,] <- c(0,cont1)
      cont1 <- cont1 + 1
    }
    else{
      lay[i,] <- c(10,cont2)
      cont2 <- cont2 + 1

    }
  }
  e1 <- ecount(net5)
  v1 <- vcount(net5)
  print(paste(
    "plotting a graph with ",as.character(v1),
    " vertices and ",e1," edges",sep = ""
  ))
  net5 <- set.vertex.attribute(net5,"pos",
                               index = 1:vcount(net12),value = 1)
  net5 <- set.vertex.attribute(net5,"pos",
                               index = (vcount(net12) + 1):(vcount(net22) +
                                                              vcount(net12)),
                               value = - 1)

  plot(
    net5,layout = layout.norm(lay),rescale = TRUE,
    vertex.label.dist = V(net5)$pos,
    vertex.label.degree = pi,edge.curved = E(net5)$edge.curved,...
  )
}






