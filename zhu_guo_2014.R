# include necessary libraries
library(FNN)
library(ggplot2)
library(dplyr)
library(dplyrExtras)
library(readr)
library(reshape2)

# define parameters
k <- 10   # number of nearest neighbors

# # generate dummy origin and destination points
# n <- 100    # number of dummy points
# set.seed(281179)
# o <- data.frame(id = 1:n, x = rnorm(n), y = rnorm(n))
# d <- data.frame(id = 1:n, x = rnorm(n), y = rnorm(n))

# # generate dummy flow matrix
# f <- data.frame(id = 1:n, 
#                 oid = sample(o$id, n, replace = FALSE), 
#                 did = sample(d$id, n, replace = FALSE)) %>% 
#     left_join(o, by = c("oid" = "id")) %>% 
#     rename(ox = x, oy = y) %>% 
#     left_join(d, by = c("did" = "id")) %>% 
#     rename(dx = x, dy = y)


# load data from csv files
centroids <- read_csv("msoa_popweightedcentroids.csv") %>% 
    rename(code = Code, 
           name = Name, 
           east = East, 
           north = North)
flows <- read_csv("wu03ew_v1.csv") %>% 
    select(o = `Area of usual residence`, 
           d = `Area of workplace`, 
           weight = `All categories: Method of travel to work`) %>% 
    filter(weight >= 1000,                # remove flows with weight < 20 and 
           o != d)                       # internal flows

# extract locations
usedCentroids <- c(flows$o, flows$d)
l <- centroids %>% 
    filter(code %in% unique(usedCentroids)) %>%   # limit to required centroids
    mutate(id = row_number()) %>% 
    select(id, 
           code, 
           name, 
           x = east, 
           y = north)

# generate flow matrix
f <- flows %>% 
    left_join(l, by = c("o" = "code")) %>% 
    rename(ox = x, oy = y, oid = id, oname = name) %>% 
    left_join(l, by = c("d" = "code")) %>% 
    rename(dx = x, dy = y, did = id, dname = name) %>% 
    filter(!is.na(oid),                 # remove flows with missing origin and
           !is.na(did)) %>%             # missing destination
    arrange(oid) %>% 
    mutate(fid = row_number()) %>% 
    select(id = fid, 
           oid, 
           oname, 
           did, 
           dname, 
           ox, oy, 
           dx, dy, 
           weight)

# draw plot of flows
xquiet <- scale_x_continuous("", breaks=NULL)
yquiet <- scale_y_continuous("", breaks=NULL)
quiet <- list(xquiet, yquiet)
ggplot(f, aes(ox, oy)) + 
    geom_segment(aes(x = ox, y = oy, xend = dx, yend = dy)) + 
    quiet + 
    coord_equal()

# calculate k nearest neighbors for point locations
knn <- get.knnx(l %>% 
                  arrange(id) %>% 
                  select(x, y) %>% 
                  as.matrix(), 
               l %>% 
                  arrange(id) %>% 
                  select(x, y) %>% 
                  as.matrix(), 
               k = k, 
               algorithm = "cover_tree")$nn.index

# calculate dissimilarity matrices for point locations
diss <- matrix(0, nrow = nrow(l), ncol = nrow(l))
for(i in 1:nrow(l)) {
    for(j in 1:nrow(l)) {
        diss[i, j] <- sum(knn[i, ] %in% knn[j, ])
    }
}

# calculate shared nearest neighbor flow dissimilarity matrix
fdiss <- matrix(0, nrow = nrow(f), ncol = nrow(f))
for(i in 1:nrow(f)) {
    p <- f[i, ]
    for(j in 1:nrow(f)) {
        q <- f[j, ]
        fdiss[i, j] <- 1 - (diss[p$oid, q$oid] / k) * (diss[p$did, q$did] / k)
    }
}

# generate ordered list of contiguous flow pairs
fp <- melt(fdiss) %>% 
    rename(p = Var1, q = Var2, diss = value) %>% 
    filter(p < q) %>%                   # remove diag and lower triangle
    arrange(diss) %>% 
    tbl_df()

# inintialize individual flow clusters
f <- f %>% 
    mutate(cluster = row_number()) %>%  
    tbl_df()

# cluster contiguous flow pairs
for(i in 1:nrow(fp)) {
    # which clusters does current pair (p,q) belong to?
    cx <- f %>% 
        filter(id == as.integer(fp[i, "p"])) %>% 
        select(cluster) %>% 
        as.integer()
    cy <- f %>% 
        filter(id == as.integer(fp[i, "q"])) %>% 
        select(cluster) %>% 
        as.integer()
    # calculate both clusters' origin and destination centroids
    ocx <- f %>% 
        filter(cluster == cx) %>% 
        group_by(cluster) %>% 
        summarise(ox = mean(ox), 
                  oy = mean(oy))
    dcx <- f %>% 
        filter(cluster == cx) %>% 
        group_by(cluster) %>% 
        summarise(dx = mean(dx), 
                  dy = mean(dy))
    ocy <- f %>% 
        filter(cluster == cy) %>% 
        group_by(cluster) %>% 
        summarise(ox = mean(ox), 
                  oy = mean(oy))
    dcy <- f %>% 
        filter(cluster == cy) %>% 
        group_by(cluster) %>% 
        summarise(dx = mean(dx), 
                  dy = mean(dy))
    # find nearest original points for origin and destination cluster centroids 
    ocxnn <- get.knnx(l %>% 
                         select(x, y) %>% 
                         as.matrix(), 
                     ocx %>% 
                         select(ox, oy) %>% 
                         as.matrix(), 
                     k = 1, 
                     algorithm = "cover_tree")$nn.index[1, 1]
    dcxnn <- get.knnx(l %>% 
                          select(x, y) %>% 
                          as.matrix(), 
                      dcx %>% 
                          select(dx, dy) %>% 
                          as.matrix(), 
                      k = 1, 
                      algorithm = "cover_tree")$nn.index[1, 1]
    ocynn <- get.knnx(l %>% 
                          select(x, y) %>% 
                          as.matrix(), 
                      ocy %>% 
                          select(ox, oy) %>% 
                          as.matrix(), 
                      k = 1, 
                      algorithm = "cover_tree")$nn.index[1, 1]
    dcynn <- get.knnx(l %>% 
                          select(x, y) %>% 
                          as.matrix(), 
                      dcy %>% 
                          select(dx, dy) %>% 
                          as.matrix(), 
                      k = 1, 
                      algorithm = "cover_tree")$nn.index[1, 1]
    # calculate flow dissimilarity of median flows (fdiss')
    odissp <- sum(knn[ocxnn, ] %in% knn[ocynn, ])
    ddissp <- sum(knn[dcxnn, ] %in% knn[dcynn, ])
    fdissp <- 1 - (odissp / k) * (ddissp / k)
    # cluster if dissimilarity of median flows < 1
    if(fdissp < 1) {
        f <- mutate_if(f, 
                       id == fp[i, ]$q, 
                       cluster = cx)
    }
}

# calculate clustered flows' origin and destination centroids
oc <- f %>% 
    group_by(cluster) %>% 
    summarise(ox = mean(ox), 
              oy = mean(oy))
dc <- f %>% 
    group_by(cluster) %>% 
    summarise(dx = mean(dx), 
              dy = mean(dy), 
              weight = n())

# generate clustered flows
cf <- left_join(oc, dc, by = "cluster")
cat("generated", nrow(cf), "clustered flows")

# draw plot of clustered flows
ggplot() + 
    geom_segment(data = f, aes(x = ox, y = oy, xend = dx, yend = dy), color = "grey") + 
    geom_segment(data = cf, aes(x = ox, y = oy, xend = dx, yend = dy, 
                                size = weight, col = factor(cluster))) +  
    quiet + 
    coord_equal()