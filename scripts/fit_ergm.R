#workspace
#devtools::install_github("tedhchen/multilayer.ergm", build_vignettes = FALSE) #not available via conda-forge
#devtools::install_github("tedhchen/tc.ergmterms", build_vignettes = FALSE) #not available via conda-forge
library(Rglpk)
library(igraph)
library(tc.ergmterms)
library(multilayer.ergm)

#MODEL0
gm = read_graph("mention.xml", format = "graphml")
gr = read_graph("complete.xml", format = "graphml")
ax = read.table("attributex.csv", sep = ";", header = TRUE, quote = "")

ax$layer.mem[ax$realm == "external" | ax$realm == "regional"] = 1
ax$layer.mem[ax$realm == "national"] = 2

ax$closeness = scales::rescale(ax$closeness, to = c(0, 1))
ax$betweenness = scales::rescale(ax$betweenness, to = c(0, 1))
ax$boomerang[ax$coalition == "forest_watchdog"] = "transboundary_coalition"
ax$boomerang[ax$type == "IDA"] = "donor government"
ax$boomerang[ax$type == "IGO"| ax$type == "DFI"] = "intergovernmental organisation"
ax$boomerang[ax$type == "MRY" | ax$type == "POC"] = "state_apparatus"
ax$boomerang[is.na(ax$boomerang)] = "other"
ax = ax[order(ax$actor, ax$layer.mem),]
rownames(ax) = ax$actor

mm = matrix(0, nrow(ax), nrow(ax))
dimnames(mm) = list(ax$actor, ax$actor)
mx = as_adjacency_matrix(gm, attr = "weight", sparse = FALSE)
rw = rownames(mx)[rownames(mx) %in% rownames(mm)]
cl = colnames(mx)[colnames(mx) %in% colnames(mm)]
mm[rw, cl] = mx[rw, cl]

rt = as_adjacency_matrix(gr, attr = "weight", sparse = FALSE)
mm[mm > 1] = 1
rt[rt > 1] = 1

edge_density(graph_from_adjacency_matrix(rt[ax$actor[ax$layer == 1], ax$actor[ax$layer == 1]], mode = "undirected"))
edge_density(graph_from_incidence_matrix(rt[ax$actor[ax$layer == 2], ax$actor[ax$layer == 1]]))
edge_density(graph_from_adjacency_matrix(rt[ax$actor[ax$layer == 2], ax$actor[ax$layer == 2]], mode = "undirected"))

edge_density(graph_from_adjacency_matrix(mm[ax$actor[ax$layer == 1], ax$actor[ax$layer == 1]], mode = "directed"))
edge_density(graph_from_incidence_matrix(mm[ax$actor[ax$layer == 2], ax$actor[ax$layer == 1]]))
edge_density(graph_from_adjacency_matrix(mm[ax$actor[ax$layer == 2], ax$actor[ax$layer == 2]], mode = "directed"))

#network
nw = network(mm, directed = TRUE, loops = FALSE, matrix.type = "adjacency")
nw %v% "type" = as.character(ax$type)
nw %v% "coalition" = as.character(ax$coalition)
nw %v% "boomerang" = as.character(ax$boomerang)
nw %v% "closeness" = as.numeric(ax$closeness)
nw %v% "betweenness" = as.numeric(ax$betweenness)
nw %v% "layer.mem" = as.numeric(ax$layer.mem)
check.multilayer(nw)

#nodemix
node = as.data.frame(summary(nw ~ nodemix(c("layer.mem", "boomerang"), levels2 = TRUE)))
colnames(node) = "count"
rownames(node) = substring(rownames(node), 25)
node$index = 1:nrow(node)
print(node[order(rownames(node)),])

#1.transboundary_coalition.1.donor government                         90     4
#1.transboundary_coalition.1.intergovernmental organisation          160    11
#1.transboundary_coalition.2.state_apparatus                          26    39

#2.transboundary_coalition.1.donor government                          9     7
#2.transboundary_coalition.1.intergovernmental organisation           20    14
#2.transboundary_coalition.2.state_apparatus                          27    42

#1.donor government.1.intergovernmental organisation                 129     8
#1.intergovernmental organisation.1.donor government                 126     2

#1.donor government.2.state_apparatus                                 86    36
#1.intergovernmental organisation.2.state_apparatus                   54    37

#2.state_apparatus.1.donor government                                 77     6
#2.state_apparatus.1.intergovernmental organisation                   46    13
#2.state_apparatus.1.transboundary_coalition                           3    27
#2.state_apparatus.2.transboundary_coalition                           3    48

#model
time = Sys.time()
m = ergm(nw
~ edges
+ nodematch("coalition", diff = TRUE, levels = c(1, 2, 4, 5, 6))
+ nodeicov("closeness")
+ nodeicov("betweenness")
+ nodeocov("closeness")
+ nodeocov("betweenness")
+ mutual
+ gwidegree(decay = 0.5, fixed = TRUE)#, attr = "type", levels = c(2, 3, 7, 8, 9))
+ gwodegree(decay = 0.5, fixed = TRUE)#, attr = "type", levels = c(2, 3, 7, 8, 9))
+ gwdsp(decay = 0.5, fixed = TRUE)
+ gwesp_scc_senderattr(decay = 0.5, type = 2, sender_attr = "layer.mem", value = 2)
+ nodemix(c("layer.mem", "boomerang"), levels2 = c(4, 11, 39,  7, 14, 42,  8, 2,  6, 13, 27, 48,  36, 37))
, eval.loglik = TRUE
, check.degeneracy = TRUE
, verbose = 0
, control = control.ergm(seed = 0, MCMC.interval = 10000, MCMC.burnin = 100000, MCMC.samplesize = 10000, parallel = 5, parallel.type = "PSOCK", MCMLE.termination = "Hummel", MCMLE.effectiveSize = NULL))
print(Sys.time() - time)
saveRDS(m, "\\\\ad.helsinki.fi/home/a/ajmalkam/Documents/congoERGM_0.rds")

#MODEL1
gm = read_graph("mention.xml", format = "graphml")
gr = read_graph("complete.xml", format = "graphml")
ax = read.table("attributex.csv", sep = ";", header = TRUE, quote = "")
ax$layer.mem[ax$realm == "external" | ax$realm == "regional"] = 1
ax$layer.mem[ax$realm == "national"] = 2

ax$closeness = scales::rescale(ax$closeness, to = c(0, 1))
ax$betweenness = scales::rescale(ax$betweenness, to = c(0, 1))
ax$boomerang[ax$coalition == "forest_watchdog"] = "transboundary_coalition"
ax$boomerang[ax$type == "IDA" | ax$type == "IGO"| ax$type == "DFI"] = "powerful_third_party"
ax$boomerang[ax$type == "MRY" | ax$type == "POC"] = "state_apparatus"
ax$boomerang[is.na(ax$boomerang)] = "other"
ax = ax[order(ax$actor, ax$layer.mem),]
rownames(ax) = ax$actor

mm = matrix(0, nrow(ax), nrow(ax))
dimnames(mm) = list(ax$actor, ax$actor)
mx = as_adjacency_matrix(gm, attr = "weight", sparse = FALSE)
rw = rownames(mx)[rownames(mx) %in% rownames(mm)]
cl = colnames(mx)[colnames(mx) %in% colnames(mm)]
mm[rw, cl] = mx[rw, cl]

rt = as_adjacency_matrix(gr, attr = "weight", sparse = FALSE)
mm[mm > 1] = 1
rt[rt > 1] = 1

edge_density(graph_from_adjacency_matrix(rt[ax$actor[ax$layer == 1], ax$actor[ax$layer == 1]], mode = "undirected"))
edge_density(graph_from_incidence_matrix(rt[ax$actor[ax$layer == 2], ax$actor[ax$layer == 1]]))
edge_density(graph_from_adjacency_matrix(rt[ax$actor[ax$layer == 2], ax$actor[ax$layer == 2]], mode = "undirected"))

edge_density(graph_from_adjacency_matrix(mm[ax$actor[ax$layer == 1], ax$actor[ax$layer == 1]], mode = "directed"))
edge_density(graph_from_incidence_matrix(mm[ax$actor[ax$layer == 2], ax$actor[ax$layer == 1]]))
edge_density(graph_from_adjacency_matrix(mm[ax$actor[ax$layer == 2], ax$actor[ax$layer == 2]], mode = "directed"))

#network
nw = network(mm, directed = TRUE, loops = FALSE, matrix.type = "adjacency")
nw %v% "type" = as.character(ax$type)
nw %v% "coalition" = as.character(ax$coalition)
nw %v% "boomerang" = as.character(ax$boomerang)
nw %v% "closeness" = as.numeric(ax$closeness)
nw %v% "betweenness" = as.numeric(ax$betweenness)
nw %v% "layer.mem" = as.numeric(ax$layer.mem)
check.multilayer(nw)

#nodemix
node = as.data.frame(summary(nw ~ nodemix(c("layer.mem", "boomerang"), levels2 = TRUE)))
colnames(node) = "count"
rownames(node) = substring(rownames(node), 25)
node$index = 1:nrow(node)
print(node[order(rownames(node)),])

#1.transboundary_coalition.1.powerful_third_party      250     9
#1.transboundary_coalition.2.state_apparatus            26    27
#
#2.transboundary_coalition.1.powerful_third_party       29    12
#2.transboundary_coalition.2.state_apparatus            27    30
#
#2.state_apparatus.1.powerful_third_party              123    11
#2.state_apparatus.1.transboundary_coalition             3    17
#2.state_apparatus.2.transboundary_coalition             3    35
#
#1.powerful_third_party.2.state_apparatus              140    26

#model
time = Sys.time()
m = ergm(nw
~ edges
+ nodematch("coalition", diff = TRUE, levels = c(1, 2, 4, 5, 6))
+ nodeicov("closeness")
+ nodeicov("betweenness")
+ nodeocov("closeness")
+ nodeocov("betweenness")
+ mutual
+ gwidegree(decay = 0.5, fixed = TRUE)#, attr = "type", levels = c(2, 3, 7, 8, 9)) #0.2?????????????????????????
+ gwodegree(decay = 0.5, fixed = TRUE)#, attr = "type", levels = c(2, 3, 7, 8, 9))
+ gwdsp(decay = 0.5, fixed = TRUE)
+ gwesp_scc_senderattr(decay = 0.5, type = 2, sender_attr = "layer.mem", value = 2)
+ nodemix(c("layer.mem", "boomerang"), levels2 = c(9, 27,  12, 30,  11, 17, 35,  26))
, eval.loglik = TRUE
, check.degeneracy = TRUE
, verbose = 0
, control = control.ergm(seed = 0, MCMC.interval = 10000, MCMC.burnin = 100000, MCMC.samplesize = 10000, parallel = 5, parallel.type = "PSOCK", MCMLE.termination = "Hummel", MCMLE.effectiveSize = NULL))
print(Sys.time() - time)
saveRDS(m, "\\\\ad.helsinki.fi/home/a/ajmalkam/Documents/congoERGM_1.rds")

#print
m0 = readRDS("/m/triton/scratch/work/malkama5/congoCoalition/congoERGM_0.rds")
m1 = readRDS("/m/triton/scratch/work/malkama5/congoCoalition/congoERGM_1.rds")
texreg::htmlreg(l = list(m0, m1), file = "congoERGM.html", stars = c(0.01, 0.05, 0.1), digits = 2, single.row = TRUE)

