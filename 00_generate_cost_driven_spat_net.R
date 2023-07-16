# Generate cost-driven spatial network (following Louf et al. 2013) ------------

# This R script generates a cost-driven spatial network following
# <div class="csl-entry">Louf, R., Jensen, P., &#38; Barthelemy, M. (2013).
# Emergence of hierarchy in cost-driven growth of spatial networks. 
# <i>Proceedings of the National Academy of Sciences of the United States of 
#America</i>, <i>110</i>(22), 8824???8829. https://doi.org/10.1073/pnas.1222441110</div>

# 0 Load libraries ----
library(spatstat.random)
library(tidyverse)
library(igraph)

# 1. Generate a set of nodes in 2D space ------
set.seed(1) # set random seed
n = 400 # Set number of "cities"

# 2. Randomly assign "populations" to each of the nodes using a power law -----
# Pm(X) = mu/x^(mu + 1)
# Parameter values used in the Louf, Jensen and Barthelemy (2012) paper
mu = 1.1
a = 1.1
k = 1
beta = 1.1 #1.1
# eta = Not specified here because k assumed 1 (see paper)

pop <- function(x){ # Generate city populations using a power law with exponent mu + 1
  mu/x^(mu + 1)
}

# Make network data frame
dnet <- rpoint(n) %>% as.data.frame() %>%
  mutate(pop = pop(runif(n, min = 0, max = 1)),
         # add node id
         id = 1,
         id.val = cumsum(id),
         id = paste0('n', cumsum(id))) #%>% 
#arrange(-pop)

# Important notes:
# L is the typical edge length of the whole system
# beta is k/eta, or the relative importance of the cost with regards to the benefits

# 3. Randomly select a node as the root ----
root <- sample_n(dnet, size = 1)
dnet %>% nrow()
message('The randomly chosen root is ', root$id,
        ' located at x = ',round(root$x, 2),
        ', y = ',round(root$y, 2), '. It has a population of ' , root$pop)

# Then establish a link between the last node and the node that maximizes:
# R'ij = k*Mi*Mj/dij^(a-1) - beta*dij

# 4. Connect them according to a Beta*/B value -----

network <- data.frame()

# Make network data frame
for (j in 1:(nrow(dnet)-2)){
  if(j == 1){
    i <- root$id[1]
  }
  
  if(j<= n-1 & j>=3){
    root <- root %>% rbind(dnet %>% filter(!id %in% root$id) %>% sample_n(., size = 1))
  }
  if(j>=3){
    i <- root$id[j]
  }
  
  
  popi <- dnet %>% filter(id == i) %>% select(pop)
  csi <- dnet %>% filter(id == i) %>% select(x, y)
  # Calc euclidean distances
  if(i == root$id[1] & j == 1){
    ds <-(sqrt(rowSums(((dnet %>% filter(id !=i) %>% select(x, y)) - c(csi))^2)))  # calc Rnotij
    rnotij <- (popi$pop*((dnet%>%filter(id != i))$pop)) / (ds^(a-1)) - beta * ds 
    # calc Rnotij
    
    network <- network %>% rbind(data.frame(from = (dnet %>% filter(id!=i) %>% 
                                                      filter(rnotij == rnotij %>% max()))$id,
                                            to = i,
                                            length = ds[rnotij == rnotij%>%max()] ))
    root <- root %>% rbind(dnet %>% filter(id == network$from))
  }
  if(j>=3){
    ds <-(sqrt(rowSums(((dnet %>% filter(id !=i, id %in% network$from | id %in% network$to #, pop > popi$pop
    ) %>% select(x, y)) - c(csi))^2)))  # calc Rnotij
    rnotij <- (popi$pop*((dnet%>%filter(id != i, id %in% network$from | id %in% network$to #, pop > popi$pop
    ))$pop)) / (ds^(a-1)) - beta * ds
    # calc Rnotij 
    
    network <- network %>% rbind(data.frame(from = i, to = (dnet %>% filter(id!=i, id %in% network$from | id %in% network$to #pop>popi$pop
    ) %>% 
      filter(rnotij == rnotij %>% max()))$id,
    length = ds[rnotij == rnotij%>%max()] ))
    network
  }
}

# Check it looks ok   
network

# Plot it in 2-D space
plot(graph_from_data_frame(network,directed = FALSE), 
     layout = as.matrix(dnet %>% select(x, y)), 
     vertex.size = 3,# arrow.size = 0.5
     vertex.label = NA) 


# Calculate Beta*/B ------

# Total length
sum(network$length) # Actual
(Ltot = mean(network$length)*n) # Approximation
# Typical length
mean(network$length)

# Assuming typical length is mean of length (seems correct based on the paper)
bnot <- function(){
  res<-k * (dnet$pop %>% mean())^2 *( (n/(mean(network$length)^2))^a/2)
  return(res)
}

message(paste('B* = ' , bnot()), 'and B = ', beta, '. So B*/B = ', beta/bnot())

# Potential future work: ----
# Simulate dynamics at 1 time scale and record all activity in the network
# Select inter-observation time (tau), percent of observed nodes and gini coefficient-----
# Reconstruct-----
# Compare reconstruction to network -----

# Troubleshooting -----
## Some things to check if you run into issues: ----
# What does the whole curve look like?
# plot(y = seq(0.0001, 1, by = 0.001), x = log(pop(seq(0.0001, 1, by = 0.001)) ))
# What about a sample of n populations?
# runif(n, min = 0, max = 1)
# hist(log(pop(runif(n, min = 0, max = 1)))) # Looks good

