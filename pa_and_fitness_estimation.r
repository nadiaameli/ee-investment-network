# Accompanying code for paper titled “The hidden structure of energy efficiency finance”
#
# This code is used to generate Preferential Attachment and Node Fitness estimations using the PAFit library in R
# Please see the paper and accompanying supplementary material for further details on the method
#
# This generic code accepts an input data file in a txt format
# The data is organised in 3 columns: Source Node, Target Node and Timestep
#

set.seed(1)
library(PAFit)

# read data from the file and construct the network
real_net1 <- graph_from_file("inputdata.txt")
summary(real_net1)

# generate summary statistics
stats_real_net1 <- get_statistics(real_net1,binning=FALSE, deg_threshold=1)
summary(stats_real_net1)

# initiate joint estimation of Preferential Attachment and Node Fitness
result_real_net1 <- joint_estimate(real_net1,stats_real_net1)
summary(result_real_net1)

# display node degrees
result_real_net1$estimate_result$k

# display attachment function corresponding to each k
result_real_net1$estimate_result$A

# display alpha parameter for Preferential Attachment when fitted to A(k) = k ^ alpha
result_real_net1$estimate_result$alpha

# display confidence intervals for attachment function
result_real_net1$estimate_result$ci

# display node fitnesses
result_real_net1$estimate_result$f

# display objective function value
result_real_net1$estimate_result$objective_value

# plot attachment function
plot(result_real_net1, stats_real_net1)
lines(stats_real_net1$center_k, stats_real_net1$center_k,col="red")

# plot fitness distribution
plot(result_real_net1, stats_real_net1, plot="f")

# Estimate the Preferential Attachment function in isolation
result_pa_only <- only_A_estimate(real_net1, stats_real_net1)
summary(result_pa_only)

# display node degrees
result_pa_only$estimate_result$k

# display attachment function corresponding to each k
result_pa_only$estimate_result$A

# display alpha parameter for Preferential Attachment when fitted to A(k) = k ^ alpha
result_pa_only$estimate_result$alpha

# display confidence intervals for attachment function
result_pa_only$estimate_result$ci

# display objective function value
result_pa_only$estimate_result$objective_value

# plot attachment function
plot(result_pa_only, stats_real_net1)
lines(stats_real_net1$center_k, stats_real_net1$center_k,col="red")

# Estimate Node Fitness using Bianconi-Barabasi model formulation
result_fit_only_BB <- only_F_estimate(real_net1,stats_real_net1, model_A="Linear")
summary(result_fit_only_BB)

# display node fitnesses
result_fit_only_BB$estimate_result$f

# display objective function value
result_fit_only_BB$estimate_result$objective_value

# plot fitness distribution
plot(result_fit_only_BB, stats_real_net1, plot="f")

# Estimate Node Fitness using Caldarelli model formulation
result_fit_only_C <- only_F_estimate(real_net1,stats_real_net1, model_A="Constant")
summary(result_fit_only_C)

# display node fitnesses
result_fit_only_C$estimate_result$f

# display objective function value
result_fit_only_C$estimate_result$objective_value

# plot fitness distribution
plot(result_fit_only_C, stats_real_net1, plot="f")
