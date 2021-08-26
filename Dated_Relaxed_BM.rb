#######################
# Reading in the Data #
#######################
# Create the full list of taxa (including all fossils and extant species #
taxa <- readTaxonData("vasetaxafull.tsv")

# this file contains only the taxa for which sequence data are available #
data <- readContinuousCharacterData("fullpca30nk.nex")

## helpers
n_taxa <- taxa.size()

moves = VectorMoves()

##########################################################################################
# Joint Fossilized Birth-Death Process prior on the topology and fossil occurrence times #
##########################################################################################


# Define exponential priors on the birth rate and death rate #
lambda <- ln(10)
speciation_rate ~ dnExponential(lambda)
extinction_rate ~ dnExponential(lambda)

# Specify a scale move on the speciation_rate parameter #
# This move will be applied with 3 different tuning values (lambda) to help improve mixing # 
moves.append( mvScale(speciation_rate, lambda=0.01, weight=1) )
moves.append( mvScale(speciation_rate, lambda=0.1,  weight=1) )
moves.append( mvScale(speciation_rate, lambda=1.0,  weight=1) )

# Specify a sliding-window move on the extinction_rate parameter #
# This move will be applied with 3 different window widths (delta) to help improve mixing # 
moves.append( mvScale(extinction_rate, lambda=0.01, weight=1) )
moves.append( mvScale(extinction_rate, lambda=0.1,  weight=1) )
moves.append( mvScale(extinction_rate, lambda=1,    weight=1) )

# Create deterministic nodes for the diversification and turnover rates so that they can be monitored #
diversification := speciation_rate - extinction_rate
turnover := extinction_rate/speciation_rate

# Fix the probability of sampling parameter (rho) to 1, #
# because all extant bears are represented in this analysis #
rho <- 1.0

# Assume an exponential prior on the rate of sampling fossils (psi) #
psi ~ dnExponential(1000) 

# Specify a scale move on the psi parameter #
# This move will be applied with 3 different tuning values (lambda) to help improve mixing # 
moves.append( mvScale(psi, lambda=0.01, weight=1) )
moves.append( mvScale(psi, lambda=0.1,  weight=1) )
moves.append( mvScale(psi, lambda=1,    weight=1) )

# The FBD is conditioned on a starting time for the process, which is the origin time #
# Specify a uniform prior on the origin #
origin_time ~ dnUniform(425, 600)
moves.append( mvSlide(origin_time, delta = 0.01, weight =2))
moves.append( mvSlide(origin_time, delta = 0.1, weight =2))
moves.append( mvSlide(origin_time, delta = 1, weight =2))

### Define the tree-prior distribution as the fossilized birth-death process ###
fbd_tree ~ dnFBDRP(origin=origin_time, lambda=speciation_rate, mu=extinction_rate, psi=psi, rho=rho, taxa=taxa)

# Specify moves on the tree and node times #
# These moves update the tree topology 
moves.append( mvFNPR(fbd_tree, weight=15.0) )
moves.append( mvCollapseExpandFossilBranch(fbd_tree, origin_time, weight=6.0) )

# These moves update the node ages #
# Because we are conditioning on the origin time, we must also sample the root node age #
moves.append( mvNodeTimeSlideUniform(fbd_tree, weight=40.0) )
moves.append( mvRootTimeSlideUniform(fbd_tree, origin_time, weight=5.0) )

##########################
# Specify the rate model #
##########################
ntips <- fbd_tree.ntips()
nbranches <- 2 * ntips - 2

sigma2_root ~ dnLoguniform(1e-5, 1)
moves.append( mvScale(sigma2_root,lambda=0.01, weight=1.0) )
moves.append( mvScale(sigma2_root,lambda=0.1, weight=1.0) )
moves.append( mvScale(sigma2_root,lambda=1, weight=1.0) )

# specify the prior on the number of rate shifts
expected_number_of_shifts <- 5
rate_shift_probability    <- expected_number_of_shifts / nbranches

# specify the prior on the magnitude of rate shifts
sd = 0.578
rate_shift_distribution = dnLognormal(-sd^2/2, sd)

# specify the branch-specific rates
for(i in nbranches:1) {

    # draw the rate multiplier from a mixture distribution
    branch_rate_multiplier[i] ~ dnReversibleJumpMixture(1, rate_shift_distribution, Probability(1 - rate_shift_probability) )

    # compute the rate for the branch
    if ( fbd_tree.isRoot( fbd_tree.parent(i) ) ) {
       branch_rates[i] := sigma2_root * branch_rate_multiplier[i]
    } else {
       branch_rates[i] := branch_rates[fbd_tree.parent(i)] * branch_rate_multiplier[i]
    }

    # keep track of whether the branch has a rate shift
    branch_rate_shift[i] := ifelse( branch_rate_multiplier[i] == 1, 0, 1 )

    # use reversible-jump to move between models with and without
    # shifts on the branch
    moves.append( mvRJSwitch(branch_rate_multiplier[i], weight=1) )

    # include proposals on the rate mutliplier (when it is not 1)
    moves.append( mvScale(branch_rate_multiplier[i], weight=1) )

}

# keep track of the number of rate shifts
num_rate_changes := sum( branch_rate_shift )

##########################
# Specify the BM process #
##########################

X ~ dnPhyloBrownianREML(fbd_tree, branchRates=branch_rates^0.5)
X.clamp(data)


########
# MCMC #
########

# initialize the model object #
mymodel = model(X)

monitors = VectorMonitors()

# Create a vector of monitors #
monitors.append( mnModel(filename="output/fossilvase_rbmbig5.log", printgen=10) )
monitors.append( mnExtNewick(filename="output/fossilvase_rbmbig5.trees", isNodeParameter=TRUE, printgen=10, separator=TAB, tree=fbd_tree, branch_rates) )
monitors.append( mnScreen(printgen=25, origin_time, num_rate_changes) )

# Initialize the MCMC object #
mymcmc = mcmc(mymodel, monitors, moves, nruns=4, combine="mixed")

# Run the MCMC #
mymcmc.burnin(generations=10000, tuningInterval=500)
mymcmc.run(generations=100000)

# Read in the tree trace and construct the maximum clade credibility (MCC) tree #
trace = readTreeTrace("output/fossilvase_rbmbig5.trees")

# Summarize tree trace and save MCC tree to file
mccTree(trace, file="output/fossilbig5.tre" )
