

### Read in the trees
T <- readTrees("output/fossilbigf.tre")[1]
ntips     <- T.ntips()
nbranches <- 2 * ntips - 2

### Read in the character data
data <- readContinuousCharacterData("fullpca30nk.nex")

# Create some vector for the moves and monitors of this analysis
moves    = VectorMoves()
monitors = VectorMonitors()

##########################
# Specify the tree model #
##########################

tree <- T

##########################
# Specify the rate model #
##########################

# specify the rate at the root
sigma2_root ~ dnLoguniform(1e-3, 1)
moves.append( mvScale(sigma2_root, weight=1.0) )


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
    if ( tree.isRoot( tree.parent(i) ) ) {
       branch_rates[i] := sigma2_root * branch_rate_multiplier[i]
    } else {
       branch_rates[i] := branch_rates[tree.parent(i)] * branch_rate_multiplier[i]
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

X ~ dnPhyloBrownianREML(tree, branchRates=branch_rates^0.5)
X.clamp(data)



#############
# The Model #
#############

mymodel = model(sigma2_root)

### set up the monitors that will output parameter values to file and screen
monitors.append( mnModel(filename="output/fossilrbm.log", printgen=10) )
monitors.append( mnScreen(printgen=1000, sigma2_root, num_rate_changes) )
monitors.append( mnExtNewick(filename="output/fossilrbm.trees", isNodeParameter=TRUE, printgen=10, separator=TAB, tree=tree, branch_rates) )



################
# The Analysis #
################

### workspace mcmc ###
mymcmc = mcmc(mymodel, monitors, moves, nruns=4, combine="mixed")


### run the MCMC ###
mymcmc.burnin(generations=10000, tuningInterval=500)
mymcmc.run(generations=100000)

## create the annotated tree
treetrace = readTreeTrace("output/fossilrbm.trees")
map_tree = mapTree(treetrace,"output/fossilrbm.tre")



