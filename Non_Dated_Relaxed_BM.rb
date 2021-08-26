contData<- readContinuousCharacterData("30pca.nex")
#extracting the names of the taxa
names <- contData.names()
n_taxa <- contData.ntaxa()
taxa <- contData.taxa()
# set the lambda value
lambda <- ln(10)

# parameterise exponential distributin with lmabda
speciation ~ dnExponential(lambda) #speciation rate
extinction <- 0

## this is a pure-birth model ##
# speciation rate is an exponentially distributed stochastic variable
# parameterised with lambda

moves    = VectorMoves()
monitors = VectorMonitors()

# assigning moves to the parameter vector;
# this is to tell the MCMC what to vary how in order to sample from the
# posterior distribution of phylogenetic trees
moves.append( mvScale(speciation, lambda=0.01, weight=2))
moves.append( mvScale(speciation, lambda=0.1, weight=2))
moves.append( mvScale(speciation, lambda=1.0, weight=2))

sampling_fraction <- 1  # sample size / population size = sampling fraction

rootAge ~ dnUniform(4.25,6)

# more moves --- this time on the rootAge
moves.append( mvSlide(rootAge, delta = 0.01, weight =2))
moves.append( mvSlide(rootAge, delta = 0.1, weight =2))
moves.append( mvSlide(rootAge, delta = 1.0, weight =2))

psi ~ dnBirthDeath(lambda=speciation, mu=extinction, rho=sampling_fraction, rootAge=rootAge, samplingStrategy="uniform", condition="nTaxa", taxa=taxa)


moves.append( mvNarrow(psi, weight=n_taxa/2.0) )
moves.append( mvFNPR(psi, weight=n_taxa/10.0) )
moves.append( mvNodeTimeSlideUniform(psi, weight=n_taxa) )
moves.append( mvNNI(psi, weight=20.0))


ntips <- psi.ntips()
nbranches <- 2 * ntips - 2


#########################
# Specify the rate model #
##########################

# specify the rate at the root
sigma2_root ~ dnLoguniform(1e-3, 1)
moves.append( mvScale(sigma2_root, weight=1.0) )


# specify the prior on the number of rate shifts
expected_number_of_shifts <- 9
rate_shift_probability    <- expected_number_of_shifts / nbranches

# specify the prior on the magnitude of rate shifts
sd = 0.578
rate_shift_distribution = dnLognormal(-sd^2/2, sd)

# specify the branch-specific rates
for(i in nbranches:1) {

    # draw the rate multiplier from a mixture distribution
    branch_rate_multiplier[i] ~ dnReversibleJumpMixture(1, rate_shift_distribution, Probability(1 - rate_shift_probability) )

    # compute the rate for the branch
    if ( psi.isRoot( psi.parent(i) ) ) {
       branch_rates[i] := sigma2_root * branch_rate_multiplier[i]
    } else {
       branch_rates[i] := branch_rates[psi.parent(i)] * branch_rate_multiplier[i]
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

traits ~ dnPhyloBrownianREML(psi, branchRates=branch_rates^0.5)
traits.clamp(contData)



#############
# The Model #
#############

mymodel = model(traits)

## set up the monitors that will output parameter values to file and screen
monitors.append( mnModel(filename="writeup/relaxed_BM30.log", printgen=10) )
monitors.append( mnScreen(printgen=1000, sigma2_root, num_rate_changes) )
monitors.append( mnExtNewick(filename="writeup/relaxed_BM30.trees", isNodeParameter=TRUE, printgen=10, separator=TAB, tree=psi,branch_rates))

################
# The Analysis #
################

### workspace mcmc ###
mymcmc = mcmc(mymodel, monitors, moves, nruns=4, combine="mixed")

### run the MCMC ###
mymcmc.burnin(generations=10000, tuningInterval=500)
mymcmc.run(generations=100000)

## create the annotated tree
#treetrace = readTreeTrace("fffrelaxed_BM30.trees", treetype="clock")

#mcc_tree = mccTree(treetrace,"fffrelaxed_BM_MCC30.nex", ccp=TRUE)






