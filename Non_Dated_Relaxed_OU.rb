contData<- readContinuousCharacterData("30pca.nex")

# print out the data
contData

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

rootAge ~ dnUniform(4.25,6.00)

# more moves --- this time on the rootAge
moves.append( mvSlide(rootAge, delta = 0.01, weight =2))
moves.append( mvSlide(rootAge, delta = 0.1, weight =2))
moves.append( mvSlide(rootAge, delta = 1.0, weight =2))

psi ~ dnBirthDeath(lambda=speciation, mu=extinction, rho=sampling_fraction, rootAge=rootAge, samplingStrategy="uniform", condition="nTaxa", taxa=taxa)


moves.append(mvNarrow(psi, weight=n_taxa/2.0) )
moves.append(mvFNPR(psi, weight=n_taxa/10.0) )
moves.append(mvNodeTimeSlideUniform(psi, weight=n_taxa) )
moves.append(mvNNI(psi, weight=20.0))


ntips     <- psi.ntips()
nbranches <- 2 * ntips - 2



##########################
# Specify the rate model #
##########################

# specify the rate parameter
sigma2 ~ dnLoguniform(1e-3, 1)
moves.append( mvScale(sigma2, weight=1.0) )

# specify the strength parameter
alpha ~ dnExponential(1)
moves.append( mvScale(alpha, weight=1.0) )

# specify theta at the root of the tree
theta_root ~ dnUniform(-5, 5)
moves.append( mvSlide(theta_root, weight=1.0) )

# specify the prior on the number of optimum shifts
expected_number_of_shifts <- 9
shift_probability    <- expected_number_of_shifts / nbranches

# specify the prior on the magnitude of optimum shifts
shift_distribution = dnNormal(0, 0.587)

# specify the branch-specific thetas
for(i in nbranches:1) {

    # draw the theta shift from a mixture distribution
    branch_deltas[i] ~ dnReversibleJumpMixture(0, shift_distribution, Probability(1 - shift_probability) )

    # compute the theta for the branch
    if ( psi.isRoot( psi.parent(i) ) ) {
       branch_thetas[i] := theta_root + branch_deltas[i]
    } else {
       branch_thetas[i] := branch_thetas[psi.parent(i)] + branch_deltas[i]
    }

    # keep track of whether the branch has a shift
    branch_theta_shift[i] := ifelse( branch_deltas[i] == 0, 0, 1 )

    # use reversible-jump to move between models with and without
    # shifts on the branch
    moves.append( mvRJSwitch(branch_deltas[i], weight=1) )

    # include proposals on the shift (when it is not 1)
    moves.append( mvScale(branch_deltas[i], weight=1) )

}

# keep track of the number of theta shifts
num_theta_changes := sum( branch_theta_shift )

##########################
# Specify the BM process #
##########################

X ~ dnPhyloOrnsteinUhlenbeckREML(psi, alpha, branch_thetas, sigma2^0.5, rootStates=theta_root)
X.clamp(contData)

#############
# The Model #
#############

mymodel = model(X)

### set up the monitors that will output parameter values to file and screen
monitors.append(mnModel(filename="writeup/r_OU30.log", printgen=10) )
monitors.append(mnScreen(printgen=1000, sigma2, num_theta_changes) )
monitors.append(mnExtNewick(filename="writeup/r_OU30.trees", isNodeParameter=TRUE, printgen=10, separator=TAB, tree=psi, branch_thetas))


################
# The Analysis #
################

### workspace mcmc ###
mymcmc = mcmc(mymodel, monitors, moves)


### run the MCMC ###
mymcmc.burnin(generations=10000, tuningInterval=500)
mymcmc.run(generations=100000)

#### annotate the tree with the average theta per branch
#treetrace = readTreeTrace("writeup/r_OU30.trees")
#mcc_tree = mccTree(treetrace,"f5/relaxed_OU30.tre")
#map_tree = mapTree(treetrace,"f5/relaxed_OU30.tre")




