## RevBayes code to be run on v1.0.13
## Authored by Benjamin Tenmann [01/07/2020]
## Written for reconstructing phylogenies
## Inferres from continuous data
## Edited by Rahim Dina 


contData<- readContinuousCharacterData("30pca.nex")

# print out the data
contData

#extracting the names of the taxa
names <- contData.names()
n_taxa <- contData.size()
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
moves.append(mvScale(speciation, lambda=0.01, weight=2))
moves.append(mvScale(speciation, lambda=0.1, weight=2))
moves.append(mvScale(speciation, lambda=1.0, weight=2))

sampling_fraction <- 1  # sample size / population size = sampling fraction

rootAge ~ dnUniform(3.75, 6.0)  # uniformly distributed root age

# more moves --- this time on the rootAge
moves.append(mvSlide(rootAge, delta = 0.01, weight =2))
moves.append(mvSlide(rootAge, delta = 0.1, weight =2))
moves.append(mvSlide(rootAge, delta = 1.0, weight =2))

#instantiate a Birth-Death prior for tree topology
psi ~ dnBDP(lambda=speciation, mu=extinction, rho=sampling_fraction, rootAge=rootAge, condition = "nTaxa", taxa=names)


#create our node height and topological rearrangement MCMC moves


moves.append(mvNarrow(psi, weight=n_taxa))
moves.append(mvFNPR(psi, weight=n_taxa/4)) 
moves.append(mvNodeTimeSlideUniform(psi, weight=n_taxa))
moves.append(mvSubtreeScale(psi, weight=n_taxa/5.0)) 

##########################
# Specify the rate model #
##########################

sigma2 ~ dnLoguniform(1e-3, 1)
moves.append( mvScale(sigma2, weight=1.0) )

alpha ~ dnExponential(1)
moves.append( mvScale(alpha, weight=1.0) )


theta ~ dnUniform(-5, 5)
moves.append( mvSlide(theta, weight=1.0) )

##########################
# Specify the BM process #
##########################

X ~ dnPhyloOrnsteinUhlenbeckREML(psi, alpha, theta, sigma2^0.5, rootStates=theta)
X.clamp(contData)


bmv = model(X) #link sigma param w/ BM model
monitors.append(mnScreen(printgen=5000))
monitors.append(mnFile(filename="writeup/OUv30.log", printgen=50, separator = TAB,speciation, rootAge, sigma2, alpha, theta))
monitors.append(mnExtNewick(filename="writeup/OUv30.trees", printgen=10,separator = TAB, psi, speciation, rootAge, sigma2, alpha, theta))

# set up MCMC
chain = mcmc(bmv, monitors, moves,nruns=4,combine="mixed")
# do a burn-in; i.e. throw away a certain number of samples at the start of the MCMC
chain.burnin(generations=10000,tuningInterval=500)

# run the MCMc
chain.run(100000)


# read in the tree-trace
#treetrace = readTreeTrace(file = "f4/OUv100.trees",treetype="clock")


# create a maximum a posteriori tree from the tree-trace
#map = mapTree( file="f4/OUv100.nex", treetrace, ccp = TRUE )

#pow_p = powerPosterior(bmv, moves, monitors, "output/model1.out", cats=50)
#pow_p.burnin(generations=10000,tuningInterval=1000)
#pow_p.run(generations=1000)
#ss = steppingStoneSampler(file="output/model1.out", powerColumnName="power", likelihoodColumnName="likelihood")
#ss.marginal()


