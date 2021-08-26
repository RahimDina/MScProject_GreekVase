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

rootAge ~ dnUniform(4.25,6.00)  # uniformly distributed root age

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

sigma2 ~ dnLoguniform(1e-3, 2.5)  # branch-rate prior

# branch-rate moves
moves.append(mvScale(sigma2, lambda=0.01, weight=1) )
moves.append(mvScale(sigma2, lambda=0.1,  weight=1) )
moves.append(mvScale(sigma2, lambda=1.0,  weight=1) ) 

# specify that we are going calculate BM likelihood using the REML PIC algorithm (see Felsenstein 1973)
traits ~ dnPhyloBrownianREML(psi, branchRates=sigma2^0.5,nSites=contData.nchar())

traits.clamp(contData) #match traits to tips
bmv = model(traits) #link sigma param w/ BM model
monitors.append(mnScreen(printgen=500))
monitors.append(mnFile(filename="writeup/BMv30.log", printgen=10, separator = TAB,speciation, rootAge, sigma2))
monitors.append(mnFile(filename="writeup/BMv30.trees", printgen=10,separator = TAB, psi, rootAge, sigma2))

# set up MCMC
chain = mcmc(bmv, monitors, moves,nruns=4,combine="mixed")
# do a burn-in; i.e. throw away a certain number of samples at the start of the MCMC
chain.burnin(generations=10000,tuningInterval=500)

# run the MCMc
chain.run(100000)

# read in the tree-trace
#treetrace = readTreeTrace(file = "f10/BMv100.trees",treetype="clock")

# create a maximum a posteriori tree from the tree-trace
#map = mapTree( file="f100/BMv100.nex", treetrace, ccp = TRUE )

