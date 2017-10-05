#Pseudocode for MC simulations
#Farah
#11.01.16

#Import count data for replicates of a specific chemostat RATEseq condition

#Estimate parameters from real data to infer NB probability distribution
  #Calculate mean and variance (dispersion) of each gene and store values (using a loop with iterations=number of simulations chosen [1000])
    #Should this be estimated for each gene and for each timepoint individually by using replicates? This is a problem if we have less than triplicate samples.
    #Alternatively, dispersion can be calculated like in DEseq, by using genes with similar count values.
    #Alternatively, dispersion can be calculated across timepoints for the same gene with the assumption that dispersion does not change over time, not sure this assumption holds however.

  #For each gene, randomly sample count number(?) from gene-specific NB using saved gene estimates.
  #I will choose the number of timepoints to be sampled (sequentially dropping the number of timepoints that are sampled).
  #Potentially, this can be accomplished as follows:
  #rnegbin(n, mu = n, theta = stop("'theta' must be specified")) [from the package MASS]
  #variance is calculated mu + mu^2/theta 
  #One issue is that I cannot input a dispersion value with this method.
  #Any other methods that would allow me to use dispersion?
  #What about rnbinom?

  #This should generate a fake matrix of counts
    #Loop through matrix and apply Dan's method of estimating degradation rate constants on simulated data.
    #Store estimated parameters.
    
#Finally, use some method on the simulated dataset to decide which is the smallest number of timepoints that can be used before degradation constant estimations are compromised.

#End

##################################
#Lingering Question:
  #What type of input should be used and what type of simulated data should be generated, i.e counts or normalized counts? (I am leaning towards normalized counts).
  #Does labeling + RNA-seq also have a NB distribution?
  #How should I quantify whether a timepoint number is compromising the estimation of degradation rate constants?
