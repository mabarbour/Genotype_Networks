Plant genetic variation shapes the structure and function of herbivore-parasitoid networks (a work in progress)
========================================================

Matthew A. Barbour, Jordi Bascompte, and Greg W. Crutsinger

### Objective

This paper examines how plant genetic variation influences herbivore-parasitoid interaction networks and how these networks shape selection pressures acting on herbivores.


### Research Questions

**(1)** How does plant genotype influence gall-parasitoid interaction networks?

**(2)** How do different interaction networks influence the form of natural selection acting on herbivore traits?

*Possible (3)* How does the loss of plant genotypes influence gall-parasitoid interaction networks?

### Methods

#### Q1
Build a meta-web of herbivore-parasitoid interactions by pooling across all genotypes (try a both qualitative and quantitative version, hopefully both will yield the same results)

Run an appropriate modularity algorithm to identify modules on the meta-web.

***Figure out how to quantify the frequency of representation of the different modules.***

Run a chi-square test based on the frequencies of each module occurring on each genotype to see whether there is an association.

### Q2
Build a meta-web of herbivore-parasitoid interactions by pooling across all genotypes (try both qualitative and quantitative version, hopefully both will yield the same results)

Evaluate the structure of the meta-web (nested or modular or both).

Compute the appropriate network measure for each replicate of each genotype and see whether I can detect any correlations between network structure and genetic distance or traits.

Could also use the method proposed in Poisot 2012 (Dissimilarity of Interaction Networks), but it doesn't examine how network structure changes (this may or may not be important, depending on whether I see any noticeable changes in my graphs of the different networks)

Could also see how the functional role of species (Guimera and Amaral Functional Cartography paper in Nature) change across genotypes.  Test this using a chi-square analysis based on genotypes and the seven different functional roles (I'm thinking now I won't see a big difference...)

***What I think may be a promising method***: Look for dissimilarity in module participation, by creating a matrix of all possible interactions (e.g H1,P1 H1,P2 H2,P1 etc.) for a particular module on a tree (or genotype), then calculate a bray-curtis dissimilarity within modules across the trees/genotypes.  Then see if this is associated with any particular gall and/or plant traits.**this is a more likely area where we will see gall trait variation influencing the network**

***Conversation with Rudolf***: Look for mantel correlation between dissimilarity in herbivore traits (e.g. gall size, multilocularity, etc.) and dissimilarity in network topology (e.g. see Bersier et al. 2002, for quantitative versions of network measures).  And see whether this correlation varies across genotypes.  Or simply look to see whether the quantitative versions of these network metrics vary across genotypes (similar to Tylianakis et al. 2007, but with different genotypes instead of habitats)

### Q3
Option 1: Use partial mantel tests and path analysis to see how genetic distance and trait dissimilarity influence dissimilarity in module frequency (either all modules combined, or modules individually).  Consider at first using all plant traits and then looking at individual traits.

Also, look for patterns between genotypes and herbivore traits, as well as herbivore traits and network topology.  Could use a mantel test for genetic distance and herbivore trait dissimilarity and measures of network topology.







