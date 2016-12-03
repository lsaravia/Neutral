# Model description


## Classical Neutral Model Description 

This model is a stochastic cellular automata (CA) or also called interactive particle system [1]. In these kind of models space is discretized into a grid and only one individual can occupy a particular position. Each position represents an area fixed by the investigator to represent the real system. 

As this is a neutral model all individuals have the same parameters, besides they should belong to different species [4]. The only difference between species is that they have a different frequency in the metacommunity and also can have different abundances in the local community. The size of the community is given by *J = dimX* x *dimY*, where *dimX* and *dimY* are the dimension of the grid. Thus we have a maximum number of individuals in a given area. 

There are four processes included in the model: death, local dispersal, and migration, starting with an empty site the following events can happen:

(1) When the grid is not full, with probability *m* an individual of a species *i* can migrate from metacommunity at a rate proportional to its frequency $X_i$ in the metacommunity.

(2) When the grid is not full individuals give birth with rate *BirthRate* to a new individual that disperse to the neighborhood with a dispersal kernel. If the place is already occupied the disperser is lost.

(3) Individuals die a rate $\mu$

(4) When an individual dies it is replaced by a migrant from metacommunity with probability *m* and with probability *1-m* or by an individual of the neighborhood. The neighborhood is established using the dispersal kernel with average distance *d*. Once the grid is full it keeps full, because when an individual dies is immediately replace by another. This is called the zero-sum assumption. 

The dispersal kernels available are:

1. Uniform: the probability is uniform up to a maximun distance. 

2. Exponential: $p(x) = \lambda e^{-\lambda x}$ with $mean=1/\lambda$ where $x\ge 0$.

3. Power law: $p(x) =  \frac{\alpha -1}{x_{min}} \left(\frac{x}{x_{min}} \right)^{-\alpha}$ with $mean =\frac{\alpha-1}{\alpha-2}x_{min}$ where $\alpha > 1$ and $x \ge x_{min}$. In all cases I used $x_{min} = 1$.

## Non saturated neutral model

The event (4) of the previous model is modified:

(4) When an individual dies it is not replaced by another so there are empty sites and the rules (1) and (2) are applied.


## Hierarchical competition model

In events (1), (2) and (4) if the migrant or disperser find and individual from other specie with number greater than itself the individual is replaced with rate *ReplacementRate*. In case the specie's number is less than the actual specie there is no replacement. 

## Biomass calculation

The output of the model is the species distribution but calculation of biomass can be turned on with the parameter *bioCalc*. I assume that the community follows Damuth power law for density-size relationship N = aM^b [1] with b calculated from the metacommunity abundances and a minimum and maximum biomass values as parameters. The average biomass is calculated using the inverse function M=(N/a)^(-4/3). I assume that species-size distribution follows a universal function, thus I calculated the biomass for each individual using a log-normal distribution with the previously calculated average [2].


## Command Line arguments 

ipsNeutral file.par file.inp file.sed/file.set

- All files are plain text. 
- The first two files are required.

- **file.par** have the simulation's parameters
- **file.inp** have the species's parameters

The other files **sed** and **inp** set initial conditions, they can not be used both at the same time.

## File structure

### file.par

It has the simulation's parameters and select the outputs the model will calculate.

1. *gr*		: Graphic output S/N

2. *nRuns*		: Number of runs with the same set of parameters

3. *nEvals*		: Number of steps in a run (End time)

4. *rndSeed*	: Seed to random number generator if 0 takes time.

5. *inter*		: Step interval for doing output graphic/data/spatial state

6. *init*		: Step for start data collection or graphics 

7. *modType*	: 1=Non saturated neutral model, 2=Neutral saturated (zero-sum) model, 
	 	3=Hierarchical non saturated, 4=Hierarchical saturated.

8.	*de*		: Output species density, richness, and Shannon index (natural logarithm). 

9.	*sa*		: Save state of the model in sed format each [inter] steps from [init] time  

10.	*baseName*	: Base file name for output.

11. *idrPal*	: Color palette for graphics (Not working)

13.	*mfDim* 	: Multifractal spectrum calculation (S/N), to specifiy the q uses a file named q.sed.
		
14. *minBox* 	: Min box size for Multifractal spectrum

15. *maxBox*	: Max box size for Multifractal spectrum.

16.	*deltaBox*	: The program uses powers of 2 intervals, but if you want to restric the number of intervals you can use this parameter.
		
8.	*pomac*		: (S/N) Takes parameters line by line from a file named pomacFile.lin and does a special output coded in the function IPSNeutral::PrintPomac.

9. *pomacFile*	: name fo the file to do multiple simulations

9.	*minProp*	: If >0 calculates H and richness for the species that have a proportion in density greater than *minProp*, if = 0 doesn't calculate additional H & richness. 

9.	*bioCalc*	: Calculates biomass based in an inverse Damuth rule, M=aN^b  between a range of biomasses. The biomass is added to any other outputs requested.

10. *bioMax*	: Maximun biomass

11. *bioMin*	: Minimun biomass 

12. *clusters*  : [S|A] S: Calculates the species with max cluster size, max cluster size, total mumber of clusters, total of species with max cluster, spanning species 
						A: Calculates the sizes of clusters of all species, the first record is the spanning species, if there is no spanning species is 0, 0
						
### File.sed

File used for initial conditions or to save the state of the system

The structure line by line is:

1.	DimX	DimY
2.	SP
3.	A matrix with dimension DimX x DimY with the position of species. If 0 there is no species.


### File.set

File to seed a number of individuals with random positions in the landscape. Each line have the following structure:

	specie age quantity 

Age is not used in this model so it can take any value.
	

### FILE.**inp**

Model parameters and frequency of species in the metacommunity.
With the following structure:

	DimX	DimY
	numSpecies
	specieParameters
	MetacommunityFrequency (numSpecies lines)
	
The structure of the line specieParameters is:

	sp	BirthRate	MortalityRate	DispersalDistance	ColonizationProb	[ReplacementRate]
	
	sp=0 represents the parameters of all species 
	DispersalDistance is the parameter of the dispersal kernel, only with uniform distribution
	represent the true maximal dispersal distance.

The structure of the lines MetacommunityFrequency is:

	sp	0	0	0	Freq
	
	sp  = 1.. numSpecies 
	Freq = Metacommunity Frequency 


	
## Dispersal functions

There are 3 different kinds of dispersal functions set at compile time using pre-processor macro definitions in the makefile, for the following formulas: dd = DispersalDistance

EXP_DISP defines an exponential distribution f(x) = dd*exp(-dd*x) with mean = 1/dd

UNIFORM_DISP defines a uniform distribution between 0 and dd with mean =dd/2

POWER_DISP defines an inverse power distribution f(x) = (dd-1)*pow(x,-dd) with mean = dd-1/(dd-2)

## Hierarchical competition and neutral model

The model can be compiled for simulating the transition between a competition Tilman's like model [3] and a neutral model [4], modifying the parameter *ReplacementRate*. This parameter only works for hierarchical models: type 3 and 4. 
The pre-processor variable HIERARCHICAL_CONT if defined permits the use of ReplacementRate, if not defined the model assign ReplacementRate=1 for Hierarchical models, this slightly faster because uses less random number calls. 
  

## References

1. White EP, Ernest SKM, Kerkhoff AJ, Enquist BJ (2007) Relationships between body size and abundance in ecology. Trends Ecol Evol 22: 323–330. doi:10.1016/j.tree.2007.03.007.

2. Giometto A, Altermatt F, Carrara F, Maritan A, Rinaldo A (2013) Scaling body size fluctuations. Proc Natl Acad Sci: 201301552. doi:10.1073/pnas.1301552110.

3. Tilman D (1994) Competition and biodiversity in spatially structured habitats. Ecology 75: 2–16.

4. Hubbell SP (2001) The unified neutral theory of biodiversity and biogeography. Princeton University Press. 

5. Durrett R, Levin SA (1994) Stochastic spatial models: a user’s guide to ecological aplications. Philosophical transactions of the Royal Society of London Series B 343: 329–350.

6. Marco DE, Montemurro MA, Cannas SA (2011) Comparing short and long-distance dispersal: modelling and field case studies. Ecography 34: 671–682. doi:10.1111/j.1600-0587.2010.06477.x.
