# Model description

## Biomass calculation

I assume that the community follows Damuth power law for density-size relationship N = aM^b [1] with b=-3/4. Each time step I calculate the average biomass using the actual densities with the inverse function M=(N/a)^(-4/3). The second assumption is that species-size distribution follows a universal function, thus I calculated the biomass for each individual using a lognormal distribution with the previously calculated average [2].


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

8.	*de*		: Output species density, richness, and Shannon index. 

9.	*sa*		: Save state of the model in sed format each [inter] steps from [init] time  

10.	*baseName*	: Base file name for output.

11. *idrPal*	: Color palette for graphics (Not working)

13.	*mfDim* 	: Multifractal spectrum calculation (S/N), to specifiy the q uses a file named q.sed.
		
14. *minBox* 	: Min box size for Multifractal spectrum

15. *maxBox*	: Max box size for Multifractal spectrum.

16.	*deltaBox*	: The program uses powers of 2 intervals, but if you want to restric the number of intervals you can use this parameter.
		
8.	*pomac*		: (S/N) Takes parameters line by line from a file named pomac.lin and does a special output coded in the function IPSNeutral::PrintPomac.

9.	*bioCalc*	: Calculates biomass based in an inverse Damuth rule, M=aN^b  between a range of biomasses. The biomass is added to any other outputs requested.

10. *bioMax*	: Maximun biomass

11. *bioMin*	: Minimun biomass 

### File.sed

File used for initial conditions or to save the state of the system

The structure line by line is:

1.	DimX	DimY
2.	SP
3.	A matrix with dimension DimX x DimY with the position of species. If 0 there is no species.


### File.set

File to seed a number of individuals with random positions in the landscape. Each line have the following structure:

	specie age quantity 

Age is not used in this model.
	

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

The model can be compiled for simulating the transition between a competition Tilman's like model [3] and a neutral model [4], modifying the parameter ReplacementRate. This parameter only works for hierarchical models: type 2 and 3. You have to define the following 

The pre-processor variable HIERARCHICAL_CONT if defined permits the use of ReplacementRate, if not defined the model assign ReplacementRate=1 for Hierarchical models, this slightly faster because uses less random numbers. 
  

## References

1. White EP, Ernest SKM, Kerkhoff AJ, Enquist BJ (2007) Relationships between body size and abundance in ecology. Trends Ecol Evol 22: 323–330. doi:10.1016/j.tree.2007.03.007.

2. Giometto A, Altermatt F, Carrara F, Maritan A, Rinaldo A (2013) Scaling body size fluctuations. Proc Natl Acad Sci: 201301552. doi:10.1073/pnas.1301552110.

3. Tilman D (1994) Competition and biodiversity in spatially structured habitats. Ecology 75: 2–16.

4. Hubbell SP (2001) The unified neutral theory of biodiversity and biogeography. Princeton University Press. 