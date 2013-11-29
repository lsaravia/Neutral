# Model description

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

9.	*bioCalc*	: Calculates biomass based in an inverse Damuth rule, M=aN^(-4/3)  between a range of biomasses. The biomass is added to any other outputs requested.

10. *bioMax*	: Maximun biomass

11. *bioMin*	: Minimun biomass 

### File.sed

sed: archivo de condiciones iniciales -> Posiciones espaciales, especie, edad, Habitat.

DimX, DimY	: dimensiones del lattice

La estructura es la siguiente 

	DimX	DimY
	SP
	Matriz de numeros de DIMX x DIMY con la posicion de las especies


### File.set

archivo de condiciones iniciales -> siembra al azar individuos de cierta especie y edad
esta compuesto de lineas con la siguiente estructura

	especie edad cantidad [posicion x maxima] [habitat]
	

### FILE.**inp**

Parametros del modelo y frecuencia de especies en la comunidad
Con la siguiente estructura de lineas

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

There are 3 diferent kinds of dispersal functions set at compile time using a macro definition in the makefile, for the following formulas: dd = DispersalDistance

EXP_DISP defines an exponential distribution f(x) = dd*exp(-dd*x) with mean = 1/dd

UNIFORM_DISP defines a uniform distribution between 0 and dd with mean =dd/2

POWER_DISP defines an inverse power distribution f(x) = (dd-1)*pow(x,-dd) with mean = dd-1/(dd-2)
  

