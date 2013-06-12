# PARAMETROS de Linea de Commandos

ipsNeutral file.par file.inp file.sed/file.set

All files are plain text. The first two files are required.

**file.par** have the simulation's parameters
**file.inp** have the species's parameters

The other files **sed** and **inp** are for set initial conditions, they can not be used  both at the same time.

# file.par

It has the simulation's parameters and kind of outputs the model will calculate.

1. *gr*		: Graphic output S/N

2. *nRuns*		: número de corridas con este set de parametros

3. *nEvals*		: número de pasos de tiempo (tiempo maximo)

4. *rndSeed*		: semilla para el generador de nros al azar, si es 0 toma el time

5. *inter*		: intervalo en pasos para mostrar grafico y/o 
		  hacer salidad de datos y/o tomar informacion para promedios

6. *init*		: numero de pasos a partir del cual inicia la toma de datos o muestra en pantalla

7. *modType*	: 1=Modelo neutral no saturado, 2=Modelo neutral saturado, 
		 	3=Modelo jerarquico no saturado, 4=Modelo jerarquico saturado.

8	*de*		: imprime densidad de especies, Riqueza y diversidad de Shanonn. 

9.	*sa*		: Guarda el estado del modelo en formato sed a partir de [inicio] 
				y con la periodicidad especificada en [inter]

10.	*baseName*	: Nombre de Archivo para salvar archivos sed y demas salidas.

11. *idrPal*	: Paleta de colores idrisi utilizada en salida grafica (NO IMPLEMENTADO)

12.	*patchStat*	: Calculo de estadisticas de parche (NO IMPLEMENTADO)

13.	*mfDim* 	: Calculo de espectro multifractal Dq por box-counting toma un 
		archivo q.sed para los q del espectro completo.
		
14: *minBox* 	:  tamaño de box minimo
15: *maxBox*	:  tamaño de box máximo, este valor se modifica para que no sea mayor
		que la mitad del tamaño de la grilla DimX
16.	*deltaBox*	: cantidad de intervalos, tambien se modifica para que sean potencias
		de 2
		
17. *moranI*	: Calculo de indice de correlacion I de Moran (No Implementado)

18.	*pomac*		: S/N Toma los parametros de un archivo llamado pomac.lin y hace una 
	salida especial resumida


# ARCHIVO sed

sed: archivo de condiciones iniciales -> Posiciones espaciales, especie, edad, Habitat.

DimX, DimY	: dimensiones del lattice

La estructura es la siguiente 

	DimX	DimY
	SP
	Matriz de numeros de DIMX x DIMY con la posicion de las especies


# ARCHIVO set

archivo de condiciones iniciales -> siembra al azar individuos de cierta especie y edad
esta compuesto de lineas con la siguiente estructura

	especie edad cantidad [posicion x maxima] [habitat]
	

# FILE **inp**

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


	
# Dispersal functions

There are 3 diferent kinds of dispersal functions set at compile time using a macro definition in the makefile, for the following formulas: dd = DispersalDistance

EXP_DISP defines an exponential distribution f(x) = dd*exp(-dd*x) with mean = 1/dd

UNIFORM_DISP defines a uniform distribution between 0 and dd with mean =dd/2

POWER_DISP defines an inverse power distribution f(x) = (dd-1)*pow(x,-dd) with mean = dd-1/(dd-2)
  
The mean is not a good measure to the size of the clusters that form in the model.
