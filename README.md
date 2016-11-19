
# Multiespecific neutral/hierachical stochastic spatial model

This is C++ code for neutral/hierarchical competition models. In these all species have the same parameters except the probability of colonization from metacommunity. More details are given in the file ModelDescription.md

There are actually 5 similar models built in

1. Non-saturated neutral (Etienne 2007)
2. Saturated Neutral 
3. Non-saturated Hierarchical (Tilman 1994)
4. Saturated Hierarchical 
5. Continuum model between 2 and 4 (or 1 and 3)

And three dispersal kernels

1. Uniform
2. Exponential
3. Inverse power 

It was compiled using Gnu C++ (g++ v5.3.0) to compile it. For a crude and slow graphical output you will need the X11 SDL development libraries and the GRX graphics library http://grx.gnu.de/grx246um.htm

You also need the code from: 
	
1. <https://github.com/lsaravia/mfsba> if you want the multifractal spectra output otherwise you have to comment the calls to the function *MultifractalSBA* and the includes to *mf.h* header.
2. <https://github.com/lsaravia/Clusters> to calculate cluster statistics. 

To compile under Windows there is a separate makefileWin.mak, you might have to modify the folders of the *mfsba* and *Clusters* sources.


## Citation

If you use the model please cite:

1. Saravia, L.A. (2015). A new method to analyse species abundances in space using generalized dimensions. Methods Ecol. Evol., 6, 1298–1310

2. Saravia L.A. (2014) Neutral: A sofwtware to simulate spatially explicit neutral/hierarchical models.  <http://dx.doi.org/10.6084/m9.figshare.969692>.

## Bibliography

Tilman, D., 1994. Competition and biodiversity in spatially structured habitats. Ecology, 75(1), pp.2-16.

Etienne, R.S., Alonso, D. & McKane, A.J., 2007. The zero-sum assumption in neutral biodiversity theory. Journal of Theoretical Biology, 248(3), pp.522-536.

## License

	Copyright 2011 Leonardo A. Saravia
 
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
