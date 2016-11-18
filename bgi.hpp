#ifndef __BGI_HPP
#define __BGI_HPP

void IGraph2(int dimX,int dimY,char * idrPal=NULL);
void PPix1(int x,int y, int color);
void PPix2(int x,int y, int color);

void IGraph(int dimX,int dimY,char * idrPal=NULL);
void PPix(int x,int y, int color, int type=0);
void EGraph();

void GLabel(char * buff,int color);
void BLabel();

void SetRGBdefaults();
void SetRGBpalette(int color, int red, int green, int blue);
void ReadIdrisiPalette(char * fileName);

#endif

