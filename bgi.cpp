#include <stdlib.h>

#ifdef GRAPHICS

#ifdef __GNUC__

#include <grx20.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "bgi.hpp"


static int Maxx,Maxy,Lenx,Leny;
GrContext *w3,*w1,*w2;
static char IdrPalette[40];
static GrColor egacolors[256];

void IGraph(int dimX, int dimY,char * idrPal)
	{

	GrSetMode(GR_default_graphics);
//	GrSetMode(GR_width_height_color_graphics,640,480,256);
//	GrSetMode(GR_width_height_color_graphics,800,600,256);
//	GrSetMode(GR_biggest_graphics);



    if( idrPal!=NULL )
		strcpy(IdrPalette,idrPal);

	int  x = GrSizeX();
	int  y = GrSizeY();
	Maxx = GrMaxX();
	Maxy = GrMaxY();

	GrSetRGBcolorMode();
	SetRGBdefaults();
	
  	GrBox(0,0,x-1,y-20,GrWhite());

	w1 = GrCreateSubContext(1,1,x-2,y-21,NULL,w1);

	w2 = GrCreateSubContext(0,y-19,x-1,y,NULL,NULL);

	GrSetContext(w1);
	GrSetUserWindow(0,0,dimX,dimY);
	//ReadIdrisiPalette(IdrPalette);
	GrSaveContext(w1);

	GrSetContext(w2);
	//ReadIdrisiPalette(IdrPalette);
	GrSaveContext(w2);


//	Maxx = GrMaxX();
//	Maxy = GrMaxY();
	Lenx = 1; //Maxx/dimX;
	Leny = 1; //Maxy/dimY;

	}

void IGraph2(int dimX, int dimY, char * idrPal)
	{

//	GrSetMode(GR_default_graphics);
//	GrSetMode(GR_width_height_color_graphics,640,480,256);
	GrSetMode(GR_width_height_color_graphics,800,600,256);

	int  x = GrSizeX();
	int  y = GrSizeY();
	Maxx = GrMaxX();
	Maxy = GrMaxY();

    if( idrPal!=NULL )
    	strcpy(IdrPalette,idrPal);

	//GrSetRGBcolorMode();
  	GrBox(0,0,x-1,y-20,GrWhite());

    int x0 = x/2-1;
    int x1 = x0+1;

	w1 = GrCreateSubContext(1,1,x0,y-21,NULL,NULL);
    GrSetContext(w1);
	GrSetUserWindow(0,0,dimX,dimY);
    ReadIdrisiPalette(IdrPalette);
	GrSaveContext(w1);
    
    w3 = GrCreateSubContext(x1,1,x-2,y-21,NULL,NULL);
    GrSetContext(w3);
	GrSetUserWindow(0,0,dimX,dimY);
    ReadIdrisiPalette(IdrPalette);
	GrSaveContext(w3);

	w2 = GrCreateSubContext(0,y-19,x-1,y,NULL,NULL);
    GrSetContext(w2);
    ReadIdrisiPalette(IdrPalette);
	GrSaveContext(w2);


//	Maxx = GrMaxX();
//	Maxy = GrMaxY();
	Lenx = 1; //Maxx/dimX;
	Leny = 1; //Maxy/dimY;
	}

void PPix1(int x, int y, int color)
{
	GrSetContext(w1);
	if( color==0)
		GrUsrFilledBox( x, y,x+Lenx,y+Leny, GrWhite() );
	else
		GrUsrFilledBox( x, y,x+Lenx,y+Leny, color+16 );
}

void PPix2(int x, int y, int color)
{
	GrSetContext(w3);
	if( color==0)
		GrUsrFilledBox( x, y,x+Lenx,y+Leny, GrWhite() );
	else
		GrUsrFilledBox( x, y,x+Lenx,y+Leny, color+16 );
}

void PPix(int x, int y, int color, int type)
	{
//   	int x1,y1,x2,y2;
	GrSetContext(w1);
//   	GrGetUserWindow(&x1,&y1,&x2,&y2);
	if( color==0)
		GrUsrFilledBox( x, y,x+Lenx,y+Leny, GrWhite() );
	else
		GrUsrFilledBox( x, y,x+Lenx,y+Leny, egacolors[color] );
/*
	int posx = Lenx * x;
	int posy = Leny * y;

	if( type==0 )
		{
		if( color==0)
			GrFilledBox( posx, posy, posx+Lenx, posy+Leny, GrWhite() );
		else
			GrFilledBox( posx, posy, posx+Lenx, posy+Leny, color+16);
		}
	else
		{
		if( color==0)
			GrFilledCircle( posx+Lenx/2, posy+Leny/2, Lenx/2, GrWhite() );
		else
			GrFilledCircle( posx+Lenx/2, posy+Leny/2, Lenx/2, color+16 );
		}
*/
	}

void GLabel(char * buff,int color)
	{
	GrSetContext(w2);
	int y = (GrSizeY() - 20) ;
	int x = 10 + ((color-1)%12) * 50;
	GrTextXY(x,y , buff,GrWhite(), 0);

//	GrTextXY(x,y , buff,GrWhite(), GrBlack());

//	if( color==0)
//		GrFilledBox( x+24  , y, x+38, y+16, GrWhite() );
//	else
		GrFilledBox( x+24  , y, x+38, y+16, egacolors[color] );
	}

void BLabel()
{
	GrSetContext(w2);
	int y = (GrSizeY() - 20);
	int x;
	int xMax = GrSizeX();
    int color=0;
	for(x=0;x<xMax;x+=4)
    {
		GrFilledBox( x , y, x+3, y+16, egacolors[color] );
    	color++;
        if(color>256)
        	break;
	}
}

void EGraph()
	{
	GrSetMode(GR_default_text);
	}


static unsigned char _dac_normal[256][3] = {
	{ 0,   0,   0   },      /* black */
	{ 0,   0,   170 },      /* blue */
	{ 191, 191, 222 },      /* Celeston */
	{ 0,   170, 0   },      /* green */	

	{ 85,  255, 85  },      /* light green */
	{ 139, 0,   170 },      /* magenta */
	{ 170, 85,  0   },      /* brown */
	{ 170, 21,  0 },        /* Rojo */	

	{ 237, 144,  131},      /* Naranja */
	{ 85,  85,  85  },      /* dark gray */
	{ 85,  85,  255 },      /* light blue */
	{ 85,  255, 255 },      /* light cyan */
	
	{ 255, 85,  85  },      /* light red */
	{ 255, 85,  255 },      /* light magenta */
	{ 255, 255, 85  },      /* yellow */
	{ 153, 153, 40 },       /* white */
	
	{ 0, 0, 0},     { 5, 5, 5},     { 8, 8, 8},     {11,11,11},
	{14,14,14},     {17,17,17},     {20,20,20},     {24,24,24},
	{28,28,28},     {32,32,32},     {36,36,36},     {40,40,40},
	{45,45,45},     {50,50,50},     {56,56,56},     {63,63,63},
	{ 0, 0,63},     {16, 0,63},     {31, 0,63},     {47, 0,63},
	{63, 0,63},     {63, 0,47},     {63, 0,31},     {63, 0,16},
	{63, 0, 0},     {63,16, 0},     {63,31, 0},     {63,47, 0},
	{63,63, 0},     {47,63, 0},     {31,63, 0},     {16,63, 0},
	{ 0,63, 0},     { 0,63,16},     { 0,63,31},     { 0,63,47},
	{ 0,63,63},     { 0,47,63},     { 0,31,63},     { 0,16,63},
	{31,31,63},     {39,31,63},     {47,31,63},     {55,31,63},
	{63,31,63},     {63,31,55},     {63,31,47},     {63,31,39},
	{63,31,31},     {63,39,31},     {63,47,31},     {63,55,31},
	{63,63,31},     {55,63,31},     {47,63,31},     {39,63,31},
	{31,63,31},     {31,63,39},     {31,63,47},     {31,63,55},
	{31,63,63},     {31,55,63},     {31,47,63},     {31,39,63},
	{45,45,63},     {49,45,63},     {54,45,63},     {58,45,63},
	{63,45,63},     {63,45,58},     {63,45,54},     {63,45,49},
	{63,45,45},     {63,49,45},     {63,54,45},     {63,58,45},
	{63,63,45},     {58,63,45},     {54,63,45},     {49,63,45},
	{45,63,45},     {45,63,49},     {45,63,54},     {45,63,58},
	{45,63,63},     {45,58,63},     {45,54,63},     {45,49,63},
	{ 0, 0,28},     { 7, 0,28},     {14, 0,28},     {21, 0,28},
	{28, 0,28},     {28, 0,21},     {28, 0,14},     {28, 0, 7},
	{28, 0, 0},     {28, 7, 0},     {28,14, 0},     {28,21, 0},
	{28,28, 0},     {21,28, 0},     {14,28, 0},     { 7,28, 0},
	{ 0,28, 0},     { 0,28, 7},     { 0,28,14},     { 0,28,21},
	{ 0,28,28},     { 0,21,28},     { 0,14,28},     { 0, 7,28},
	{14,14,28},     {17,14,28},     {21,14,28},     {24,14,28},
	{28,14,28},     {28,14,24},     {28,14,21},     {28,14,17},
	{28,14,14},     {28,17,14},     {28,21,14},     {28,24,14},
	{28,28,14},     {24,28,14},     {21,28,14},     {17,28,14},
	{14,28,14},     {14,28,17},     {14,28,21},     {14,28,24},
	{14,28,28},     {14,24,28},     {14,21,28},     {14,17,28},
	{20,20,28},     {22,20,28},     {24,20,28},     {26,20,28},
	{28,20,28},     {28,20,26},     {28,20,24},     {28,20,22},
	{28,20,20},     {28,22,20},     {28,24,20},     {28,26,20},
	{28,28,20},     {26,28,20},     {24,28,20},     {22,28,20},
	{20,28,20},     {20,28,22},     {20,28,24},     {20,28,26},
	{20,28,28},     {20,26,28},     {20,24,28},     {20,22,28},
	{ 0, 0,16},     { 4, 0,16},     { 8, 0,16},     {12, 0,16},
	{16, 0,16},     {16, 0,12},     {16, 0, 8},     {16, 0, 4},
	{16, 0, 0},     {16, 4, 0},     {16, 8, 0},     {16,12, 0},
	{16,16, 0},     {12,16, 0},     { 8,16, 0},     { 4,16, 0},
	{ 0,16, 0},     { 0,16, 4},     { 0,16, 8},     { 0,16,12},
	{ 0,16,16},     { 0,12,16},     { 0, 8,16},     { 0, 4,16},
	{ 8, 8,16},     {10, 8,16},     {12, 8,16},     {14, 8,16},
	{16, 8,16},     {16, 8,14},     {16, 8,12},     {16, 8,10},
	{16, 8, 8},     {16,10, 8},     {16,12, 8},     {16,14, 8},
	{16,16, 8},     {14,16, 8},     {12,16, 8},     {10,16, 8},
	{ 8,16, 8},     { 8,16,10},     { 8,16,12},     { 8,16,14},
	{ 8,16,16},     { 8,14,16},     { 8,12,16},     { 8,10,16},
	{11,11,16},     {12,11,16},     {13,11,16},     {15,11,16},
	{16,11,16},     {16,11,15},     {16,11,13},     {16,11,12},
	{16,11,11},     {16,12,11},     {16,13,11},     {16,15,11},
	{16,16,11},     {15,16,11},     {13,16,11},     {12,16,11},
	{11,16,11},     {11,16,12},     {11,16,13},     {11,16,15},
	{11,16,16},     {11,15,16},     {11,13,16},     {11,12,16},
	{21,21,21},     {21,21,63},     {21,63,21},     {21,63,63},
	{63,21,21},     {63,21,63},     {63,63,21},     {63,63,63}
};

void SetRGBdefaults()
{
  int i;

  for (i = 0; i < 256; ++i)
    egacolors[i]=GrAllocColor(_dac_normal[i][0],_dac_normal[i][1],_dac_normal[i][2]);
    //SetRGBpalette(i,_dac_normal[i][0],_dac_normal[i][1],_dac_normal[i][2]);
}

#define SIX(a)  ((a)&63)

void SetRGBpalette(int color, int red, int green, int blue)
	{
  	GrSetColor( color, SIX(red)<<2, SIX(green)<<2, SIX(blue)<<2);
	}

void ReadIdrisiPalette(char * fileName)
	{
	unsigned int k=0;
	char buf[20];
	FILE *in;
	unsigned char lee[3];

	in = fopen(fileName,"rb");
   	if( in==NULL )
   	{
		GrSetMode(GR_default_text);
		fprintf(stderr, "Cannot open palette file.\n");
		exit(1);
   	}
       
	fread(buf,18,1,in);
	int i=0;
	while (!feof(in))
	{
		fread(&lee,3,1,in);
		//_dac_normal[k][0] = lee[0];
		//_dac_normal[k][1] = lee[1];
		//_dac_normal[k][2] = lee[2];
        k = GrAllocCell();
        if( k==GrNOCOLOR )
        	break;

//	  	GrSetColor( k, lee[0], lee[1], lee[2]);
		egacolors[i++] = GrAllocColor( lee[0], lee[1], lee[2]);

	}

	fclose(in);
	//GrRefreshColors();
	}

#endif

#else
void IGraph2(int dimX,int dimY,char * idrPal=NULL){}
void PPix1(int x,int y, int color){}
void PPix2(int x,int y, int color){}

void IGraph(int dimX,int dimY,char * idrPal=NULL){}
void PPix(int x,int y, int color, int type=0){}
void EGraph(){}

void GLabel(char * buff,int color){}
void BLabel(){}

void SetRGBdefaults(){}
void SetRGBpalette(int color, int red, int green, int blue){}
void ReadIdrisiPalette(char * fileName){}
#endif


