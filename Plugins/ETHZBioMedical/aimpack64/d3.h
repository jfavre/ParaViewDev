/*	  
**  D3.h
*/	  

#ifndef _d3_h
#define _d3_h
#include "d2.h"

#ifdef __cplusplus
extern "C" {
#endif

/*	  
**  Constants
*/	  

#define    D3idNORM	    0
#define    D3idSUP	    (1<<24)
#define    D3idSUB	    (2<<24)

/*	  
**  Macros D3
*/	  

     /* Gives the file offset, given a 3D position (x,y,z) and 
        the dimension of the data set.                         */
#define D3Offset(pos,dim) \
      ( ((pos).z*(dim).y + (pos).y)*(dim).x + (pos).x ) 

#define D3Position(off,dim,pos){\
      (pos).x = (off)%(dim).x;\
      (pos).y = ((off)/(dim).x)%(dim).y;\
      (pos).z = (off)/((dim).x*(dim).y);\
}

#define D3Cell(dat,dim,pos)\
      ( *( (dat) + D3Offset(pos,dim) ) )

#define D3Bit8Offset(pos,dim)\
  ((((pos).z>>1)*((dim).y+1>>1)+((pos).y>>1))*((dim).x+1>>1)+((pos).x>>1))
#define D3Bit8Cell(dat,dim,pos)\
  ( ( *( (dat) + D3BinOffset(pos,dim) ) ) &\
    (((pos).z&1)<<2) +(((pos).y&1)<<1) + ((pos).x&1) )

#define D3BinIsSet(dat,offs)    ( *(dat + (offs>>3)) & (1<<(offs&7)) ) 
#define D3BinIsOff(dat,offs) 	( !D3BinIsSet(dat,offs) )
#define D3BinSet(dat,offs)	*(dat + (offs>>3)) |= (1<<(offs&7)) 
#define D3BinClear(dat,offs)	*(dat + (offs>>3)) &= (~(1<<(offs&7))) 


#define D3Round(v,u) { (u).x = D1Round((v).x); (u).y = D1Round((v).y); (u).z = D1Round((v).z); }
#define D3fTrunc(v) { D1fTrunc((v).x); D1fTrunc((v).y); D1fTrunc((v).z); }

#define D3Min(v)  Min3( (v).x, (v).y, (v).z )
#define D3Max(v)  Max3( (v).x, (v).y, (v).z )

#define D3MaxFabs(v)  Max3( fabs((v).x), fabs((v).y), fabs((v).z) )

#define D3PrintVector(v) printf("( %6.3f , %6.3f , %6.3f )\n",(v).x,(v).y,(v).z); 
#define D3PrintfVector(v,str,mess){ \
		printf("%s:\n",mess); \
		printf(str,(v).x); printf(str,(v).y); printf(str,(v).z); \
		printf("\n");}

#define D3iPrintf(v,mess){ \
		printf("%10s: ",mess); \
		printf("%5Ld %5Ld %5Ld",(v).x,(v).y,(v).z); \
		printf("\n");}

#define D3p(v,mess) D3iPrintf(v,mess)

#define D3fPrintf(v,mess){ \
		printf("%10s: ",mess); \
		printf("%6.3f %6.3f %6.3f\n",(v).x,(v).y,(v).z);}

#define D3pf(v,mess) D3fPrintf(v,mess)

//#define D3iPrintLog(v,name)\
//  printf("%-30s   %7Ld %7Ld %7Ld\n",name,(v).x,(v).y,(v).z);

// works on windows, macro above does not work correctly on windows platform
// unclear why?
#define D3iPrintLog(v,name)\
  printf("%-30s   %7Ld %7Ld %7Ld\n",name,(int)(v).x,(int)(v).y,(int)(v).z);

#define D3fPrintLog(v,name)\
  printf("%-30s   %7.4f %7.4f %7.4f\n",name,(v).x,(v).y,(v).z);

#define D3ComponentProc(v,func,u) {\
  (u).x = (func)((v).x); \
  (u).y = (func)((v).y); \
  (u).z = (func)((v).z); \
}

#define D3Set(v,a,b,c)	    {(v).x=(a); (v).y=(b); (v).z=(c);} 
#define D3SetAll(v,a)	    (v).x = (v).y = (v).z = (a) 

#define D3Convert(v,u,type) {(u).x = (type) (v).x;  \
                             (u).y = (type) (v).y;  \
                             (u).z = (type) (v).z;  }

#define D3SetVector(v,a,b,c) {(v).x=(a); (v).y=(b); (v).z=(c);} 
#define D3SetfVector(v,a,b,c) {(v).x=(float)(a); (v).y=(float)(b); (v).z=(float)(c);} 
#define	D3SetNormVector(v,a,b,c) {(v).z  = 1.0/sqrt((a)*(a)+(b)*(b)+(c)*(c));\
			         (v).x  = (a)*(v).z;  	\
			         (v).y  = (b)*(v).z;	\
			         (v).z *= (c); 		} 

#define D3Prod(v)		((v).x * (v).y * (v).z)
#define D3Sum(v)		((v).x + (v).y + (v).z)

#define D3VProd(a,b,vprod)	{ (vprod).x = (a).y*(b).z - (a).z*(b).y;\
				  (vprod).y = (a).z*(b).x - (a).x*(b).z;\
				  (vprod).z = (a).x*(b).y - (a).y*(b).x;}

#define D3SProd(a,b) 		( (a).x*(b).x + (a).y*(b).y + (a).z*(b).z )
#define D3AbsSqr(a)		D3SProd(a,a)
#define D3Abs2(a)		D3SProd(a,a)
#define D3Abs(a)		sqrt( (double) D3SProd(a,a))
#define D3Length(a)		D3Abs(a)
#define D3Length2(a)		D3Abs2(a)
#define D3LengthCB(a)		( abs((a).x) + abs((a).y) + abs((a).z) )
    /* CB: city block: sum component distance:  */
#define D3LengthMC(a)		Max3(abs((a).x),abs((a).y),abs((a).z))
    /* MC: max component: chess board distance */

#define D3Distance2(a,b)\
  (Sqr((a).x - (b).x) + Sqr((a).y - (b).y) + Sqr((a).z - (b).z))
#define D3Distance(a,b) sqrt(D3Distance2(a,b))

#define	D3Norm(v)	  { D1fTMP = D3Length(v);\
			    (v).x /= D1fTMP;	 \
			    (v).y /= D1fTMP;	 \
			    (v).z /= D1fTMP;	 } 

#define D3NormVector(v) D3Norm(v);

#define D3DirectionProd(v)\
( (((v).x==0)?1:(v).x) * (((v).y==0)?1:(v).y) * (((v).z==0)?1:(v).z) )

#define D3HalfSpaceZPos(v)\
(  (v).z >  0 ||\
 ( (v).z == 0 && ( (v).y > 0 || ((v).y == 0 && (v).x > 0) ) ) )

#define D3Angle(v,u)	   (acos((D3SProd(v,u)/D3Length(v)/D3Length(u))))
    /* between 0 pi */

#define D3AbsAngle(v,u)	   (acos(fabs((D3SProd(v,u)/D3Length(v)/D3Length(u)))))
    /* between 0 pi/2 */

#define D3PolarCoord(v,polar){\
  (polar).r     = D3Length(v);\
  (polar).phi   = atan2((v).y,(v).x);\
  (polar).theta = atan2(D2Length(v),(v).z);\
}

#define D3LinComb(a,v,b,u,lincmb) {(lincmb).x = (a)*(v).x + (b)*(u).x;  \
                                   (lincmb).y = (a)*(v).y + (b)*(u).y;  \
                                   (lincmb).z = (a)*(v).z + (b)*(u).z;  }

/* Comparison.    NB! D3Gt() !=  !D3Le(). Look at the corresp. octants */  

#define D3Eq(v,u) ( ((v).x == (u).x) && ((v).y == (u).y) && ((v).z == (u).z) )
#define D3Gt(v,u) ( ((v).x >  (u).x) && ((v).y >  (u).y) && ((v).z >  (u).z) )
#define D3Ge(v,u) ( ((v).x >= (u).x) && ((v).y >= (u).y) && ((v).z >= (u).z) )
#define D3Le(v,u) ( ((v).x <= (u).x) && ((v).y <= (u).y) && ((v).z <= (u).z) )
#define D3Lt(v,u) ( ((v).x <  (u).x) && ((v).y <  (u).y) && ((v).z <  (u).z) )

#define D3EqAnyComponent(v,u)(\
       ((v).x == (u).x) || ((v).y == (u).y) || ((v).z == (u).z) )

#define D3MaxComponents(v,u,w) {\
  (w).x = Max2((v).x,(u).x); \
  (w).y = Max2((v).y,(u).y); \
  (w).z = Max2((v).z,(u).z); \
}

#define D3MinComponents(v,u,w) {\
  (w).x = Min2((v).x,(u).x); \
  (w).y = Min2((v).y,(u).y); \
  (w).z = Min2((v).z,(u).z); \
}

#define D3Dec(v,dec)      {  (v).x -= (dec).x ;  \
                             (v).y -= (dec).y ;  \
                             (v).z -= (dec).z ;  }

#define D3Inc(v,inc)      {  (v).x += (inc).x ;  \
                             (v).y += (inc).y ;  \
                             (v).z += (inc).z ;  }

#define D3Scale(v,scale)  {  (v).x *= (scale).x;  \
                             (v).y *= (scale).y;  \
                             (v).z *= (scale).z;  }

#define D3CScale(v,factor) { (v).x *= (factor);  \
                             (v).y *= (factor);  \
                             (v).z *= (factor);  }


#define D3CShrink(v,factor){ (v).x /= (factor);  \
                             (v).y /= (factor);  \
                             (v).z /= (factor);  }

#define D3Sub(v,u,diff)     {(diff).x = (v).x - (u).x;  \
                             (diff).y = (v).y - (u).y;  \
                             (diff).z = (v).z - (u).z;  }

#define D3Add(v,u,sum)      {(sum).x = (v).x + (u).x;  \
                             (sum).y = (v).y + (u).y;  \
                             (sum).z = (v).z + (u).z;  }

#define D3Mul(v,u,prod)      {(prod).x = (v).x * (u).x;  \
                              (prod).y = (v).y * (u).y;  \
                              (prod).z = (v).z * (u).z;  }

#define D3MulComponents(v)   ( (v).x * (v).y * (v).z )


#define D3Div(v,u,quote)    {(quote).x = (v).x / (u).x;  \
                             (quote).y = (v).y / (u).y;  \
                             (quote).z = (v).z / (u).z;  }

#define D3Inv(v,u)	    {(u).x = 1.0 / (v).x;  \
                             (u).y = 1.0 / (v).y;  \
                             (u).z = 1.0 / (v).z;  }

#define D3Mod(v,u,rest)     {(rest).x = (v).x % (u).x;  \
                             (rest).y = (v).y % (u).y;  \
                             (rest).z = (v).z % (u).z;  }

#define D3ModExt(v,u,rest)  {(rest).x = ModExt((v).x,(u).x);  \
                             (rest).y = ModExt((v).y,(u).y);  \
                             (rest).z = ModExt((v).z,(u).z);  }

#define D3CMul(v,const,u)   {(u).x = (v).x*(const); \
                             (u).y = (v).y*(const); \
                             (u).z = (v).z*(const); }

#define D3FMul(v,fract,u)   {(u).x = ((v).x*(fract).num)/(fract).den; \
			     (u).y = ((v).y*(fract).num)/(fract).den; \
			     (u).z = ((v).z*(fract).num)/(fract).den; }

#define D3Chs(v,u) D3CMul(v,-1,u)

#define D3CDiv(v,const,u)   {(u).x = (v).x/(const);\
                             (u).y = (v).y/(const);\
                             (u).z = (v).z/(const);}

#define D3FDiv(v,fract,u)   {(u).x = ((v).x*(fract).den)/(fract).num; \
			     (u).y = ((v).y*(fract).den)/(fract).num; \
			     (u).z = ((v).z*(fract).den)/(fract).num; }

#define D3AddCMul(v,c,u,sum)      {(sum.x) = (v).x + (c)*(u).x;  \
                             	   (sum.y) = (v).y + (c)*(u).y;  \
                             	   (sum.z) = (v).z + (c)*(u).z;  }

#define D3AddCDiv(v,u,c,sum)      {(sum).x = (v).x + (u).x/(c);  \
                             	   (sum).y = (v).y + (u).y/(c);  \
                             	   (sum).z = (v).z + (u).z/(c);  }
      /* nb! arguments not same as D3AddCMul */


#define D3SubCMul(v,c,u,diff)      {(diff).x = (v).x - (c)*(u).x;  \
                              	    (diff).y = (v).y - (c)*(u).y;  \
                             	    (diff).z = (v).z - (c)*(u).z;  }

#define D3DecCMul(v,c,u)           {(v).x -= (c)*(u).x;  \
                              	    (v).y -= (c)*(u).y;  \
                             	    (v).z -= (c)*(u).z;  }

#define D3IncCMul(v,c,u)           {(v).x += (c)*(u).x;  \
                              	    (v).y += (c)*(u).y;  \
                             	    (v).z += (c)*(u).z;  }

#define	D3ModLength(v,len,u) {	D1fTMP = (len)/D3Length(v);	  \
				D3CMul(v,D1fTMP,u); }

#define	D3MulPow2(v, p2, u) {\
  (u).x = ( (p2).x >= 0 ) ? (v).x * (1<<(p2).x) : (v).x / (1<<(-(p2).x));\
  (u).y = ( (p2).y >= 0 ) ? (v).y * (1<<(p2).y) : (v).y / (1<<(-(p2).y));\
  (u).z = ( (p2).z >= 0 ) ? (v).z * (1<<(p2).z) : (v).z / (1<<(-(p2).z));\
}

#define	D3DivPow2(v, p2, u) {\
  (u).x = ( (p2).x >= 0 ) ? (v).x / (1<<(p2).x) : (v).x * (1<<(-(p2).x));\
  (u).y = ( (p2).y >= 0 ) ? (v).y / (1<<(p2).y) : (v).y * (1<<(-(p2).y));\
  (u).z = ( (p2).z >= 0 ) ? (v).z / (1<<(p2).z) : (v).z * (1<<(-(p2).z));\
}

#define	D3AdjustToPow2(v, u){\
  (u).x = (Fract(Log2((v).x)) != 0.0) ? 1<<(Trunc(Log2((v).x)+1)) : (v).x;\
  (u).y = (Fract(Log2((v).y)) != 0.0) ? 1<<(Trunc(Log2((v).y)+1)) : (v).y;\
  (u).z = (Fract(Log2((v).z)) != 0.0) ? 1<<(Trunc(Log2((v).z)+1)) : (v).z;\
}

#define D3PlaneProj(v,n,vn,p)    {(p).x = (v).x - (n).x*(vn);  \
                                  (p).y = (v).y - (n).y*(vn);  \
                                  (p).z = (v).z - (n).z*(vn);  }

#define D3PlaneProjection(v,n,p){\
    (p).x = (v).x - (n).x*D3SProd(v,n);\
    (p).y = (v).y - (n).y*D3SProd(v,n);\
    (p).z = (v).z - (n).z*D3SProd(v,n);}

#define D3MoveOrgin(v,orgin)     {(v).x = (v).x - (orgin).x;  \
                                  (v).y = (v).y - (orgin).y;  \
                                  (v).z = (v).z - (orgin).z;  }

/* nb! D3IsInside: max inclusive */ 

#define	D3IsInside(p,min,max) \
      ( (p).x >= (min).x && (p).y >= (min).y && (p).z >= (min).z && \
	(p).x <= (max).x && (p).y <= (max).y && (p).z <= (max).z    )

#define	D3InVol(p,min,max) \
      ( (p).x >= (min).x && (p).y >= (min).y && (p).z >= (min).z && \
	(p).x <  (max).x && (p).y <  (max).y && (p).z <  (max).z    )

#define	D3InVolOff(p,dim,off) \
      ( (p).z >= (off).z && (p).z < ((dim).z-(off).z) &&\
        (p).y >= (off).y && (p).y < ((dim).y-(off).y) &&\
        (p).x >= (off).x && (p).x < ((dim).x-(off).x) )

#define	D3IsValid(p,im) \
      ( (p).z >= (im).off.z && (p).z < ((im).dim.z-(im).off.z) &&\
        (p).y >= (im).off.y && (p).y < ((im).dim.y-(im).off.y) &&\
        (p).x >= (im).off.x && (p).x < ((im).dim.x-(im).off.x) )

#define D3StrSetNumerical(str,p){\
  (p).x = (str);\
  while ( strspn((p).x,NUMERICAL_CHAR_SET) == 0 ) (p).x++;\
  (p).y =  p.x + strspn((p).x,NUMERICAL_CHAR_SET);\
  while ( strspn((p).y,NUMERICAL_CHAR_SET) == 0 ) (p).y++;\
  (p).z =  p.y + strspn((p).y,NUMERICAL_CHAR_SET);\
  while ( strspn((p).z,NUMERICAL_CHAR_SET) == 0 ) (p).z++;\
}

#define D3Atoi(ap,ip){\
  (ip).x = atoi(ap.x);\
  (ip).y = atoi(ap.y);\
  (ip).z = atoi(ap.z);\
}

#define D3Atof(ap,fp){\
  (fp).x = atof(ap.x);\
  (fp).y = atof(ap.y);\
  (fp).z = atof(ap.z);\
}

/*	  
**  Macros D3v
*/	  

#define D3vMinIndex(v)\
((v)[0] < (v)[1] ? ( (v)[0] < (v)[2] ? 0 : 2 ) : ( (v)[1] < (v)[2] ? 1 : 2 ))

#define D3vMaxIndex(v)\
((v)[0] > (v)[1] ? ( (v)[0] > (v)[2] ? 0 : 2 ) : ( (v)[1] > (v)[2] ? 1 : 2 ))

#define D3vSetUnitVector(v,i){\
  (v)[i]         = 1;\
  (v)[((i)+1)%3] = 0;\
  (v)[((i)+2)%3] = 0;\
}

/*	  
**  Macros D3aim
*/	  

#define D3AnySetCell(im,pos,value)\
  switch( (im).type ) {\
    case D1Tfloat:\
      D3Cell((float *)(im).dat,(im).dim,pos) = (float) value;\
      break;\
    case D1Tshort:\
      D3Cell((short *)(im).dat,(im).dim,pos) = (short) value;\
      break;\
    case D1Tchar:\
      D3Cell((char *)(im).dat,(im).dim,pos) = (char) value;\
      break;\
    default:\
      break;\
  }

#define D3AnyGetCell(im,pos,value)\
  switch( (im).type ) {\
    case D1Tfloat:\
      value = (float) D3Cell((float *)(im).dat,(im).dim,pos);\
      break;\
    case D1Tshort:\
      value = (short) D3Cell((short *)(im).dat,(im).dim,pos);\
      break;\
    case D1Tchar:\
      value = (char) D3Cell((char *)(im).dat,(im).dim,pos);\
      break;\
    default:\
      break;\
  }

#define D3IsSuperior(im)\
  ( !D3Eq((im).supdim,D3iNULL) )

#define D3PrintLogInfoAny(im,name){\
  printf("!-------------------------------------------------------------------------------\n");\
  printf("%-30s%-50s\n","Volume",name);\
  printf("%-30s%20.1f\n","AIM Version", (float)(im).version/10.);\
  printf("!\n");\
  D3iPrintLog((im).dim,"dim");\
  D3iPrintLog((im).off,"off");\
  D3iPrintLog((im).pos,"pos");\
  D3fPrintLog((im).el_size_mm,"element size in mm");\
  printf("%-30s   %7.4f %7.4f %7.4f\n","physical size in mm",\
    ((im).el_size_mm.x*(im).dim.x),\
    ((im).el_size_mm.y*(im).dim.y),\
    ((im).el_size_mm.z*(im).dim.z));\
  if( D3IsSuperior(im) ){\
    printf("!\n");\
    D3iPrintLog((im).supdim,"supdim");\
    D3iPrintLog((im).suppos,"suppos");\
    D3iPrintLog((im).subdim,"subdim");\
    D3iPrintLog((im).testoff,"testoff");\
  }\
  printf("!\n");\
  switch( (im).type ) {\
    case D1Tshort:\
      printf("%-30s%20s\n","Type of data","short (16 bit)");\
      break;\
    case D1Tchar:\
      printf("%-30s%20s\n","Type of data","char   (8 bit)");\
      break;\
    default:\
      printf("%-30s%20d\n","Size of data in bytes",D1TSize((im).type));\
  }\
  printf("!-------------------------------------------------------------------------------\n");\
}

#define D3ip(im,name) D3PrintLogInfoAny(im,name)


#define D3INITIMAGE(im,val) \
  for( D1iTMP1=0; D1iTMP1<(im).dim.y*(im).dim.x*(im).dim.z; D1iTMP1++ ) \
    *((im).dat + D1iTMP1) = (val);


#define D3ImageToAnyCopy(im,any_im) {\
  (any_im).dat = (void *) (im).dat;\
  (any_im).dim = (im).dim;\
  (any_im).pos = (im).pos;\
  (any_im).off = (im).off;\
  (any_im).el_size_mm = (im).el_size_mm;\
}

#define D3ImageGeometryCopy(im,out_im) {\
  (out_im).pos = (im).pos;\
  (out_im).dim = (im).dim;\
  (out_im).off = (im).off;\
  (out_im).suppos  = (im).suppos;\
  (out_im).supdim  = (im).supdim;\
  (out_im).subdim  = (im).subdim;\
  (out_im).testoff = (im).testoff;\
  (out_im).el_size_mm = (im).el_size_mm;\
}

#define D3InitAnyImage(any_im,typ) {\
  (any_im).version    = 030;\
  (any_im).dat	      = 0;\
  (any_im).id	      = 0;\
  (any_im).ref	      = 0;\
  (any_im).type	      = typ;\
  (any_im).proc_log   = 0;\
  (any_im).assoc.dat  = 0;\
  (any_im).assoc.nr   = 0;\
  (any_im).assoc.size = 0;\
  D3InitAnyImageGeometry(any_im);\
}

#define D3InitAnyImageGeometry(any_im) {\
  (any_im).pos	      = D3iNULL;\
  (any_im).dim	      = D3iNULL;\
  (any_im).off	      = D3iNULL;\
  (any_im).suppos     = D3iNULL;\
  (any_im).supdim     = D3iNULL;\
  (any_im).subdim     = D3iNULL;\
  (any_im).testoff    = D3iNULL;\
  (any_im).el_size_mm = D3fNULL;\
}

#define D3InitAnyFileImage020(any_im,typ) {\
  (any_im).version    = 020;\
  (any_im).id	      = 0;\
  (any_im).ref	      = 0;\
  (any_im).type	      = typ;\
  (any_im).assoc.nr   = 0;\
  (any_im).assoc.size = 0;\
  D3InitAnyFileImageGeometry020(any_im);\
}

#define D3InitAnyFileImageGeometry020(any_im) {\
  (any_im).pos	      = D3iNULL32;\
  (any_im).dim	      = D3iNULL32;\
  (any_im).off	      = D3iNULL32;\
  (any_im).suppos     = D3iNULL32;\
  (any_im).supdim     = D3iNULL32;\
  (any_im).subdim     = D3iNULL32;\
  (any_im).testoff    = D3iNULL32;\
  (any_im).el_size_mm = D3fNULL;\
}


#define D3AnyClearSup(any_im) {\
  (any_im).suppos     = D3iNULL;\
  (any_im).supdim     = D3iNULL;\
  (any_im).subdim     = D3iNULL;\
  (any_im).testoff    = D3iNULL;\
}

#define D3InitAny(any_im,typ) D3InitAnyImage(any_im,typ)

#define D3ReInitAny(any_im,typ){\
  D3FreeAnyImage(any_im);\
  D3InitAnyImage(any_im,typ);\
}


#define D3InitAnyImage011(any_im,typ) {\
  (any_im).dat	      = 0;\
  (any_im).id	      = 0;\
  (any_im).ref	      = 0;\
  (any_im).type	      = typ;\
  (any_im).proc_log   = 0;\
  (any_im).assoc.dat  = 0;\
  (any_im).assoc.size = 0;\
  (any_im).dim	      = D3iNULL32;\
  (any_im).off	      = D3iNULL32;\
  (any_im).subdim     = D3iONE32;\
  (any_im).pos	      = D3iNULL32;\
  (any_im).el_size_mm = D3fNULL;\
  (any_im).version    = 011;\
}

#define D3InitAny011(any_im,typ) D3InitAnyImage011(any_im,typ)

#define D3InitAny010(any_im,typ) {\
  (any_im).dat	      = 0;\
  (any_im).id	      = 0;\
  (any_im).ref	      = 0;\
  (any_im).type	      = typ;\
  (any_im).proc_log   = 0;\
  (any_im).assoc.dat  = 0;\
  (any_im).assoc.size = 0;\
  (any_im).dim	      = D3iNULL32;\
  (any_im).off	      = D3iNULL32;\
  (any_im).subdim     = D3iONE32;\
  (any_im).pos	      = D3iNULL32;\
  (any_im).el_size_mm = 0;\
}

#define D3FreeAnyImage(any_im) {\
  FREE((any_im).dat);\
  Free((any_im).proc_log);\
  Free((any_im).assoc.dat);\
}

#define D3InitSupImage(sup_im,im,sdim,tdim,type){\
\
  D3InitAny(sup_im,type);\
\
  (sup_im).subdim   = sdim;\
\
  (sup_im).dim.x = ((im).dim.x - 2*(im).off.x - ((tdim).x-(sdim).x))/(sdim).x;\
  (sup_im).dim.y = ((im).dim.y - 2*(im).off.y - ((tdim).y-(sdim).y))/(sdim).y;\
  (sup_im).dim.z = ((im).dim.z - 2*(im).off.z - ((tdim).z-(sdim).z))/(sdim).z;\
\
  (sup_im).pos.x = (im).pos.x + (im).off.x + ((tdim).x-(sdim).x)/2;\
  (sup_im).pos.y = (im).pos.y + (im).off.y + ((tdim).y-(sdim).y)/2;\
  (sup_im).pos.z = (im).pos.z + (im).off.z + ((tdim).z-(sdim).z)/2;\
\
  (sup_im).el_size_mm = (im).el_size_mm;\
}

#define D3InitAnySup(im,sup_im,type){\
  D3InitAny(sup_im,type);\
  D3ImageGeometryCopy(im,sup_im);\
  (sup_im).id     = D3idSUP;\
  (sup_im).supdim = D3iONE;\
  (sup_im).suppos = (im).off;\
  D3SubCMul((im).dim,2,(im).off,(sup_im).subdim);\
}

#define D3InitAnySub(sup_im,sub_im,type){\
  D3InitAny(sub_im,type);\
  D3ImageGeometryCopy(sup_im,sub_im);\
  (sup_im).id     = D3idSUB;\
}


#define D3INITVOL(vol)	\
{			\
  (vol).dat = 0;	\
  D3Set(vol.dim,0,0,0);	\
  D3Set(vol.off,0,0,0);	\
  D3Set(vol.pos,0,0,0);	\
} 

#define D3TypeName(im) D1TName((im).type)

  /* #define D3Volume(im) ((im).dim.x*(im).dim.y*(im).dim.z) */

#define D3Volume(im)\
  ((im).id == D3idSUP ?\
     ((im).supdim.x*(im).supdim.y*(im).supdim.z) :\
   (im).id == D3idSUB ?\
     ((im).supdim.x*(im).supdim.y*(im).supdim.z*\
      (im).subdim.x*(im).subdim.y*(im).subdim.z) :\
     ((im).dim.x*(im).dim.y*(im).dim.z) )
      

#define D3PhysVoxVolume(im)\
	((im).el_size_mm.x*(im).el_size_mm.y*(im).el_size_mm.z)

#define D3PhysVolume(im) (D3Volume(im)*D3PhysVoxVolume(im))

#define D3ElementSize(im) D1TSize((im).type)

#define D3MemorySize(im) (\
  (im).type == D3Tbit8 ?\
    ((((im).dim.x+1)>>1)*(((im).dim.y+1)>>1)*(((im).dim.z+1)>>1)* \
    D3ElementSize(im) + 1) :\
  (im).version == 020 && ((im).type == D1TcharCmp || (im).type == D1TbinCmp || (im).type == D1TcharCmp2 ) ?\
    ((int64) (* ((int32 *) (im).dat))) :\
  (im).type == D1TcharCmp || (im).type == D1TbinCmp || (im).type == D1TcharCmp2 ?\
    ((int64) (* ((int64 *) (im).dat))) :\
  (im).type == D1TZBencoded ?\
  		((int64) (* ((int64 *) (im).dat))) :\
  (D3Volume(im)*D3ElementSize(im)) )

#define D3MemorySizeUncompress(im) (\
  (D3Volume(im)*D3ElementSize(im)) )

#ifdef USING_VLM_MEMORY
/* it could be a small volume at a position very far off --> check individually */
#define D3Try020FileAim(im) {\
  if (D3MemorySize(im) < VLM_2GB \
    && D3Max((im).dim) < VLM_2GB \
    && D3Max((im).pos) < VLM_2GB \
    && D3Max((im).off) < VLM_2GB \
    ) (im).version = 020;\
  else {\
    printf("!%% Writing AimVersion020 requested... ");\
    printf("    not possible because of data-size or dimensions > 2 GB.\n");\
    printf("!%% Writing current AimVersion...\n ");\
  }\
}
#endif

#define D3MallocAny(im)\
{\
    MALLOC((im).dat, D3MemorySize(im), char);\
}

#define D3ReMallocAny(im)\
{\
  FREE((im).dat);\
  D3MallocAny(im);\
}

#define D3ReInitAnyImageMalloc(im, template_im, new_type)\
{\
  if( (im).dat != 0 &&\
      D3MemorySize(im) == (D3Prod((template_im).dim)*D1TSize(new_type)) ){\
    D3ImageGeometryCopy(template_im, im);\
    (im).type = new_type;\
  } else {\
    D3ImageGeometryCopy(template_im, im);\
    (im).type = new_type;\
    D3ReMallocAny(im);\
  }\
}

#define D3ReInitAnyImageCalloc(im, template_im, new_type)\
{\
  if( (im).dat != 0 &&\
      D3MemorySize(im) == (D3Prod((template_im).dim)*D1TSize(new_type)) ){\
    D3ImageGeometryCopy(template_im, im);\
    (im).type = new_type;\
    MEMSET((im).dat, 0, D3MemorySize(im));\
  } else {\
    D3ImageGeometryCopy(template_im, im);\
    (im).type = new_type;\
    D3ReCallocAny(im);\
  }\
}

#define D3CallocAny(im)\
{\
  D3MallocAny(im);\
  if ((im).dat != 0) MEMSET((im).dat, 0, D3MemorySize(im));\
}

#define D3ReCallocAny(im)\
{\
  FREE((im).dat);\
  D3CallocAny(im);\
}

/*	  
**  Macros D3Image
*/	  

#define D3Malloc(im,typ)						       \
{\
    MALLOC((im).dat, (im).dim.x*(im).dim.y*(im).dim.z, typ);\
}

#define D3ReMalloc(im,typ)						       \
{\
    REMALLOC((im).dat, (im).dim.x*(im).dim.y*(im).dim.z, typ);\
}

#define D3Calloc(im,typ)						       \
{\
    CALLOC((im).dat, (im).dim.x*(im).dim.y*(im).dim.z, typ);\
}

#define D3ReCalloc(im,typ)						       \
{\
    RECALLOC((im).dat, (im).dim.x*(im).dim.y*(im).dim.z, typ);\
}


#define D3CopyImage(im_in,im_out) {\
  switch( (im_in).type ) {\
    case D1Tfloat:\
      MEMCPY((im_out).dat, (im_in).dat, \
		    (im_in).dim.x*(im_in).dim.y*(im_in).dim.z*sizeof(D1float));\
      break;\
    case D1Tshort:\
      MEMCPY((im_out).dat, (im_in).dat, \
		    (im_in).dim.x*(im_in).dim.y*(im_in).dim.z*sizeof(D1short));\
      break;\
    case D1Tchar:\
      MEMCPY((im_out).dat, (im_in).dat, \
		    (im_in).dim.x*(im_in).dim.y*(im_in).dim.z*sizeof(D1char));\
      break;\
  }\
}

#define D3CopyAnyImage(im_in,im_out) {\
  if ((im_out).dat) FREE((im_out).dat);\
  if((im_out).assoc.dat != 0 && (im_out).assoc.nr != 0 ){\
    Free((im_out).assoc.dat);\
    (im_out).assoc.nr = 0;\
  }\
  if ((im_out).proc_log) Free((im_out).proc_log);\
  (im_out) = (im_in);\
  if(D3MemorySize(im_in)){\
    D3MallocAny(im_out);\
    MEMCPY((im_out).dat,(im_in).dat,D3MemorySize(im_in));\
  }\
  if ((im_in).proc_log){\
    Malloc((im_out).proc_log,strlen((im_in).proc_log)+1,char);\
    aimpack_memcpy((im_out).proc_log,(im_in).proc_log,strlen((im_in).proc_log)+1);\
  }\
  if((im_in).assoc.dat != 0 && (im_in).assoc.nr != 0 ){\
    Malloc((im_in).assoc.dat, (im_in).assoc.nr*(im_in).assoc.size, char);\
    aimpack_memcpy((im_out).assoc.dat,(im_in).assoc.dat,\
	 (im_in).assoc.nr*(im_in).assoc.size);\
  }\
}

#define D3CopyAssocData(im_in,im_out) {\
  if((im_in).assoc.dat != 0 && (im_in).assoc.nr != 0 ){\
    Malloc((im_out).assoc.dat, (im_in).assoc.nr*(im_in).assoc.size, char);\
    aimpack_memcpy((im_out).assoc.dat,(im_in).assoc.dat,\
	 (im_in).assoc.nr*(im_in).assoc.size);\
  }\
}

/*	  
**  D3 Types
*/	  

typedef struct { 
	D1char			x,y,z;
} D3char;

typedef struct { 
	D1char			*x,*y,*z;
} D3charP;

typedef struct { 
	D1short			x,y,z;
} D3short;

#ifndef _d3int
#define _d3int

typedef struct { 
	int64			x,y,z;
}D3int,D3iVector,D3iValue,D3iPoint,D3int64,D3int_64,D3iVector_64,D3iValue_64,D3iPoint_64;

typedef struct { 
	int32			x,y,z;
} D3int32,D3int_32,D3iVector_32,D3iValue_32,D3iPoint_32;

#endif	    /* _d3int */

#ifndef _d3float
#define _d3float

typedef struct { 
	float			x,y,z;
} D3float,D3fVector,D3fValue,D3fPoint;

#endif	    /* _d3float */

typedef struct { 
	double		x,y,z;
} D3double,D3dVector,D3dValue,D3dPoint;

typedef struct { 
	float			r,phi,theta;
} D3fPolar;

typedef struct { 
	D3float			n;	    /*  normal  to plane       */
	D3float			ex,ey;	    /*  unit vectors in plane  */
	float			scale;	    /*  coordinate scale       */
} D3coord_plane;

typedef struct {
	D2fVector 		dx,dy,dz;
} D3D2fVectorInc;			    /* vector increment in mil.c */

typedef struct {
	D1float	 		dx,dy,dz;
} D3fValueInc;				    /* scalar increment in mil.c */

typedef struct {
      int32         nr;
      float	  d_phi;
      D3float     *vec;      
} D3DirectionLUT;

typedef struct {
      int32         nr;
      D3float     *v;      
      float	  d_phi;
} D3f_directions;

typedef D3fVector	        D3Principal[3];	

#ifndef _d3f_principal
#define _d3f_principal

typedef D3fVector	        D3f_principal[3];	

#endif	    /* _d3f_principal */


typedef D1float			D3QuadraticForm[3][3];
typedef D1float			D3f_quad_form[3][3];

typedef struct {
	D3fVector		*dir;	    /* norm direction */
	D1float			*val;	    /* value in corresp. direction */
	int32			nr;
} D3DirectionalDistr;

typedef D3DirectionalDistr  D3f_dir_distribution;

typedef struct {
	D3fVector		*dir;	    /* norm direction */
	double			*val;	    /* value in corresp. direction */
	int32			nr;
} D3f_dir_ddistribution;  /*  double  */

typedef struct {
	int32	    nr;
	int64	    *offs;	    
	float	    *val;	    
} D3f_offset_LUT;		

typedef struct {
	int32	    nr;
	int64	    *offs;	    
} D3_offset_LUT;		

typedef struct {
	int32		  nr;
	D3_offset_LUT	  *set;
} D3_offset_LUTs;		

/**Check if this makes sense with MAXNBLUT**/
# define MAXNBLUT 3375
typedef struct {
	int	    nr;
	int64	    offs[MAXNBLUT];	    
	D3char	    pos[MAXNBLUT];
} D3_nb_offset_LUT;	

typedef struct {
	D3iPoint_64	pos;	   /* global position 		            */
	D3iValue_64	sep;	   /* sample/separation dist 	     	    */
	D3iValue_64	tsize;	   /* t-vol size size.k = sep.k + 2n, n>=0  */
				   /* offs.k = (size.k-sep.k)/2 = n	    */
	D3iValue_64	dim;	   /* dim.k in terms of sep.k		    */
} D3SupVol,D3SupImage,D3Spatial,D3SupVol_64,D3SupImage_64,D3Spatial_64;

typedef struct {
	D3iPoint_32	pos;	   /* global position 		            */
	D3iValue_32	sep;	   /* sample/separation dist 	     	    */
	D3iValue_32	tsize;	   /* t-vol size size.k = sep.k + 2n, n>=0  */
				   /* offs.k = (size.k-sep.k)/2 = n	    */
	D3iValue_32	dim;	   /* dim.k in terms of sep.k		    */
} D3SupVol_32,D3SupImage_32,D3Spatial_32;

typedef struct { 
	D3iPoint_64		pos;
	D3iValue_64		dim;
} D3Subvol,D3SubImage,D3Subvol_64,D3SubImage_64;

typedef struct { 
	D3iPoint_32		pos;
	D3iValue_32		dim;
} D3Subvol_32,D3SubImage_32;

/** D3AnyFileImage -> version without pointers to write on disk**/
/**New aim Version 030**/
typedef struct {
	char	    version;	    /* version 3.0: 030    */
	char	    *proc_log;
	void	    *dat;
	int32	    id;
	int32	    ref;
	int32	    type;	    /* type of data	    */
	D3int64	    pos;	    /* in pixels (absolut)  */
	D3int64	    dim;	    /* in pixels	    */
	D3int64	    off;	    /* in pixels	    */
	D3int64	    supdim;	    /* in subdims ( = 0 -> no sup image) */
	D3int64	    suppos;	    /* in pixels (relative) */
	D3int64	    subdim;	    /* in pixels	    */
	D3int64	    testoff;	    /* in pixels	    */
	D3float	    el_size_mm;
	D1AnyData   assoc;	    /* associated data	*/
} D3AnyImage,D3aim,D3AnyImage030;

typedef struct {
	char	    version;	    /* version 3.0: 030    */
	int32	    id;
	int32	    ref;
	int32	    type;	    /* type of data	    */
	D3int64	    pos;	    /* in pixels (absolut)  */
	D3int64	    dim;	    /* in pixels	    */
	D3int64	    off;	    /* in pixels	    */
	D3int64	    supdim;	    /* in subdims ( = 0 -> no sup image) */
	D3int64	    suppos;	    /* in pixels (relative) */
	D3int64	    subdim;	    /* in pixels	    */
	D3int64	    testoff;	    /* in pixels	    */
	D3int64	    el_size_nano;   /* in NANO meters !!!   */
	D1AnyFileData assoc;	    /* associated data	    */
} D3AnyFileImage,D3AnyFileImage030;

typedef struct {
	char	    version;	    /* version 2.0: 020    */
	char	    *proc_log;
	void	    *dat;
	int32	    id;
	int32	    ref;
	int32	    type;	    /* type of data	    */
	D3int32	    pos;	    /* in pixels (absolut)  */
	D3int32	    dim;	    /* in pixels	    */
	D3int32	    off;	    /* in pixels	    */
	D3int32	    supdim;	    /* in subdims ( = 0 -> no sup image) */
	D3int32	    suppos;	    /* in pixels (relative) */
	D3int32	    subdim;	    /* in pixels	    */
	D3int32	    testoff;	    /* in pixels	    */
	D3float	    el_size_mm;
	D1AnyData   assoc;	    /* associated data	*/
    /**D1AnyData is the same for all version!**/
} D3AnyImage020;

typedef struct {
	char	    version;	    /* version 2.0: 020    */
    int32	    proc_log;	    /*  length of 32bit pointer */
	int32	    dat;	    /*  length of 32bit pointer */
	int32	    id;
	int32	    ref;
	int32	    type;	    /* type of data	    */
	D3int32	    pos;	    /* in pixels (absolut)  */
	D3int32	    dim;	    /* in pixels	    */
	D3int32	    off;	    /* in pixels	    */
	D3int32	    supdim;	    /* in subdims ( = 0 -> no sup image) */
	D3int32	    suppos;	    /* in pixels (relative) */
	D3int32	    subdim;	    /* in pixels	    */
	D3int32	    testoff;	    /* in pixels	    */
	D3float	    el_size_mm;     /* in millimeters	    */
	D1AnyFileData020 assoc;	    /* associated data	    */
} D3AnyFileImage020;

typedef struct {
	char	    *proc_log;
	void	    *dat;
	int32	    id;
	int32	    ref;
	int32	    type;	    /* type of data	*/
	D3int32	    dim;	    /* in subdims	*/
	D3int32	    off;	    /* in subdims	*/
	D3int32	    subdim;	    /* in pixels	*/
	D3int32	    pos;	    /* in pixels	*/
	D3float	    el_size_mm;
	D1AnyData   assoc;	    /* associated data	*/
	char	    version;	    /* version 1.1: 011 */
} D3AnyImage011;

typedef struct {
	int32	    proc_log;
	int32	    dat;
	int32	    id;
	int32	    ref;
	int32	    type;	    /* type of data	*/
	D3int32	    dim;	    /* in subdims	*/
	D3int32	    off;	    /* in subdims	*/
	D3int32	    subdim;	    /* in pixels	*/
	D3int32	    pos;	    /* in pixels	*/
	D3float	    el_size_mm;
	D1AnyFileData020   assoc;	    /* associated data	*/
	char	    version;	    /* version 1.1: 011 */
} D3AnyFileImage011;

typedef struct {
	char	    *proc_log;
	void	    *dat;
	int32	    id;
	int32	    ref;
	int32	    type;	    /* type of data	*/
	D3iValue_32    dim;	    /* in subdims	*/
	D3iValue_32    off;	    /* in subdims	*/
	D3iValue_32    subdim;	    /* in pixels	*/
	D3iPoint_32    pos;	    /* in pixels	*/
	D1float	    el_size_mm;
	D1AnyData   assoc;	    /* associated data	*/
} D3AnyImage010;

typedef struct {
	int32	    proc_log;
	int32	    dat;
	int32	    id;
	int32	    ref;
	int32	    type;	    /* type of data	*/
	D3iValue_32    dim;	    /* in subdims	*/
	D3iValue_32    off;	    /* in subdims	*/
	D3iValue_32    subdim;	    /* in pixels	*/
	D3iPoint_32    pos;	    /* in pixels	*/
	D1float	    el_size_mm;
	D1AnyFileData020   assoc;	    /* associated data	*/
} D3AnyFileImage010;

typedef struct {
	D3fVector		n;
	D1float			len;
	D1float			rad;
} D3Cylinder;

typedef struct {
	D1float			rad;
} D3Sphere;


/*	 
**  Global constant variables
*/

  static D3int64	D3iNULL = {0,0,0};
  static D3int32	D3iNULL32 = {0,0,0};
  static D3int64	D3iONE  = {1,1,1}; 
  static D3int32	D3iONE32  = {1,1,1};
  static D3char		D3cNULL = {0,0,0};
  static D3char		D3cONE  = {1,1,1}; 

  static D3float	D3fNULL = {0.,0.,0.};
  static D3float	D3fONE  = {1.,1.,1.}; 

/*	  
**  Global temporary variables, should only be used in macros.
*/	  

  static  D3float	D3fTMP;
  static  D3int64		D3iTMP;


#ifdef __cplusplus
}  /* Close scope of 'extern "C"' declaration which encloses file. */
#endif

#endif	    /* _d3_h */
/* DON'T ADD ANYTHING AFTER THIS #endif */
