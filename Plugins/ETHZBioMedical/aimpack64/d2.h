#ifndef _d2_h
#define _d2_h

#ifdef __cplusplus
extern "C" {
#endif


#ifdef __VMS
#define	PrintSourceInfo  printf("at line %d in file %s\n", __LINE__, __FILE__);
#   ifdef __INITIAL_POINTER_SIZE
#	    include "vlm_memory.h"
#	    define USING_VLM_MEMORY
#       define aimpack_malloc _malloc64
#       define aimpack_calloc _calloc64
#       define aimpack_realloc _realloc64
#       define aimpack_memchr _memchr64
#       define aimpack_memcpy _memcpy64
#       define aimpack_memmove _memmove64
#       define aimpack_memset _memset64 
#       define aimpack_strcat _strcat64
#       define aimpack_strchr _strchr64
#       define aimpack_strcpy _strcpy64
#       define aimpack_strncat _strncat64
#       define aimpack_strncpy _strncpy64
#       define aimpack_strpbrk _strpbrk64
#       define aimpack_strrchr _strrchr64
#       define aimpack_strstr _strstr64
#       define aimpack_strtok _strtok64
#   endif
#else
#define	PrintSourceInfo  ;
#define aimpack_malloc malloc
#define aimpack_calloc calloc
#define aimpack_realloc realloc
#define aimpack_strcpy strcpy
#define aimpack_strstr strstr
#define aimpack_memchr memchr
#define aimpack_memcpy memcpy
#define aimpack_memmove memmove
#define aimpack_memset memset
#define aimpack_strcat strcat
#define aimpack_strchr strchr
#define aimpack_strcpy strcpy
#define aimpack_strncat strncat
#define aimpack_strncpy strncpy
#define aimpack_strpbrk strpbrk
#define aimpack_strrchr strrchr
#define aimpack_strstr strstr
#define aimpack_strtok strtok
#endif

#ifdef sun
#define __int64 long long
#define __int32 int
#endif

#ifdef linux 
#define __int64 long long 
#define __int32 int 
#endif

#ifndef NULL
#define NULL 0
#endif

#define	PrintUp	  printf("%s",SC_UP); 

#ifndef PI
#define PI	3.1415926535898
#endif

#ifdef VTK
#ifndef Pi
#define Pi	3.1415926535898
#endif
#endif

#ifndef True
#define True	1
#endif

#ifndef False
#define False	0
#endif

#ifndef On
#define On	1
#endif

#ifndef Off
#define Off	0
#endif

#ifndef Undef
#define Undef (-1)
#endif

#define StrTrueFalse(i)\
(((i) ==  0  )   ? "False" : "True")

#define StrTrueFalseUndef(i)\
(((i) ==  0  )   ? "False" :\
 ((i) == -1  )   ? "Undef" :\
                   "True")
#define StrOnOff(i)\
(((i) ==  0  )   ? "Off" : "On")

#define StrYesNo(i)\
(((i) ==  0  )   ? "No" : "Yes")

#define NUMERICAL_CHAR_SET  "1234567890.eE+-"
#define ALPHABETIC_CHAR_SET \
			"abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"

#define D1fMAX  ((float)9.999e10)

#define pow2(a)	  ( (a)*(a) )
#define pow3(a)	  ( (a)*(a)*(a) )
#define sqr(a)	  pow2(a)
#define Sqr(a)	  pow2(a)
#define Sqrt(a)	  sqrt((double)(a))
#define Trunc(a)  ((long) (a))
#define Round(a)  ( ((a)>=0.) ? Trunc((a) + 0.5) : Trunc((a) - 0.5) )
#define RoundP(a) Trunc((a) + 0.5)	/* positive values only */
#define RoundPos(a)  RoundP(a)
#define RoundN(a,n)  (Round((a)*pow(10,n))/pow(10,n))
#define Fract(a)       ( (a) - Trunc(a) )
#define FractRound(a)  ( (a) - Round(a) )
#define RoundUpPos(a)  (Fract((a)) ? Trunc((a))+1 : Trunc((a)) )
#define RoundDown(a)   ( ((a)>=0.) ? Trunc((a)) : -RoundUpPos(-(a)) )
#define Log(a) ( log((double)(a)) )
#define Log2(a) ( log(a)/log(2) )
/* mod for negative a's by circular shift */ 
#define ModExt(a,b) ( ((a)<0) ? ((b) + (a)%(b))%(b) : (a)%(b) )
#define xModExt(a,b) ( ((1 + abs((a)/(b)))*(b) + (a))%(b)  ) /* old */
#ifndef Min2
#define Min2(a,b) ( ((a)<(b))?(a):(b) )
#define Max2(a,b) ( ((a)>(b))?(a):(b) )
#define Min3(a,b,c) ( ((a)<(b)) ? Min2(a,c) : Min2(b,c) )
#define Max3(a,b,c) ( ((a)>(b)) ? Max2(a,c) : Max2(b,c) )
#endif

/* Bit operations */

#define IsSetBit(data,bit_mask) ((data) & (bit_mask))
#define IsOffBit(data,bit_mask) (!IsSetBit(data,bit_mask))
#define SetBit(data,bit_mask) ((data) |= (bit_mask))
#define ClearBit(data,bit_mask) ((data) &= (~(bit_mask)))
#define ToggleBit(data,bit_mask) if( IsSetBit(data,bit_mask) )\
   ClearBit(data,bit_mask); else SetBit(data,bit_mask);

#define Sort3(v,i)\
{\
  if( v[1] > v[0] ){\
      i[0] = 1;\
      i[1] = 0;\
  } else {\
      i[0] = 0;\
      i[1] = 1;\
  }\
  if( v[2] > v[i[1]] ){\
      i[2] = i[1];\
      i[1] = 2;\
  } else {\
      i[2] = 2;\
  }\
  if( v[i[1]] > v[i[0]] ){\
      i[0] = i[1];\
      i[1] = 3 - i[0] - i[2];\
  }\
}

#define StrEq(s1,s2) ( strcmp(s1, s2) == 0 )
#define StrEqAbbr(abbr,key) ( strncmp(abbr, key, strlen(abbr)) == 0 )
#define StrCpyAlloc(out,in) {\
  (out) = (char *) aimpack_malloc(strlen(in)+1);\
  aimpack_strcpy(out,in);\
}
#define StrToUpper(str)\
for( intTMP=0; intTMP<strlen(str); intTMP++ ) {\
  (str)[intTMP] = (char) toupper((str)[intTMP]);\
}
/* True <-- yes,true,on,1 */
#define StrBooleanIsTrue(str) (aimpack_strchr("YyTtOo1",(str)[0]) != 0 )   
#define StrReplaceSuffix(str,old_suffix,new_suffix)\
    aimpack_strcpy(aimpack_strstr(str, old_suffix), new_suffix);

#define Free(dat)\
{\
  if (dat) {\
    free(dat);\
    dat = 0;\
  }\
}

#define Malloc(dat,nr,typ)\
{\
  (dat) = (typ *) aimpack_malloc((nr)*sizeof(typ));\
  if ( (dat) == 0 ) {\
    if( ((nr)*sizeof(typ) < 1024) ) {\
      printf("Couldn't allocate (Malloc) %Ld Bytes of memory\n",(nr)*sizeof(typ));\
    }\
    if( ((nr)*sizeof(typ) >= 1024) && ((nr)*sizeof(typ) < 1048576) ) {\
      printf("Couldn't allocate (Malloc) %5.1f kBytes of memory\n",((float)(nr)*sizeof(typ))/1024.0);\
    }\
    if( ((nr)*sizeof(typ) >= 1048576) ) {\
      printf("Couldn't allocate (Malloc) %5.1f MBytes of memory\n",((float)(nr)*sizeof(typ))/1048576.0);\
    }\
    PrintSourceInfo;\
	exit(0);\
  }\
}

#define MallocMax(dat,nr,max_nr,buff_nr,typ)\
{\
  if( (nr) > (max_nr) + (buff_nr) || (max_nr) == 0 || (dat) == 0 ){\
    if( dat ) free(dat);\
    (max_nr) = (nr);\
    (dat) = (typ *) aimpack_malloc(((nr)+(buff_nr))*sizeof(typ));\
    if ( (dat) == 0 ) {\
      if( (((nr)+(buff_nr))*sizeof(typ) < 1024) ) {\
        printf("Couldn't allocate (MallocMax) %Ld Bytes of memory\n",((nr)+(buff_nr))*sizeof(typ));\
      }\
      if( (((nr)+(buff_nr))*sizeof(typ) >= 1024) && (((nr)+(buff_nr))*sizeof(typ) < 1048576) ) {\
        printf("Couldn't allocate (MallocMax) %5.1f kBytes of memory\n",((float)((nr)+(buff_nr))*sizeof(typ))/1024.0);\
      }\
      if( (((nr)+(buff_nr))*sizeof(typ) >= 1048576) ) {\
        printf("Couldn't allocate (MallocMax) %5.1f MBytes of memory\n",((float)((nr)+(buff_nr))*sizeof(typ))/1048576.0);\
      }\
      PrintSourceInfo;\
 	  exit(0);\
    }\
  }\
}

#define MallocVoid(dat,nr,typ)\
{\
  (dat) = (void *) aimpack_malloc((nr)*sizeof(typ));\
  if ( (dat) == 0 ) {\
    if( ((nr)*sizeof(typ) < 1024) ) {\
      printf("Couldn't allocate (MallocVoid) %Ld Bytes of memory\n",(nr)*sizeof(typ));\
    }\
    if( ((nr)*sizeof(typ) >= 1024) && ((nr)*sizeof(typ) < 1048576) ) {\
      printf("Couldn't allocate (MallocVoid) %5.1f kBytes of memory\n",((float)(nr)*sizeof(typ))/1024.0);\
    }\
    if( ((nr)*sizeof(typ) >= 1048576) ) {\
      printf("Couldn't allocate (MallocVoid) %5.1f MBytes of memory\n",((float)(nr)*sizeof(typ))/1048576.0);\
    }\
    PrintSourceInfo;\
	exit(0);\
  }\
}

#define VoidMalloc(dat,size)\
{\
  (dat) = (void *) aimpack_malloc(size);\
  if ( (dat) == 0 ) {\
    if( (size < 1024) ) {\
      printf("Couldn't allocate (VoidMalloc) %Ld Bytes of memory\n",size);\
    }\
    if( (size >= 1024) && (size < 1048576) ) {\
      printf("Couldn't allocate (VoidMalloc) %5.1f kBytes of memory\n",(float)(size)/1024.0);\
    }\
    if( (size >= 1048576) ) {\
      printf("Couldn't allocate (VoidMalloc) %5.1f MBytes of memory\n",(float)(size)/1048576.0);\
    }\
    PrintSourceInfo;\
    exit(0);\
  }\
}

#define ReMalloc(dat,nr,typ)\
{\
  if (dat) free(dat);\
  (dat) = (typ *) aimpack_malloc((nr)*sizeof(typ));\
  if ( (dat) == 0 ) {\
    if( ((nr)*sizeof(typ) < 1024) ) {\
      printf("Couldn't allocate (ReMalloc) %Ld Bytes of memory\n",(nr)*sizeof(typ));\
    }\
    if( ((nr)*sizeof(typ) >= 1024) && ((nr)*sizeof(typ) < 1048576) ) {\
      printf("Couldn't allocate (ReMalloc) %5.1f kBytes of memory\n",((float)(nr)*sizeof(typ))/1024.0);\
    }\
    if( ((nr)*sizeof(typ) >= 1048576) ) {\
      printf("Couldn't allocate (ReMalloc) %5.1f MBytes of memory\n",((float)(nr)*sizeof(typ))/1048576.0);\
    }\
    PrintSourceInfo;\
    exit(0);\
  }\
}

#define Calloc(dat,nr,typ)\
{\
  (dat) = (typ *) aimpack_calloc((nr),sizeof(typ));\
  if ( (dat) == 0 ) {\
    if( ((nr)*sizeof(typ) < 1024) ) {\
      printf("Couldn't allocate (Calloc) %Ld Bytes of memory\n",(nr)*sizeof(typ));\
    }\
    if( ((nr)*sizeof(typ) >= 1024) && ((nr)*sizeof(typ) < 1048576) ) {\
      printf("Couldn't allocate (Calloc) %5.1f kBytes of memory\n",((float)(nr)*sizeof(typ))/1024.0);\
    }\
    if( ((nr)*sizeof(typ) >= 1048576) ) {\
      printf("Couldn't allocate (Calloc) %5.1f MBytes of memory\n",((float)(nr)*sizeof(typ))/1048576.0);\
    }\
    PrintSourceInfo;\
    exit(0);\
  }\
}

#define ReCalloc(dat,nr,typ)\
{\
  if (dat) free(dat);\
  (dat) = (typ *) aimpack_calloc((nr),sizeof(typ));\
  if ( (dat) == 0 ) {\
    if( ((nr)*sizeof(typ) < 1024) ) {\
      printf("Couldn't allocate (ReCalloc) %Ld Bytes of memory\n",(nr)*sizeof(typ));\
    }\
    if( ((nr)*sizeof(typ) >= 1024) && ((nr)*sizeof(typ) < 1048576) ) {\
      printf("Couldn't allocate (ReCalloc) %5.1f kBytes of memory\n",((float)(nr)*sizeof(typ))/1024.0);\
    }\
    if( ((nr)*sizeof(typ) >= 1048576) ) {\
      printf("Couldn't allocate (ReCalloc) %5.1f MBytes of memory\n",((float)(nr)*sizeof(typ))/1048576.0);\
    }\
    PrintSourceInfo;\
    exit(0);\
  }\
}

#define MemoryManagmentInit(mem_man_func_ptr)\
{\
  voidP_fp	MemoryManageAlloc=NULL;\
  MemoryManageAlloc = (voidP_fp) mem_man_func_ptr;\
}

#define Limit(val,min,max)  \
        ( ((val) >= (max)) ? (max) : ( ((val)<=(min)) ? (min) : (val) ) ) 

#define FEqRel(a,b,e)						     \
   ( (b) > 0 ?  ( (a) <= (b)*(1. + (e)) ) && ( (a) >= (b)*(1.- (e)) )   \
	     :  ( (a) >= (b)*(1. + (e)) ) && ( (a) <= (b)*(1.- (e)) ) ) 

#define feq_rel(a,b,e) FEqRel(a,b,e)

#define FEq(a,b,e)					      \
( (b) > 0 ?  ( (a) <= ((b)+(e)) ) && ( (a) >= ((b)-(e)) )     \
          :  ( (a) >= ((b)+(e)) ) && ( (a) <= ((b)-(e)) ) )   

#define feq_abs(a,b,e) FEq(a,b,e)

typedef char Str160[161],Str80[81],Str40[41],Str20[21],Str10[11];

typedef	void   (*void_fp)();
typedef	float  (*float_fp)();
typedef	double (*double_fp)();
typedef	void *(*voidP_fp)();

/*-----------------------------------------------------------------------------+
| Types used to read CT header files.					       |
+-----------------------------------------------------------------------------*/

#ifndef _CTFileHeader
#define _CTFileHeader
typedef union {

    struct {
	short	samples;
	short	scans;
	short	firstchannel;
	short	lastchannel;
	long		scandistum;
	short	sampletime[8];
	short	anzproj;
	short	pulsrate[8];
	short	lnI0[8];
	short	day;
	char		month[4];
	short	year;
	short	hour;
	short	min;
	short	sec;
	char		name[28];
	short	measnr;
	short	age;
	char		sex[2];
	short	height;
	short	weight;
	char		diagnose[22];
	char		therapy[10];
	char		site[12];
	char		group[10];
	char		remarks[68];
	short	backpro_shift;
	short	convol_shift;
	char		ctrl_filename[16];
	float		scout_begin_pos;
	short	scout_scale;
	short	scout_slices;
	short	object_size;
	float		zinc;
	long		zpos;
	float		refline;
	short	dummy1[10];
	short	sampletime2[8];
	long		dummy2;
	short	pulsrate2[8];
	long		dummy3;
	short	lnI02[8];
	float		scandistmm;
	short	dummy4;
	char		cast[4];
	short	scanner_id;
	short	dummy5;
	char		filename[8];
    } fh;

    struct {
	short	header[256];
    } ih;

} CTFileHeader, *CTFileHeaderP, **CTFileHeaderPP;
#endif  /* _CTFileHeader */

#ifndef __RmBoneHeader
#define __RmBoneHeader
typedef struct _RmBoneHeader {
	long	 x;		/* start x position in the xy - planes 	*/
	long	 y;		/* start y position in the xy - planes	*/
	long	 dimx;		/* number of voxels in the x direction	*/
	long	 dimy;		/* number of voxels in the y direction	*/
	long	 dimz;		/* number of voxels in the z direction	*/
	long	 support;	/* number of voxels in all directions,	*/
				/* which can't be used (supposed zero)	*/
	char	 info[256];	/* informations on the image processing	*/
	char	 first[32];	/* name of the first original slice	*/
	char	 last [32];	/* name of the last original slice	*/
	float	 pixsize_in_mm;	/* results out of the expressions in	*/
				/* the CT header: scandistum / samples 	*/
} RmBoneHeader, *RmBoneHeaderP;
#endif /* __RmBoneHeader */	   

/* Global temporary variables, should only be used in macros */

  static long			intTMP;

/* Global constants, special characters SC */

  static char		SC_BOLD[]      = {27,'[','1','m'}; 
  static char		SC_UNDERLINE[] = {27,'[','4','m'}; 
  static char		SC_BLINK[]     = {27,'[','5','m'}; 
  static char		SC_REVERSE[]   = {27,'[','7','m'}; 
  static char		SC_OFF[]       = {27,'[','0','m'}; 

  static char		SC_UP[]	       = {27,'[','1','A'}; 

  static char		SC_BELL	       = 7; 

/*-----------------------------------------------------------------------------+
|     D1 Macros and Types						       |
+-----------------------------------------------------------------------------*/

#define D1fMAX	    ((float)9.999e10)
#define D1fMIN	    ((float)-9.999e10)
#define D1cMAX	    ((1<<7)-1)
#define D1cMIN	    (-((1<<7)))
#define D1sMAX	    ((1<<15)-1)
#define D1sMIN	    (-((1<<15)))
#define D1sMAX	    ((1<<15)-1)
#define D1sMAX_LUT  ((1<<14)-1)

#define D1EQ(a,b,e)	 FEq(a,b,e)

#define D1fTrunc(a)   { (a) += 100.; (a) -= 100.; } /* truncate round errors  */
#define D1Round(a)    Round(a)			  

/*	  
**  D1C Macros (Complex) old, use D1Z Macros 
*/	  

#define D1CRe(z)		    ((z).re)
#define D1CIm(z)		    ((z).im)
#define D1CAmplitudeSqr(z)	    ( Sqr((z).re) + Sqr((z).im) )
#define D1CAmplitude(z)		    sqrt((double) D1CAmplitudeSqr(z))

/*	  
**  D1Z Macros (Complex) New definitions 
*/	  

#define D1ZRe(z) D1CRe(z)

#define D1ZIm(z) D1CIm(z)

#define D1ZAmplitudeSqr(z) D1CAmplitudeSqr(z)

#define D1ZAmplitude(z) D1CAmplitude(z)

/* nb! not recursive e.g. D1ZMul(a,b,a) is not correct */

#define D1ZMul(a,b,z){\
  (z).re = (a).re*(b).re - (a).im*(b).im;\
  (z).im = (a).re*(b).im + (a).im*(b).re;\
}

#define D1ZConj(in,out){\
  (out).re =  (in).re;\
  (out).im = -(in).im;\
}

#define D1ZCMul(in,const,out){\
  (out).re = (in).re*(const);\
  (out).im = (in).im*(const);\
}

#define D1ZCDiv(in,const,out){\
  (out).re = (in).re/(const);\
  (out).im = (in).im/(const);\
}
 



/*	  
**  D1F Macros (Fraction)
*/	  

#define D1FNumerator(f)		    ((f).num)
#define D1FDenominator(f)	    ((f).den)

/*	  
**  Type identifier constants
*/	  

#define D1Tundef	    	0
#define D1Tchar		    	((1<<16) + sizeof(char))
#define D1Tpixel	   		((13<<16) + sizeof(char))
#define D1Tshort	    	((2<<16) + sizeof(short))
#define D1Tint		    	((3<<16) + sizeof(int64))
#define D1Tint32			((3<<16) + sizeof(int32))
#define D1Tfloat	    	((4<<16) + sizeof(float))
#define D1Tcomplex	    	((5<<16) + sizeof(D1complex))
#define D1Tfraction	    	((7<<16) + sizeof(D1fraction))
#define D1TcharP	    	((9<<16) + sizeof(D1charP))
#define D1TcharCmp	    	((8<<16) + sizeof(D1charCmp))
#define D1TcharCmp2	    	((8<<16) + sizeof(D1charCmp2))
#define D3Tbit8		    	((6<<16) + sizeof(D1byte))
#define D3Tchar		    	((6<<16) + sizeof(D3char))
#define D3Tint		   		((10<<16) + sizeof(D3int64))
#define D3Tint32		   	((10<<16) + sizeof(D3int32))
#define D3Tint64		   	((10<<16) + sizeof(D3int64))
#define D3Tfloat	    	((6<<16) + sizeof(D3float))
#define D3Taim		   		((11<<16) + sizeof(D3aim))
#define D1Tboolean	   		((12<<16) + sizeof(D1boolean))
#define D1Tint_tri	   		((14<<16) + sizeof(D1int64_tri))
#define D1Tint_tet	   		((15<<16) + sizeof(D1int64_tet))
#define D1Tint_hex	   		((15<<16) + sizeof(D1int64_hex))
#define D1Tint64_tri	   	((14<<16) + sizeof(D1int64_tri))
#define D1Tint64_tet	   	((15<<16) + sizeof(D1int64_tet))
#define D1Tint64_hex	   	((15<<16) + sizeof(D1int64_hex))
#define D1Tint32_tri	   	((14<<16) + sizeof(D1int32_tri))
#define D1Tint32_tet	   	((15<<16) + sizeof(D1int32_tet))
#define D1Tint32_hex	   	((15<<16) + sizeof(D1int32_hex))
#define D2Tint		   		((16<<16) + sizeof(D2int64))
#define D2Tint64		   	((16<<16) + sizeof(D2int64))
#define D2Tint32		   	((16<<16) + sizeof(D2int32))
#define D2Tfloat	   		((17<<16) + sizeof(D2float))
#define D1Trgb		   		((18<<16) + sizeof(D1rgb))
/* #define D3PTsp_db_type	   	((19<<16) + sizeof(D3P_Sp_db_type)) */
#define D1Tdouble	   		((20<<16) + sizeof(double))
#define D1TbinCmp	   		((21<<16) + sizeof(char))
#define D1Tuchar	   		((22<<16) + sizeof(unsigned char))
#define D1Tushort	   		((23<<16) + sizeof(unsigned short))
#define D1Ti3efloat	   		((26<<16) + sizeof(float))

#define D1TZBencoded	   	((100<<16)+ sizeof(D1byte)) /* XAimpack */

#define D1Tsize_mask	    0x0000FFFF			  /* first two bytes */
#define D1TSize(t)	    ((t) & D1Tsize_mask )

#define D3Type_aim	      "aim"			  /* type name */
#define D1Type_char	      "char"
#define D1Type_short	      "short"
#define D1Type_float	      "float"
#define D1Type_charCmp	      "charCmp"
#define D1Type_run_length     "run_length"
#define D1Type_binCmp	      "binCmp"
#define D1Type_bin_run_length "bin_run_length"
#define D1Type_charCmp2	      "charCmp2"
#define D1Type_bit8	      "bit8"

#define DnTtype(name) (\
  StrEqAbbr(name,D1Type_char)		? D1Tchar :\
  StrEqAbbr(name,D1Type_short)		? D1Tshort :\
  StrEqAbbr(name,D1Type_float)		? D1Tfloat :\
  StrEqAbbr(name,D1Type_i3efloat)	? D1Ti3efloat :\
  StrEqAbbr(name,D1Type_bit8)		? D3Tbit8 :\
  StrEqAbbr(name,D1Type_charCmp)	? D1TcharCmp:\
  StrEqAbbr(name,D1Type_run_length)	? D1TcharCmp:\
  StrEqAbbr(name,D1Type_binCmp)		? D1TbinCmp:\
  StrEqAbbr(name,D1Type_bin_run_length) ? D1TbinCmp:\
  StrEqAbbr(name,D1Type_charCmp2)	? D1TcharCmp2:\
                                          D1Tundef)

#define D1TName(t)\
(((t) == D1Tchar  )	  ? "Char" :\
 ((t) == D1Tpixel )	  ? "Pixel" :\
 ((t) == D1Tshort )	  ? "Short" :\
 ((t) == D1Tint	  )	  ? "Int" :\
 ((t) == D1Tfloat )	  ? "Float" :\
 ((t) == D1Ti3efloat )	  ? "I3EFloat (little endian)" :\
 ((t) == D1Tdouble)	  ? "Double" :\
 ((t) == D1Trgb   )	  ? "RGB" :\
 ((t) == D1TcharCmp )	  ? "CharCmp" :\
 ((t) == D1TcharCmp2)	  ? "CharCmp2" :\
 ((t) == D1TbinCmp)	  ? "BinCmp" :\
 ((t) == D3Tbit8  )	  ? "Bit8" :\
 /*((t) == D3PTsp_db_type)  ? "SP_DB_rec" :*/\
 		        "Spec")

/*	  
**  D1 Types
*/	  
 
typedef char		D1char,*D1charP;
typedef char		D1str10[11];
typedef char		D1str20[21];
typedef char		D1str40[41];
typedef char		D1str80[81];
typedef unsigned char	D1uchar,D1Bin,D1byte,D1pixel;

typedef short		D1short;
typedef unsigned short	D1ushort;

typedef long		D1long;
typedef unsigned long   D1ulong;

#ifdef __GNUC__
typedef signed long long int64,D1int,D1boolean,*int64P;
typedef signed long long D1int64, *D1int64P;
typedef signed long int int32,*int32P,D1int32, *D1int32P,D1boolean32;
/** a long is suposed to be 32 bit**/
typedef unsigned long long D1uint64,D1uint;
typedef unsigned long int D1uint32;
#else
typedef signed __int64 int64,D1int,D1boolean,*int64P;
typedef signed __int64 D1int64, *D1int64P;
typedef signed __int32 int32,*int32P,D1int32, *D1int32P,D1boolean32;
/** a long is suposed to be 32 bit**/
typedef unsigned __int64 D1uint64,D1uint;
typedef unsigned __int32 D1uint32;
#endif

typedef short int16;

typedef float		D1float;    /* als double def um round errors zu det. */

typedef double		xxxD1float;    /* als double def um round errors zu det. */
typedef double		D1double;	


typedef struct  {
	int64 idx[3];
} D1int_tri,D1int64_tri;

typedef struct  {
	int32 idx[3];
} D1int32_tri;

typedef struct  {
	int64 idx[4];
} D1int_tet,D1int64_tet;

typedef struct  {
	int32 idx[4];
} D1int32_tet;

typedef struct  {
	int64 idx[8];
} D1int_hex, D1int64_hex;

typedef struct  {
	int32 idx[8];
} D1int32_hex;

typedef struct {
	float	    re,im;
} D1complex;

typedef struct {
	int32	    num,den;
} D1fraction;

typedef struct {
	unsigned char r,g,b;	/* red,green,blue */
} D1rgb;

typedef struct {
	int64	    n;
	double	    sum;
	double	    sqr_sum;
	double	    max_val;
	double	    min_val;
} D1Stat,D1Stat64;

typedef struct {
	int32	    n;
	double	    sum;
	double	    sqr_sum;
	double	    max_val;
	double	    min_val;
} D1Stat32;

typedef struct {
	D1char	    val;
	D1uchar	    len;
} D1charCmp;  /* compressed */

typedef struct {
	D1char	    val;
	D1char	    dif;
	D1uchar	    len;
} D1charCmp2;  /* compressed level 2: differeces */

typedef struct {
	int32	    id;		    /*  identifier	      */
	void	    *dat;
	int32	    nr;		    /*  number of elements    */
	int32	    size;	    /*  size of elements      */
	int32	    type;	    /*  type of element data  */
} D1AnyData;

typedef struct {
	int32	    id;		    /*  identifier	      */
	int32	    nr;		    /*  number of elements    */
	int32	    size;	    /*  size of elements      */
	int32	    type;	    /*  type of element data  */
} D1AnyFileData,D1AnyFileData030;

typedef struct {
	int32	    id;		    /*  identifier	      */
	int32	    dat;	    /*  pointer 32bit long */
	int32	    nr;		    /*  number of elements    */
	int32	    size;	    /*  size of elements      */
	int32	    type;	    /*  type of element data  */
} D1AnyFileData020;

/* D1AnyParam Types  */

typedef struct {
	int32	    id;		    /*  identifier	      */
	int32	    valid;
	void	    *dat;
	int32	    nr;		    /*  number of elements    */
	int32	    type;	    /*  type of element data  */
	char	    *name;
} D1AnyParam;

/* D1AnyList Types  */

typedef struct  {
	int32		type;
	void		*dat;
	int64		nr;
	int64		max_nr;	     /* for efficient memory allocation  */	
	int64		buff_nr;     /* for efficient memory allocation  */	
} D1AnyList, D1AnyList64;  

typedef struct  {
	int32		type;
	void		*dat;
	int32		nr;
	int32		max_nr;	     /* for efficient memory allocation  */	
	int32		buff_nr;     /* for efficient memory allocation  */	
} D1AnyList32;

#define D1AnyListInit(lst,typ) {\
  (lst).type	      = typ;\
  (lst).dat	      = 0;\
  (lst).nr	      = 0;\
}

#define D1AnyListDynamInit(lst,typ,b_nr) {\
  D1AnyListInit(lst,typ)\
  (lst).buff_nr	      = b_nr;\
  (lst).max_nr	      = 0;\
}

#ifdef USING_VLM_MEMORY

#define MALLOC MallocVLMDat
#define REMALLOC ReMallocVLMDat
#define CALLOC CallocVLMDat
#define RECALLOC ReCallocVLMDat
#define FREE FreeVLMDat
#define MEMSET VLM_memset
#define MEMCPY VLM_memcpy

#define D1AnyListFree(lst)  FreeVLMDat((lst).dat)

#define D1AnyListMemSize(lst) ((lst).nr*D1TSize((lst).type))

/*  Must be Freed by FreeVLMDat !!! */
/*  char for memory allocation */
#define D1AnyListMalloc(lst)\
{\
  MallocVLMDat((lst).dat,D1AnyListMemSize(lst),char);\
}

#define D1AnyListMallocMax(lst, new_nr)\
{\
  (lst).nr = (new_nr);\
  if( (lst).nr > (lst).max_nr + (lst).buff_nr || (lst).dat == 0){\
    if( (lst).dat ) FreeVLMDat( (lst).dat );\
    (lst).max_nr = (lst).nr;\
    (lst).nr += (lst).buff_nr;\
    D1AnyListMalloc(lst);\
    (lst).nr = (lst).max_nr;\
  }\
}

#define D1AnyListMallocMaxInfo(lst, new_nr)\
{\
  (lst).nr = (new_nr);\
  if( (lst).nr > (lst).max_nr + (lst).buff_nr || (lst).dat == 0){\
    printf("D1AnyListMallocMax: MORE MEM NEEDED  nr = %Ld\n\n",new_nr);\
    if( (lst).dat ) FreeVLMDat( (lst).dat );\
    (lst).max_nr = (lst).nr;\
    (lst).nr += (lst).buff_nr;\
    D1AnyListMalloc(lst);\
    (lst).nr = (lst).max_nr;\
  }\
}

#else

#define MALLOC Malloc
#define REMALLOC ReMalloc
#define CALLOC Calloc
#define RECALLOC ReCalloc
#define FREE Free
#define MEMSET aimpack_memset
#define MEMCPY aimpack_memcpy

#define D1AnyListFree(lst)  Free((lst).dat)

#define D1AnyListMemSize(lst) ((lst).nr*D1TSize((lst).type))

#define D1AnyListMalloc(lst)\
{\
    VoidMalloc((lst).dat, D1AnyListMemSize(lst));\
}

#define D1AnyListMallocMax(lst, new_nr)\
{\
  (lst).nr = (new_nr);\
  if( (lst).nr > (lst).max_nr + (lst).buff_nr || (lst).dat == 0){\
    if( (lst).dat ) free( (lst).dat );\
    (lst).max_nr = (lst).nr;\
    (lst).nr += (lst).buff_nr;\
    D1AnyListMalloc(lst);\
    (lst).nr = (lst).max_nr;\
  }\
}

#define D1AnyListMallocMaxInfo(lst, new_nr)\
{\
  (lst).nr = (new_nr);\
  if( (lst).nr > (lst).max_nr + (lst).buff_nr || (lst).dat == 0){\
    printf("D1AnyListMallocMax: MORE MEM NEEDED  nr = %Ld\n\n",new_nr);\
    if( (lst).dat ) free( (lst).dat );\
    (lst).max_nr = (lst).nr;\
    (lst).nr += (lst).buff_nr;\
    D1AnyListMalloc(lst);\
    (lst).nr = (lst).max_nr;\
  }\
}
#endif

/*	  
**  D2 Macros
*/	  

#define D2Round(v,u) { (u).x = D1Round(v.x); (u).y = D1Round(v.y); }
 
#define D2fTrunc(v) { D1fTrunc(v.x); D1fTrunc(v.y); }

#define D2Set(v,a,b)	    {(v).x=(a); (v).y=(b);} 

#define D2Eq(v,u) ( ((v).x == (u).x) && ((v).y == (u).y) )

#define D2Min(v)  Min2( (v).x, (v).y )
#define D2Max(v)  Max2( (v).x, (v).y )


#define D2MaxComponents(v,u,w)   {(w).x = Max2((v).x,(u).x); \
                             	  (w).y = Max2((v).y,(u).y); }

#define D2MinComponents(v,u,w)   {(w).x = Min2((v).x,(u).x); \
                             	  (w).y = Min2((v).y,(u).y); }


#define D2Copy(v,u)	    {(u).x = (v).x;\
                             (u).y = (v).y;}

#define D2Inc(v,inc)        {(v).x += (inc).x ;  \
                             (v).y += (inc).y ;  }

#define D2Dec(v,dec)        {(v).x -= (dec).x ;  \
                             (v).y -= (dec).y ;  }

#define D2Sub(v,u,diff)		{diff.x	     = (v).x - (u).x;  \
				 diff.y	     = (v).y - (u).y;  }

#define D2Add(v,u,sum)		{sum.x	     = (v).x + (u).x;  \
				 sum.y	     = (v).y + (u).y;  }

#define D2AddCMul(v,c,u,sum)      {(sum.x) = (v).x + (c)*(u).x;  \
                             	   (sum.y) = (v).y + (c)*(u).y;  }

#define D2SubCMul(v,c,u,sum)      {(sum.x) = (v).x - (c)*(u).x;  \
                             	   (sum.y) = (v).y - (c)*(u).y;  }

#define D2AddCDiv(v,u,c,sum)      {(sum).x = (v).x + (u).x/(c);  \
                             	   (sum).y = (v).y + (u).y/(c);  }
      /* nb! parameters not same as D2AddCMul */

#define D2Div(v,u,quotient)	{quotient.x  = (v).x / (u).x;  \
				 quotient.y  = (v).y / (u).y;  }

#define D2Mod(v,u,remainder)    {remainder.x = (v).x % (u).x;  \
				 remainder.y = (v).y % (u).y;  }

#define D2Mul(v,u,prod)		{prod.x	     = (v).x * (u).x;  \
				 prod.y	     = (v).y * (u).y;  }

#define D2CMul(v,c,u)	    {(u).x = (v).x*(c);\
                             (u).y = (v).y*(c);}

#define D2CDiv(v,c,u)	    {(u).x = (v).x/(c);\
                             (u).y = (v).y/(c);}


#define D2LinComb(a,v,b,u,lincmb) {lincmb.x = (a)*(v).x + (b)*(u).x;  \
                                   lincmb.y = (a)*(v).y + (b)*(u).y;  }
   
#define D2SProd(a,b)	    ( (a).x*(b).x + (a).y*(b).y  )

#define D2AbsSqr(a)	    D2SProd(a,a)
#define D2Abs2(a)	    D2SProd(a,a)
#define D2Abs(a)	    sqrt( (double) D2SProd(a,a))
#define D2Length(a)	    D2Abs(a)
#define D2Length2(a)	    D2Abs2(a)

#define D2Distance2(a,b)\
  (Sqr((a).x - (b).x) + Sqr((a).y - (b).y))
#define D2Distance(a,b) sqrt(D2Distance2(a,b))

#define	D2Norm(v)	    { D1fTMP = D2Length(v) ;\
			      (v).x /= D1fTMP;	  			  \
			      (v).y /= D1fTMP;				  }
   
#define	D2NormVector(v)	    D2Norm(v)

#define D2InRectangle(p,offs,size) \
      ( (p).x >= (offs).x 	       && (p).y >= (offs).y 		 && \
        (p).x < ((offs).x + (size).x)  && (p).y < ((offs).y + (size).y)    )


#define D2InRectangleTop(p,offs,size) \
      ( (p).x >= (offs).x 	       && (p).y >= (offs).y 		 && \
        (p).x <= ((offs).x + (size).x) && (p).y <= ((offs).y + ((size).y/2)) )

#define D2RectInRect(sub,sup)\
 ( (sub).pos.x >= (sup).pos.x  &&  (sub).pos.y >= (sup).pos.y &&\
  ((sub).pos.x  + (sub).dim.x) <= ((sup).pos.x +  (sup).dim.x) &&\
  ((sub).pos.y  + (sub).dim.y) <= ((sup).pos.y +  (sup).dim.y) )

#define D2InRect(p,rect)\
 ( (p).x >= (rect).pos.x  &&  (p).y >= (rect).pos.y &&\
   (p).x < ((rect).pos.x +  (rect).dim.x) &&\
   (p).y < ((rect).pos.y +  (rect).dim.y) )

#define D2iPrintf(v,mess){ \
		printf("%10s: ",mess); \
		printf("%5Ld %5Ld",(v).x,(v).y); \
		printf("\n");}
    

#define D2StrSetNumerical(str,p){\
  (p).x = (str);\
  while ( strspn((p).x,NUMERICAL_CHAR_SET) == 0 ) (p).x++;\
  (p).y =  p.x + strspn((p).x,NUMERICAL_CHAR_SET);\
  while ( strspn((p).y,NUMERICAL_CHAR_SET) == 0 ) (p).y++;\
}

#define D2Atoi(ap,ip){\
  (ip).x = atoi(ap.x);\
  (ip).y = atoi(ap.y);\
}

#define D2Atof(ap,fp){\
  (fp).x = atof(ap.x);\
  (fp).y = atof(ap.y);\
}

#define D2SetImage(im,val) \
  for( D1iTMP1=0; D1iTMP1<(im).dim.y*(im).dim.x; D1iTMP1++ ) \
    *((im).dat + D1iTMP1) = (val);

#define D2InitImage(im) (im).dat = 0;

#define D2MallocAny(im)\
{\
  switch( im.type ) {\
    case D1Tshort:\
      MallocVoid((im).dat, (im).dim.x*(im).dim.y, D1short);\
      break;\
    case D1Tchar:\
      MallocVoid((im).dat, (im).dim.x*(im).dim.y, D1char);\
      break;\
    case D1Tpixel:\
      MallocVoid((im).dat, (im).dim.x*(im).dim.y, D1pixel);\
      break;\
    case D1Tfloat:\
      MallocVoid((im).dat, (im).dim.x*(im).dim.y, D1float);\
      break;\
    default:\
      printf("Type not supported (D2MallocAny)\n");\
      PrintSourceInfo;\
      exit(0);\
  }\
}



#define D2Malloc(im,typ)\
{\
    Malloc((im).dat, (im).dim.x*(im).dim.y, typ);\
}

#define D2Calloc(im,typ)\
{\
    Calloc((im).dat, (im).dim.x*(im).dim.y, typ);\
}

#define D2ReMalloc(im,typ)\
{\
    ReMalloc((im).dat, (im).dim.x*(im).dim.y, typ);\
}

#define D2ReCalloc(im,typ)\
{\
    ReCalloc((im).dat, (im).dim.x*(im).dim.y, typ);\
}

#define DnFree(var)\
{\
  if ((var).dat == 0) {\
    printf("Warning: Attempt to free() zero pointer\n");\
    PrintSourceInfo;\
  } else {\
    free((var).dat);\
    (var).dat=0;\
  }\
}

#define DnFreeAny(var)\
{\
  if ((var).dat == 0) {\
    printf("Warning: Attempt to free() zero pointer (data) \n");\
    PrintSourceInfo;\
  } else {\
    free((var).dat);\
    (var).dat=0;\
  }\
  if ((var).proc_log) {\
    free((var).proc_log);\
    (var).proc_log=0;\
  }\
}

#define DnProcLogFree(im) Free((im).proc_log)

#define DnProcLogInit(im,str){\
  ReMalloc((im).proc_log, strlen(str)+1, char);\
  aimpack_memcpy((im).proc_log, str, strlen(str)+1);\
}

#define DnProcLogCopy(im1,im2){\
  if( (im1).proc_log != 0 ){\
  	ReMalloc((im2).proc_log, strlen((im1).proc_log)+1, char);\
    aimpack_memcpy((im2).proc_log, (im1).proc_log, strlen((im1).proc_log)+1);\
  } else {\
    (im2).proc_log = 0;\
  }\
}

#define DnProcLogAppend(im,str){\
  if( ((im).proc_log) != 0 ) {\
    Malloc(D1cpTMP, strlen((im).proc_log)+strlen(str)+1, char );\
    aimpack_memcpy(D1cpTMP, (im).proc_log, strlen((im).proc_log) );\
    aimpack_memcpy(D1cpTMP+strlen((im).proc_log), str, strlen(str)+1 );\
    if( (im).proc_log !=0 ) free( (im).proc_log );\
    (im).proc_log = D1cpTMP;\
  } else {\
    printf("Can't append string to undefined proc log (=0). The string was:\n%s\n",str);\
    PrintSourceInfo;\
  }\
}


#define DnProcLogAddString(im,str){\
  sprintf(D1strTMP,"%s\n",str);\
  DnProcLogAppend(im,D1strTMP);\
}

#define DnProcLogAddPar(im, format_str, par){\
  sprintf(D1strTMP, format_str, par);\
  DnProcLogAppend(im,D1strTMP);\
}

#define DnProcLogAddFloatPar(im,name,val){\
  sprintf(D1strTMP,"%-42s%11.5f\n",name,val);\
  DnProcLogAppend(im,D1strTMP);\
}

#define DnProcLogAddD3floatPar(im,name,par){\
  sprintf(D1strTMP,"%-42s%11.5f%11.5f%11.5f\n",name,par.x,par.y,par.z);\
  DnProcLogAppend(im,D1strTMP);\
}

#define DnProcLogAddD3intPar(im,name,par){\
  sprintf(D1strTMP,"%-42s%11Ld%11Ld%11Ld\n",name,par.x,par.y,par.z);\
  DnProcLogAppend(im,D1strTMP);\
}

#define DnProcLogAddD2intPar(im,name,par){\
  sprintf(D1strTMP,"%-42s%11Ld%11Ld\n",name,par.x,par.y);\
  DnProcLogAppend(im,D1strTMP);\
}

#define DnProcLogAddD2floatPar(im,name,par){\
  sprintf(D1strTMP,"%-42s%11.5f%11.5f\n",name,par.x,par.y);\
  DnProcLogAppend(im,D1strTMP);\
}

#define DnProcLogAddIntPar(im,name,val){\
  sprintf(D1strTMP,"%-42s%11Ld\n",name,val);\
  DnProcLogAppend(im,D1strTMP);\
}

#define DnProcLogAddBooleanYNPar(im,name,val){\
  sprintf(D1strTMP,"%-42s%11s\n",name,val?"yes":"no");\
  DnProcLogAppend(im,D1strTMP);\
}

#define DnProcLogAddBooleanTFPar(im,name,val){\
  sprintf(D1strTMP,"%-42s%11s\n",name,val?"true":"false");\
  DnProcLogAppend(im,D1strTMP);\
}

#define DnProcLogAddStringPar(im,name,val){\
  sprintf(D1strTMP,"%-30s%-50s\n",name,val);\
  DnProcLogAppend(im,D1strTMP);\
}

#define DnProcLogAddComment(im,str){\
  sprintf(D1strTMP,"! %s\n",str);\
  DnProcLogAppend(im,D1strTMP);\
}

#define DnProcLogAddSeparator(im){\
  DnProcLogAppend(im,"!-------------------------------------------------------------------------------\n");\
}

#define DnProcLogPrint(im){\
  for( D1iTMP1=0; D1iTMP1<strlen((im).proc_log); D1iTMP1++ )\
    printf("%c",(im).proc_log[D1iTMP1]);\
  printf("\n");\
}

/*	  
**  D2 Types
*/	  

typedef struct { 
	D1int64			x,y;
    } D2int64,D2iVector_64,D2iValue_64,D2iPoint_64,D2int,D2iVector,D2iValue,D2iPoint;

typedef struct { 
	D1int32			x,y;
	} D2int32,D2int_32,D2iVector_32,D2iValue_32,D2iPoint_32;

typedef struct { 
	D1float			x,y;
	} D2float,D2fVector,D2fValue,D2fPoint;

typedef struct { 
	D1double			x,y;
	} D2double,D2dVector,D2dValue,D2dPoint;

typedef struct { 
	D1char			*x,*y;
} D2charP;

typedef struct { 
	D2iPoint_64		pos;
	D2iValue_64		size;
	} D2SubImage_64,D2SubImage;

typedef struct { 
	D2iPoint_32		pos;
	D2iValue_32		size;
	} D2SubImage_32;

typedef struct { 
	D2int64			pos;
	D2int64			dim;
} D2iRect_64,D2iRect;

typedef struct { 
	D2int32			pos;
	D2int32			dim;
} D2iRect_32;

typedef struct {
	D1short			*dat;
	D2iPoint_64		pos;
	D2iValue_64		dim;
	D2iValue_64		off;
	} D2Image_64,D2ShortImage_64,D2Image,D2ShortImage;

typedef struct {
	D1short	    *dat; 
	D2iPoint_32		pos;
	D2iValue_32		dim;
	D2iValue_32		off;
	} D2Image_32,D2ShortImage_32;

typedef struct {
	D1int64			*dat;
	D2iPoint_64		pos;
	D2iValue_64		dim;
	D2iValue_64		off;
	} D2IntImage_64,D2IndexImage_64,D2IntImage,D2IndexImage;

typedef struct {
	D1int32	    *dat;
	D2iPoint_32		pos;
	D2iValue_32		dim;
	D2iValue_32		off;
	} D2IntImage_32,D2IndexImage_32;

typedef struct {
	D1short			*dat; 
	D2iPoint_64		pos;
	D2iValue_64		dim;
	D2iValue_64		off;
	CTFileHeader	*header;	
	} D2CTSlice_64,D2CTSlice;

typedef struct {
	D1short		*dat;
	D2iPoint_32		pos;
	D2iValue_32		dim;
	D2iValue_32		off;
	CTFileHeader	*header;	
	} D2CTSlice_32;

typedef struct {
	D1Bin			*dat;
	D2iPoint_64		pos;
	D2iValue_64		dim;
	D2iValue_64		off;
	} D2BinImage_64,D2UCharImage_64,D2BinImage,D2UCharImage;

typedef struct {
	D1Bin		*dat;
	D2iPoint_32		pos;
	D2iValue_32		dim;
	D2iValue_32		off;
	} D2BinImage_32,D2UCharImage_32;

typedef struct {
	D1float			*dat;
	D2iPoint_64		pos;
	D2iValue_64		dim;
	D2iValue_64		off;
	} D2FloatImage_64,D2FloatImage;

typedef struct {
	D1float		*dat;
	D2iPoint_32		pos;
	D2iValue_32		dim;
	D2iValue_32		off;
	} D2FloatImage_32;

typedef struct {
	D1char			*dat;
	D2iPoint_64		pos;
	D2iValue_64		dim;
	D2iValue_64		off;
	} D2CharImage_64,D2CharImage;

typedef struct {
	D1char			*dat;
	D2iPoint_32		pos;
	D2iValue_32		dim;
	D2iValue_32		off;
	} D2CharImage_32;


/* Global temporary variables, should only be used in macros */

  static  D1char	D1cTMP;
  static  D1char	*D1cTMPp;
  static  D1float	D1fTMP;
  static  int64 D1iTMP1;
  static  int64	D1iTMP2;
  static  int64	D1iTMP3;

/* Global constants */

  static D2float		D2fNULL = {0.,0.}; 
  static D2int64		D2iNULL = {0,0}; 
  static D2int32		D2iNULL32 = {0,0}; 
  static D2int64		D2iONE	= {1,1}; 
  static D2int32		D2iONE32	= {1,1}; 
  static D2int64		D2iTWO	= {2,2}; 
  static D2int32		D2iTWO32	= {2,2}; 

/* Global temporary variables, should only be used in macros */

  static D2float		D2fTMP;
/**carfully check the use of this counter and its size because it is used in
the malloc function!**/
  static D1char			*D1cpTMP;
  static D1char			D1strTMP[256];

#ifdef __cplusplus
}  /* Close scope of 'extern "C"' declaration which encloses file. */
#endif

#endif	    /* _d2_h */
/* DON'T ADD ANYTHING AFTER THIS #endif */
