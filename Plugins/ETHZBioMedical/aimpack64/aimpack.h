/*==============================================================================

aimpack.h (see aimpack.c for documentation)

Steve Boyd
Aug 10, 2001
IBT/ETH, Moussonstrasse 18, CH-8044, Zuerich

==============================================================================*/

#ifndef _AIMPACK_H
#define _AIMPACK_H

#ifdef __cplusplus
extern "C" {
#endif

#ifdef AIMPACK_SILENT
void silent(char* str,...) {}
#define Aimpack_printf silent
#else
#define Aimpack_printf printf
#endif


#include "d2.h"
#include "d3.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/*------------------------------------------------------------------------------
Special definitions that are system dependent.
------------------------------------------------------------------------------*/
/* VMS -----------------------------------------------------------------------*/
#ifdef __VMS
#define AIMPACK_EXPORT

#define PLATFORM_LITTLE_ENDIAN_VMS
#define PLATFORM_32BIT    0
#define PLATFORM_2GB 2147483647

#include "dp_time.h"
#include <file.h>
#include <cvtdef.h> 		/*included for VAX to IEEE float conversion */
int CVT$FTOF( void *dat,int a,void *out,int b,int c); /* tested by Andres */
#include <lib$routines.h>   /*allocating more than 2GB memory */
#include <ssdef.h>          /*allocating more than 2GB memory */

#define PRId64 "Ld"
#define PRIu64 "Lu"

#define PLATFORM_FTOF_VAX(dat) {\
    CVT$FTOF((unsigned int *) (dat), CVT$K_VAX_F, (unsigned int *) (dat), CVT$K_IEEE_S,0);\
}
#define PLATFORM_FTOF_IEEE(dat){\
    CVT$FTOF((void *) (dat), CVT$K_IEEE_S, (void *) (dat), CVT$K_VAX_F,0);\
}
#if __IEEE_FLOAT == 1
#define PLATFORM_VAX_TO_I3E_FLOAT(image030){ (image030)->type = D1Ti3efloat;}
#define PLATFORM_I3E_TO_VAX_FLOAT(image030){ (image030)->type = D1Tfloat; }
#else
#define PLATFORM_VAX_TO_I3E_FLOAT(image030){\
    AimpackConvertAimFloat2I3EFloat(image030);\
}
#define PLATFORM_I3E_TO_VAX_FLOAT(image030){\
    AimpackConvertAimI3EFloat2Float(image030);\
} 
#endif

#define PLATFORM_FILE_POINTER FILE *fp;
#define PLATFORM_FILE_OPEN_WRITE (fp = fopen(filename,"wb","rfm=fix",alq,bls,mrs)) == NULL
#define PLATFORM_FILE_OPEN_READ (fp = fopen(filename,"rb","rfm=fix")) == NULL
#define PLATFORM_FILE_CLOSE fclose(fp);
#define PLATFORM_FILE_WRITE_C_MB fwrite(c_mb,rec_size,1,fp) != 1
#define PLATFORM_FILE_WRITE_BUF  fwrite(buf,rec_size,1,fp) != 1
#define PLATFORM_FILE_READ_C_MB  fread(c_mb,rec_size,1,fp) != 1
#define PLATFORM_FILE_READ_BUF   fread(buf,rec_size,1,fp) != 1

#define PLATFORM_SET_REC_SIZE {\
	int			f_desc;\
	f_desc = fileno(fp);\
	fstat(f_desc, &st);\
	rec_size  = (int) st.st_fab_mrs;\
}

#define PLATFORM_STRUCT_STAT struct stat		st;

#define PLATFORM_DP_TIMES      DP_times		time;
#define PLATFORM_DP_TIMES_INIT  _DP_TimeInit(&time);
#define PLATFORM_DP_TIMES_GET   _DP_TimeGet(&time);
#define PLATFORM_DP_TIMES_COMPUTE     *rate_mb_s = (float) AimpackMemBlockSize(mblist)/time.tot.elap/(1<<20);
#endif


/* SUN -----------------------------------------------------------------------*/
#ifdef sun
#define AIMPACK_EXPORT

#define PLATFORM_BIG_ENDIAN
#define PLATFORM_32BIT    0
#define PLATFORM_2GB 2147483647

#include <inttypes.h> 		/*file used for PRId64 and other types definition on 
							  VMS and Solaris */

#define PLATFORM_FTOF_VAX(dat) {}
#define PLATFORM_FTOF_IEEE(dat){} 
#define PLATFORM_VAX_TO_I3E_FLOAT(image030){ (image030)->type = D1Ti3efloat;}
#define PLATFORM_I3E_TO_VAX_FLOAT(image030){ (image030)->type = D1Tfloat; }

#define PLATFORM_FILE_POINTER FILE *fp;
#define PLATFORM_FILE_OPEN_WRITE (fp = fopen(filename,"wb")) == NULL
#define PLATFORM_FILE_OPEN_READ (fp = fopen(filename,"rb")) == NULL
#define PLATFORM_FILE_CLOSE fclose(fp);
#define PLATFORM_FILE_WRITE_C_MB fwrite(c_mb,rec_size,1,fp) != 1
#define PLATFORM_FILE_WRITE_BUF  fwrite(buf,rec_size,1,fp) != 1
#define PLATFORM_FILE_READ_C_MB  fread(c_mb,rec_size,1,fp) != 1
#define PLATFORM_FILE_READ_BUF   fread(buf,rec_size,1,fp) != 1

#define PLATFORM_SET_REC_SIZE {\
	stat(filename,&st);\
	rec_size  = (int) Aimpack_RECORD_SIZE;\
}

#define PLATFORM_STRUCT_STAT struct stat		st;

#define PLATFORM_DP_TIMES      
#define PLATFORM_DP_TIMES_INIT 
#define PLATFORM_DP_TIMES_GET  
#define PLATFORM_DP_TIMES_COMPUTE     *rate_mb_s = (float) 0.0;
#endif

/* LINUX ---------------------------------------------------------------------*/
#ifdef linux

#define AIMPACK_EXPORT

#define PLATFORM_LITTLE_ENDIAN_WIN32
#define PLATFORM_2GB 2147483647
#if __WORDSIZE == 64
#  define PLATFORM_32BIT    0
#else 
#  define PLATFORM_32BIT    1
#endif

/* include file used for PRId64 and other types defininion on Linux */
#define __STDC_FORMAT_MACROS 1
#include <inttypes.h> 
 
#define PLATFORM_FTOF_VAX(dat) {}
#define PLATFORM_FTOF_IEEE(dat){} 
#define PLATFORM_VAX_TO_I3E_FLOAT(image030){ (image030)->type = D1Ti3efloat;}
#define PLATFORM_I3E_TO_VAX_FLOAT(image030){ (image030)->type = D1Tfloat; }

#define PLATFORM_FILE_POINTER FILE *fp;
#define PLATFORM_FILE_OPEN_WRITE (fp = fopen(filename,"wb")) == NULL
#define PLATFORM_FILE_OPEN_READ (fp = fopen(filename,"rb")) == NULL
#define PLATFORM_FILE_CLOSE fclose(fp);
#define PLATFORM_FILE_WRITE_C_MB fwrite(c_mb,rec_size,1,fp) != 1
#define PLATFORM_FILE_WRITE_BUF  fwrite(buf,rec_size,1,fp) != 1
#define PLATFORM_FILE_READ_C_MB  fread(c_mb,rec_size,1,fp) != 1
#define PLATFORM_FILE_READ_BUF   fread(buf,rec_size,1,fp) != 1

#define PLATFORM_SET_REC_SIZE {\
	stat(filename,&st);\
	rec_size  = (int) Aimpack_RECORD_SIZE;\
}

#define PLATFORM_STRUCT_STAT struct stat st;

#define PLATFORM_DP_TIMES      
#define PLATFORM_DP_TIMES_INIT 
#define PLATFORM_DP_TIMES_GET  
#define PLATFORM_DP_TIMES_COMPUTE     *rate_mb_s = (float) 0.0;
#endif

/* SGI -----------------------------------------------------------------------*/
#ifdef sgi
#define AIMPACK_EXPORT

#define PLATFORM_BIG_ENDIAN
#define PLATFORM_32BIT    1
#define PLATFORM_2GB 2147483647

#define PLATFORM_FTOF_VAX(dat) {}
#define PLATFORM_FTOF_IEEE(dat){} 
#define PLATFORM_VAX_TO_I3E_FLOAT(image030){ (image030)->type = D1Ti3efloat;}
#define PLATFORM_I3E_TO_VAX_FLOAT(image030){ (image030)->type = D1Tfloat; }

#define PLATFORM_FILE_POINTER FILE *fp;
#define PLATFORM_FILE_OPEN_WRITE (fp = fopen(filename,"wb")) == NULL
#define PLATFORM_FILE_OPEN_READ (fp = fopen(filename,"rb")) == NULL
#define PLATFORM_FILE_CLOSE fclose(fp);
#define PLATFORM_FILE_WRITE_C_MB fwrite(c_mb,rec_size,1,fp) != 1
#define PLATFORM_FILE_WRITE_BUF  fwrite(buf,rec_size,1,fp) != 1
#define PLATFORM_FILE_READ_C_MB  fread(c_mb,rec_size,1,fp) != 1
#define PLATFORM_FILE_READ_BUF   fread(buf,rec_size,1,fp) != 1

#define PLATFORM_SET_REC_SIZE {\
	stat(filename,&st);\
	rec_size  = (int) Aimpack_RECORD_SIZE;\
}

#define PLATFORM_STRUCT_STAT struct stat		st;

#define PLATFORM_DP_TIMES      
#define PLATFORM_DP_TIMES_INIT 
#define PLATFORM_DP_TIMES_GET  
#define PLATFORM_DP_TIMES_COMPUTE     *rate_mb_s = (float) 0.0;
#endif

/* CRAY ----------------------------------------------------------------------*/
#ifdef cray
#define AIMPACK_EXPORT

#define PLATFORM_BIG_ENDIAN
#define PLATFORM_32BIT    0
#define PLATFORM_2GB 2147483647

#include <ffio.h>
#include <fcntl.h>

#define PLATFORM_FTOF_VAX(dat) {}
#define PLATFORM_FTOF_IEEE(dat){} 
#define PLATFORM_VAX_TO_I3E_FLOAT(image030){ (image030)->type = D1Ti3efloat;}
#define PLATFORM_I3E_TO_VAX_FLOAT(image030){ (image030)->type = D1Tfloat; }

#define PLATFORM_FILE_POINTER int fp;
#define PLATFORM_FILE_OPEN_WRITE (fp = ffopen(filename,O_WRONLY | O_CREAT, 0644)) == -1
#define PLATFORM_FILE_OPEN_READ (fp = ffopen(filename,O_RDONLY)) == -1
#define PLATFORM_FILE_CLOSE ffclose(fp);
#define PLATFORM_FILE_WRITE_C_MB ffwrite(fp,c_mb,rec_size) != rec_size
#define PLATFORM_FILE_WRITE_BUF  ffwrite(fp,buf,rec_size) != rec_size
#define PLATFORM_FILE_READ_C_MB  ffread(fp,c_mb,rec_size) != rec_size
#define PLATFORM_FILE_READ_BUF   ffread(fp,buf,rec_size) != rec_size

#define PLATFORM_SET_REC_SIZE {\
	stat(filename,&st);\
	rec_size  = (int) Aimpack_RECORD_SIZE;\
}

#define PLATFORM_STRUCT_STAT struct stat		st;

#define PLATFORM_DP_TIMES      
#define PLATFORM_DP_TIMES_INIT 
#define PLATFORM_DP_TIMES_GET  
#define PLATFORM_DP_TIMES_COMPUTE     *rate_mb_s = (float) 0.0;
#endif

/* WIN32 ---------------------------------------------------------------------*/
#ifdef WIN32
/* For generating a DLL in windows */
#include "aimpack64DllConfig.h"
#ifdef AIMPACK64_DLL
#define AIMPACK_EXPORT __declspec( dllexport )	/* both .cxx and .h needed */
#else
#define AIMPACK_EXPORT __declspec( dllimport )	/* only .h need.  Implementation .cxx from dll */
#endif

#define PLATFORM_LITTLE_ENDIAN_WIN32
#define PLATFORM_32BIT    1
#define PLATFORM_2GB 2147483647

#define PRId64 "ld"
#define PRIu64 "lu"

#define PLATFORM_FTOF_VAX(dat) {}
#define PLATFORM_FTOF_IEEE(dat){} 
#define PLATFORM_VAX_TO_I3E_FLOAT(image030){ (image030)->type = D1Ti3efloat;}
#define PLATFORM_I3E_TO_VAX_FLOAT(image030){ (image030)->type = D1Tfloat; }

#define PLATFORM_FILE_POINTER FILE *fp
#define PLATFORM_FILE_OPEN_WRITE (fp = fopen(filename,"wb")) == NULL
#define PLATFORM_FILE_OPEN_READ (fp = fopen(filename,"rb")) == NULL
#define PLATFORM_FILE_CLOSE fclose(fp);
#define PLATFORM_FILE_WRITE_C_MB fwrite(c_mb,rec_size,1,fp) != 1
#define PLATFORM_FILE_WRITE_BUF  fwrite(buf,rec_size,1,fp) != 1
#define PLATFORM_FILE_READ_C_MB  fread(c_mb,rec_size,1,fp) != 1
#define PLATFORM_FILE_READ_BUF   fread(buf,rec_size,1,fp) != 1

#define PLATFORM_SET_REC_SIZE {\
	stat(filename,&st);\
	rec_size  = (int) Aimpack_RECORD_SIZE;\
}

#define PLATFORM_STRUCT_STAT struct stat st

#define PLATFORM_DP_TIMES      
#define PLATFORM_DP_TIMES_INIT 
#define PLATFORM_DP_TIMES_GET  
#define PLATFORM_DP_TIMES_COMPUTE     *rate_mb_s = (float) 0.0;
#endif


/*------------------------------------------------------------------------------
Lists defined for storing the four blocks that comprise a single AIM file.
------------------------------------------------------------------------------*/
/*  Single memory block  */
typedef  struct {
	char		      *adr;	/*  pointer to memory block        */
	int64		      size;	/*  size of memory block in bytes  */ 
} AimpackMemBlock;

/*  List of memory blocks  */
typedef  struct {
	int32		      nr;	/*  nr of memory blocks in list    */
	AimpackMemBlock        *mb;	/*  vector of AimpackMemBlock's       */
} AimpackMemBlockList;

/*------------------------------------------------------------------------------
Set the block size for reading and writing the AIM file.
------------------------------------------------------------------------------*/
#define Aimpack_RECORD_SIZE        512    /* = n512, n <= 128   */

/*  New string included in pre-pre-header in file first ! */
static char  version030_string[] = "AIMDATA_V030   ";  /*  15 char plus \0 */


/*------------------------------------------------------------------------------
Platform dependencies:
Big endian are: AIX, HP-UX, IBM mainframe, Macintosh, and Solaris.
Little endian are: AXP/VMS, Digital UNIX, Intel ABI, OS/2, VAX/VMS, Windows.
NOTE: Windows32 float order is different than other little endian machines, hence
there is a special definition for windows32 platforms that fixes floats, but leaves
integers alone.
------------------------------------------------------------------------------*/
#ifdef PLATFORM_BIG_ENDIAN
#define little_to_big_endian_SwapBytes(bytes,nbytes) AimpackRmSwapByteOrder((char *)&(bytes),nbytes);
#define little_to_big_endian_ConvertFloat(data,forward) AimpackRmConvertFloatVMS2Sun(&data,forward);
#define little_to_big_endian_SwapImageHeader010(image,forward) {\
	little_to_big_endian_SwapBytes((image)->id,sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->ref,sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->type,sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->dim.x,sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->dim.y,sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->dim.z,sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->off.x,sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->off.y,sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->off.z,sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->subdim.x,sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->subdim.y,sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->subdim.z,sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->pos.x,sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->pos.y,sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->pos.z,sizeof(int32)); \
	little_to_big_endian_ConvertFloat((image)->el_size_mm,forward); \
	little_to_big_endian_SwapBytes((image)->assoc.id, sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->assoc.nr, sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->assoc.size, sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->assoc.type, sizeof(int32)); \
}
#define little_to_big_endian_SwapImageHeader011(image,forward) {\
	little_to_big_endian_SwapBytes((image)->id,sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->ref,sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->type, sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->dim.x, sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->dim.y, sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->dim.z, sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->off.x, sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->off.y, sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->off.z, sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->subdim.x, sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->subdim.y, sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->subdim.z, sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->pos.x, sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->pos.y, sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->pos.z, sizeof(int32)); \
	little_to_big_endian_ConvertFloat((image)->el_size_mm.x,forward); \
	little_to_big_endian_ConvertFloat((image)->el_size_mm.y,forward); \
	little_to_big_endian_ConvertFloat((image)->el_size_mm.z,forward); \
	little_to_big_endian_SwapBytes((image)->assoc.id, sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->assoc.nr, sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->assoc.size, sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->assoc.type, sizeof(int32)); \
}
#define little_to_big_endian_SwapImageHeader020(image,forward) {\
	little_to_big_endian_SwapBytes((image)->id, sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->ref, sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->type, sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->pos.x, sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->pos.y, sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->pos.z, sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->dim.x, sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->dim.y, sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->dim.z, sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->off.x, sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->off.y, sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->off.z, sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->supdim.x, sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->supdim.y, sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->supdim.z, sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->suppos.x, sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->suppos.y, sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->suppos.z, sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->subdim.x, sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->subdim.y, sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->subdim.z, sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->testoff.x, sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->testoff.y, sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->testoff.z, sizeof(int32)); \
	little_to_big_endian_ConvertFloat((image)->el_size_mm.x,forward); \
	little_to_big_endian_ConvertFloat((image)->el_size_mm.y,forward); \
	little_to_big_endian_ConvertFloat((image)->el_size_mm.z,forward); \
	little_to_big_endian_SwapBytes((image)->assoc.id, sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->assoc.nr, sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->assoc.size, sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->assoc.type, sizeof(int32)); \
}

#define little_to_big_endian_SwapImageHeader030(image,forward) {\
	little_to_big_endian_SwapBytes((image)->id, sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->ref, sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->type, sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->pos.x, sizeof(int64)); \
	little_to_big_endian_SwapBytes((image)->pos.y, sizeof(int64)); \
	little_to_big_endian_SwapBytes((image)->pos.z, sizeof(int64)); \
	little_to_big_endian_SwapBytes((image)->dim.x, sizeof(int64)); \
	little_to_big_endian_SwapBytes((image)->dim.y, sizeof(int64)); \
	little_to_big_endian_SwapBytes((image)->dim.z, sizeof(int64)); \
	little_to_big_endian_SwapBytes((image)->off.x, sizeof(int64)); \
	little_to_big_endian_SwapBytes((image)->off.y, sizeof(int64)); \
	little_to_big_endian_SwapBytes((image)->off.z, sizeof(int64)); \
	little_to_big_endian_SwapBytes((image)->supdim.x, sizeof(int64)); \
	little_to_big_endian_SwapBytes((image)->supdim.y, sizeof(int64)); \
	little_to_big_endian_SwapBytes((image)->supdim.z, sizeof(int64)); \
	little_to_big_endian_SwapBytes((image)->suppos.x, sizeof(int64)); \
	little_to_big_endian_SwapBytes((image)->suppos.y, sizeof(int64)); \
	little_to_big_endian_SwapBytes((image)->suppos.z, sizeof(int64)); \
	little_to_big_endian_SwapBytes((image)->subdim.x, sizeof(int64)); \
	little_to_big_endian_SwapBytes((image)->subdim.y, sizeof(int64)); \
	little_to_big_endian_SwapBytes((image)->subdim.z, sizeof(int64)); \
	little_to_big_endian_SwapBytes((image)->testoff.x, sizeof(int64)); \
	little_to_big_endian_SwapBytes((image)->testoff.y, sizeof(int64)); \
	little_to_big_endian_SwapBytes((image)->testoff.z, sizeof(int64)); \
	little_to_big_endian_ConvertFloat((image)->el_size_mm.x,forward); \
	little_to_big_endian_ConvertFloat((image)->el_size_mm.y,forward); \
	little_to_big_endian_ConvertFloat((image)->el_size_mm.z,forward); \
	little_to_big_endian_SwapBytes((image)->assoc.id, sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->assoc.nr, sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->assoc.size, sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->assoc.type, sizeof(int32)); \
}

#define little_to_big_endian_SwapImageFileHeader030(image,forward) {\
	little_to_big_endian_SwapBytes((image)->id, sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->ref, sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->type, sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->pos.x, sizeof(int64)); \
	little_to_big_endian_SwapBytes((image)->pos.y, sizeof(int64)); \
	little_to_big_endian_SwapBytes((image)->pos.z, sizeof(int64)); \
	little_to_big_endian_SwapBytes((image)->dim.x, sizeof(int64)); \
	little_to_big_endian_SwapBytes((image)->dim.y, sizeof(int64)); \
	little_to_big_endian_SwapBytes((image)->dim.z, sizeof(int64)); \
	little_to_big_endian_SwapBytes((image)->off.x, sizeof(int64)); \
	little_to_big_endian_SwapBytes((image)->off.y, sizeof(int64)); \
	little_to_big_endian_SwapBytes((image)->off.z, sizeof(int64)); \
	little_to_big_endian_SwapBytes((image)->supdim.x, sizeof(int64)); \
	little_to_big_endian_SwapBytes((image)->supdim.y, sizeof(int64)); \
	little_to_big_endian_SwapBytes((image)->supdim.z, sizeof(int64)); \
	little_to_big_endian_SwapBytes((image)->suppos.x, sizeof(int64)); \
	little_to_big_endian_SwapBytes((image)->suppos.y, sizeof(int64)); \
	little_to_big_endian_SwapBytes((image)->suppos.z, sizeof(int64)); \
	little_to_big_endian_SwapBytes((image)->subdim.x, sizeof(int64)); \
	little_to_big_endian_SwapBytes((image)->subdim.y, sizeof(int64)); \
	little_to_big_endian_SwapBytes((image)->subdim.z, sizeof(int64)); \
	little_to_big_endian_SwapBytes((image)->testoff.x, sizeof(int64)); \
	little_to_big_endian_SwapBytes((image)->testoff.y, sizeof(int64)); \
	little_to_big_endian_SwapBytes((image)->testoff.z, sizeof(int64)); \
	little_to_big_endian_SwapBytes((image)->el_size_nano.x,sizeof(int64)); \
	little_to_big_endian_SwapBytes((image)->el_size_nano.y,sizeof(int64)); \
	little_to_big_endian_SwapBytes((image)->el_size_nano.z,sizeof(int64)); \
	little_to_big_endian_SwapBytes((image)->assoc.id, sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->assoc.nr, sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->assoc.size, sizeof(int32)); \
	little_to_big_endian_SwapBytes((image)->assoc.type, sizeof(int32)); \
}

#define little_to_big_endian_SwapImageData(image) {\
	if ((image)->type == D1Tshort) { \
	int    i,size; \
	short  *dat; \
	dat  = (short *) (image)->dat; \
	size = (image)->dim.x*(image)->dim.y*(image)->dim.z; \
	for (i=0;i<size;i++,dat++) \
	little_to_big_endian_SwapBytes(*dat,sizeof(short)); \
	} \
}
#endif

#ifdef PLATFORM_LITTLE_ENDIAN_VMS
#define little_to_big_endian_SwapBytes(bytes,nbytes) {}
#define little_to_big_endian_ConvertFloat(data,forward) {}
#if __IEEE_FLOAT == 1
# define little_to_big_endian_SwapImageHeader010(image,forward) {\
    if(forward) {\
        PLATFORM_FTOF_VAX(&((image)->el_size_mm));\
    }\
    if(!forward) {\
        PLATFORM_FTOF_IEEE(&((image)->el_size_mm));\
    }\
}
# define little_to_big_endian_SwapImageHeader011(image,forward) {\
    if(forward) {\
        PLATFORM_FTOF_VAX(&((image)->el_size_mm.x));\
        PLATFORM_FTOF_VAX(&((image)->el_size_mm.y));\
        PLATFORM_FTOF_VAX(&((image)->el_size_mm.z));\
    }\
    if(!forward) {\
        PLATFORM_FTOF_IEEE(&((image)->el_size_mm.x));\
        PLATFORM_FTOF_IEEE(&((image)->el_size_mm.y));\
        PLATFORM_FTOF_IEEE(&((image)->el_size_mm.z));\
    }\
}
# define little_to_big_endian_SwapImageHeader020(image,forward) {\
    if(forward) {\
        PLATFORM_FTOF_VAX(&((image)->el_size_mm.x));\
        PLATFORM_FTOF_VAX(&((image)->el_size_mm.y));\
        PLATFORM_FTOF_VAX(&((image)->el_size_mm.z));\
    }\
    if(!forward) {\
        PLATFORM_FTOF_IEEE(&((image)->el_size_mm.x));\
        PLATFORM_FTOF_IEEE(&((image)->el_size_mm.y));\
        PLATFORM_FTOF_IEEE(&((image)->el_size_mm.z));\
    }\
}
# define little_to_big_endian_SwapImageHeader030(image,forward) {\
    if(forward) {\
        PLATFORM_FTOF_VAX(&((image)->el_size_mm.x));\
        PLATFORM_FTOF_VAX(&((image)->el_size_mm.y));\
        PLATFORM_FTOF_VAX(&((image)->el_size_mm.z));\
    }\
    if(!forward) {\
        PLATFORM_FTOF_IEEE(&((image)->el_size_mm.x));\
        PLATFORM_FTOF_IEEE(&((image)->el_size_mm.y));\
        PLATFORM_FTOF_IEEE(&((image)->el_size_mm.z));\
    }\
}
#else
# define little_to_big_endian_SwapImageHeader010(image,forward) {}
# define little_to_big_endian_SwapImageHeader011(image,forward) {}
# define little_to_big_endian_SwapImageHeader020(image,forward) {}
# define little_to_big_endian_SwapImageHeader030(image,forward) {}
#endif
# define little_to_big_endian_SwapImageFileHeader030(image,forward) {}
# define little_to_big_endian_SwapImageData(image) {}
#endif

#ifdef PLATFORM_LITTLE_ENDIAN_WIN32
#define little_to_big_endian_SwapBytes(bytes,nbytes) {}
#define little_to_big_endian_ConvertFloat(data,forward) AimpackSkbConvertFloatVMS2Win32(&data,forward);
#define little_to_big_endian_SwapImageHeader010(image,forward) {\
	little_to_big_endian_ConvertFloat((image)->el_size_mm,forward); \
}
#define little_to_big_endian_SwapImageHeader011(image,forward) {\
	little_to_big_endian_ConvertFloat((image)->el_size_mm.x,forward); \
	little_to_big_endian_ConvertFloat((image)->el_size_mm.y,forward); \
	little_to_big_endian_ConvertFloat((image)->el_size_mm.z,forward); \
}
#define little_to_big_endian_SwapImageHeader020(image,forward) {\
	little_to_big_endian_ConvertFloat((image)->el_size_mm.x,forward); \
	little_to_big_endian_ConvertFloat((image)->el_size_mm.y,forward); \
	little_to_big_endian_ConvertFloat((image)->el_size_mm.z,forward); \
}
#define little_to_big_endian_SwapImageHeader030(image,forward) {\
	little_to_big_endian_ConvertFloat((image)->el_size_mm.x,forward); \
	little_to_big_endian_ConvertFloat((image)->el_size_mm.y,forward); \
	little_to_big_endian_ConvertFloat((image)->el_size_mm.z,forward); \
}

# define little_to_big_endian_SwapImageFileHeader030(image,forward) {}
#define little_to_big_endian_SwapImageData(image) {}
#endif


/*
 * ANSI function declarations 
*/

/* 
 * WRITE FUNCTIONS: 
 * Write an AIM file to disk based on data stored in D3AnyImage
 * structure.  AimpackWriteImage is the base function that calls all the 
 * other necessary functions for writing the AIM to file.  
*/
AIMPACK_EXPORT void AimpackWriteImage(D3AnyImage030 *, char *, int *);
AIMPACK_EXPORT void AimpackWriteImageRate(D3AnyImage030 *, char *, int *, float *);
AIMPACK_EXPORT void AimpackWriteImageProcLog(D3AnyImage030 *, char *, int *);
AIMPACK_EXPORT void AimpackWriteImageData(D3AnyImage030 *, char *, int *);
AIMPACK_EXPORT void AimpackWriteImage020(D3AnyImage030 *, char *,int *);

/* 
 * READ FUNCTIONS: 
 * Read the AIM file from disk.  AimpackReadImage or AimpackReadImageUncompress
 * are the base functions that call all the necessary other functions for 
 * reading an AIM file from disk. 
*/
AIMPACK_EXPORT void AimpackReadImage(D3AnyImage030 *, char *, int *);
AIMPACK_EXPORT void AimpackReadImageRate(D3AnyImage030 *, char *, int *, float *);
AIMPACK_EXPORT void AimpackReadImageUncompress(D3AnyImage030 *, char *, int *);
AIMPACK_EXPORT void AimpackReadImageUncompressRate(D3AnyImage030 *, char *, int *, float *);
AIMPACK_EXPORT void AimpackReadImageInfo(D3AnyImage030 *, char *, int *);
AIMPACK_EXPORT void AimpackReadImage020(D3AnyImage030 *, char *, int *);
AIMPACK_EXPORT void AimpackReadImage011(D3AnyImage011 *, char *, int *);
AIMPACK_EXPORT void AimpackReadImage010(D3AnyImage010 *, char *, int *);

/* 
 * SET MEM BLOCK: 
 * Initialize the four memory blocks prior to writing with AimpackWriteMemBlock.
 * This is a pointer structure containing the address and size (in bytes) of 
 * the four memory blocks.  
*/
AIMPACK_EXPORT int AimpackMemBlockSize(AimpackMemBlockList);
AIMPACK_EXPORT void AimpackSetImageMemBlock(D3AnyImage *, D3AnyFileImage *, AimpackMemBlockList *, int32);
AIMPACK_EXPORT void AimpackSetImageMemBlock020(D3AnyImage*,D3AnyFileImage020*,AimpackMemBlockList*,int32);

/* SET FILE IMAGE FROM IMAGE */
AIMPACK_EXPORT void AimpackSetFileImage(D3AnyImage, D3AnyFileImage *);
AIMPACK_EXPORT void AimpackSetFileImage020(D3AnyImage, D3AnyFileImage020 *);

/* SET IMAGE FROM FILE IMAGE */
AIMPACK_EXPORT void AimpackSetImageFromFileImage(D3AnyFileImage030,D3AnyImage030*);
AIMPACK_EXPORT void AimpackSetImageFromFileImage020(D3AnyFileImage020,D3AnyImage030 *);


/* 
 * RESTORE MEMORY BLOCKS: 
 * Functions called after AIM data from file has been read into a 4 part 
 * memory block with AimpackReadMemBlock. 
*/
AIMPACK_EXPORT void AimpackRestoreImage030MemBlock(AimpackMemBlockList *mblist, D3AnyImage030 *image);
AIMPACK_EXPORT void AimpackRestoreImage020MemBlock(AimpackMemBlockList *,D3AnyImage030 *);
AIMPACK_EXPORT void AimpackRestoreImage011MemBlock(AimpackMemBlockList *,D3AnyImage030 *);
AIMPACK_EXPORT void AimpackRestoreImage010MemBlock(AimpackMemBlockList *,D3AnyImage030 *);

AIMPACK_EXPORT void AimpackRestoreImage030InfMemBlo(AimpackMemBlockList *,D3AnyImage030 *);
AIMPACK_EXPORT void AimpackRestoreImage020InfMemBlo(AimpackMemBlockList *,D3AnyImage030 *);
AIMPACK_EXPORT void AimpackRestoreImage011InfMemBlo(AimpackMemBlockList *,D3AnyImage030 *);
AIMPACK_EXPORT void AimpackRestoreImage010InfMemBlo(AimpackMemBlockList *,D3AnyImage030 *);

/* RESTORE:  
 * Legacy, functions called after ReadImage01x()
*/
AIMPACK_EXPORT void AimpackRestoreImage011(AimpackMemBlockList *, D3AnyImage011 *);
AIMPACK_EXPORT void AimpackRestoreImage010(AimpackMemBlockList *, D3AnyImage010 *);

/* 
 * CONVERT: 
 * Legacy. New using FileImage instead of Images 
*/
AIMPACK_EXPORT void AimpackConvertFile011ToFile020(D3AnyFileImage011, D3AnyFileImage020 *);
AIMPACK_EXPORT void AimpackConvertFile010ToFile011(D3AnyFileImage010, D3AnyFileImage011 *);

/* 
 * WRITE MEM BLOCKS: 
 * All image data is divided into four memory blocks and then written to disk. 
 * AimpackSetImageMemBlock must be called first to initialize the four memory 
 * blocks. Prior to writing the 4 blocks, a pre-header is written containing 
 * 5 integers: sizeof(pre-header), sizeof(memblock[0]), sizeof(memblock[1]),...
*/
AIMPACK_EXPORT void AimpackWriteMemBlock(AimpackMemBlockList, char *, int *);
AIMPACK_EXPORT void AimpackWriteMemBlockPart(AimpackMemBlockList, char *, int, int *);
AIMPACK_EXPORT void AimpackWriteMemBlock32bit(AimpackMemBlockList, char *,int *);

/* 
 * READ MEM BLOCKS:
 * All image data is read from the 4 memory blocks and stored into memory 
 * before calling AimpackRestoreImage020MemBlock to set data into the structure
 * D3AnyImage. The pre-header must be read first which contains the number 
 * of bytes of the four main memory blocks. 
*/
AIMPACK_EXPORT void AimpackReadMemBlock(AimpackMemBlockList *, char *, int *);
AIMPACK_EXPORT void AimpackReadMemBlockParts(AimpackMemBlockList *, int, char *, int *);

/* 
 * PROCEDURES: 
 * Run-length encoding.  Encoding is done before AimpackWriteImage and decoding
 *  is done after AimpackReadImage (or included in AimpackReadImageUncompress 
*/
AIMPACK_EXPORT void AimpackCompress(D3AnyImage030, D3AnyImage *, int);
AIMPACK_EXPORT void AimpackUncompress(D3AnyImage030 *, char);

/* 
 * SWAP BYTES: 
 * Used to convert between big/little endian. 
*/
AIMPACK_EXPORT void AimpackRmSwapByteOrder(char *, int);
AIMPACK_EXPORT void AimpackRmConvertFloatVMS2Sun(float *, int);
AIMPACK_EXPORT void AimpackSkbConvertFloatVMS2Win32(float *, int);

/* CONVERT AIM DATA IEEE <-> FLOAT */
AIMPACK_EXPORT void AimpackConvertAimFloat2I3EFloat(D3AnyImage030 *);
AIMPACK_EXPORT void AimpackConvertAimI3EFloat2Float(D3AnyImage030 *);


#ifdef __cplusplus
}
#endif

#endif


