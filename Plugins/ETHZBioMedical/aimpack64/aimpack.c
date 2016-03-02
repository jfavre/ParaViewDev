/*==VERSION 1.1=================================================================
aimpack.c

This module was written to put all the necessary files for reading and writing
Scanco AIM files into one program and header: aimpack.c, aimpack.h

The original files used to make this software include d3_io.c, th_io.c, th_io.h,
d3_compress.c, dp_time.c, dp_time.h.  If there are some procedures that were not
included in this package, then refer to the original files.  Note, however, that
many of the procedure names have been changed so that this module could be run in
parallel with the original procedures.  Therefore, it might take a bit of hunting
to find the original files.

Notes on usage:
1. Write an AIM
D3InitAny(), AimpackCompress(), AimpackWriteImage()
2. Read an AIM
D3InitAnyImage(), AimpackReadImageUncompress() or Aimpack
ReadImage()
3. Check the AIM header
D3PrintLogInfoAny()

Features and changes: 
1. Compatible for reading and writing on the following platforms:
Read/Write on Sun, SGI (not tested on SGI as of Dec., 2001)
Read/Write on VMS, Win32
Read/Write on Cray (64 bit operating system - can only read AIMs written on the cray!!!).
2. The original files used to include procedures that calculated IO time
for the reading and writing functions.  These functions have all been
removed due to the platform dependency hassles.
3. March 1, 2002.  I updated the prototypes in .h and the function declarations
in .cxx to include AIMPACK_EXPORT.  This is necessary for compiling on the PC
so that a DLL can be created.  On the PC, if using the skbdll.lib and 
skbdll.dll for your application (which includes aimpack), then you should
not define SKBDLL.  If you are not using skbdll.lib and skbdll.dll, then
your project should have SKBDLL defined.  If compiling on any other
platform, AIMPACK_EXPORT should have a blank definition (ie. #define AIMPACK_EXPORT).

WARNING: The float conversion on the SUN, SGI, and Win32 are good for the
range of floats that we expect (i.e., 0.034 or 1.00, etc).  The byte swap
on the floats and shorts may not have the bitwise change properly.

Float representation with four bytes, byte order according to platform:
VMS	  [0 1 2 3]
SUN   [3 2 1 0]
WIN32 [2 3 0 1]

Acknowledgements:
1. Most of this software was written by Tor Hildebrand for use on VMS.
2. Compatability for SUN systems by Ralph Mueller.
3. Compatability for Cray and Win32 by Steve Boyd.

Steve Boyd
Aug 10, 2001
IBT/ETH, Moussonstrasse 18, CH-8044, Zuerich

Version 1.1: (Revision 1) December 7, 2001. Steven Boyd
I added the compatability for Win32.  Also, I put all the #ifdef
statements in the aimpack.h file to make it easier to understand.
Version 1.2: (Revision 1) March 1, 2002. Steven Boyd
I had to include AIMPACK_EXPORT in the .h and .cxx file.  See comment 3
above in features and changes.  Note that d2.h and d3.h were also 
slightly modified to initialize voidP_fp.  Also removed prototypes 
from within functions (these shouldn't be there).

==============================================================================*/    

/*-----------------------------------------------------------------------------+
| Libraries
+-----------------------------------------------------------------------------*/    
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>

#include "d2.h"
#include "d3.h"
#include "aimpack.h"

/******************************************************************************/
/* WRITE AIM FUNCTIONS     */
/******************************************************************************/

/*	  
**  AimpackWriteImage() - write 030 aim
**  Calls all necessary function to write the contents of D3AnyImage structure.
**  First it initializes the memory block, then writes the content of the
**  memory block to disk.
*/

AIMPACK_EXPORT void AimpackWriteImage(D3AnyImage030 *image, char *filename, int *status)
{
	AimpackMemBlockList    mblist;
	D3AnyImage		ima;
	D3AnyFileImage	fileima;       

	/* if version of image is 020 the function tries to write a 020 aim 
	** else an 030 aim is written
	*/
	if(image->version == 020) {    
		Aimpack_printf("!>  Writing AimVersion020 requested... \n");
		AimpackWriteImage020(image,filename,status);
	} else { 
		/* its an 030 aim */
		
		/* Convert VAX float to IEEE float if it is a float image */
		if (image->type == D1Tfloat) {
			Aimpack_printf("!>  Converting Image from VAX float to IEEE float...\n");
			PLATFORM_VAX_TO_I3E_FLOAT(image);
		}

		ima = *image; /* struct copy */

		ima.version = 030;
		AimpackSetFileImage(ima,&fileima);
		AimpackSetImageMemBlock(&ima,&fileima,&mblist,1);
		AimpackWriteMemBlock(mblist,filename,status);
		free(mblist.mb);

		if (image->type == D1Ti3efloat) {
			Aimpack_printf("!>  After writing on disk, converting Image in memory back to VAX float...\n");
			PLATFORM_I3E_TO_VAX_FLOAT(image);
		}

	}
	Aimpack_printf("!>  Writing done!\n");
}

/*	  
**  AimpackWriteImageRate() - write 030 aim
**  Calls all necessary function to write the contents of D3AnyImage structure.
**  First it initializes the memory block, then writes the content of the
**  memory block to disk.
*/

AIMPACK_EXPORT void AimpackWriteImageRate(
    D3AnyImage030 *image, 
    char *filename, 
    int *status,
    float *rate_mb_s)
{
    AimpackMemBlockList    mblist;
    D3AnyImage		ima = *image;	/* struct copy */
    D3AnyFileImage	fileima;

    PLATFORM_DP_TIMES;
    PLATFORM_DP_TIMES_INIT;

    AimpackWriteImage(image,filename,status);

    PLATFORM_DP_TIMES_GET;    

  if( status ){ 
    AimpackSetFileImage(ima,&fileima);
    AimpackSetImageMemBlock(&ima,&fileima,&mblist,1);  
    PLATFORM_DP_TIMES_COMPUTE;
    free(mblist.mb);
  } else
    *rate_mb_s = 0.;

}

/*	
**  legacy 
**  AimpackWriteImage020()
**  writes a 020 version aim if image is small enough and image.type is not 
**  D1Tfloat
**  if it is not possible to write a 020 aim a 030 version is written
*/

AIMPACK_EXPORT void AimpackWriteImage020(D3AnyImage030 *image, char *filename, int *status)
{
	AimpackMemBlockList    mblist;  
	D3AnyImage		ima = *image;	/* struct copy */
	D3AnyFileImage020	fileima;

	/** Image bigger than 2GB write a 030 **/
	if ( D3MemorySize(ima) > PLATFORM_2GB) {
		/* implicitly no single dim is then bigger than signed int */
		Aimpack_printf("!>  Cannot write as AimVersion020, file bigger than 2 GB\n");
		Aimpack_printf("!>  Writing as current AIM 030 instead...\n");
		image->version = 030;
		AimpackWriteImage(image,filename,status);
		return;
	}

	/** Image type D1Tfloat: write a 030  **/
	if ( ima.type == D1Tfloat) {
		Aimpack_printf("!>  Writing current AIM Version 030 as IEEE floats.\n");
		image->version = 030;
		AimpackWriteImage(image,filename,status);
		return;
	}

	/** write a 020 if none of the above condition was true **/
	Aimpack_printf("!>  approved.\n");

	/* transfer D3AnyImage into D3AnyFileImage020 */
	AimpackSetFileImage020(ima,&fileima); 

	/* write size and adr for imagestruct, proc_log, data and assocdata
	** into the proper fields of the MemBlockList
	*/
	AimpackSetImageMemBlock020(&ima,&fileima,&mblist,1);

	/* write data to disk with help of the MemBlockList */
	AimpackWriteMemBlock32bit(mblist,filename,status);
	free(mblist.mb);
	Aimpack_printf("!>  Aim Version 020 written.\n");
}


/*	  
**  AimpackWriteImageProcLog()
**  Updates the image processing log.
*/	  

AIMPACK_EXPORT void AimpackWriteImageProcLog(D3AnyImage030 *image, char *filename, int *status)
{
	FILE*		fp;

	if ((fp = fopen(filename,"w")) == 0) {
		*status = 0;
		return;
	} else {
		*status = 1;
	}

	fprintf(fp,"%-30s   %7" PRId64 "  %7" PRId64 "  %7" PRId64 " \n","!> dim", image->dim.x,image->dim.y,image->dim.z);
    fprintf(fp,"%-30s   %7" PRId64 "  %7" PRId64 "  %7" PRId64 " \n","!> off", image->off.x,image->off.y,image->off.z);
	fprintf(fp,"%-30s   %7" PRId64 "  %7" PRId64 "  %7" PRId64 " \n","!> pos", image->pos.x,image->pos.y,image->pos.z);
	fprintf(fp,"%-30s   %7.4f %7.4f %7.4f\n","!> voxel size in mm", 
		image->el_size_mm.x,image->el_size_mm.y,image->el_size_mm.z);

	/*fprintf(fp,"!> %-28s%-10s %" PRId64 "  byte/voxel\n","Type of data",D3TypeName(*image),D1TSize(image->type));
	*/
	if( D3MemorySize(*image) >= (1<<20) )
		fprintf(fp,"!> %-28s%-10.1f Mbyte\n","Total memory size", (float)D3MemorySize(*image)/(1<<20));
	else   if( D3MemorySize(*image) >= (1<<10) )
		fprintf(fp,"!> %-28s%-10.1f Kbyte\n","Total memory size", (float)D3MemorySize(*image)/(1<<10));
	else
		fprintf(fp,"!> %-28s%-10" PRId64 "  byte\n","Total memory size", D3MemorySize(*image));

	if( fprintf(fp,"%s",image->proc_log) != (signed)strlen(image->proc_log) ) {
		*status = 0;
		fclose(fp);
		return;
	}

	fclose(fp);
	*status = 1;

}


/*	  
**  AimpackWriteImageData()
**  Updates the image data.
*/	  

AIMPACK_EXPORT void AimpackWriteImageData(D3AnyImage030 *image030, char *filename, int *status)
{
	AimpackMemBlockList    mblist;
	D3AnyImage030		   ima = *image030;	/* struct copy */
	D3AnyFileImage  	   fileima;

	AimpackSetFileImage(ima,&fileima);
	AimpackSetImageMemBlock(&ima,&fileima,&mblist,1);
	AimpackWriteMemBlockPart(mblist,filename,2,status);
	free(mblist.mb);
}


/******************************************************************************/
/* READ AIM FUNCTIONS   */
/******************************************************************************/

/*	  
 *  AimpackReadImage()
 */	  

AIMPACK_EXPORT void AimpackReadImage(D3AnyImage030 *image030, char *filename, int *status)
{
	AimpackMemBlockList    mblist;

	Aimpack_printf("!>  Reading memory blocks...\n"); 
	AimpackReadMemBlock(&mblist,filename,status); 
	Aimpack_printf("!>  Reading memory blocks done!\n"); 
	 
	if( *status ) {

		D3FreeAnyImage(*image030);
		
		/* check which version aim is and restore it with the proper function */
		if( mblist.mb[0].size == sizeof(D3AnyFileImage030) ){
			Aimpack_printf("!>  Restoring D3AnyImage030...\n"); 
			AimpackRestoreImage030MemBlock(&mblist,image030);  

			if (image030->type == D1Ti3efloat) {
				Aimpack_printf("!>  Converting IEEE float to VAX float...\n");
				PLATFORM_I3E_TO_VAX_FLOAT(image030);
			}
			Aimpack_printf("!>  Reading done!\n"); 
			return;
		}


		if( mblist.mb[0].size == sizeof(D3AnyFileImage020) ){
			Aimpack_printf("!>  Restoring D3AnyImage020...\n");
			AimpackRestoreImage020MemBlock(&mblist,image030);
			image030->version = 020;
			Aimpack_printf("!>  Restoring done!\n");  
			return;
		}

		if( mblist.mb[0].size == sizeof(D3AnyFileImage011) ){
			Aimpack_printf("!>  Restoring D3AnyImage011...\n");
			AimpackRestoreImage011MemBlock(&mblist,image030);
			image030->version = 020;
			Aimpack_printf("!>  Restoring done!\n");
			return;
		} 

		if( mblist.mb[0].size == sizeof(D3AnyFileImage010) ){
			Aimpack_printf("!>  Restoring D3AnyImage010...\n");
			AimpackRestoreImage010MemBlock(&mblist,image030);
			image030->version = 020;
			Aimpack_printf("!>  Restoring done!\n");
			return;
		} 

		/* If none of the above */
		Aimpack_printf("!>  AimpackReadImage: AIM version not recognized.\n");
		*status=0;
	}
	else {
		printf("!>  Error in AimpackReadMemBlock - returned status 0! \n");
	}

}

/*	  
**  AimpackReadImageRate()
**  Aim030 is distinguished by the header size, not the aim.version entry !
**  aim.version determines the Uncompression !
*/	  


AIMPACK_EXPORT void AimpackReadImageRate(
    D3AnyImage030 *image030, 
    char *filename, 
    int *status,
    float *rate_mb_s)
{

    AimpackMemBlockList    mblist;
    D3AnyImage		ima;
    D3AnyFileImage	fileima;

    PLATFORM_DP_TIMES;
    PLATFORM_DP_TIMES_INIT;

  AimpackReadImage(image030,filename,status);

    PLATFORM_DP_TIMES_GET; 

  if( status ){ 
    ima = *image030;
    AimpackSetFileImage(ima,&fileima);
    AimpackSetImageMemBlock(&ima,&fileima,&mblist,1);  
    PLATFORM_DP_TIMES_COMPUTE;
    free(mblist.mb);
  }
  else
    *rate_mb_s = 0.;

}

/*	  
**  AimpackReadImageUncompress()
**  This address the uncompression procedure after reading the AIM image.
**  If the AIM was not compressed, the AimpackUncompress procedure will do nothing.
*/	  

AIMPACK_EXPORT void AimpackReadImageUncompress(D3AnyImage030 *image, char *filename, int *status)
{ 
	AimpackReadImage(image,filename,status);

	if( *status ) {
		AimpackUncompress(image,0);
	}
}

/*	  
**  AimpackReadImageUncompressRate()
*/	  

AIMPACK_EXPORT void AimpackReadImageUncompressRate(
    D3AnyImage030 *image030, 
    char *filename, 
    int *status,
    float *rate_mb_s)
{

    AimpackMemBlockList    mblist;
    D3AnyImage		ima;
    D3AnyFileImage	fileima;

    PLATFORM_DP_TIMES;
    PLATFORM_DP_TIMES_INIT;

  AimpackReadImageUncompress(image030,filename,status);

    PLATFORM_DP_TIMES_GET;

  if( status ){ 
    ima = *image030;
    AimpackSetFileImage(ima,&fileima);
    AimpackSetImageMemBlock(&ima,&fileima,&mblist,1);  
    PLATFORM_DP_TIMES_COMPUTE;
    free(mblist.mb);
  }
  else
    *rate_mb_s = 0.;

}


/*	  
**  AimpackReadImageInfo()
**  Just reads enough of the AIM file (first two memory blocks) to process
**  the header.
*/	  

AIMPACK_EXPORT void AimpackReadImageInfo(D3AnyImage030 *image030, char *filename, int *status)
{
	AimpackMemBlockList    mblist;

	AimpackReadMemBlockParts(&mblist,2,filename,status);

	if( *status ) {

		D3FreeAnyImage(*image030);

		if( mblist.mb[0].size == sizeof(D3AnyFileImage030) ){
			AimpackRestoreImage030InfMemBlo(&mblist,image030);  
			return;
		}

		if( mblist.mb[0].size == sizeof(D3AnyFileImage020) ){
			AimpackRestoreImage020InfMemBlo(&mblist,image030);  
			return;
		}

		if( mblist.mb[0].size == sizeof(D3AnyFileImage011) ){
			AimpackRestoreImage011InfMemBlo(&mblist,image030);
			return;
		} 

		if( mblist.mb[0].size == sizeof(D3AnyFileImage010) ){
			AimpackRestoreImage010InfMemBlo(&mblist,image030);
			return;
		} 

		/*  If none of the above */
		Aimpack_printf("!>  AimpackReadImageInfo: aim type not recognized.\n");
		*status=0;

	}
}


/*	  
**  legacy
**  AimpackReadImage020()
*/	  

AIMPACK_EXPORT void AimpackReadImage020(D3AnyImage030 *image030, char *filename, int *status)
{
	AimpackMemBlockList    mblist;

	AimpackReadMemBlock(&mblist,filename,status);  

	if( *status ) {
		/* frees .dat, .proc_log and .assoc.dat fields */
		D3FreeAnyImage(*image030);

		if( mblist.mb[0].size == sizeof(D3AnyFileImage020) ){
			AimpackRestoreImage020MemBlock(&mblist,image030);  
			return;
		}
		if( mblist.mb[0].size == sizeof(D3AnyFileImage011) ){
			AimpackRestoreImage011MemBlock(&mblist,image030);
			return;
		} 

		if( mblist.mb[0].size == sizeof(D3AnyFileImage010) ){
			AimpackRestoreImage010MemBlock(&mblist,image030);
			return;
		} 
	}
}


/*	  
**  AimpackReadImage011()
**  Legacy.
*/	  

AIMPACK_EXPORT void AimpackReadImage011(D3AnyImage011 *image, char *filename, int *status)
{
	AimpackMemBlockList        mblist;

	AimpackReadMemBlock(&mblist,filename,status);

	if( *status ){
		D3FreeAnyImage(*image);
		AimpackRestoreImage011(&mblist,image);
	}

}

/*	  
**  AimpackReadImage010()
**  Legacy.
*/	  

AIMPACK_EXPORT void AimpackReadImage010(D3AnyImage010 *image, char *filename, int *status)
{
	AimpackMemBlockList        mblist;

	AimpackReadMemBlock(&mblist,filename,status);  

	if( *status ){
		D3FreeAnyImage(*image);
		AimpackRestoreImage010(&mblist,image); 
	}
}


/******************************************************************************/
/* GENERATE MEMBLOCKLIST STRUCT FUNCTIONS   */
/******************************************************************************/

/*	  
**  AimpackMemBlockSize()
*/	  

AIMPACK_EXPORT int AimpackMemBlockSize(
            AimpackMemBlockList            mblist)
{
    int64	  i,tot=0;

  for( i=0; i<mblist.nr; i++ ) {    
      tot += mblist.mb[i].size;
  }
  return(tot);

}

/*	  
**  AimpackSetImageMemBlock()
**  New from Version 3.0 on -> fileimage w/o pointers gets set to MemBlock
*/	  

AIMPACK_EXPORT void AimpackSetImageMemBlock(D3AnyImage *image, D3AnyFileImage *fileimage, AimpackMemBlockList *mblist, int32 swap)
{

	mblist->nr  = 4;
	Malloc(mblist->mb, mblist->nr, AimpackMemBlock);

	/* Assign mblock (the output structure) */

	/* mblist->mb[0] is the image structure: fileimage; it does not contain any
	 * pointer fields anymore!*/
	mblist->mb[0].adr    =  (char*) fileimage;\
	mblist->mb[0].size   =  sizeof(D3AnyFileImage);\

	/* mblist->mb[1] is the processing log string: image->proc_log */
	mblist->mb[1].adr    =  (char*) image->proc_log;
	mblist->mb[1].size   =  (image->proc_log != 0) ? strlen(image->proc_log) + 1: 0 ;

	/* mblist->mb[2] is the image data: image->dat */
	mblist->mb[2].adr    =  (char*) image->dat;
	mblist->mb[2].size   =  (image->dat != 0) ? D3MemorySize(*image) : 0;
	
	/* mblist->mb[3] is the image associated data: image->assoc */
	mblist->mb[3].adr    =  (char*) image->assoc.dat;
	mblist->mb[3].size   =  (image->assoc.dat != 0) ?
		image->assoc.nr*image->assoc.size : 0;

	/* the only numbers that remains unswaped until AimpackWriteMemBlock are
	 * the mb[i].size */
	if (swap) {
		little_to_big_endian_SwapImageData(image);
		/* version 030 FileImage has el_size_nano fields, the Image itselfe has
		 * el_size_mm fields instead */
		little_to_big_endian_SwapImageFileHeader030(fileimage,0);
	}
}


/*	  
**  AimpackSetImageMemBlock020()  
**  New from Version 3.0 on dummy pointers gets set to MemBlock
**  D3AnyFileImage is of type 020
*/	  

AIMPACK_EXPORT void AimpackSetImageMemBlock020(D3AnyImage *image,D3AnyFileImage020 *fileimage,AimpackMemBlockList *mblist,int32 swap)
{

	mblist->nr  = 4;
	Malloc(mblist->mb, mblist->nr, AimpackMemBlock);

	/* Assign mblock (the output structure) */

	/* mblist->mb[0] is the image structure: fileimage; the pointer fields are 
	 * set to zero!*/
	mblist->mb[0].adr    =  (char*) fileimage;\
	mblist->mb[0].size   =  sizeof(D3AnyFileImage020);\

	/* mblist->mb[1] is the processing log string: image->proc_log */
	mblist->mb[1].adr    =  (char*) image->proc_log;
	mblist->mb[1].size   =  (image->proc_log != 0) ?
		strlen(image->proc_log) + 1 : 0 ;

	/* mblist->mb[2] is the image data: image->dat */
	mblist->mb[2].adr    =  (char*) image->dat;
	mblist->mb[2].size   =  (image->dat != 0) ? D3MemorySize(*image) : 0;

	/* mblist->mb[3] is the image associated data: image->assoc */
	mblist->mb[3].adr    =  (char*) image->assoc.dat;
	mblist->mb[3].size   =  (image->assoc.dat != 0) ?
		image->assoc.nr*image->assoc.size : 0;

	if (swap) {
		little_to_big_endian_SwapImageData(image);
		little_to_big_endian_SwapImageHeader020(fileimage,0);
	}
}

/******************************************************************************/
/* GENERATE FILE IMAGE STRUCTS OF IMAGE STRUCTS AND VICE VERSA
/******************************************************************************/

/*	 
**  AimpackSetFileImage() Fill in FileImage without pointers from image
*/	 

AIMPACK_EXPORT void AimpackSetFileImage(D3AnyImage image, D3AnyFileImage *fileimage)
{

	int *dummy;

	/* Clear all fields of 4-byte aligned fileimage->version */
	dummy = (int *)(&(fileimage->version));
	*dummy= 0;

	fileimage->version  = image.version;
	fileimage->id	    = image.id;
	fileimage->ref      = image.ref;
	fileimage->type     = image.type;
	fileimage->pos      = image.pos;
	fileimage->dim      = image.dim;
	fileimage->off      = image.off;
	fileimage->supdim   = image.supdim;
	fileimage->suppos   = image.suppos;
	fileimage->subdim   = image.subdim;
	fileimage->testoff  = image.testoff;
	fileimage->el_size_nano.x   = Round(image.el_size_mm.x*1.E6);
	fileimage->el_size_nano.y   = Round(image.el_size_mm.y*1.E6);
	fileimage->el_size_nano.z   = Round(image.el_size_mm.z*1.E6);
	fileimage->assoc.id	    = image.assoc.id;
	fileimage->assoc.nr	    = image.assoc.nr;
	fileimage->assoc.size	    = image.assoc.size;
	fileimage->assoc.type	    = image.assoc.type;

}


/*	 
**  AimpackSetFileImage020() Fill in FileImage with dummy pointers from image
**  ! the integer fields in D3AnyImage are implicitly casted from 64 bit to 32
**  bit for the D3AnyFileImage struct !
*/	 

AIMPACK_EXPORT void AimpackSetFileImage020(D3AnyImage image, D3AnyFileImage020  *fileimage)
{
	int *dummy;

	/*  Clear all fields of 4-byte aligned fileimage->version */
	dummy = (int *)(&(fileimage->version));
	*dummy= 0;

	fileimage->version  = 020;
	fileimage->id	      = image.id;
	fileimage->ref      = image.ref;
	fileimage->type     = image.type;
	fileimage->pos.x          = image.pos.x;  /* implicit cast from int64 to int */
	fileimage->pos.y          = image.pos.y;
	fileimage->pos.z          = image.pos.z;
	fileimage->dim.x          = image.dim.x;
	fileimage->dim.y          = image.dim.y;
	fileimage->dim.z          = image.dim.z;
	fileimage->off.x          = image.off.x;
	fileimage->off.y          = image.off.y;
	fileimage->off.z          = image.off.z;
	fileimage->supdim.x       = image.supdim.x;
	fileimage->supdim.y       = image.supdim.y;
	fileimage->supdim.z       = image.supdim.z;
	fileimage->suppos.x       = image.suppos.x;
	fileimage->suppos.y       = image.suppos.y;
	fileimage->suppos.z       = image.suppos.z;
	fileimage->subdim.x       = image.subdim.x;
	fileimage->subdim.y       = image.subdim.y;
	fileimage->subdim.z       = image.subdim.z;
	fileimage->testoff.x      = image.testoff.x;
	fileimage->testoff.y      = image.testoff.y;
	fileimage->testoff.z      = image.testoff.z;
	fileimage->el_size_mm.x	    = image.el_size_mm.x   ;
	fileimage->el_size_mm.y	    = image.el_size_mm.y   ;
	fileimage->el_size_mm.z	    = image.el_size_mm.z   ;
	fileimage->assoc.id	    = image.assoc.id;
	fileimage->assoc.nr	    = image.assoc.nr;
	fileimage->assoc.size	    = image.assoc.size;
	fileimage->assoc.type	    = image.assoc.type;

	/*  setdummy pointers */
	fileimage->dat = 0;
	fileimage->proc_log = 0;
	fileimage->assoc.dat = 0;

}


/*	 
**  AimpackSetImageFromFileImage() 
**  Fill in Image from FileImage pointers 
*/	 

AIMPACK_EXPORT void AimpackSetImageFromFileImage(D3AnyFileImage030 fileimage,D3AnyImage030 *image)
{
	/* First swap the header in order to not get a segmentation fault */
	little_to_big_endian_SwapImageFileHeader030(&fileimage,1);

	image->version      = fileimage.version;
	image->id	        = fileimage.id;
	image->ref	        = fileimage.ref;
	image->type	        = fileimage.type;
	image->pos	        = fileimage.pos;
	image->dim	        = fileimage.dim;
	image->off	        = fileimage.off;
	image->supdim	    = fileimage.supdim;
	image->suppos	    = fileimage.suppos;
	image->subdim	    = fileimage.subdim;
	image->testoff      = fileimage.testoff;
	image->el_size_mm.x	= fileimage.el_size_nano.x/1.E6;
	image->el_size_mm.y	= fileimage.el_size_nano.y/1.E6;
	image->el_size_mm.z	= fileimage.el_size_nano.z/1.E6;
	image->assoc.id	    = fileimage.assoc.id;
	image->assoc.nr	    = fileimage.assoc.nr;
	image->assoc.size	= fileimage.assoc.size;
	image->assoc.type	= fileimage.assoc.type;
}

/*	 
**  AimpackSetImageFromFileImage020() 
**  Fill in Image from FileImage pointers 
*/	 

AIMPACK_EXPORT void AimpackSetImageFromFileImage020(D3AnyFileImage020 fileimage, D3AnyImage030 *image)
{
	
	/* First swap the header in order to not get a segmentation fault */
	little_to_big_endian_SwapImageHeader020(&fileimage,1);
		
	image->version      = fileimage.version;  /*  keep version 020 */ 
	/* so uncompress recognizes it ! A.Laib 26.4.1002 */
	image->id	      = fileimage.id;
	image->ref	      = fileimage.ref;
	image->type	      = fileimage.type;
	image->pos.x          = (int64)fileimage.pos.x;
	image->pos.y          = (int64)fileimage.pos.y;
	image->pos.z          = (int64)fileimage.pos.z;
	image->dim.x          = (int64)fileimage.dim.x;
	image->dim.y          = (int64)fileimage.dim.y;
	image->dim.z          = (int64)fileimage.dim.z;
	image->off.x          = (int64)fileimage.off.x;
	image->off.y          = (int64)fileimage.off.y;
	image->off.z          = (int64)fileimage.off.z;
	image->supdim.x       = (int64)fileimage.supdim.x;
	image->supdim.y       = (int64)fileimage.supdim.y;
	image->supdim.z       = (int64)fileimage.supdim.z;
	image->suppos.x       = (int64)fileimage.suppos.x;
	image->suppos.y       = (int64)fileimage.suppos.y;
	image->suppos.z       = (int64)fileimage.suppos.z;
	image->subdim.x       = (int64)fileimage.subdim.x;
	image->subdim.y       = (int64)fileimage.subdim.y;
	image->subdim.z       = (int64)fileimage.subdim.z;
	image->testoff.x      = (int64)fileimage.testoff.x;
	image->testoff.y      = (int64)fileimage.testoff.y;
	image->testoff.z      = (int64)fileimage.testoff.z;
	image->el_size_mm.x	    = fileimage.el_size_mm.x   ;
	image->el_size_mm.y	    = fileimage.el_size_mm.y   ;
	image->el_size_mm.z	    = fileimage.el_size_mm.z   ;
	image->assoc.id	    = fileimage.assoc.id	    ;
	image->assoc.nr	    = fileimage.assoc.nr	    ;
	image->assoc.size	    = fileimage.assoc.size	    ;
	image->assoc.type	    = fileimage.assoc.type	    ;
}

/******************************************************************************/
/* Generate Image struct    */
/******************************************************************************/

/* the restore functions are used to create a D3AnyFileImage030 out of a 
** memory block list mblist
*/

/*	  
**  AimpackRestoreImage030MemBlock()
*/	  

AIMPACK_EXPORT void AimpackRestoreImage030MemBlock(AimpackMemBlockList *mblist,D3AnyImage030 *image)
{
	D3AnyFileImage030              *image_p;

	/* mblist->mb[0] holds the image structure: image    */
	image_p = (D3AnyFileImage030 *) mblist->mb[0].adr;

	/* also swap data in that step*/
	AimpackSetImageFromFileImage(*image_p,image);
	free(image_p);

	/* mblist->mb[1] holds the processing log string: image->proc_log */
	image->proc_log = mblist->mb[1].adr; 

	/* mblist->mb[2] holds image data: image->data */
	image->dat = (void *) mblist->mb[2].adr; 

	/* mblist->mb[3] holds image associated data: image->assoc.dat */
	image->assoc.dat = (void *) mblist->mb[3].adr; 

	little_to_big_endian_SwapImageData(image);  
}


/*	  
**  AimpackRestoreImage020MemBlock()
**  Based on the four memory blocks, the image structure is restored.
*/	  

AIMPACK_EXPORT void AimpackRestoreImage020MemBlock(AimpackMemBlockList *mblist,D3AnyImage030 *image)
{
	D3AnyFileImage020          *image_p; 
	
	/* image_p points to the header data mb[0] */
	image_p = (D3AnyFileImage020 *) mblist->mb[0].adr;
    
	/* store the information in the D3AnyFileImage020 into a D3AnyImage030 */
	/* the information will be correctly swaped when comming out of this function */
	AimpackSetImageFromFileImage020(*image_p,image);
	free(image_p);

	/* Set pointers of D3AnyImage030 correctly */
	/* mblist->mb[1] holds the processing log string: image->proc_log */
	image->proc_log = mblist->mb[1].adr;

	/* mblist->mb[2] holds image data: image->data */
	image->dat = (void *) mblist->mb[2].adr;
	
	/* mblist->mb[3] holds image associated data: image->assoc.dat */
	image->assoc.dat = (void *) mblist->mb[3].adr; 
	little_to_big_endian_SwapImageData(image);
}


/*	  
**  AimpackRestoreImage011MemBlock()
**  Legacy.
*/	  
AIMPACK_EXPORT void AimpackRestoreImage011MemBlock(AimpackMemBlockList *mblist,D3AnyImage030 *image)
{
	D3AnyFileImage011              *image_p;
	D3AnyFileImage020              image020;

	/* mblist->mb[0] holds the image structure: image    */
	image_p = (D3AnyFileImage011 *) mblist->mb[0].adr;
	AimpackConvertFile011ToFile020(*image_p,&image020);	
	AimpackSetImageFromFileImage020(image020,image);
	free(image_p);

	/* mblist->mb[1] holds the processing log string: image->proc_log */
	image->proc_log = mblist->mb[1].adr; 

	/* mblist->mb[2] holds image data: image->data */
	image->dat = (void *) mblist->mb[2].adr; 

	/* mblist->mb[3] holds image associated data: image->assoc.dat */
	image->assoc.dat = (void *) mblist->mb[3].adr; 

	little_to_big_endian_SwapImageData(image);
}


/*	  
**  AimpackRestoreImage010MemBlock()
**  Legacy.
*/	  
/** New: - gives back a D3AnyImage020 struct
- The conversion is done here and handel with a FileImage (without 
pointers
*/

AIMPACK_EXPORT void AimpackRestoreImage010MemBlock(AimpackMemBlockList *mblist,D3AnyImage030 *image)
{
	D3AnyFileImage010              *image_p;
	D3AnyFileImage020              image020;
	D3AnyFileImage011              image011;

	/* mblist->mb[0] holds the image structure: image    */
	image_p = (D3AnyFileImage010 *) mblist->mb[0].adr;
	AimpackConvertFile010ToFile011(*image_p,&image011);
	AimpackConvertFile011ToFile020(image011,&image020);
	AimpackSetImageFromFileImage020(image020,image);
	free(image_p);

	/* mblist->mb[1] holds the processing log string: image->proc_log */
	image->proc_log = mblist->mb[1].adr; 

	/* mblist->mb[2] holds image data: image->data */
	image->dat = (void *) mblist->mb[2].adr; 

	/* mblist->mb[3] holds image associated data: image->assoc.dat */
	image->assoc.dat = (void *) mblist->mb[3].adr; 

	little_to_big_endian_SwapImageData(image);
}


/* RestoreImage01x() functions give back an 01x D3AnyImage struct, they are
** only used in ReadImage01x()
*/

AIMPACK_EXPORT void AimpackRestoreImage011(AimpackMemBlockList *mblist, D3AnyImage011 *image)
{
	D3AnyImage011              *image_p;

	image_p = (D3AnyImage011 *) mblist->mb[0].adr; 
	*image = *image_p;
	free(image_p);

	image->proc_log = mblist->mb[1].adr; 
	image->dat = (void *) mblist->mb[2].adr; 
	image->assoc.dat = (void *) mblist->mb[3].adr; 

	little_to_big_endian_SwapImageHeader011(image,1);
	little_to_big_endian_SwapImageData(image);

}

AIMPACK_EXPORT void AimpackRestoreImage010(AimpackMemBlockList *mblist, D3AnyImage010 *image)
{
	D3AnyImage010              *image_p;

	image_p = (D3AnyImage010 *) mblist->mb[0].adr;
	*image = *image_p;
	free(image_p);

	image->proc_log = mblist->mb[1].adr; 
	image->dat = (void *) mblist->mb[2].adr; 
	image->assoc.dat = (void *) mblist->mb[3].adr; 

	little_to_big_endian_SwapImageHeader010(image,1);
	little_to_big_endian_SwapImageData(image);

}


/*	  
**  AimpackRestoreImage030InfoMemBlock()
**  Restore just the first two memory blocks: image structure and 
**  processing log.
*/		  

AIMPACK_EXPORT void AimpackRestoreImage030InfMemBlo(AimpackMemBlockList *mblist, D3AnyImage030 *image)
{
	D3AnyFileImage030              *image_p;


	/* mblist->mb[0] holds the image structure: image    */
	image_p = (D3AnyFileImage030 *) mblist->mb[0].adr;

	AimpackSetImageFromFileImage(*image_p,image);
	free(image_p);

	/* mblist->mb[1] holds the processing log string: image->proc_log */
	image->proc_log = mblist->mb[1].adr; 

	image->dat = 0; 
	image->assoc.dat = 0; 

}



/*	  
**  AimpackRestoreImage020InfMemBlo()
*/	  


AIMPACK_EXPORT void AimpackRestoreImage020InfMemBlo(AimpackMemBlockList *mblist, D3AnyImage030 *image)
{
	D3AnyFileImage020          *image_p;

	/* image_p points to the header data mb[0] */
	image_p = (D3AnyFileImage020 *) mblist->mb[0].adr;

	/* store the information in the D3AnyFileImage020 into a D3AnyImage030 */
	AimpackSetImageFromFileImage020(*image_p,image);
	free(image_p);

	/* Set pointers of D3AnyImage030 correctly */
	/* mblist->mb[1] holds the processing log string: image->proc_log */
	image->proc_log = mblist->mb[1].adr; 

	/* the other memory blocks are not restored -> set pointers to zero */
	image->dat = 0;
	image->assoc.dat = 0;

}

/*	  
**  AimpackRestoreImage011InfMemBlo()
**  Legacy.
*/	  

AIMPACK_EXPORT void AimpackRestoreImage011InfMemBlo(AimpackMemBlockList *mblist, D3AnyImage030 *image)
{
	D3AnyFileImage011              *image_p;
	D3AnyFileImage020              image020;

	/* mblist->mb[0] holds the image structure: image    */
	image_p = (D3AnyFileImage011 *) mblist->mb[0].adr;

	AimpackConvertFile011ToFile020(*image_p,&image020);
	AimpackSetImageFromFileImage020(image020,image);
	free(image_p);

	/* mblist->mb[1] holds the processing log string: image->proc_log */
	image->proc_log = mblist->mb[1].adr; 

	image->dat = 0; 
	image->assoc.dat = 0; 
}

/*	  
**  AimpackRestoreImage010InfMemBlo()
**  Legacy.
*/	  

AIMPACK_EXPORT void AimpackRestoreImage010InfMemBlo(AimpackMemBlockList *mblist, D3AnyImage030 *image)
{
	D3AnyFileImage010              *image_p;
	D3AnyFileImage020              image020;
	D3AnyFileImage011              image011;

	/* mblist->mb[0] holds the image structure: image    */
	image_p = (D3AnyFileImage010 *) mblist->mb[0].adr;

	AimpackConvertFile010ToFile011(*image_p,&image011);
	AimpackConvertFile011ToFile020(image011,&image020);
	AimpackSetImageFromFileImage020(image020,image);
	free(image_p);

	/* mblist->mb[1] holds the processing log string: image->proc_log */
	image->proc_log = mblist->mb[1].adr; 
	image->dat = 0; 
	image->assoc.dat = 0; 
}


/******************************************************************************/
/* CONVERT OLD FILE IMAGE FORMAT IN NEW FILE IMAGE FORMAT     */
/******************************************************************************/

/*	  
**  AimpackConvert011To020()
**  Legacy.
*/	  
AIMPACK_EXPORT void AimpackConvertFile011ToFile020(D3AnyFileImage011 any011,D3AnyFileImage020 *any020)
{

	/* Convert */

	Aimpack_printf("!>  Converting AIM V1.1 to V2.0\n");

	D3InitAnyFileImage020(*any020,D1Tundef);
	
	any020->id           = any011.id;
	any020->ref          = any011.ref;
	any020->type         = any011.type;
	any020->pos          = any011.pos;
	any020->dim          = any011.dim;
	any020->off          = any011.off;
	any020->el_size_mm   = any011.el_size_mm;
	any020->assoc.id     = any011.assoc.id;
	any020->assoc.nr     = any011.assoc.nr;
	any020->assoc.size   = any011.assoc.size;
	any020->assoc.type   = any011.assoc.type;
	any020->version      = 020;
}  

/*	  
**  AimpackConvert010To011()
**  Legacy.
*/	  
/**Scanco version*/
AIMPACK_EXPORT void AimpackConvertFile010ToFile011(D3AnyFileImage010 any010, D3AnyFileImage011 *any011)
{

	/* Convert */

	Aimpack_printf("!>  Converting AIM V1.0 to V1.1\n");

	any011->id           = any010.id;
	any011->ref          = any010.ref;
	any011->type         = any010.type;
	any011->dim          = any010.dim;
	any011->off          = any010.off;
	any011->subdim       = any010.subdim;
	any011->pos          = any010.pos;
	any011->el_size_mm.x = any010.el_size_mm;
	any011->el_size_mm.y = any010.el_size_mm;
	any011->el_size_mm.z = any010.el_size_mm;
	any011->assoc.id     = any010.assoc.id;
	any011->assoc.nr     = any010.assoc.nr;
	any011->assoc.size   = any010.assoc.size;
	any011->assoc.type   = any010.assoc.type;
	any011->version      = 011;
}


/******************************************************************************/
/* WRITE MEMBLOCK TO DISK FUNCTIONS     */
/******************************************************************************/

/*	  
**  AimpackWriteMemBlock()
**  writes an aim file version 030 to disk 
**  with pre-header containing first a 16 byte identification string 
**  "AIMDATA_V030   "\0 followed by five 64-bit integers
**  The status flag 'status' is '1' for success and '0' for failed.
**  the data written to the file are accessed over AimpackMemBLockList
**  !The preheader size is indicated without the 16byte version string!
*/	  

AIMPACK_EXPORT void AimpackWriteMemBlock(AimpackMemBlockList mblock, char
										 *filename, int *status)
{
	AimpackMemBlock       *mb_ptr; /** points to the MemBlockList **/
	int64                 i;
	int32                 nr_mb;   /** # memory blocks **/
	int64                 mb_size;
	int64                 nr_rec;  /** # records to be written **/
	int64                 c_mb_size; /**size of data block that needs to be next written to file **/
	int32                 c_buf_size; /** # free bytes from 512 **/
	int64                 size;
	int32                 print, dbg; 
	int32                 rec_size; /** RECORD_SIZE **/
	int64                 tot_size; /** size of data **/
	int64                 *block_size;
	char		          *c_buf; /** pointer to the 1 free byte in the recordbuffer **/
	char                  *buf;  /** pointer to the beginning of the record buffer  **/
	char                  *c_mb; /** pointer to the data that needs to be next written to the file  **/ 
	char   bls[20]; /** Aimpack_RECORD_SIZE, block size**/
	char   mrs[20]; /** Aimpack_RECORD_SIZE, max record size **/
	char   alq[20]; /** # of blocks, 1block = 512bytes VMS, allocation quantity **/
	char   tmp[20];
	PLATFORM_FILE_POINTER;

	print=0;
	dbg=0;

	mb_ptr = mblock.mb;
	nr_mb  = mblock.nr; 

	Aimpack_printf("!>  Start Writing MemBlock 64 bit...\n");
	/* calculate tot_size=size of the preheader + size of the memory blocks 
	** from the mblist[i].size field
	*/
	tot_size = (nr_mb+1)*sizeof(int64)+16*sizeof(char); /*16bytes for version string AIMDATA_V030---\0 */
	for(i=0;i<nr_mb;i++) {
		tot_size += mb_ptr[i].size;
	}

	/* open file */
	nr_rec = tot_size/Aimpack_RECORD_SIZE + ((tot_size%Aimpack_RECORD_SIZE==0)?0:1);
	size   = nr_rec * Aimpack_RECORD_SIZE;  
	sprintf(bls,"bls = %d",Aimpack_RECORD_SIZE);
	sprintf(mrs,"mrs = %d",Aimpack_RECORD_SIZE);
	sprintf(alq,"alq = %" PRId64 "",size/512);          /* nr fblocks=512 byte (VMS)  */  

	if ( PLATFORM_FILE_OPEN_WRITE ){
		if (print) Aimpack_printf("!>  Can't open file %s for writing.\n",filename);
		*status=0; 
		return;
	}

	/* Allocate memory for buffer buf that saves the part of a memblock that  
	** is left over from the last 512 byte block that was written to disk 
	** (buffer buf is used if a block is not exactly 512 byte long and needs 
	** to be compounded with another block from the list
	*/
	rec_size = Aimpack_RECORD_SIZE;
	Malloc(buf, rec_size, char);
	c_buf = buf;
	c_buf_size = rec_size;

	size = (nr_mb+1)*sizeof(int64) + 16*sizeof(char);
	Malloc(block_size, size, int64);

	/* Allocate memory for the preheader and write the data into this memory 
	** space using the mblist[i].size fields
	**
	** First 16 bytes version030_string is written followed by the 64bit integer
	** referring to the preheader size itself, and the four 64 bit integer
	** referring to the header, proc_log, data and assoc_data block sizes
	** 
	** !the preheader size does not contain the 16 bytes for the version string!
	*/

	/*  NEW Version 030: Identifier String uses the first 16 bytes ! 5.6.2001 A.Laib */
	aimpack_memcpy(block_size,version030_string,16*sizeof(char) );

	block_size[2] = size - 16*sizeof(char);
	little_to_big_endian_SwapBytes(block_size[2],sizeof(int64)); /* size of preheader */
	for(i=0;i<nr_mb;i++) {
		block_size[i+3] = mb_ptr[i].size;
		little_to_big_endian_SwapBytes(block_size[i+3],sizeof(int64));
	}

	c_mb      = (char*) block_size; 
	c_mb_size = size;

	/* Write preheader to file, because it is smaller than 512 bytes it 
	** is copied first to the buffer buf in memory, the pointer c_buf is 
	** increased such that it points to the first free field in the buffer 
	** buf and the counter c_buf_size
	*/
	while (c_mb_size >= rec_size) {  
		if( PLATFORM_FILE_WRITE_C_MB ) {
			if (print) Aimpack_printf("!>  Can't write to file: %s.\n",filename);
			free(buf); 
			free(block_size);
			PLATFORM_FILE_CLOSE;
			*status=0; 
			return;
		}
		c_mb_size -= rec_size;
		c_mb      += rec_size;
	}  

	if( c_mb_size > 0 ){                                 
		aimpack_memcpy(c_buf,c_mb,c_mb_size);  
		c_buf_size = rec_size - c_mb_size;
		c_buf      = buf      + c_mb_size;
		c_mb_size  = 0;
	}

	/* Write data pointed by mblist[i].adr to the file.  
	** At the beginning buffer buf is filled with the data from the first 
	** memblock until 512 bytes are full and than written to file.
	** The following data (rest of block header, proc_log, data and assoc_data)
	**  are written in full 512 bytes block if possible otherwise buffer buf 
	** is used to compound a full 512 byte block from different memblocks.
	*/
	for(i=0;i<nr_mb;i++) {

		c_mb      = mb_ptr[i].adr;    
		c_mb_size = mb_ptr[i].size;   
		sprintf(tmp,"Working on block %" PRId64 " ",i);

		if( c_mb_size > 0 ){  /* ignore empty memory blocks */
			/* if buffer partly filled: copy current mb to the buffer  */

			if( c_buf_size < rec_size ){             /* buffer part. filled  */
				if( c_mb_size <= c_buf_size ){         /* c_mb fits completly in buffer */
					aimpack_memcpy(c_buf,c_mb,c_mb_size);  
					c_buf_size -= c_mb_size;
					c_buf      += c_mb_size;
					c_mb_size   = 0;
				} else {
					aimpack_memcpy(c_buf,c_mb,c_buf_size);       /* fill out buffer with c_mb */
					c_mb_size -= c_buf_size;
					c_mb      += c_buf_size;
					c_buf_size = 0;
				}
				if( c_buf_size == 0 ){                  /* buffer full: write to file */
					if( PLATFORM_FILE_WRITE_BUF ) {
						if (print) Aimpack_printf("!>  Can't write to file: %s. Size: %" PRId64 "  Bytes\n",filename,size);
						free(buf); 
						free(block_size);
						PLATFORM_FILE_CLOSE;
						*status=0; 
						return;
					} 
					c_buf_size = rec_size;
					c_buf      = buf;
				}
			}

			if( c_mb_size > 0 ){

				/* write memblock direct to file if possible */ 

				while (c_mb_size >= rec_size) {                   /* mid   */
					if( PLATFORM_FILE_WRITE_C_MB ) {
						if (print) Aimpack_printf("!>  Can't write to file: %s.\n",filename);
						free(buf); free(block_size);
						PLATFORM_FILE_CLOSE;
						*status=0; return;
					}
					c_mb_size -= rec_size;
					c_mb      += rec_size;
				}  

				if( c_mb_size > 0 ){                              /* tail   */
					aimpack_memcpy(c_buf,c_mb,c_mb_size);  
					c_buf_size = rec_size - c_mb_size;
					c_buf      = buf      + c_mb_size;
					c_mb_size  = 0;
				}
			}
		}
	} /* end for i */ 

	/* write buffer to file if not empty */

	if ( c_buf_size < rec_size ){
		if( PLATFORM_FILE_WRITE_BUF ) {
			if (print) Aimpack_printf("!>  Can't write to file: %s. Size: %" PRId64 "  Bytes\n",filename,size);
			free(buf); free(block_size);
			PLATFORM_FILE_CLOSE;
			*status=0; return;
		} 
	}

	free(buf);
	free(block_size);
	PLATFORM_FILE_CLOSE;
	*status = 1;
	Aimpack_printf("!>  AimVersion030 written.\n");
}


/*	  
**  AimpackWriteMemBlock32bit() 
**  functions writes an aimfile version 020 
**  with pre-header containing five 32-bit integers and no identifier string
**  The status flag 'status' is '1' for success and '0' for failed.
**  the data that are written to the file are accessed over AimpackMemBLockList
**
**  The first memory block contains the size of the header block (btyes)
**  and the 4 following blocks (i.e., 20, 140, 0, 3375, 0).
**  The remaining four blocks are then written to complete the AIM file.
**  The AIM file is written in 512 byte blocks, so if one of the list of blocks
**  doesn't complete a 512 block, then it is filled out with the next block.  If it
**  is the last of the list of blocks, then it is filled out with zeros.
*/	  

AIMPACK_EXPORT void AimpackWriteMemBlock32bit(AimpackMemBlockList mblock, char *filename, int *status)
{
	AimpackMemBlock       *mb_ptr; /** points to the MemBlockList **/
	int32                 i;
	int32                 nr_mb;   /** # memory blocks **/
	int32                 nr_rec;  /** # records to be written **/
	int32                 c_mb_size; /**size of data block that needs to be next written to file **/
	int32                 c_buf_size; /** # free bytes from 512 **/
	int32                 size, print, dbg; 
	int32                 rec_size; /** RECORD_SIZE **/
	int32                 tot_size; /** size of data **/
	int32                 *block_size;
	char		          *c_buf; /** pointer to the 1 free byte in the recordbuffer **/
	char                  *buf;  /** pointer to the beginning of the record buffer  **/
	char                  *c_mb; /** pointer to the data that needs to be next written to the file  **/ 
	char   bls[20]; /** Aimpack_RECORD_SIZE, block size**/
	char   mrs[20]; /** Aimpack_RECORD_SIZE, max record size **/
	char   alq[20]; /** # of blocks, 1block = 512bytes VMS, allocation quantity **/
	char   tmp[20];

	PLATFORM_FILE_POINTER;	
	print=0;
	dbg=0;
	Aimpack_printf("!>  Start writting MemBlock 32 bit...\n");
	mb_ptr = mblock.mb;
	nr_mb  = mblock.nr; 
	tot_size = (nr_mb+1)*sizeof(int32);

	
	for(i=0;i<nr_mb;i++) {
		if(print) Aimpack_printf("!>  memblock size to write %d = %" PRId64 "\n",i,mb_ptr[i].size);
		tot_size += mb_ptr[i].size;
	}

	nr_rec = tot_size/Aimpack_RECORD_SIZE + ((tot_size%Aimpack_RECORD_SIZE==0)?0:1);
	size   = nr_rec * Aimpack_RECORD_SIZE;

	sprintf(bls,"bls = %d",Aimpack_RECORD_SIZE);
	sprintf(mrs,"mrs = %d",Aimpack_RECORD_SIZE);
	sprintf(alq,"alq = %d",size/Aimpack_RECORD_SIZE);  /* nr fblocks=512 byte (VMS)  */  

	if (PLATFORM_FILE_OPEN_WRITE){
		if (print) Aimpack_printf("!>  Can't open file %s for writing.\n",filename);
		*status=0; return;
	}

	/* Allocate memory for buffer buf that saves the part of a memblock that 
	** is left over from the last 512 byte block that was written to disk
	*/
	rec_size = Aimpack_RECORD_SIZE;
	Malloc(buf, rec_size, char);
	c_buf = buf;
	c_buf_size = rec_size;

	/* Allocate memory for the preheader and write the data into this 
	** memory space using the mblist[i].size field
	*/
	size = (nr_mb+1)*sizeof(int32);
	Malloc(block_size, size, int32);

	block_size[0] = size;
	little_to_big_endian_SwapBytes(block_size[0],sizeof(int32));
	for(i=0;i<nr_mb;i++) {
		block_size[i+1] = mb_ptr[i].size;
		little_to_big_endian_SwapBytes(block_size[i+1],sizeof(int32));
	}

	/* Write preheader to file (five 32 bit integers = 20 byte)
	** if >= 512bytes write them down (what a strange test ;-) )
	** if < 512 bytes save them to the c_buf region in the memory and write them
	** down with the next full 512 byte block
	*/
	c_mb      = (char*) block_size; /** pointer to preheader **/
	c_mb_size = size; /** size of preheader, unswapped **/

	while (c_mb_size >= rec_size) {  
		if( PLATFORM_FILE_WRITE_C_MB ) {
			if (print) Aimpack_printf("!>  Can't write to file: %s.\n",filename);
			free(buf); 
			free(block_size);
			PLATFORM_FILE_CLOSE;
			*status=0; 
			return;
		}
		c_mb_size -= rec_size;
		c_mb      += rec_size;
	}  

	if( c_mb_size > 0 ){                               
		aimpack_memcpy(c_buf,c_mb,c_mb_size);  
		c_buf_size = rec_size - c_mb_size;
		c_buf      = buf      + c_mb_size;
		c_mb_size  = 0;
	}

	/* Write all memory blocks, pointed out by mblock, to file 
	** first c_buf is filled and than writte to file
	** memblocks are directelly written to file if 512byte big and c_buf 
	** currently empty
	** after the for loop the last part of datas is written to file  
	*/
	for(i=0;i<nr_mb;i++) {
		c_mb      = mb_ptr[i].adr;    
		c_mb_size = mb_ptr[i].size;   
		sprintf(tmp,"Working on block %d",i);

		if( c_mb_size > 0 ){  /* ignore empty memory blocks */

			/* if buffer partly filled: copy current mb to the buffer  */
			if( c_buf_size < rec_size ){  
				if( c_mb_size <= c_buf_size ){  /* current c_mb fits completlly into the buffer */
					aimpack_memcpy(c_buf,c_mb,c_mb_size);  
					c_buf_size -= c_mb_size;
					c_buf      += c_mb_size;
					c_mb_size   = 0;
				} else { /*just a part of c_mb suits into  c_buf */
					aimpack_memcpy(c_buf,c_mb,c_buf_size); 
					c_mb_size -= c_buf_size;
					c_mb      += c_buf_size;
					c_buf_size = 0;
				}

				if( c_buf_size == 0 ){  /* buffer full: write to file */
					if( PLATFORM_FILE_WRITE_BUF ) {
						if (print) Aimpack_printf("!>  Can't write to file: %s. Size: %d Bytes\n",filename,size);
						free(buf); free(block_size);
						PLATFORM_FILE_CLOSE;
						*status=0; return;
					} 
					c_buf_size = rec_size;
					c_buf      = buf;
				}
			} /**if( c_mb_size > 0 )**/

			if( c_mb_size > 0 ){
				/* write the rest of the memblock direct to file if possible */ 
				while (c_mb_size >= rec_size) {    /* full 512 byte part of memblock */
					if( PLATFORM_FILE_WRITE_C_MB ) {
						if (print) Aimpack_printf("!>  Can't write to file: %s.\n",filename);
						free(buf); 
						free(block_size);
						PLATFORM_FILE_CLOSE;
						*status=0; 
						return;
					}
					c_mb_size -= rec_size;
					c_mb      += rec_size;
				}  

				if( c_mb_size > 0 ){   /* Rest of memblock smaller than 512bytes */
					aimpack_memcpy(c_buf,c_mb,c_mb_size);  
					c_buf_size = rec_size - c_mb_size;
					c_buf      = buf      + c_mb_size;
					c_mb_size  = 0;
				}
			}
		}
	} /* end for i */ 

	/* write c_buf buffer to file if not empty */
	if ( c_buf_size < rec_size ){
		if( PLATFORM_FILE_WRITE_BUF ) {
			if (print) Aimpack_printf("!>  Can't write to file: %s. Size: %d Bytes\n",filename,size);
			free(buf); free(block_size);
			PLATFORM_FILE_CLOSE;
			*status=0; return;
		} 
	}
	free(buf);
	free(block_size);
	PLATFORM_FILE_CLOSE;
	*status = 1;
	Aimpack_printf("!>  Writting MemBlock 32 bit done!\n");
}


/*	  
**  AimpackWriteMemBlockPart
**
**  Writes a single memory block (mblock.mb[nr]) to a file with name 'filename'.
**  The status flag 'status' is '1' for success and '0' for failed.
*/	  

AIMPACK_EXPORT void AimpackWriteMemBlockPart(AimpackMemBlockList mblock, char *filename, int nr, int *status)
{
	AimpackMemBlock        *mb_ptr;
	int                 nr_rec;     /* number of 512 byte blocks */
	int64               c_mb_size;  /* size of the data still to be written */
	int                 c_buf_size; /* free space in buf */
	int64               size;       /* size of the resulting file */
	int                 rec_size;   /* 512 byte */
	int64               tot_size;   /* size of the memory block */
	char		        *c_buf;     /* pointer to the first free byte of buf */
	char                *buf;       /* 512byte memory buffer */
	char                *c_mb;      /* adr of the memblock to be written */
	int    print=0;
	char   bls[20];
	char   mrs[20];
	char   alq[20];
	PLATFORM_FILE_POINTER;

	mb_ptr = mblock.mb;  
	tot_size = mb_ptr[nr].size;     
	nr_rec = tot_size/Aimpack_RECORD_SIZE + ((tot_size%Aimpack_RECORD_SIZE==0)?0:1);
	size   = nr_rec * Aimpack_RECORD_SIZE;

	sprintf(bls,"bls = %d",Aimpack_RECORD_SIZE);
	sprintf(mrs,"mrs = %d",Aimpack_RECORD_SIZE);
	sprintf(alq,"alq = %" PRId64 " ",size/512);          /* nr fblocks=512 byte (VMS)  */  

	if ( PLATFORM_FILE_OPEN_WRITE ){
		if (print) Aimpack_printf("!>  Can't open file %s for writing.\n",filename);
		*status=0; return;
	}

	/* Allocate record buffer: buf  */

	rec_size = Aimpack_RECORD_SIZE;
	Malloc(buf, rec_size, char);
	c_buf = buf;
	c_buf_size = rec_size;


	/* Write single block to file */

	c_mb      = mb_ptr[nr].adr;    
	c_mb_size = mb_ptr[nr].size;     
	little_to_big_endian_SwapBytes(c_mb_size,sizeof(int64));

	while (c_mb_size >= rec_size) {  
		if( PLATFORM_FILE_WRITE_C_MB ) {
			if (print) Aimpack_printf("!>  Can't write to file: %s.\n",filename);
			free(buf);
			PLATFORM_FILE_CLOSE;
			*status=0; return;
		}
		c_mb_size -= rec_size;
		c_mb      += rec_size;
	}  

	if( c_mb_size > 0 ){                                  /* tail   */
		aimpack_memcpy(c_buf,c_mb,c_mb_size);  
		c_buf_size = rec_size - c_mb_size;
		c_buf      = buf      + c_mb_size;
		c_mb_size  = 0;
	}

	/* write buffer to file if not empty */

	if ( c_buf_size < rec_size )
		if( PLATFORM_FILE_WRITE_BUF ) {
			if (print) Aimpack_printf("!>  Can't write to file: %s. Size: %" PRId64 "  Bytes\n",filename,size);
			free(buf);
			PLATFORM_FILE_CLOSE;
			*status=0; return;
		} 

		free(buf); 
		PLATFORM_FILE_CLOSE;
		*status = 1;

}


/******************************************************************************/
/* READ MEMBLOCK INTO MEMORY FUNCTIONS  */
/******************************************************************************/



/*	  
**  AimpackReadMemBlock
**
**  Reads a list of memory blocks ('mblock') from a file with name 'filename'.
**  The first 20 bytes contain 5 integers.  First integer is normally 20 (block size
**  itself), followed by the size in bytes of the four following blocks in the list.
**  
**  The status flag 'status' is '1' for success and '0' for failed.
*/	  

AIMPACK_EXPORT void AimpackReadMemBlock(AimpackMemBlockList *mblock, char *filename, int *status)
{
	AimpackMemBlock	*mb_ptr;
	int64            i,j;
	int32            nr_mb; /* number of memblocks in the file */
	int64            head_mb_size; /* size of preheader */
	int64            c_mb_size;  /* size of the currently to read buffer */
	int32            c_buf_size; /* size of not yet readed data in buf */
	int64            size;
	int32            rec_size; /* record size of the file */
	int64            tot_size; 
	int64            *block_size64; /* pointer to preheader buffer 64bit */
	int32            *block_size32; /* pointer to preheader buffer 32bit */
	int			     file_32bit_flag=0; 
	int			     file_64bit_flag=0;
	char		     *c_buf; /* pointer to the first not yet read char in buf */
	char             *buf; /* buffer in memory to which the data is read */
	char             *c_mb; /* pointer to the first empty space in the currently read buffer */

	int			print= 0;
	int32 firstint = 0; /* used for swap in file_32bit_flag setting */


	PLATFORM_STRUCT_STAT;
	PLATFORM_FILE_POINTER;
	
	if ( PLATFORM_FILE_OPEN_READ ){
		if (print) Aimpack_printf("!>  Can't open file %s for reading\n",filename);
		*status=0; return;
	}


	PLATFORM_SET_REC_SIZE;

	/* Allocate record buffer: buf  */
	Malloc(buf, rec_size, char);
	c_buf = buf;
	c_buf_size = rec_size;

	/* Read first record to buffer buf 
	** The first record will contain the preheader data
	*/
	if( PLATFORM_FILE_READ_BUF ) {
		if (print) Aimpack_printf("!>  Can't read from file: %s.\n",filename);
		PLATFORM_FILE_CLOSE;
		*status=0; 
		return;
	} 
	c_buf_size  = rec_size;

	/* check if the first 2bytes = version030_string
	** and initialises the 64 resp. 32 bit flag
	*/

	if (strncmp(buf,version030_string,15)==0) {
		c_buf += 16;   /*  Skip first 16 characters */
		c_buf_size -= 16;
		file_64bit_flag = 1;
		head_mb_size = *((int64 *)c_buf);
		little_to_big_endian_SwapBytes(head_mb_size,sizeof(int64));
		if(print) Aimpack_printf("!>  64bit: head_mb_size swapped = %" PRId64 " \n",head_mb_size); 
	} else {
		/*  Aims up to Version 020 have 20-byte or 16-byte 'pre'-header */
		firstint = *((int32*)buf);
		little_to_big_endian_SwapBytes(firstint,sizeof(int32));
		if(firstint <= 20){
				file_32bit_flag = 1;
				head_mb_size = (int64) firstint;
				if(print) Aimpack_printf("!>  32bit: head_mb_size = %" PRId64 " \n",head_mb_size); 
		} else {
			Aimpack_printf("!>  File neither 32bit version nor AIM_V030. Return.\n"); 
			*status = 0;
			return;
		}
	}
	

	/* allocate buffer for preheader  */
	Malloc(block_size64, head_mb_size, int64);

	if(file_32bit_flag){ /* AIMS up to version 020 */
		block_size32 = (int32 *) block_size64;
		nr_mb        = head_mb_size/sizeof(int32) - 1;   /* do not count header */
		if(print) Aimpack_printf("!>  32bit: nr_mb = %d\n",nr_mb);
	} else{
		nr_mb        = head_mb_size/sizeof(int64) - 1;
		if(print) Aimpack_printf("!>  64bit: nr_mb = %d\n",nr_mb);
	} 
	

	/* allocate memory for MemBlockList */
	Malloc(mb_ptr, nr_mb, AimpackMemBlock);
	mblock->mb = mb_ptr;
	mblock->nr = nr_mb;
	
	/* set copy target pointer to the preheader buffer */
	c_mb_size = head_mb_size;
	if(file_32bit_flag){
		c_mb = (char*) block_size32;
	}else{
		c_mb = (char*) block_size64;
	}

	/* copy preheader data to its proper buffer block_size64/32 */         
	if( c_mb_size <= rec_size ) {     /* all in buffer */
		aimpack_memcpy(c_mb,c_buf,c_mb_size);  
		c_buf_size -= c_mb_size;
		c_buf      += c_mb_size;
		c_mb_size = 0;  
	} 
	else {               
		aimpack_memcpy(c_mb,c_buf,c_buf_size);          /* head */  
		c_mb_size -= c_buf_size;
		c_mb      += c_buf_size;
		while (c_mb_size >= rec_size) { /* read middle part direct in the preheader buffer */
			if( PLATFORM_FILE_READ_C_MB ) {
				if (print) Aimpack_printf("!>  Can't read from file: %s.\n",filename);
				PLATFORM_FILE_CLOSE;
				*status=0; return;
			}
			c_mb_size -= rec_size;
			c_mb      += rec_size;
		}  
		if( c_mb_size > 0 ){                       /* tail */
			if( PLATFORM_FILE_READ_BUF ) {
				if (print) Aimpack_printf("!>  Can't read from file: %s.\n",filename);
				PLATFORM_FILE_CLOSE;
				*status=0; return;
			}
			aimpack_memcpy(c_mb,buf,c_mb_size);  
			c_buf_size = rec_size - c_mb_size;
			c_buf      = buf      + c_mb_size;
		}
	}


	/* Test if the aimfile is bigger than 2GB and return if so on 32 bit
	** architectures
	*/
	for(i=0,tot_size=0;i<nr_mb;i++){ 
		if (file_32bit_flag) {
			little_to_big_endian_SwapBytes(block_size32[i+1],sizeof(int32));
			if(print) Aimpack_printf("!>  32bit: blocksize nr %" PRId64 " = %d\n",i+1,block_size32[i+1]);
			tot_size=tot_size+block_size32[i+1];
		} else { 
			little_to_big_endian_SwapBytes(block_size64[i+1],sizeof(int64));
			if(print) Aimpack_printf("!>  64bit: blocksize nr %" PRId64 " = %" PRId64 " \n",i+1,block_size64[i+1]);
			tot_size=tot_size+block_size64[i+1];
		}
	}


	if(PLATFORM_32BIT){
		if(tot_size>= PLATFORM_2GB) {
			Aimpack_printf("!>  File is bigger than 2GB - Read slicewise instead.Exit\n");
			free(buf); 
			free(block_size64);
			*status=0; 
			return;
		}
	}

	/* Allocate memory buffers for the memory blocks */
	for(i=0,tot_size=0;i<nr_mb;i++) {    
		if (file_32bit_flag) {
			/** bytes are already swaped in the previous loop**/
			mb_ptr[i].size = size = block_size32[i+1];
		} else {
			/** bytes are already swaped in the previous loop**/
			mb_ptr[i].size = size = block_size64[i+1];
		}     
		tot_size += size;
		if( size > 0 ) {
			/* check if memory could be allocated */
            MALLOC(mb_ptr[i].adr,size, char);
		} else { 
			mb_ptr[i].adr = 0;
		}
	}/* for */


	/* store data from the file in the prepared memory blocks fields
	** c_mb always points to the proper address mb_ptr[i].adr
	** c_mb_size indicates how big the data are
	*/
	for(i=0;i<nr_mb;i++) {
		c_mb      = mb_ptr[i].adr;
		c_mb_size = mb_ptr[i].size;

		if( c_mb_size > 0 ){             /* ignore empty memory blocks */
			if( c_mb_size <= c_buf_size ){             /* all in buffer */
				aimpack_memcpy(c_mb,c_buf,c_mb_size);
				c_buf_size -= c_mb_size;
				c_buf      += c_mb_size;
			} 
			else {  
				/* first copy head */
				aimpack_memcpy(c_mb,c_buf,c_buf_size);                       
				c_mb_size  -= c_buf_size;
				c_mb       += c_buf_size;
				c_buf_size  = 0;

				/* secondly copy middle part directly from file to buffer */ 
				while (c_mb_size >= rec_size) {               
					if( PLATFORM_FILE_READ_C_MB ) {
						if (print) Aimpack_printf("!>  Can't read from file: %s.\n",filename);
						PLATFORM_FILE_CLOSE;
						*status=0; 
						return;
					}
					c_mb_size -= rec_size;
					c_mb      += rec_size;
				} 

				/* finally copy rest to the buffer buf and therefrom the final bytes
				** that belong to this block
				*/ 
				if( c_mb_size > 0 ){                            
					if( PLATFORM_FILE_READ_BUF ) {
						if (print) Aimpack_printf("!>  Can't read from file: %s.\n",filename);
						PLATFORM_FILE_CLOSE;
						*status=0; 
						return;
					}
					aimpack_memcpy(c_mb,buf,c_mb_size);  
					c_buf_size = rec_size - c_mb_size;
					c_buf      = buf      + c_mb_size;
				}
			}
		}
	if(print) Aimpack_printf("!>  memblocks adr %" PRId64 " = %" PRId64 " \n",i,mb_ptr[i].adr );
	} /* end for i */

	free(buf); 
	free(block_size64);
	PLATFORM_FILE_CLOSE;
	*status=1; 
}    


/*	  
**  AimpackReadMemBlockParts
**
**  Reads the first 'max_nr' of memory blocks ('mblock') from a file 
**  with name 'filename'.
**  The status flag 'status' is '1' for success and '0' for failed.
**  Normally, max_nr is '2' which reads in the first two blocks (image structure
**  and processing log).
*/	  

AIMPACK_EXPORT void AimpackReadMemBlockParts(AimpackMemBlockList *mblock, int max_nr, char *filename, int *status)
{
	AimpackMemBlock        *mb_ptr;
	int64               i;
	int64               j;
	int32               nr_mb;
	int32               nr_mb_org;
	int64               head_mb_size; 
	int64               c_mb_size;
	int32               c_buf_size;
	int64               size; 
	int32               rec_size;
	int64               tot_size;
	int64               *block_size64;
	int32               *block_size32;
	int64               memory_offset;
	char		        *c_buf;
	char                *buf;
	char                *c_mb;

	short   			print=0;
	int32	    		file_32bit_flag=0;
	int32   			file_64bit_flag=0;
	int32 firstint = 0;

	PLATFORM_STRUCT_STAT;
	PLATFORM_FILE_POINTER;

	if ( PLATFORM_FILE_OPEN_READ ){
		if (print) Aimpack_printf("!>  Can't open file %s for reading\n",filename);
		*status=0; return;
	}

	PLATFORM_SET_REC_SIZE;

	/* Allocate record buffer: buf  */

	Malloc(buf, rec_size, char);
    c_buf = buf;
	c_buf_size = rec_size;

	/* Read first record to buffer */

	if( PLATFORM_FILE_READ_BUF ) {
		if (print) Aimpack_printf("!>  Can't read from file: %s.\n",filename);
		PLATFORM_FILE_CLOSE;
		*status=0; return;
	} 
	c_buf_size  = rec_size;

	/* check if the first 2bytes = version030_string
	** and initialises the 64 resp. 32 bit flag
	*/
	if (strncmp(buf,version030_string,15)==0) {
		c_buf += 16;   /*  Skip first 16 characters */
		c_buf_size -= 16;
		file_64bit_flag = 1;
		memory_offset = 16;
		head_mb_size = *((int64 *)c_buf);
		if(print) Aimpack_printf("!>  head_mb_size = %" PRId64 " \n",head_mb_size);
		little_to_big_endian_SwapBytes(head_mb_size,sizeof(int64));
		if(print) Aimpack_printf("!>  head_mb_size swapped = %" PRId64 " \n",head_mb_size); 
	} else {
		/*  Aims up to Version 020 have 20-byte or 16-byte 'pre'-header */
		firstint = *((int32*)buf);
		little_to_big_endian_SwapBytes(firstint,sizeof(int32));
		if(print) Aimpack_printf("!>  firstint swaped = %d\n",firstint);
		if(firstint <= 20){
				file_32bit_flag = 1;
				head_mb_size = (int64) firstint;
				if(print) Aimpack_printf("!>  head_mb_size (file_32bit_flag = 1) = %" PRId64 " \n",head_mb_size); 
		} else {
			Aimpack_printf("!>  File neither 32bit version nor AIM_V030. Return.\n"); 
			*status = 0;
			return;
		}
		memory_offset = 0;
	}

	/* allocate buffer for preheader  */
	Malloc(block_size64, head_mb_size, int64);
	memory_offset += head_mb_size;
	
	/* set 32 or 64 bit flag */
	if(file_32bit_flag){ /* AIMS up to version 020 */
		block_size32 = (int32 *) block_size64;
		nr_mb        = head_mb_size/sizeof(int32) - 1;   /* do not count header */
	} else{
		nr_mb        = head_mb_size/sizeof(int64) - 1;
	} 
	

	nr_mb_org = nr_mb;
	if (max_nr < nr_mb){
		nr_mb =  max_nr; 
	}

	Malloc(mb_ptr, nr_mb, AimpackMemBlock);

	mblock->mb = mb_ptr;
	mblock->nr = nr_mb;
	 
	/* set copy target pointer to the preheader buffer */
	c_mb_size = head_mb_size;
	if(file_32bit_flag){
		c_mb = (char*) block_size32;
	}else{
		c_mb = (char*) block_size64;
	}
	 

	if( c_mb_size <= rec_size ) {                   /* all in buffer */
		aimpack_memcpy(c_mb,c_buf,c_mb_size);  
		c_buf_size -= c_mb_size;
		c_buf      += c_mb_size;
		c_mb_size = 0;  
	} 
	else {               
		aimpack_memcpy(c_mb,c_buf,c_buf_size);                        /* head   */  
		c_mb_size -= c_buf_size;
		c_mb      += c_buf_size;
		while (c_mb_size >= rec_size) {                     /* mid   */
			if( PLATFORM_FILE_READ_C_MB ) {
				if (print) Aimpack_printf("!>  Can't read from file: %s.\n",filename);
				PLATFORM_FILE_CLOSE;
				*status=0; 
				return;
			}
			c_mb_size -= rec_size;
			c_mb      += rec_size;
		}  
		if( c_mb_size > 0 ){                                /* tail   */
			if( PLATFORM_FILE_READ_BUF ) {
				if (print) Aimpack_printf("!>  Can't read from file: %s.\n",filename);
				PLATFORM_FILE_CLOSE;
				*status=0; 
				return;
			}
			aimpack_memcpy(c_mb,buf,c_mb_size);  
			c_buf_size = rec_size - c_mb_size;
			c_buf      = buf      + c_mb_size;
		}
	}

	/* Allocate memory for the memory blocks */

	for(i=0,tot_size=0;i<nr_mb;i++) {    
		if (file_32bit_flag) {
			little_to_big_endian_SwapBytes(block_size32[i+1],sizeof(int32));
			mb_ptr[i].size = size = block_size32[i+1];
		} else if (file_64bit_flag) {
			little_to_big_endian_SwapBytes(block_size64[i+1],sizeof(int64));
			mb_ptr[i].size = size = block_size64[i+1];
		} 

		if (i<2) {
			memory_offset += size;
		} 
		Aimpack_printf("!>  Image Data starts at byte offset %" PRId64 " \n",memory_offset); 

		tot_size += size;
		if( size > 0 ) {
			/* check if memory could be allocated */
            MALLOC(mb_ptr[i].adr,size, char);
		} else { 
			mb_ptr[i].adr = 0;
		}
	}/* for */

	/* Read blocks */

	for(i=0;i<nr_mb;i++) {
		c_mb      = mb_ptr[i].adr;    
		c_mb_size = mb_ptr[i].size;     

		if( c_mb_size > 0 ){             /* ignore empty memory blocks */

			if( c_mb_size <= c_buf_size ){             /* all in buffer */
				aimpack_memcpy(c_mb,c_buf,c_mb_size); 
				c_buf_size -= c_mb_size;
				c_buf      += c_mb_size;
			} 
			else {  
				aimpack_memcpy(c_mb,c_buf,c_buf_size);                      /* head   */  
				c_mb_size  -= c_buf_size;
				c_mb       += c_buf_size;
				c_buf_size  = 0;

				/* read records direct to memblock if possible */ 

				while (c_mb_size >= rec_size) {                   /* mid   */
					if( PLATFORM_FILE_READ_C_MB ) {
						if (print) Aimpack_printf("!>  Can't read from file: %s.\n",filename);
						PLATFORM_FILE_CLOSE;
						*status=0; return;
					}
					c_mb_size -= rec_size;
					c_mb      += rec_size;
				}  
				if( c_mb_size > 0 ){                              /* tail   */
					if( PLATFORM_FILE_READ_BUF ) {
						if (print) Aimpack_printf("!>  Can't read from file: %s.\n",filename);
						PLATFORM_FILE_CLOSE;
						*status=0; return;
					}
					aimpack_memcpy(c_mb,buf,c_mb_size);  
					c_buf_size = rec_size - c_mb_size;
					c_buf      = buf      + c_mb_size;
				}
			}

		}
	} /* end for i */

	free(buf); 
	free(block_size64);
	PLATFORM_FILE_CLOSE;
	*status=1; 
}    

/******************************************************************************/
/* COMPRESS FUNCITONS   */
/******************************************************************************/




/*	  
**  AimpackCompress()
**  Run-length encoding of the image data.
*/	  

AIMPACK_EXPORT void AimpackCompress(D3AnyImage030 in, D3AnyImage *out, int out_type)
{
	int64		i,offs,out_offs;
	D3int64	  dim,off,out_dim,p,out_p;
	D3int64	  bit_pos;
	char	  *c_in,value;
	D1byte	  *b_out;

	D1charCmp	  *field_list;
	int64		nr_fields,field_offs;
	char	  cur_val;
	D1uchar	  cur_len;
	int64		mem_size;
	char	  *data;

	D1uchar	  *uc_field_list_len;
	char	  c_value_1,c_value_2;
	int32		c_value_2_found;
	int32		version020_flag = 0;

	D1charCmp2	  *field_list2;
	int32		undefined_dif;
	char	  cur_dif;
	char	  cur_first_val;
	char	  cur_last;

	if( in.type != D1Tchar ){
		Aimpack_printf("!>  AimpackCompress: Data not of type D1Tchar.\n");
		out->dat = 0;
		return;
	}
	Aimpack_printf("!>  AimpackCompress: start compression...\n");

	D3InitAny(*out,out_type);
	D3ImageGeometryCopy(in,*out);
	DnProcLogCopy(in,*out);	
	D3CopyAssocData(in,*out);

	/*  old compression only with 32bit size in header */
	if(in.version == 020) {
		/* printf("!> Trying compression AimVersion020..."); */
		version020_flag = 1;
		out->version = 020;
	}

	/*  Change 6.12.2000 A. Laib not to 'save' offset */
	/*  Change of Bug: 5.4.2001: old offset added to pos... A. Laib */
	D3SubCMul(out->dim,2,out->off,out->dim);
	D3Inc(out->pos, out->off);
	out->off = D3iNULL; 

	dim = in.dim;
	off = in.off;

	switch( out_type ) {

	case D3Tbit8:

		/* Take original offset and pos. 9.7.2001 A.Laib */
		D3ImageGeometryCopy(in,*out);

		D3Add(dim, D3iONE, out_dim);
		D3CDiv(out_dim, 2, out_dim);

		D3ReCallocAny(*out);

		c_in  = (D1char *) in.dat; 
		b_out = (D1byte *) out->dat; 

		for( p.z=off.z; p.z<dim.z-off.z; p.z++ ) {
			out_p.z   = p.z>>1;
			bit_pos.z = (p.z&1)<<2; 
			for( p.y=off.y; p.y<dim.y-off.y; p.y++ ) {

				out_p.y   =  p.y>>1;
				bit_pos.y = bit_pos.z + ((p.y&1)<<1); 
				offs      = (p.z*dim.y + p.y)*dim.x + off.x;
				out_offs  = (out_p.z*out_dim.y + out_p.y)*out_dim.x + (off.x>>1);

				for( p.x=off.x; p.x<dim.x-off.x; p.x++, offs++) {
					bit_pos.x = bit_pos.y + (p.x&1); 
					if( c_in[offs] != 0 ){
						b_out[out_offs] += 1<<bit_pos.x;
						value = c_in[offs];
					}
					out_offs += p.x&1;

				}	
			}
		}
		b_out[out_dim.x*out_dim.y*out_dim.z] = value;

		break;

	case D1TcharCmp:

		Free(out->dat);

		c_in  = (D1char *) in.dat; 

		cur_val = c_in[D3Offset(off,dim)];
		cur_len = 0;
		nr_fields = 1;
		for( p.z=off.z; p.z<dim.z-off.z; p.z++ ) {
			for( p.y=off.y; p.y<dim.y-off.y; p.y++ ) {

				offs      = (p.z*dim.y + p.y)*dim.x + off.x;

				for( p.x=off.x; p.x<dim.x-off.x; p.x++, offs++) {
					if( c_in[offs] == cur_val ){
						if( cur_len == 255 ){
							nr_fields++;
							cur_len = 1; 
						} else {
							cur_len++;
						}
					} else {
						cur_val = c_in[offs];
						cur_len  = 1;
						nr_fields++;
					}
				}	
			}
		}

		if (version020_flag) {
			if (nr_fields*2 + 4 < PLATFORM_2GB) {
				mem_size = nr_fields*2 + 4;
			} else {
				Aimpack_printf("!>  AIM Compression 020 not possible. Writing current compression...\n");
				version020_flag = 0;
				mem_size = nr_fields*2 + 8;
				out->version = 030;
			}
		} else {
			mem_size = nr_fields*2 + 8;
		}

		MALLOC(data, mem_size, char); 
		
		field_list = (D1charCmp *) data;

		if (version020_flag){
			*((int32 *) data) = mem_size;
			field_offs = 2;
		} else {
			*((int64 *) data) = mem_size;
			field_offs = 4; 
		}

		cur_val = c_in[D3Offset(off,dim)];
		cur_len = 0;

		for( p.z=off.z; p.z<dim.z-off.z; p.z++ ) {
			for( p.y=off.y; p.y<dim.y-off.y; p.y++ ) {
				offs = (p.z*dim.y + p.y)*dim.x + off.x;
				for( p.x=off.x; p.x<dim.x-off.x; p.x++, offs++) {
					if( c_in[offs] == cur_val ){
						if( cur_len == 255 ){
							/* write prev field */
							field_list[field_offs].val = cur_val;
							field_list[field_offs].len = cur_len;
							field_offs++;
							cur_len = 1;  /* keep old val */
						} else {
							cur_len++;
						}
					} else {
						/* write prev field */
						field_list[field_offs].val = cur_val;
						field_list[field_offs].len = cur_len;
						field_offs++;
						cur_val = c_in[offs];
						cur_len = 1;
					}
				}	
			}
		}

		/* write last field */
		field_list[field_offs].val = cur_val;
		field_list[field_offs].len = cur_len;

		out->dat = (void *) field_list;

		break;

	case D1TbinCmp:

		c_in  = (D1char *) in.dat; 

		/*  Check input  */

		c_value_1 = c_in[D3Offset(off,dim)];
		c_value_2_found = False;

		for( p.z=off.z; p.z<dim.z-off.z; p.z++ ) {
			for( p.y=off.y; p.y<dim.y-off.y; p.y++ ) {
				offs      = (p.z*dim.y + p.y)*dim.x + off.x;
				for( p.x=off.x; p.x<dim.x-off.x; p.x++, offs++) {
					if( c_in[offs] != c_value_1 ){
						if( !c_value_2_found ){
							c_value_2       =  c_in[offs];
							c_value_2_found = True;
						} else if( c_in[offs] != c_value_2 ){
							Aimpack_printf("!>  Warning compressing: Input is not binary. 3 or more values in volume.\n");
							Aimpack_printf("!>  D1TcharCmp Compression applied\n");
							AimpackCompress(in,out,D1TcharCmp);
							return;
						}
					}
				}	
			}
		}

		/*
		printf("c_value_1 = %d\n",c_value_1); 
		printf("c_value_2 = %d\n",c_value_2); 
		*/

		Free(out->dat);

		/*  Count fields  */

		cur_val = c_value_1;
		cur_len = 0;
		nr_fields = 1;
		for( p.z=off.z; p.z<dim.z-off.z; p.z++ ) {
			for( p.y=off.y; p.y<dim.y-off.y; p.y++ ) {

				offs      = (p.z*dim.y + p.y)*dim.x + off.x;

				for( p.x=off.x; p.x<dim.x-off.x; p.x++, offs++) {
					if( c_in[offs] == cur_val ){
						if( cur_len == 254 ){
							nr_fields++;
							cur_len = 1; 
						} else {
							cur_len++;
						}
					} else {
						cur_val = c_in[offs];
						cur_len  = 1;
						nr_fields++;
					}
				}	
			}
		}

		/* printf("nr_fields = %" PRId64 " \n",nr_fields); */

		if (version020_flag) {
			if ( nr_fields*1 + 4 + 2 < PLATFORM_2GB) {
				mem_size = nr_fields*1 + 4 + 2;
				/* printf("!>  approved. \n"); */
			} else {
				Aimpack_printf("!>  AIM Compression 020 not possible. Writing 030 compression...\n");
				version020_flag = 0;
				mem_size = nr_fields*1 + 8 + 2;
				out->version = 030;
			}
		} else { 
			mem_size = nr_fields*1 + 8 + 2;
		}

		MALLOC(data, mem_size, char);
		
		
		uc_field_list_len = (D1uchar *) data;
		if (version020_flag) {
			*((int32 *) data)     = mem_size;
			*((char *) (data+4)) = c_value_1;
			*((char *) (data+5)) = c_value_2;
			field_offs = 6;
		} else {
			*((int64 *) data)     = mem_size;
			*((char *) (data+8)) = c_value_1;
			*((char *) (data+9)) = c_value_2;
			field_offs = 10;
		}
		cur_val = c_value_1;
		cur_len = 0;

		for( p.z=off.z; p.z<dim.z-off.z; p.z++ ) {
			for( p.y=off.y; p.y<dim.y-off.y; p.y++ ) {

				offs = (p.z*dim.y + p.y)*dim.x + off.x;

				for( p.x=off.x; p.x<dim.x-off.x; p.x++, offs++) {
					if( c_in[offs] == cur_val ){
						if( cur_len == 254 ){
							/* write prev field */
							uc_field_list_len[field_offs] = 255;   /* means 254 and don't change value */
							field_offs++;
							cur_len = 1;  /* keep old val */
						} else {
							cur_len++;
						}
					} else {
						/* write prev field */
						uc_field_list_len[field_offs] = cur_len;
						field_offs++;
						cur_val = c_in[offs];
						cur_len = 1;
					}
				}	
			}
		}

		/* write last field */
		uc_field_list_len[field_offs] = cur_len;

		out->dat = (void *) uc_field_list_len;

		/*
		printf("field_offs = %" PRId64 " \n",field_offs); 
		printf("D3MemorySize(*out) = %" PRId64 " \n",D3MemorySize(*out)); 
		*/

		break;

#ifndef cray /* This case cannot be included on the cray (conflict with D1TcharCmp). */
	case D1TcharCmp2:

		Free(out->dat);

		c_in  = (D1char *) in.dat; 

		cur_first_val = cur_last = 0;
		cur_dif = 0;
		undefined_dif = True;
		cur_len = 0;
		nr_fields = 1;

		for( p.z=off.z; p.z<dim.z-off.z; p.z++ ) {
			for( p.y=off.y; p.y<dim.y-off.y; p.y++ ) {

				offs = (p.z*dim.y + p.y)*dim.x + off.x;

				for( p.x=off.x; p.x<dim.x-off.x; p.x++, offs++) {
					if( undefined_dif ){
						cur_dif = c_in[offs] - cur_last;
						cur_last = c_in[offs];
						cur_len++;
						undefined_dif = False;
					} else if( c_in[offs] - cur_last == cur_dif ){
						cur_last = c_in[offs];
						if( cur_len == 255 ){
							nr_fields++;
							cur_first_val = cur_last = c_in[offs];
							cur_len = 0;
							undefined_dif = True;
						} else {
							cur_len++;
						}
					} else {
						nr_fields++;
						cur_first_val = cur_last = c_in[offs];
						cur_len = 0;
						undefined_dif = True;
					}
				}	
			}
		}

		Aimpack_printf("!>  nr_fields = %" PRId64 " \n",nr_fields);
		if (version020_flag) {
			if (nr_fields*3 + 4 < PLATFORM_2GB) {
				mem_size = nr_fields*3 + 4;
			} else {
				Aimpack_printf("!>  AIM Compression 020 not possible. Writing current compression...\n");
				version020_flag = 0;
				mem_size = nr_fields*3 + 8;
				out->version = 030;
			}
		} else {
			mem_size = nr_fields*3 + 8;
		}

		MALLOC(data, mem_size, char);
		
		field_list2 = (D1charCmp2 *) data;
		if (version020_flag) {
			*((int32 *) data) = mem_size;
			field_offs = 2;
		} else {
			*((int64 *) data) = mem_size;
			field_offs = 4;
		}

		cur_first_val = cur_last = 0;
		cur_dif = 0;
		undefined_dif = True;
		cur_len = 0;
		nr_fields = 1;

		for( p.z=off.z; p.z<dim.z-off.z; p.z++ ) {
			for( p.y=off.y; p.y<dim.y-off.y; p.y++ ) {

				offs = (p.z*dim.y + p.y)*dim.x + off.x;

				for( p.x=off.x; p.x<dim.x-off.x; p.x++, offs++) {
					if( undefined_dif ){
						cur_dif = c_in[offs] - cur_last;
						cur_last = c_in[offs];
						cur_len++;
						undefined_dif = False;
					} else if( c_in[offs] - cur_last == cur_dif ){
						cur_last = c_in[offs];
						if( cur_len == 255 ){
							/* write prev field */
							field_list2[field_offs].val = cur_first_val;
							field_list2[field_offs].dif = cur_dif;
							field_list2[field_offs].len = cur_len;
							field_offs++;
							cur_first_val = cur_last = c_in[offs];
							cur_len = 0;
							undefined_dif = True;
						} else {
							cur_len++;
						}
					} else {
						/* write prev field */
						field_list2[field_offs].val = cur_first_val;
						field_list2[field_offs].dif = cur_dif;
						field_list2[field_offs].len = cur_len;
						field_offs++;
						cur_first_val = cur_last = c_in[offs];
						cur_len = 0;
						undefined_dif = True;
					}
				}	
			}
		}

		/* write last field */

		field_list2[field_offs].val = cur_first_val;
		field_list2[field_offs].dif = cur_dif;
		field_list2[field_offs].len = cur_len;

		out->dat = (void *) field_list2;

		Aimpack_printf("!>  field_offs = %" PRId64 " \n",field_offs); 
		Aimpack_printf("!>  D3MemorySize(*out) = %" PRId64 " \n",D3MemorySize(*out)); 

		for( i=0; i<10; i++ ) {   
			Aimpack_printf("!>  %" PRId64 "  val = %d len = %d\n",i,field_list2[i].val,field_list2[i].len); 
		}


		break;
#endif
	default:
		Aimpack_printf("!>  Aimpack Compress: Invalid Compression Type.\n");
		return;
	}
	Aimpack_printf("!>  AimpackCompress: Compression done!\n");

}

/*	  
**  AimpackUncompress()
**  Uncompresses the run-length encoded image data.
*/	  

AIMPACK_EXPORT void AimpackUncompress(D3AnyImage030 *in_out, char value)
{
	D3AnyImage030	  in;
	D3AnyImage030     *out;
	int64		offs;
	int64       in_offs;
	D3int64	    bit_pos;
	D3int64	    dim;
	D3int64     off;
	D3int64     in_dim;
	D3int64     p;
	D3int64     out_p;
	char	    *c_out;
	char        *out_dat;
	D1byte	    *b_in;

	D1charCmp	*field_list;
	int64		field_offs;
	char	    cur_val;
	D1uchar	    cur_len;

	D1uchar	    *uc_field_list_len;
	char	    c_value_1;
	char        c_value_2;
	int32		is_value_1; 
	int32       change_val;
	int32       version020_flag = 0;

	/*  old compression only with 32bit size in header */
	if(in_out->version == 020) {
		/*in_out->version = 030;*/  /*  update header */
		version020_flag = 1;
	}

	/*  Change 6.12.2000 A.Laib not to write offset */
	/*  Change of Bug: 5.4.2001: old offset added to pos... A. Laib */
	/*  Change of Bug: 9.7.2001: */
	/*  only done when uncompress is really performed. inside switch now  A. Laib */

	switch( in_out->type ) {

	case D3Tbit8:

		/* no offset compensation done. 9.7.2001 A.Laib */
		in  = *in_out;
		out =  in_out;

		dim = in.dim;
		off = in.off;

		/*  Change 19.12.2002 A.Laib */
		out->type = D1Tchar;
		Aimpack_printf("!>  Uncompressing   Format: D3bit8 ...\n");
		/*printf("!> Image size after D3bit8 uncompression will be: \n");*/
		/*D3PrintMemSizeInfoAny(*out);*/

		D3Add(dim, D3iONE, in_dim);
		D3CDiv(in_dim, 2, in_dim);

		CALLOC(out_dat, D3Volume(in), D1char);

		c_out = out_dat; 
		b_in  = (D1byte *) in.dat; 

		value = (value==0)? b_in[in_dim.x*in_dim.y*in_dim.z] : value;

		for( p.z=off.z; p.z<dim.z-off.z; p.z++ ) {
			out_p.z   = p.z>>1;
			bit_pos.z = (p.z&1)<<2;
			for( p.y=off.y; p.y<dim.y-off.y; p.y++ ) {
				out_p.y   =  p.y>>1;
				bit_pos.y = bit_pos.z + ((p.y&1)<<1); 
				offs      = (p.z*dim.y + p.y)*dim.x + off.x;
				in_offs   = (out_p.z*in_dim.y + out_p.y)*in_dim.x + (off.x>>1);
				for( p.x=off.x; p.x<dim.x-off.x; p.x++, offs++) {
					bit_pos.x = bit_pos.y + (p.x&1);
					c_out[offs] = (b_in[in_offs] & (1<<bit_pos.x)) != 0 ? value : 0;
					in_offs += p.x&1;
				}	

			}
		}

		free(in.dat); 

		out->type = D1Tchar;
		out->dat  = out_dat;

		Aimpack_printf("!>  Uncompressing done.\n");

		break;

	case D1TcharCmp:

		D3SubCMul(in_out->dim,2,in_out->off,in_out->dim);
		D3Inc(in_out->pos, in_out->off);
		in_out->off = D3iNULL;  

		in  = *in_out;
		out =  in_out;

		dim = in.dim;
		off = in.off;

		/*  Change 19.12.2002 A.Laib */
		out->type = D1Tchar;
		Aimpack_printf("!>  Uncompressing   Format: D1TcharCmp ...\n");
		/*printf("!> Image size after D1charCmp (run-length) uncompression will be: \n");
		D3PrintMemSizeInfoAny(*out);*/

		CALLOC(out_dat, D3Volume(in), D1char);

		c_out = out_dat; 

		if (version020_flag) {
			field_offs = 2;
		} else {
			field_offs = 4;
		}
		field_list = (D1charCmp *) in.dat;

		/*
		printf("D3MemorySize(in) = %" PRId64 " \n",D3MemorySize(in));  
		nr_fields = (D3MemorySize(in) - 4) /2;
		printf("field_offs = %" PRId64 " \n",field_offs); 
		*/

		cur_len = field_list[field_offs].len;
		cur_val =	field_list[field_offs].val;

		for( p.z=off.z; p.z<dim.z-off.z; p.z++ ) {
			for( p.y=off.y; p.y<dim.y-off.y; p.y++ ) {
				offs = (p.z*dim.y + p.y)*dim.x + off.x;
				for( p.x=off.x; p.x<dim.x-off.x; p.x++, offs++) {
					c_out[offs] = cur_val;
					cur_len--;
					if( cur_len == 0 ){
						field_offs++;
						cur_len = field_list[field_offs].len;
						cur_val =	field_list[field_offs].val;
					}
				}	
			}
		}

		free(in.dat);

		out->type = D1Tchar;
		out->dat  = out_dat;

		Aimpack_printf("!>  Uncompressing done.\n");

		break;

	case D1TbinCmp:

		D3SubCMul(in_out->dim,2,in_out->off,in_out->dim);
		D3Inc(in_out->pos, in_out->off);
		in_out->off = D3iNULL;  

		in  = *in_out;
		out =  in_out;

		dim = in.dim;
		off = in.off;

		/*  Change 19.12.2002 A.Laib */
		out->type = D1Tchar;
		Aimpack_printf("!>  Uncompressing   Format: D1TbinCmp ...\n");
		/*printf("!>  Image size after D1binCmp (run-length, binary) uncompression will be: \n");
		D3PrintMemSizeInfoAny(*out);*/

		/* if (version020_flag) printf("!>  AimVersion 020..."); */

		CALLOC(out_dat, D3Volume(in), D1char);	

		c_out = out_dat; 

		uc_field_list_len = (D1uchar *) in.dat;

		if (version020_flag){
			c_value_1 = *( ((char *) uc_field_list_len) + 4);
			c_value_2 = *( ((char *) uc_field_list_len) + 5);
			field_offs = 6;
		} else {
			c_value_1 = *( ((char *) uc_field_list_len) + 8);
			c_value_2 = *( ((char *) uc_field_list_len) + 9);
			field_offs = 10;
		}

		cur_len = uc_field_list_len[field_offs];

		if( cur_len == 255) {
			cur_len = 254;
			change_val = False;
		} else { 
			change_val = True;
		}

		cur_val    = c_value_1;
		is_value_1 = True;

		for( p.z=off.z; p.z<dim.z-off.z; p.z++ ) {
			for( p.y=off.y; p.y<dim.y-off.y; p.y++ ) {
				offs = (p.z*dim.y + p.y)*dim.x + off.x;
				for( p.x=off.x; p.x<dim.x-off.x; p.x++, offs++) {
					c_out[offs] = cur_val;
					cur_len--;
					if( cur_len == 0 ){
						if( change_val ){
							is_value_1 = !is_value_1;
							cur_val = (is_value_1)? c_value_1 : c_value_2;
						}
						field_offs++;
						cur_len = uc_field_list_len[field_offs];
						if( cur_len == 255) {
							cur_len = 254;
							change_val = False;
						} else {  
							change_val = True;
						}	   
					}
				}	
			}
		}

		free(in.dat);

		out->type = D1Tchar;
		out->dat  = out_dat;

		Aimpack_printf("!>  Uncompressing done.\n");

		break;

	case D1Tshort:
	case D1Tchar:
	case D1Tint:
	case D1Tfloat:
		Aimpack_printf("!>  Aimpack Uncompress: File not compressed.\n");
		break;
		
	default:
		Aimpack_printf("!>  Aimpack Uncompress: Invalid Compression Type.\n");
		return;

	}

}


/*-----------------------------------------------------------------------------+
| These are the functions written to make the Little Endian/ Big Endian cross
! platform compatibility.
+-----------------------------------------------------------------------------*/    

/*
** void AimpackRmSwapByteOrder()					      
** Swaps bytes in short, int or quad according to nbytes(=no of swaped bytes)
** Data needs to be smaller than 8bytes!)   
*/

AIMPACK_EXPORT void AimpackRmSwapByteOrder(char *data, int nbytes) 
{
	int	i;
	char	tmp[8];
	
	if (nbytes > 8)
	{
		Aimpack_printf("!>  Swap Byte: data to swap exceedes 8 bytes. \n ");  
	}
	aimpack_memcpy(tmp,data,nbytes);
	
		for (i=0;i<nbytes;i++) data[i] = tmp[nbytes-1-i];
	
}

/*
** void AimpackRmConvertFloatVMS2Sun(data)
** Converts a float with VMS representation to a Sun float		 
** representation if forward is true and vice versa if false.
*/

AIMPACK_EXPORT void AimpackRmConvertFloatVMS2Sun(float *data, int forward)
{
	if (!forward) *data *= 4.0;

	AimpackRmSwapByteOrder((char *)data  ,sizeof(*data)/2);
	AimpackRmSwapByteOrder((char *)data+2,sizeof(*data)/2);

	if (forward) *data /= 4.0;
}

/* void AimpackSkbConvertFloatVMS2Win32(data)
** Converts a float with VMS representation to a Windows float
** representation if forward is true and vice versa if false.
*/

AIMPACK_EXPORT void AimpackSkbConvertFloatVMS2Win32(float *data, int forward)
{
	if (!forward) *data *= 4.0;

	AimpackRmSwapByteOrder((char *)data  ,sizeof(*data)); /* first swap all four bytes */
	AimpackRmSwapByteOrder((char *)data  ,sizeof(*data)/2);
	AimpackRmSwapByteOrder((char *)data+2,sizeof(*data)/2);

	if (forward) *data /= 4.0;
}


/******************************************************************************/
/* CONVERT AIM DATA IEEE <-> FLOAT     */
/******************************************************************************/


/*	 
**  AimpackConvertAimFloat2I3EFloat();
*/	 

AIMPACK_EXPORT void AimpackConvertAimFloat2I3EFloat(D3AnyImage030 *image)
{

  int64 i;
  int64 totdim;
  float *dat;
  
  if(image->type != D1Tfloat) {
    Aimpack_printf("!>  AimpackConvertAimFloat2I3EFloat: input type not VAX float. Return.\n");
    return;
  }

  dat = (float *) image->dat; 

  /*  Changed 22.5.2002 to make use of cvt$ftof, see help rtl cvt$ cvt$f */

  totdim = image->dim.x * image->dim.y *image->dim.z;
     
  for( i=0; i<totdim; i++, dat++ ) {     
    PLATFORM_FTOF_VAX(dat);
  }
 
  image->type = D1Ti3efloat;

}


/*	 
**  AimpackConvertAimI3EFloat2Float();
*/	 

AIMPACK_EXPORT void AimpackConvertAimI3EFloat2Float(D3AnyImage030 *image)
{

  int64 i;
  int64 totdim;
  float *dat;
  
  if(image->type != D1Ti3efloat) {
    Aimpack_printf("!>  AimpackConvertAimI3EFloat2Float: input type not IEEE float. Return.\n");
    return;
  }

  dat = (float *) image->dat; 

  /*  Changed 22.5.2002 to make use of cvt$ftof, see help rtl cvt$ cvt$f */

  totdim = image->dim.x * image->dim.y *image->dim.z;
  for( i=0; i<totdim; i++, dat++) {  
    PLATFORM_FTOF_IEEE(dat);
 }


  image->type = D1Tfloat;
}
