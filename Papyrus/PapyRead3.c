/********************************************************************************/
/*                                                                              */
/* Papyrus 3 library.                                                           */
/* This library constitutes a DICOM file system which helps reading and writing */
/* DICOM files and DICOMDIR files.                                              */
/*                                                                              */
/* Copyright (C) 2004 - Service of  Medical Informatics -                       */
/* University Hospitals of Geneva (HUG), Geneva, Switzerland                    */
/*                                                                              */
/* This library is a free software; you can redistribute it and/or modify it    */
/* under the terms of the GNU Lesser General Public License as published by the */
/*  Free Software Foundation; either version 2.1 of the License, or             */
/* (at your option) any later version.                                          */
/*                                                                              */
/* This library is distributed in the hope that it will be useful,              */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of               */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                         */
/* See the GNU Lesser General Public License for more details.                  */
/*                                                                              */
/* You should have received a copy of the GNU Lesser General Public License     */
/* along with this library; if not, write to                                    */
/* the Free Software Foundation, Inc.,                                          */
/* 59 Temple Place, Suite 330,                                                  */
/* Boston, MA 02111-1307 USA                                                    */
/*                                                                              */
/* You can contact us for more information at osiris@sim.hcuge.ch               */
/* or by writing to Papyrus,                                                    */
/* Unite d'Imagerie Numerique / Service d'Informatique Medicale / HUG,          */
/* 24, Micheli-du-Crest street, 1211 Geneva 14, Switzerland.                    */
/*                                                                              */
/* The University Hopitals of Geneva, hereby disclaims all copyright interest   */
/* in the library `Papyrus' (a library for reading and writing DICOM files).    */
/*                                                                              */
/* Geneva, april 2004                                                           */
/* Antoine Geissbuhler, head of the Service of Medical Informatics,             */
/* University Hospitals of Geneva, Switzerland                                  */
/*                                                                              */
/********************************************************************************/

/********************************************************************************/
/*                                                                              */
/*	Project  : P A P Y R U S  Toolkit                                           */
/*	File     : PapyRead3.c                                                      */
/*	Function : contains all the reading functions                               */
/*	Authors  : Matthieu Funk                                                    */
/*             Christian Girard                                                 */
/*             Jean-Francois Vurlod                                             */
/*             Marianne Logean                                                  */
/*                                                                              */
/*	History  : 12.1990	version 1.0                                             */
/*             04.1991	version 1.1                                             */
/*             12.1991	version 1.2                                             */
/*             06.1993	version 2.0                                             */
/*             06.1994	version 3.0                                             */
/*             06.1995	version 3.1                                             */
/*             02.1996	version 3.3                                             */
/*             02.1999	version 3.6                                             */
/*             04.2001	version 3.7                                             */
/*             09.2001  version 3.7  on CVS                                     */
/*             10.2001  version 3.71 MAJ Dicom par CHG                          */
/*             03.2010  Fuli Wu                                                 */
/********************************************************************************/

#ifdef Mac
#pragma segment papy3
#endif

/* ------------------------- includes ---------------------------------------*/

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "setjmp.h"
//#include "jpegless.h"       /* interface for JPEG lossless decompressor */
//#include "jpeglib.h"	    /* interface for JPEG lossy decompressor */
#include "openjpeg.h"
#include "jasper.h"


#ifdef MAYO_WAVE
#include "Mayo.h"	/* interface for wavelet decompressor */
#define TO_SWAP_MAYO
#endif /* MAYO_WAVE */

#ifndef Papyrus3H 
#include "Papyrus3.h"
#endif
	  
extern PapyShort ExtractJPEGlossy8 (PapyShort inFileNb, PapyUChar *ioImage8P, PapyULong inPixelStart, PapyULong *inOffsetTableP, int inImageNb, int inDepth, int mode);
extern PapyShort ExtractJPEGlossy12 (PapyShort inFileNb, PapyUChar *ioImage8P, PapyULong inPixelStart, PapyULong *inOffsetTableP, int inImageNb, int inDepth, int mode);
extern PapyShort ExtractJPEGlossy16 (PapyShort inFileNb, PapyUChar *ioImage8P, PapyULong inPixelStart, PapyULong *inOffsetTableP, int inImageNb, int inDepth, int mode);
extern void PapyrusLockFunction( int lock);
extern short UseOpenJpeg;

static char **globalElementPtrs = 0L;
static int JasperInitialized = 0;


/********************************************************************************/
/*										*/
/*	Papy3GetElement : gets the value(s) of the specified element		*/
/* 	return : the value(s) of the element					*/
/*										*/
/********************************************************************************/

UValue_T * CALLINGCONV
Papy3GetElement (SElement *inGrOrModP, int inElement, PapyULong *outNbValueP, int *outElemTypeP)

/*SElement 	*inGrOrModP;		     ptr on the group or the module */
/*PapyShort	inElement;	   the position of the element in the group */
/*PapyULong *outNbValueP;		       the number of values to read */
/*PapyShort	*outElemTypeP;		    what is the type of the element */
{
  SElement *theElemP;	       /* work pointer on the elements of the group */
  UValue_T *theValueP;			    /* the value we are looking for */
  
  
  if (inGrOrModP == NULL) return NULL;
  
  theElemP  = inGrOrModP;
  theElemP += inElement;			   /* points on the desired element */
  
  *outElemTypeP = theElemP->vr;	 /* is it a short a long or an ASCII char ? */
  
  if (theElemP->nb_val > 0L)		    /* there is an introduced value */
  {
    *outNbValueP = theElemP->nb_val;
    theValueP    = theElemP->value;
  } /* then */
  
  else
  {
     *outNbValueP = 0L;
     theValueP    = NULL;
  } /* else */
  
  return theValueP;
  
} /* endof Papy3GetElement */

/*
sample error callback expecting a FILE* client object
*/
static void error_callback(const char *msg, void *a)
{
//	printf( "%s\r\r", msg);
}
/**
sample warning callback expecting a FILE* client object
*/
static void warning_callback(const char *msg, void *a)
{
//	printf( "%s\r\r", msg);
}

/**
sample debug callback expecting no client object
*/
static void info_callback(const char *msg, void *a)
{
//	printf( "%s\r\r", msg);
}

static inline int int_ceildivpow2(int a, int b) {
	return (a + (1 << b) - 1) >> b;
}


bool read_JPEG2000_file (void* raw, char *inputdata, size_t inputlength)
{
	opj_dparameters_t parameters;  /* decompression parameters */
	opj_event_mgr_t event_mgr;      /* event manager */
	opj_image_t *image = nil;
	opj_dinfo_t* dinfo;             /* handle to a decompressor */
	opj_cio_t *cio;
	unsigned char *src = (unsigned char*)inputdata; 
	int file_length = inputlength;

	/* configure the event callbacks (not required) */
	memset(&event_mgr, 0, sizeof(opj_event_mgr_t));
	event_mgr.error_handler = error_callback;
	event_mgr.warning_handler = warning_callback;
	event_mgr.info_handler = info_callback;

  /* set decoding parameters to default values */
  opj_set_default_decoder_parameters(&parameters);
 
  // default blindly copied
  parameters.cp_layer=0;
  parameters.cp_reduce=0;

  /* JPEG-2000 codestream */
  parameters.decod_format = 0;
  parameters.cod_format = 1;

  /* get a decoder handle */
  dinfo = opj_create_decompress(CODEC_J2K);

  /* catch events using our callbacks and give a local context */
  opj_set_event_mgr((opj_common_ptr)dinfo, &event_mgr, NULL);

  /* setup the decoder decoding parameters using user parameters */
  opj_setup_decoder(dinfo, &parameters);
  
  /* open a byte stream */
  cio = opj_cio_open((opj_common_ptr)dinfo, src, file_length);
  
  /* decode the stream and fill the image structure */
  image = opj_decode(dinfo, cio);
  if(!image)
  {
    opj_destroy_decompress(dinfo);
    opj_cio_close(cio);
    return false;
  }
      
  /* close the byte stream */
  opj_cio_close(cio);

  /* free the memory containing the code-stream */

   // Copy buffer
   for (int compno = 0; compno < image->numcomps; compno++)
   {
      opj_image_comp_t *comp = &image->comps[compno];

      int w = image->comps[compno].w;
      int wr = int_ceildivpow2(image->comps[compno].w, image->comps[compno].factor);
	    int numcomps = image->numcomps;
	   
      int hr = int_ceildivpow2(image->comps[compno].h, image->comps[compno].factor);
	   
	   if( wr == w && numcomps == 1)
	   {
		   if (comp->prec <= 8)
		   {
			   uint8_t *data8 = (uint8_t*)raw + compno;
			   int *data = image->comps[compno].data;
			   int i = wr * hr;
			   while( i -- > 0)
				   *data8++ = (uint8_t) *data++;
		   }
		   else if (comp->prec <= 16)
		   {
			   uint16_t *data16 = (uint16_t*)raw + compno;
			   int *data = image->comps[compno].data;
			   int i = wr * hr;
				while( i -- > 0)
					*data16++ = (uint16_t) *data++;
		   }
		   else
		   {
			  //printf( "****** 32-bit jpeg encoded is NOT supported\r");
		   }
	   }
	   else
	   {
			if (comp->prec <= 8)
			{
			 uint8_t *data8 = (uint8_t*)raw + compno;
			 for (int i = 0; i < wr * hr; i++)
			 {
				*data8 = (uint8_t) (image->comps[compno].data[i / wr * w + i % wr]);
				data8 += numcomps;
			 }
			}
			else if (comp->prec <= 16)
			{
			 uint16_t *data16 = (uint16_t*)raw + compno;
			 for (int i = 0; i < wr * hr; i++)
			 {
				*data16 = (uint16_t) (image->comps[compno].data[i / wr * w + i % wr]);
				data16 += numcomps;
			 }
			}
			else
			{
				//printf( "****** 32-bit jpeg encoded is NOT supported\r");
			}
	   }
   }


  /* free remaining structures */
  if(dinfo) {
    opj_destroy_decompress(dinfo);
  }

  /* free image data structure */
  if( image)
	opj_image_destroy(image);

  return true;
}

PapyShort ExtractJPEG2000 (PapyShort inFileNb, PapyUChar *ioImage8P, PapyULong inPixelStart, PapyULong *inOffsetTableP, int inImageNb, int inDepth, long offsetSize)
{
	int				fmtid;
	PapyUChar		theTmpBuf [256];
	PapyUChar		*theTmpBufP;
	PapyULong		i, thePos, theLength, theULong, x;
	PapyShort		theErr;
	PapyUShort		theUShort1, theUShort2;
	PapyUChar		*theCompressedP;
	long			ok = FALSE;
	
	Papy3FSeek (gPapyFile [inFileNb], SEEK_SET, (PapyLong) (inPixelStart + inOffsetTableP [inImageNb - 1]));
	
	theLength = 0;
	ok = FALSE;
	while (!ok)
	{
		/* read 8 chars from the file */
		i 	      = 8L;
		thePos      = 0L;
		theTmpBufP  = (unsigned char *) &theTmpBuf [0];
		if ((theErr = (PapyShort) Papy3FRead (gPapyFile [inFileNb], &i, 1L, theTmpBufP)) < 0)
		{
			return -1;
		} /* if */

		thePos = 0L;
		theUShort1 = Extract2Bytes (inFileNb, theTmpBufP, &thePos);
		theUShort2 = Extract2Bytes (inFileNb, theTmpBufP, &thePos);
		theULong = Extract4Bytes (inFileNb, theTmpBufP, &thePos);
		theLength += theULong;

		/* offset table found ? */
		if ((theUShort1 == 0xFFFE) && (theUShort2 == 0xE000))
		{
			Papy3FSeek (gPapyFile [inFileNb], SEEK_CUR, theULong);
		} /* if */
		else if ((theUShort1 == 0xFFFE) && (theUShort2 == 0xE0DD)) ok = TRUE;

	} /* while */
	
	long allocatedLength = theLength;
	
	theCompressedP = malloc( theLength);

	Papy3FSeek (gPapyFile [inFileNb], SEEK_SET, (PapyLong) (inPixelStart + inOffsetTableP [inImageNb - 1]));
	
	theLength = 0;
	ok = FALSE;
	while (!ok)
	{
		/* read 8 chars from the file */
		i 	      = 8L;
		thePos      = 0L;
		theTmpBufP  = (unsigned char *) &theTmpBuf [0];
		if ((theErr = (PapyShort) Papy3FRead (gPapyFile [inFileNb], &i, 1L, theTmpBufP)) < 0)
		{
			return -1;
		} /* if */

		thePos = 0L;
		theUShort1 = Extract2Bytes (inFileNb, theTmpBufP, &thePos);
		theUShort2 = Extract2Bytes (inFileNb, theTmpBufP, &thePos);
		theULong = Extract4Bytes (inFileNb, theTmpBufP, &thePos);
		
		/* offset table found ? */
		if ((theUShort1 == 0xFFFE) && (theUShort2 == 0xE000))
		{
			if ((theErr = (PapyShort) Papy3FRead (gPapyFile [inFileNb], &theULong, 1L, theCompressedP + theLength)) < 0)
			{
				Papy3FClose (&gPapyFile [inFileNb]);
				free( theCompressedP) ;
				RETURN (theErr);
			} /* if */
			
			theLength += theULong;
			
			//if( theLength > allocatedLength)	printf( "****** theLength > allocatedLength\r\r");
			
			//Papy3FSeek (gPapyFile [inFileNb], SEEK_CUR, theULong);
		} /* if */
		else if ((theUShort1 == 0xFFFE) && (theUShort2 == 0xE0DD)) ok = TRUE;
		
	} /* while */
	
	int succeed = 0;
		
	if( UseOpenJpeg == 1)
	{
		PapyrusLockFunction( 0);
		succeed = read_JPEG2000_file( ioImage8P, (char*) theCompressedP, theLength);
		//if( succeed == 0)	printf( "**** OpenJPEG 2000 failed to open this file. Will try Jasper2000.\r");
		PapyrusLockFunction( 1);
	}
	
  if( succeed == 0)
	{
		jas_image_t *jasImage;
		jas_matrix_t *pixels[4];
		char *fmtname;
		
		if( JasperInitialized == 0)
		{
			JasperInitialized = 1;
			jas_init();
		}
		
		PapyrusLockFunction( 0);
		
		jas_stream_t *jasStream = jas_stream_memopen((char *)theCompressedP, theLength);
		
		if ((fmtid = jas_image_getfmt(jasStream)) < 0)
		{
			PapyrusLockFunction( 1);
			RETURN( -32);
		}
			// Decode the image. 
		if (!(jasImage = jas_image_decode(jasStream, fmtid, 0)))
		{
			PapyrusLockFunction( 1);
			RETURN( -35);
		}
		
		// Close the image file. 
		jas_stream_close(jasStream);
		int numcmpts = jas_image_numcmpts(jasImage);
		int width = jas_image_cmptwidth(jasImage, 0);
		int height = jas_image_cmptheight(jasImage, 0);
		int depth = jas_image_cmptprec(jasImage, 0);
		fmtname = jas_image_fmttostr(fmtid);
		
		int bitDepth = 0;
		if (depth == 8)
			bitDepth = 1;
		else if (depth <= 16)
			bitDepth = 2;
		else if (depth > 16)
			bitDepth = 4;
		
		unsigned char *newPixelData = ioImage8P;
		
		if( gArrPhotoInterpret [inFileNb] == MONOCHROME1 || gArrPhotoInterpret [inFileNb] == MONOCHROME2) numcmpts = 1;
		
		for (i=0; i < numcmpts; i++)
			pixels[ i] = jas_matrix_create( height, width);
		
		if( numcmpts == 1)
		{
			if (depth > 8)
			{
				jas_image_readcmpt(jasImage, 0, 0, 0, width, height, pixels[0]);
				
				unsigned short *px = (unsigned short*) newPixelData;
				
				int_fast32_t	*ptr = &(pixels[0])->rows_[0][0];
				x = width * height;
				while( x-- > 0) *px++ = *ptr++;
			}
			else
			{
				jas_image_readcmpt(jasImage, 0, 0, 0, width, height, pixels[0]);
				
				char *px = (char*) newPixelData;
				
				int_fast32_t	*ptr = &(pixels[0])->rows_[0][0];
				x = width * height;
				while( x-- > 0) *px++ =	*ptr++;
			}
		}
		else
		{
			for( i = 0 ; i < numcmpts; i++)
				jas_image_readcmpt(jasImage, i, 0, 0, width, height, pixels[ i]);
			
			char *px = (char*) newPixelData;
			
			int_fast32_t	*ptr1 = &(pixels[0])->rows_[0][0];
			int_fast32_t	*ptr2 = &(pixels[1])->rows_[0][0];
			int_fast32_t	*ptr3 = &(pixels[2])->rows_[0][0];
			
			x = width * height;
			while( x-- > 0)
			{
				*px++ =	*ptr1++;
				*px++ =	*ptr2++;
				*px++ =	*ptr3++;
			}
		}
		
		for (i=0; i < numcmpts; i++)
			jas_matrix_destroy( pixels[ i]);
		
		jas_image_destroy(jasImage);
		
		PapyrusLockFunction( 1);		
	}
	
	free( theCompressedP);
	return 0;
}

/********************************************************************************/
/*									 	*/
/*	ExtractWavelet : gets and decode Wavelet pixel data			*/
/* 	return : the image							*/
/*										*/
/********************************************************************************/

#ifdef MAYO_WAVE
PapyShort
ExtractWavelet (PapyShort inFileNb, PapyUChar *ioImage8P, PapyULong inPixelStart,
		PapyULong *inOffsetTableP, int inImageNb, int inDepth)
{
		 
  PapyUShort		theGroup, theElement;
  MayoCompressedImage	*theCompressedP; 
  MayoRawImage		*theRawP ;  
  int			theJs, theIs;
  PapyUChar	        theTmpBuf [256];
  PapyUChar		*theTmpBufP;
  PapyUChar		*tmpBufPtr2;
  PapyULong		i, thePos, theSize, theLength;
  PapyShort		theErr;
  PapyUShort		*theImage16P;
  PapyUChar		*theValTempP, *theValFinalP;
  PapyUChar		theHigh, theLow;


  Papy3FSeek (gPapyFile [inFileNb], SEEK_SET, (PapyLong) (inPixelStart + inOffsetTableP [inImageNb - 1]));
  

  theTmpBufP = (PapyUChar *) &theTmpBuf [0];
  i = 8L; 					
  if ((theErr = (PapyShort) Papy3FRead (gPapyFile [inFileNb], &i, 1L, theTmpBufP)) < 0)
  {
    Papy3FClose (&gPapyFile [inFileNb]);
    RETURN (theErr);
  } 
    
  thePos     = 0L;
  theGroup   = Extract2Bytes (theTmpBufP, &thePos);
  theElement = Extract2Bytes (theTmpBufP, &thePos);
    
  /* Pixel data fragment not found when expected */
  if ((theGroup != 0xFFFE) || (theElement != 0xE000)) RETURN (papBadArgument);

/* Load the compressed file into memory */
  /*
  theCompressedP = MayoReadCompressed (gPapyFile [inFileNb]) ; 
  if ( theCompressedP == NULL ) { 
    exit(MayoGetError()) ; } 
*/
  tmpBufPtr2 = (PapyUChar *) &theTmpBuf [0];
  theCompressedP = (MayoCompressedImage *) emalloc3(sizeof(MayoCompressedImage)) ; 
  if ( theCompressedP == NULL ) 
  { 
    return (-1);
  } /* if */

  i = 8L;
  if ((theErr = (PapyShort) Papy3FRead (gPapyFile [inFileNb], &i, 1L,theTmpBufP)) < 0)
  {
    Papy3FClose (&gPapyFile [inFileNb]);
    efree3((void **) &theCompressedP) ; 
    RETURN (theErr);
  } /* if */
    
  thePos 		  = 0L;
  theCompressedP->length  = Extract4Bytes (inFileNb, theTmpBufP, &thePos);
  theCompressedP->version = Extract4Bytes (inFileNb, theTmpBufP, &thePos);
    

  /* Allocate memory for the image data */ 
  theCompressedP->buf = (unsigned char *) emalloc3 (theCompressedP->length) ; 
  if (theCompressedP->buf == NULL) 
  {   
    efree3 ((void **) &theCompressedP); 
    return(-1); 
  } /* if */
    
  /* Read the image data */ 
  theLength = (PapyULong) theCompressedP->length;
  if ((theErr = (PapyShort) Papy3FRead (gPapyFile [inFileNb], &theLength, 1L, theCompressedP->buf)) < 0)
  {
    Papy3FClose (&gPapyFile [inFileNb]);
    efree3((void **) &theCompressedP->buf) ; 
    efree3((void **) &theCompressedP) ;
    RETURN (theErr);
  } /* if */
  theCompressedP->length = (int) theLength;

  /* Run the decompressor */ 
  theRawP = MayoDecompress (theCompressedP); 
  if (theRawP == NULL) 
  {   
    return (MayoGetError()); 
  } /* if */

  /* Copy decompressed image */
  theSize = (PapyULong) theRawP->xsize * theRawP->ysize * theRawP->bytesperpixel;
  if (inDepth == 8) memcpy (ioImage8P, theRawP->buf, theSize);
  else
  {
    theImage16P = (PapyUShort *) ioImage8P;
    memcpy (ioImage8P, theRawP->buf, theSize);
  } /* else */

  /* Swap bytes if it is a 16-bit image*/
#ifdef TO_SWAP_MAYO

  theValTempP  = (PapyUChar *) ioImage8P;
  theValFinalP = (PapyUChar *) ioImage8P;

  if (inDepth == 16)
  {
    for (theJs = 0; theJs < theRawP->xsize; theJs++) 
    {
      for (theIs = 0; theIs < theRawP->ysize; theIs++) 
      {
	theLow  	= *theValTempP;
	theValTempP++;
	theHigh		= *theValTempP;
	theValTempP++;
	*theValFinalP 	= theHigh;
	theValFinalP++;
	*theValFinalP 	= theLow;
	theValFinalP++;
      } /* for */
    } /* for */
  } /* if ...inDepth = 16 */
#endif /* TO_SWAP_MAYO */

/* Free allocated memory  
  MayoFreeCompressed(theCompressedP) ; 
  MayoFreeRaw(theRawP) ; */

  if (theRawP != NULL ) 
  { 
    if (theRawP->buf != NULL ) 
    { 
      efree3 ((void **) &(theRawP->buf));
    } /* if */
    efree3 ((void **) &theRawP);
  } /* if */ 
  
  if (theCompressedP != NULL) 
  { 
    if (theCompressedP->buf != NULL)
    {
      efree3 ((void **) &(theCompressedP->buf));
    } /* if */ 
    efree3 ((void **) &theCompressedP);
  } /* if */ 

  return (0);

} /* endof ExtractWavelet */
#endif /* MAYO_WAVE */



/********************************************************************************/
/*									 	*/
/*     				DecodeRLESegment				*/
/*									 	*/
/********************************************************************************/

void
DecodeRLESegment (PapyShort inFileNb, PapyUShort *ioImageP, PapyUChar *inRleP, 
                  PapyULong inLength, int inSegtot, int inSegNb)
/* decode a RLE segment                                                         */
/* ioImageP  	: pointer on real image (8 or 16 bits)                          */
/* inRleP    	: pointer on rle buffer (8bits)                                 */
/* inLength 	: length of rle buffer                                          */
/* inSegtot 	: total number of segments (1, 2 or 3)	                        */
/* inSegNb  	: number of current segment (1, 2 or 3) (only if inSegtot > 2)	*/
{
  PapyLong		j, theIndj;
  PapyUChar		*thePixP;
  PapyUChar		theVal;
  char                  theCode;
  PapyShort		i, theIMax;


  /* *** single segment *** */
  /* ********************** */
  
  if (inSegtot == 1) 
  { 
    /* convert rle into real image */
    thePixP = (PapyUChar *) ioImageP;
    theIndj = 0L;
    for (j = 0L; j < (int) inLength;) 
    {
      theCode = (char) inRleP [j];
      j++; /* yes, I know but do not move it */
      
      /* sequence of different bytes */
      if (theCode == 0) 
      {
        if (j < (int) (inLength - 1)) thePixP [theIndj++] = inRleP [j++];
      } /* if */
      
      /* repetition of the same byte */
      else if ((theCode <= -1) && (theCode >= -127)) 
      {
        theVal = inRleP [j++];
        theIMax = -theCode;
        for (i = 0; i <= theIMax; i++) 
          thePixP [theIndj++] = theVal;
      } /* if */
      
      else /* if ((theCode > 0) && (theCode <= 127)) */
      {
        for (i = 0; i < (theCode + 1); i++) 
          thePixP [theIndj++] = inRleP [j++];
      } /* if */
	      
	  
    } /* for */
  } /* if ...single segment */
  
  /* *** two segments *** */
  /* ******************** */
  
  else if (inSegtot == 2) 
  {
    /* we assume it is a 16 bit image	*/
    /* convert rle into real image	*/
    thePixP = (PapyUChar *) ioImageP;
    theIndj = 0L;
    if (inSegNb == 2) theIndj++;
    for (j = 0L; j < (int)inLength; ) 
    {
      theCode = (char) inRleP [j];
      j++; /* yes, I know but do not move it */
      /* sequence of different bytes */
      if (theCode == 0) 
      {
        if (j < (int) (inLength - 1)) thePixP [theIndj] = inRleP [j++];
        theIndj = theIndj + 2;
      } /* if */
	  
      /* repetition of the same byte */
      else if ((theCode <= -1) && (theCode >= -127)) 
      {
        theVal  = inRleP [j++];
        theIMax = -theCode;
        for (i = 0; i <= theIMax; i++) 
        {
          thePixP [theIndj] = theVal;
          theIndj = theIndj + 2;
        } /* for */
      } /* if */
      
      else  /* if ((theCode > 0) && (theCode <= 127)) */
      {
        for (i = 0; i < (theCode + 1); i++) 
        {
          thePixP [theIndj] = inRleP [j++];
          theIndj = theIndj + 2;
        } /* for */
      } /* if */
	   
      
    } /* for */
  } /* if ...two segments */
  
  /* *** three segments *** */
  /* ******************** */
  
  else if (inSegtot == 3) 
  {
    /* this must be a RGB or YBR image */
    /* so convert each channel at a time */
    thePixP = (PapyUChar *) ioImageP;
    
    /* computes the offset in the resulting pixmap */
    /* assuming that each plane is 8 bits depth    */
    theIndj = 0L;
    theIndj += ((PapyLong) gx0028Rows [inFileNb] * (PapyLong) gx0028Columns [inFileNb]) 
    	        * (PapyLong) (inSegNb - 1);
    for (j = 0L; j < (int)inLength; ) 
    {
      theCode = (char) inRleP [j];
      j++; /* yes, I know but do not move it */
      /* sequence of different bytes */
      if (theCode == 0) 
      {
        if (j < (int)(inLength - 1)) thePixP [theIndj++] = inRleP [j++];
      }/* if */
	  
      /* repetition of the same byte */
      else if ((theCode <= -1) && (theCode >= -127)) 
      {
        theVal = inRleP [j++];
        theIMax = -theCode;
        for (i = 0; i <= theIMax; i++) thePixP [theIndj++] = theVal;
      } /* if */
      
      else /* if ((theCode > 0) && (theCode <= 127)) */
      {
        for (i = 0; i < (theCode + 1); i++) thePixP [theIndj++] = inRleP [j++];
      } /* if */
	  
      
	  
    } /* for */    
  } /* if ...three segments */
  
} /* endof DecodeRLESegment */



/********************************************************************************/
/*									 	*/
/*	ExtractRLE : gets and decode a RLE pixel data element			*/
/* 	return : the image							*/
/*										*/
/********************************************************************************/

PapyShort
ExtractRLE (PapyShort inFileNb, PapyUShort *ioImage16P, PapyULong inPixelStart,
	    PapyULong *inOffsetTableP, int inImageNb)
{
  PapyUChar	theTmpBuf [256];
  PapyUChar	*theTmpBufP;
  PapyUShort	theGroup, theElement;
  PapyShort	theErr;
  PapyULong 	theNbOfSegments, i, thePos, theLength;
  PapyUChar	*theRleP;
  long		theOffset1, theOffset2, theOffset3, theRleLen;
    
    
  /* for each image						*/
  /* FFFE E000 length RLE_header RLE_segment1 RLE_segment2 ...	*/
  /* length is 4 bytes, in the case of a single image		*/
  Papy3FSeek (gPapyFile [inFileNb], SEEK_SET, (PapyLong) (inPixelStart + inOffsetTableP [inImageNb - 1]));
  
  /* read 8 chars from the file */
  theTmpBufP = (PapyUChar *) &theTmpBuf [0];
  i = 8L; 					/* grNb, elemNb & elemLength */
  if ((theErr = (PapyShort) Papy3FRead (gPapyFile [inFileNb], &i, 1L, theTmpBufP)) < 0)
  {
    Papy3FClose (&gPapyFile [inFileNb]);
    RETURN (theErr);
  } /* if */
    
  thePos     = 0L;
  theGroup   = Extract2Bytes (inFileNb, theTmpBufP, &thePos);
  theElement = Extract2Bytes (inFileNb, theTmpBufP, &thePos);
  theLength  = Extract4Bytes (inFileNb, theTmpBufP, &thePos);
    
  /* Pixel data fragment not found when expected */
  if ((theGroup != 0xFFFE) || (theElement != 0xE000)) RETURN (papBadArgument);
  
  /* read 4 chars from the file = number of segments */
  theTmpBufP = (PapyUChar *) &theTmpBuf [0];
  i = 4L;
  thePos = 0L;					/* grNb, elemNb & elemLength */
  if ((theErr = (PapyShort) Papy3FRead (gPapyFile [inFileNb], &i, 1L, theTmpBufP)) < 0)
  {
    Papy3FClose (&gPapyFile [inFileNb]);
    RETURN (theErr);
  } /* if */
  theNbOfSegments = Extract4Bytes (inFileNb, theTmpBufP, &thePos);
  if (theNbOfSegments > 3L) RETURN (papWrongValue); /* we allow to read 8, 16 and 32 bit images */
    
  /* read theOffset1, theOffset2, theOffset3 and skip 48 bytes */
  theTmpBufP = (PapyUChar *) &theTmpBuf [0];
  i          = 12L;
  thePos     = 0L;	
  
  /* grNb, elemNb & elemLength */
  if ((theErr = (PapyShort) Papy3FRead (gPapyFile [inFileNb], &i, 1L, theTmpBufP)) < 0)
  {
    Papy3FClose (&gPapyFile [inFileNb]);
    RETURN (theErr);
  } /* if */
  theOffset1 = Extract4Bytes (inFileNb, theTmpBufP, &thePos);
  theOffset2 = Extract4Bytes (inFileNb, theTmpBufP, &thePos);
  theOffset3 = Extract4Bytes (inFileNb, theTmpBufP, &thePos);
  Papy3FSeek (gPapyFile [inFileNb], SEEK_CUR, (PapyLong) 48L);
    
  if (theNbOfSegments == 1) 
  {
    /* read rle image */
    theRleLen = theLength - 64L;
    theRleP = (PapyUChar *) emalloc3 ((PapyULong) (theRleLen + 10L));
    /* extract the image from the file */
    theErr = Papy3FRead (gPapyFile [inFileNb], (PapyULong *) &theRleLen, 1L, (void *) theRleP);

    DecodeRLESegment (inFileNb, ioImage16P, theRleP, theRleLen, theNbOfSegments, 1);
    /* delete rle image */
    efree3 ((void **) &theRleP);
  }/* if ...single segment */
    
  else if (theNbOfSegments == 2) 
  {
    /* deal with first segment */
    theRleLen = theOffset2 - 64L;
    theRleP   = (PapyUChar *) emalloc3 ((PapyULong) (theRleLen + 10L));
    /* extract the image from the file */
    theErr = Papy3FRead (gPapyFile [inFileNb], (PapyULong *) &theRleLen, 1L, (void *) theRleP);
    DecodeRLESegment (inFileNb, ioImage16P, theRleP, theRleLen, theNbOfSegments, 2);
    /* delete rle image */
    efree3 ((void **) &theRleP);
      
    /* deal with second segment */
    theRleLen = theLength - theOffset2;
    theRleP   = (PapyUChar *) emalloc3 ((PapyULong) (theRleLen + 10L));
    /* extract the image from the file */
    theErr = Papy3FRead (gPapyFile [inFileNb], (PapyULong *) &theRleLen, 1L, (void *) theRleP);
    DecodeRLESegment (inFileNb, ioImage16P, theRleP, theRleLen, theNbOfSegments, 1);
    /* delete rle image */
    efree3 ((void **) &theRleP);
      
  }/* if ...two segments */

  else if (theNbOfSegments == 3) 
  {
    /* deal with first segment */
    theRleLen = theOffset2 - 64L;
    theRleP   = (PapyUChar *) emalloc3 ((PapyULong) (theRleLen + 10L));
    /* extract the image from the file */
    theErr = Papy3FRead (gPapyFile [inFileNb], (PapyULong *) &theRleLen, 1L, (void *) theRleP);
    DecodeRLESegment (inFileNb, ioImage16P, theRleP, theRleLen, theNbOfSegments, 1);
    /* delete rle image */
    efree3 ((void **) &theRleP);
      
    /* deal with second segment */
    theRleLen = theOffset3 - theOffset2;
    theRleP   = (PapyUChar *) emalloc3 ((PapyULong) (theRleLen + 10L));
    /* extract the image from the file */
    theErr = Papy3FRead (gPapyFile [inFileNb], (PapyULong *) &theRleLen, 1L, (void *) theRleP);
    DecodeRLESegment (inFileNb, ioImage16P, theRleP, theRleLen, theNbOfSegments, 2);
    /* delete rle image */
    efree3 ((void **) &theRleP);
      
    /* deal with third segment */
    theRleLen = theLength - theOffset3;
    theRleP   = (PapyUChar *) emalloc3 ((PapyULong) (theRleLen + 10L));
    /* extract the image from the file */
    theErr = Papy3FRead (gPapyFile [inFileNb], (PapyULong *) &theRleLen, 1L, (void *) theRleP);
    DecodeRLESegment (inFileNb, ioImage16P, theRleP, theRleLen, theNbOfSegments, 3);
    /* delete rle image */
    efree3 ((void **) &theRleP);
      
  } /* if ...three segments */
  
  return 0;
    
} /* endof ExtractRLE */


static inline unsigned short readUint16(const unsigned char *data)
{
  return (((unsigned short)(*data) << 8) | ((unsigned short)(*(data+1))));
}

static unsigned char scanJpegDataForBitDepth(
  const unsigned char *data,
  const long fragmentLength)
{
  long offset = 0;
  while(offset+4 < fragmentLength)
  {
	unsigned short val = readUint16(data+offset);
    switch( val)
    {
      case 0xffc0: // SOF_0: JPEG baseline
        return data[offset+4];
        /* break; */
      case 0xffc1: // SOF_1: JPEG extended sequential DCT
        return data[offset+4];
        /* break; */
      case 0xffc2: // SOF_2: JPEG progressive DCT
        return data[offset+4];
        /* break; */
      case 0xffc3 : // SOF_3: JPEG lossless sequential
        return data[offset+4];
        /* break; */
      case 0xffc5: // SOF_5: differential (hierarchical) extended sequential, Huffman
        return data[offset+4];
        /* break; */
      case 0xffc6: // SOF_6: differential (hierarchical) progressive, Huffman
        return data[offset+4];
        /* break; */
      case 0xffc7: // SOF_7: differential (hierarchical) lossless, Huffman
        return data[offset+4];
        /* break; */
      case 0xffc8: // Reserved for JPEG extentions
        offset += readUint16(data+offset+2)+2;
        break;
      case 0xffc9: // SOF_9: extended sequential, arithmetic
        return data[offset+4];
        /* break; */
      case 0xffca: // SOF_10: progressive, arithmetic
        return data[offset+4];
        /* break; */
      case 0xffcb: // SOF_11: lossless, arithmetic
        return data[offset+4];
        /* break; */
      case 0xffcd: // SOF_13: differential (hierarchical) extended sequential, arithmetic
        return data[offset+4];
        /* break; */
      case 0xffce: // SOF_14: differential (hierarchical) progressive, arithmetic
        return data[offset+4];
        /* break; */
      case 0xffcf: // SOF_15: differential (hierarchical) lossless, arithmetic
        return data[offset+4];
        /* break; */
      case 0xffc4: // DHT
        offset += readUint16(data+offset+2)+2;
        break;
      case 0xffcc: // DAC
        offset += readUint16(data+offset+2)+2;
        break;
      case 0xffd0: // RST m
      case 0xffd1:
      case 0xffd2:
      case 0xffd3:
      case 0xffd4:
      case 0xffd5:
      case 0xffd6:
      case 0xffd7:
        offset +=2;
        break;
      case 0xffd8: // SOI
        offset +=2;
        break;
      case 0xffd9: // EOI
        offset +=2;
        break;
      case 0xffda: // SOS
        offset += readUint16(data+offset+2)+2;
        break;
      case 0xffdb: // DQT
        offset += readUint16(data+offset+2)+2;
        break;
      case 0xffdc: // DNL
        offset += readUint16(data+offset+2)+2;
        break;
      case 0xffdd: // DRI
        offset += readUint16(data+offset+2)+2;
        break;
      case 0xffde: // DHP
        offset += readUint16(data+offset+2)+2;
        break;
      case 0xffdf: // EXP
        offset += readUint16(data+offset+2)+2;
        break;
      case 0xffe0: // APPn
      case 0xffe1:
      case 0xffe2:
      case 0xffe3:
      case 0xffe4:
      case 0xffe5:
      case 0xffe6:
      case 0xffe7:
      case 0xffe8:
      case 0xffe9:
      case 0xffea:
      case 0xffeb:
      case 0xffec:
      case 0xffed:
      case 0xffee:
      case 0xffef:
        offset += readUint16(data+offset+2)+2;
        break;
      case 0xfff0: // JPGn
      case 0xfff1:
      case 0xfff2:
      case 0xfff3:
      case 0xfff4:
      case 0xfff5:
      case 0xfff6:
      case 0xfff7:
      case 0xfff8:
      case 0xfff9:
      case 0xfffa:
      case 0xfffb:
      case 0xfffc:
      case 0xfffd:
        offset += readUint16(data+offset+2)+2;
        break;
      case 0xfffe: // COM
        offset += readUint16(data+offset+2)+2;
        break;
      case 0xff01: // TEM
        break;
      default:
        if ((data[offset]==0xff) && (data[offset+1]>2) && (data[offset+1] <= 0xbf)) // RES reserved markers
        {
          offset += 2;
        }
        else return 0; // syntax error, stop parsing
        break;
    }
  } // while
  return 0; // no SOF marker found
}


/********************************************************************************/
/*									 	*/
/*	Papy3GetPixelData : gets the specified image or icon and put it either	*/
/* 	in the passed module or group. The moduleId parameter should contain 	*/
/*	the value IconImage if one wants to extract an icon or ImagePixel if 	*/
/*	one wants to extract the image itsself, wether a group or a module has	*/
/*	been passed to the routine (this is important).				*/
/*	BEWARE : in case of extracting the pixel data to a module, you should 	*/
/*	have gotten the module before calling this routine.			*/
/*		 in case of extracting the pixel data to a group, you should 	*/
/*	have read the group 0x0028 and the group 0x7FE0 before calling this 	*/
/* 	routine.								*/
/* 	return : the image, or NULL if something went wrong			*/
/*										*/
/********************************************************************************/

PapyUShort * CALLINGCONV
Papy3GetPixelData (PapyShort inFileNb, int inImageNb, SElement *inGrOrModP, int inModuleId)
{
  PapyUChar	 *theBufP, theTmpBuf [256], *theTmpBufP;
  PapyUChar	 *theCharP, theChar0, theChar1;
  PapyUShort	 *theUShortP, theUShort1, theUShort2;
  PapyShort	 theErr;
  int		 theFrameCount = 1, theLoop, ok, theIsModule;
  PAPY_FILE	 theFp;
  PapyULong	 theBytesToRead, i, theULong, thePos, *theOffsetTableP=NULL;
  PapyULong	 theRefPoint, thePixelStart;
  SElement 	 *theElemP;	/* work pointer on the element of the module */
  
  
  /* some usefull tests */
  if (inImageNb > gArrNbImages [inFileNb] || inModuleId > END_MODULE) return NULL;
  
  /* test to learn if the routine was passed a module or a group in parameter */
  if (inGrOrModP->group == 0x0028) theIsModule = TRUE;
  else theIsModule = FALSE;
  

  
  /* get the file pointer from the file number */
  theFp = gPapyFile [inFileNb];

  /* position the file pointer to the pixel data to read */
  switch (inModuleId)
  {
    case IconImage :
      /* only allow to get an icon from a PAPYRUS 3 file */
      if (gIsPapyFile [inFileNb] != PAPYRUS3) return NULL;
      
      /* it is one of the pointer sequence module, so go to the given ptr sequence */
      if (Papy3FSeek (gPapyFile [inFileNb], (int) SEEK_SET, (PapyLong) gOffsetToPtrSeq [inFileNb] + 8L) != 0)
        return NULL;
      
      /* look for the given item of the ptr seq */
      for (i = 1L; i < (PapyULong) inImageNb; i++)
      {
        theBytesToRead = Papy3ExtractItemLength (inFileNb);
        if (Papy3FSeek (gPapyFile [inFileNb], (int) SEEK_CUR, (PapyLong) theBytesToRead) != 0)
          return NULL;
      } /* for */
      
      /* then points to the first element of the item */
      if (Papy3FSeek (gPapyFile [inFileNb], (int) SEEK_CUR, (PapyLong) 8L) != 0) 
        return NULL;
        
      /* look now for the right group, i.e. image */
      if ((theErr = Papy3GotoGroupNb (inFileNb, 0x7FE0)) < 0) return NULL;
      /* ... then the right element */
      theErr = Papy3GotoElemNb (inFileNb, 0x7FE0, 0x0010, &theBytesToRead);
        
      /* jump over the description of the element */
      if (gArrTransfSyntax [inFileNb] == LITTLE_ENDIAN_IMPL)
      {
        if (Papy3FSeek (gPapyFile [inFileNb], (int) SEEK_CUR, (PapyLong) 8L) != 0) 
          return NULL;
      } /* if */
      else if (gArrTransfSyntax [inFileNb] == LITTLE_ENDIAN_EXPL || gArrTransfSyntax [inFileNb] == BIG_ENDIAN_EXPL) 
      {
        if (Papy3FSeek (gPapyFile [inFileNb], (int) SEEK_CUR, (PapyLong) 12L) != 0) 
          return NULL;
      } /* else */

      /* position to the right element knowing if it is a group or a module */
      if (theIsModule)
        theElemP = inGrOrModP + papPixelDataII;
      else
        theElemP = inGrOrModP + papPixelDataGr;
      break;
    
    case ImagePixel :
      /* go to the begining of the specified image */
      if (Papy3FSeek (gPapyFile [inFileNb], (int) SEEK_SET, (PapyLong) *(gRefPixelOffset [inFileNb] + inImageNb - 1)) != 0)
        return NULL;

      /* position to the right element knowing if it is a group or a module */
      if (theIsModule)
        theElemP = inGrOrModP + papPixelData;
      else
        theElemP = inGrOrModP + papPixelDataGr;
      break;
    
    default :
      return NULL;
      break;
  } /* switch */
    
  
  /* get the size of the pixel data */
  if (inModuleId == IconImage) 
    theBytesToRead = theElemP->length;
  else
    theBytesToRead = (PapyULong) gx0028Rows [inFileNb] * (PapyULong) gx0028Columns [inFileNb] * 
    		     (PapyULong) (((gx0028BitsAllocated [inFileNb] - 1) / 8) + 1L);

    
	if( gx0028BitsAllocated [inFileNb] == 16 && gx0028BitsStored [inFileNb] == 8 && gArrPhotoInterpret [inFileNb] == RGB)   //ANTOINE: Sometimes these DICOMs are really strange???
	{
		gx0028BitsAllocated[inFileNb] = 8;
		
		theBytesToRead = (PapyULong) gx0028Rows [inFileNb] * (PapyULong) gx0028Columns [inFileNb] * 
    		     (PapyULong) (((gx0028BitsAllocated [inFileNb] - 1) / 8) + 1L);
	}
	

  /* if it is a RGB or a YBR_FULL image, multiply the bytes to read by 3 */
  if (inModuleId == ImagePixel && 
      (gArrPhotoInterpret [inFileNb] == RGB ||
	   gArrPhotoInterpret [inFileNb] == UNKNOWN_COLOR  ||
       gArrPhotoInterpret [inFileNb] == YBR_FULL  ||
	   gArrPhotoInterpret [inFileNb] == YBR_ICT  ||
	   gArrPhotoInterpret [inFileNb] == YUV_RCT  ||
	   gArrPhotoInterpret [inFileNb] == YBR_RCT)) theBytesToRead *= 3L;
  else if (inModuleId == ImagePixel && 
           (gArrPhotoInterpret [inFileNb] == YBR_FULL_422 ||
            gArrPhotoInterpret [inFileNb] == YBR_PARTIAL_422)) theBytesToRead *= 3L;
  
  goImageSize[ inFileNb] = theBytesToRead;
  /* allocate the memory for the pixel data */
  theBufP = (PapyUChar *) emalloc3 ((PapyULong) theBytesToRead);
  
  
  /* image reading depending on the image encoding */
  
  /* first test if the images is not encoded */
  if (inModuleId == IconImage || 
      (gArrCompression [inFileNb]     == NONE &&
       (gArrPhotoInterpret [inFileNb] == MONOCHROME1 	||
        gArrPhotoInterpret [inFileNb] == MONOCHROME2 	||
        gArrPhotoInterpret [inFileNb] == PALETTE     	||
        gArrPhotoInterpret [inFileNb] == RGB         	||
        gArrPhotoInterpret [inFileNb] == YBR_FULL  	||
		gArrPhotoInterpret [inFileNb] == UNKNOWN_COLOR  ||
        gArrPhotoInterpret [inFileNb] == YBR_FULL_422	||
        gArrPhotoInterpret [inFileNb] == YBR_RCT  	||
        gArrPhotoInterpret [inFileNb] == YBR_ICT	||
		gArrPhotoInterpret [inFileNb] == YUV_RCT	||
        gArrPhotoInterpret [inFileNb] == YBR_PARTIAL_422)))
  {    
    /* if it is a DICOM file then jump to the right image */
    if (gIsPapyFile [inFileNb] == DICOM10 || gIsPapyFile [inFileNb] == DICOM_NOT10)
      theErr = Papy3FSeek (theFp, SEEK_CUR, (PapyLong) (theBytesToRead * (inImageNb - 1)));
    
    /* read theBytesToRead bytes from the file */
    if ((theErr = (PapyShort) Papy3FRead (theFp, &theBytesToRead, 1L, theBufP)) < 0)
    {
      theErr = Papy3FClose (&theFp);
      efree3 ((void **) &theBufP);
      return NULL;
    } /* if */
    
    /* swap the bytes if necessary */
    if (inModuleId == ImagePixel && gx0028BitsAllocated [inFileNb] > 8)
    {
	    if( gx0028BitsAllocated [inFileNb] > 16)
	    {
		    register PapyULong	*theULongP = (PapyULong *) theBufP;
		    register long			ii;
  		  
		    ii = theBytesToRead / 4;
  		  
		    #if __BIG_ENDIAN__
		    {
			  if( gArrTransfSyntax [inFileNb] != BIG_ENDIAN_EXPL)
			  {
			    while( ii-- > 0)
			    {
			      *theULongP++ = OSSwapLittleToHostInt32( *theULongP);
			    }
			  }
		    }
		    #else
		    if( gArrTransfSyntax [inFileNb] != BIG_ENDIAN_EXPL)
		    {
		    }
		    else
		    {
			    while( ii-- > 0)
			    {
				  *theULongP++ = OSSwapBigToHostInt32( *theULongP);
			    }
		    }
		    #endif
	    }
	    else
	    {
		    register PapyUShort	 *theUShortP = (PapyUShort *) theBufP;
		    register long			ii;
  		  
		    ii = theBytesToRead / 2;
			  
		    #if __BIG_ENDIAN__
		    {
			  if( gArrTransfSyntax [inFileNb] != BIG_ENDIAN_EXPL)
			  {
			    while( ii-- > 0)
			    {
			      *theUShortP++ = OSSwapLittleToHostInt16( *theUShortP);
			    }
			  }
		    }
		    #else
		    if( gArrTransfSyntax [inFileNb] != BIG_ENDIAN_EXPL)
		    {
  		  
		    }
		    else
		    {
			  while( ii-- > 0)
			    {
				  *theUShortP++ = OSSwapBigToHostInt16( *theUShortP);
			    }
		    }
		    #endif
	    }
  	  
    } /* if ...more than 8 bits depth image */
    
  } /* if ...module IconImage or photometric interpretation is monochrome/palette/rgb */
  
  /* *** not IconImage module and the pixels are compressed *** */
  else
  {
    /* if the image conforms to the DICOM standard there should be an offset table */
    if (!(gArrTransfSyntax [inFileNb] == LITTLE_ENDIAN_IMPL && 
    	  gArrCompression [inFileNb]  != NONE))
    {
      /* check to see if there is an offset table, as expected */
      /* so read 8 chars from the file */
      theTmpBufP = (unsigned char *) &theTmpBuf [0];
      i = 8L; 					/* grNb, elemNb & elemLength */
      if ((theErr = (PapyShort) Papy3FRead (theFp, &i, 1L, theTmpBufP)) < 0)
      {
        theErr = Papy3FClose (&theFp);
        return NULL;
      } /* if */
    
      thePos     = 0L;
		  theGroup   = Extract2Bytes (inFileNb, theTmpBufP, &thePos);
		  theElement = Extract2Bytes (inFileNb, theTmpBufP, &thePos);

      /* test if the values are correct */
      if (theUShort1 != 0xFFFE || theUShort2 != 0xE000)
        return NULL;
    
      /* offset table size */
      /* extract the element length according to the little-endian syntax */
      theULong = Extract4Bytes (inFileNb, theTmpBufP, &thePos);
    
      {
        ok = FALSE;
        theFrameCount = 0;
      
      
	    if( gCachedFramesMap[ inFileNb] == 0)
	    {
        /* initialize a file reference point */
        Papy3FTell (theFp, (PapyLong *) &theRefPoint);
      	Papy3FSeek (theFp, SEEK_SET, theRefPoint + theULong);
  			
			  /* allocate memory for the offset table */
			  #define MAX_NUMBER_OF_FRAMES 100000L
			  theOffsetTableP = (PapyULong *) emalloc3 ((PapyULong) (100000L * sizeof (PapyULong)));
  		  
			  while (!ok)
			  {
			    /* read fragment information : 0xFFFE, 0xE000, length */
			    Papy3FTell (theFp, (PapyLong *) &thePixelStart);
  			  
			    /* read 8 chars from the file */
			    i 	      = 8L;
			    thePos      = 0L;
			    theTmpBufP  = (unsigned char *) &theTmpBuf [0];
			    if ((theErr = (PapyShort) Papy3FRead (theFp, &i, 1L, theTmpBufP)) < 0)
			    {
				  theErr = Papy3FClose (&theFp);
				  efree3 ((void **) &theOffsetTableP);
				  return NULL;
			    } /* if */
  			
			    thePos = 0L;
			    theUShort1 = Extract2Bytes (inFileNb, theTmpBufP, &thePos);
			    theUShort2 = Extract2Bytes (inFileNb, theTmpBufP, &thePos);
			    theULong   = Extract4Bytes (inFileNb, theTmpBufP, &thePos);
  			
			    /* offset table found ? */
			    if ((theUShort1 == 0xFFFE) && (theUShort2 == 0xE000))
			    {
				  theOffsetTableP [theFrameCount] = thePixelStart - theRefPoint;
  				
				  if( theFrameCount >= MAX_NUMBER_OF_FRAMES)
					  fprintf(stdout, "*********** MAJOR MEMORY BUG : theFrameCount > MAX_NUMBER_OF_FRAMES -> THIS APP WILL CRASH !\r");
				  else theFrameCount ++;
  				
				  Papy3FSeek (theFp, SEEK_CUR, theULong);
			    } /* if */
			    else if ((theUShort1 == 0xFFFE) && (theUShort2 == 0xE0DD)) ok = TRUE;
  			  
			  } /* while */
  	      
        /* position the file pointer on the first image */
        Papy3FSeek (theFp, SEEK_SET, theRefPoint);
      
       	gCachedFramesMap[ inFileNb] = theOffsetTableP;
      }
		  else theOffsetTableP = gCachedFramesMap[ inFileNb];

      } /* else */
    
    } /* if ...not a Papyrus compressed image */
    else
    {
      /* allocate room to store the offset table */
      theOffsetTableP = (PapyULong *) emalloc3 ((PapyULong) (sizeof (PapyULong)));
      
      /* there is no offset to the JPEG encoded image */
      theOffsetTableP [0] = 0L;
    } /* else ...Papyrus compressed image */
  
    /* get the position of the first pixel */
    Papy3FTell (theFp, (PapyLong *) &thePixelStart);            
  
    /* in case of a PAPYRUS file, there should be only one frame. */
    /* The positioning of the file pointer to the right image has already been performed */
    if (gIsPapyFile [inFileNb] == PAPYRUS3) inImageNb = 1;
    
    
    /*  *** different ways of reading depending on the compression algorithm *** */

  
    /********************************************************************/
    /*******************     Lossless JPEG     **************************/
    /********************************************************************/
    if (gArrCompression [inFileNb] == JPEG_LOSSLESS)
    {
      if (gArrTransfSyntax [inFileNb] == LITTLE_ENDIAN_EXPL)
        theErr = ExtractJPEGlosslessDicom (inFileNb, theBufP, thePixelStart, theOffsetTableP, inImageNb);
      else /* little-endian-explicit VR */
        theErr = ExtractJPEGlosslessPap (inFileNb, theBufP, thePixelStart, theElemP->length);
    } /* if ...JPEG lossless */

  
    /********************************************************************/
    /*******************     Lossy JPEG     *****************************/
    /********************************************************************/
    else if (gArrCompression [inFileNb] == JPEG_LOSSY)
    {	
      theErr = ExtractJPEGlossy (inFileNb, theBufP, thePixelStart, theOffsetTableP, inImageNb, 
                              (int) gx0028BitsAllocated [inFileNb]);
		
    } /* if ...JPEG lossy */

#ifdef MAYO_WAVE    
    /********************************************************************/
    /*******************     MAYO WAVELET   *****************************/
    /********************************************************************/
    else if (gArrCompression [inFileNb] == MAYO_WAVELET)
    {
      theErr = ExtractWavelet ((PapyShort) inFileNb, (PapyUChar *) theBufP, thePixelStart,
                            theOffsetTableP, inImageNb, (int) gx0028BitsAllocated [inFileNb]);

    } /* if ...Mayo Wavelet */
#endif
  
    /********************************************************************/
    /*******************     RLE     ************************************/
    /********************************************************************/
    else if (gArrCompression [inFileNb] == RLE)
    {
      theErr = ExtractRLE (inFileNb, (PapyUShort *) theBufP, thePixelStart, theOffsetTableP, inImageNb);
    } /* if ...Run Length Encoding */

  
    /********************************************************************/
    /*******************     unknown     ********************************/
    /********************************************************************/
    else
    {
      /* black image, that is better than nothing ... */
      for (i = 0L; i < theBytesToRead; i++) theBufP [i] = 0;
    } /* if ...nothing known */
    
  } /* else ...not icon image or compressed pixel data */
  
  
  /* allocate room in the element in order to put the pixel data in the module */
  theElemP->value  = (UValue_T *) emalloc3 ((PapyULong) sizeof (UValue_T));
  theElemP->nb_val = 1L;
  
  /* extract the pixel data depending on the value representation */
  if (inModuleId == IconImage || gArrCompression [inFileNb] == NONE)
  {
    theElemP->vr = OW;
    theElemP->value->ow = (PapyUShort *) theBufP;
  } /* if ...icon image or uncompressed file */
  else
  {
    theElemP->vr = OB;
    theElemP->value->a = (char *) theBufP;
  } /* else ...compressed image */
  
  if (theOffsetTableP != NULL) efree3 ((void **) &theOffsetTableP);
    
  return (PapyUShort *) theBufP;
  
} /* endof Papy3GetPixelData */


/********************************************************************************/
/*										*/
/*	Extract2Bytes : extract a 2-Bytes value (USS, SS or AT) from the buf and*/
/*	increment pos accordingly.						*/
/* 	return : the extracted value						*/
/*										*/
/********************************************************************************/

PapyUShort
Extract2Bytes (unsigned char *inBufP, PapyULong *ioPosP)

/*unsigned char *inBufP;				 the buffer to read from */
/*PapyULong 	*ioPosP;			      the position in the buffer */
{
  PapyUShort 		theUShort;
  unsigned char		*theCharP;


  /* points to the right place in the buffer */
  theCharP  = inBufP;
  theCharP += *ioPosP;
  /* updates the current position in the read buffer */
  *ioPosP += 2;
    
  /* extract the element according to the little-endian syntax */
  theUShort  = (PapyUShort) (*(theCharP + 1));
  theUShort  = theUShort << 8;
  theUShort |= (PapyUShort) *theCharP;
    
  return theUShort;

} /* endof Extract2Bytes */



/********************************************************************************/
/*										*/
/*	Extract4Bytes : extract a 4-Bytes value (UL, SL or FL) of the buf and 	*/
/*	increment pos accordingly.						*/
/* 	return : the extracted value					 	*/
/*										*/
/********************************************************************************/

PapyULong
Extract4Bytes (unsigned char *inBufP, PapyULong *ioPosP)

/*unsigned char *inBufP;				 the buffer to read from */
/*PapyULong 	*ioPosP;			      the position in the buffer */
{
  unsigned char	*theCharP;
  PapyULong	theULong = 0L, theTmpULong;
    
    
  /* points to the right place in the buffer */
  theCharP  = inBufP;
  theCharP += *ioPosP;
  /* updates the current position in the read buffer */
  *ioPosP += 4;
    
  /* extract the element according to the little-endian syntax */
  theTmpULong  = (PapyULong) (*(theCharP + 3));
  theTmpULong  = theTmpULong << 24;
  theULong    |= theTmpULong;
  theTmpULong  = (PapyULong) (*(theCharP + 2));
  theTmpULong  = theTmpULong << 16;
  theULong    |= theTmpULong;
  theTmpULong  = (PapyULong) (*(theCharP + 1));
  theTmpULong  = theTmpULong << 8;
  theULong    |= theTmpULong;
  theTmpULong  = (PapyULong) *theCharP;
  theULong    |= theTmpULong;
    
  return theULong;
    
} /* endof Extract4Bytes */



/********************************************************************************/
/*										*/
/*	Extract8Bytes : extract a 8-Bytes value (FD) of the buf and 		*/
/*	increment pos accordingly.						*/
/* 	return : the extracted value					 	*/
/*										*/
/********************************************************************************/

PapyFloatDouble
Extract8Bytes (unsigned char *inBufP, PapyULong *ioPosP)

/*unsigned char *inBufP;				 the buffer to read from */
/*PapyULong 	*ioPosP;			      the position in the buffer */
{
  unsigned char		*theCharP, theDoubleArr [8], i;
  PapyFloatDouble	*theFloatDoubleP;
    
    
  /* points to the right place in the buffer */
  theCharP  = inBufP;
  theCharP += *ioPosP;
  /* updates the current position in the read buffer */
  *ioPosP  += 8;
    
  /* extract the element according to the little-endian syntax */
  for (i = 0; i < 4; i++)
  {
    theDoubleArr [2 * i]       = *theCharP;
    theDoubleArr [(2 * i) + 1] = *(theCharP + 1);
    theCharP += 2;
  } /* for ...extraction of the value */
    
  theFloatDoubleP = (PapyFloatDouble *) &theDoubleArr;
    
  return *theFloatDoubleP;
    
} /* endof Extract8Bytes */



/********************************************************************************/
/*										*/
/*	ExtractString : extract a string from the buffer and put it in the 	*/
/*	given element. Increment pos accordingly.				*/
/*										*/
/********************************************************************************/

void
ExtractString (SElement *ioElemP, unsigned char *inBufP, PapyULong *ioBufPosP, 
	       PapyULong inElemLength)
{
  char			*theStringP, *theP, *theCharValP, *theCharWrkP;
  unsigned char		*theTmpP;
  int			ii, j, theStringLength;
  
		  				   /* 1 for the string terminator */
  theStringP = (char *) emalloc3 ((PapyULong) (inElemLength + 1));
  theP = theStringP;
  theTmpP = inBufP;
  /* extract the element from the buffer */
  for (ii = 0L; ii < (int) inElemLength; ii++, (*ioBufPosP)++)
    *(theP++) = theTmpP [*ioBufPosP];
    
  theStringP [ii] = '\0';
    
  theCharValP = theStringP;
 
  theStringLength = strlen (theCharValP); 
          
  ioElemP->nb_val = 1L;     /* number of strings */
  theCharWrkP = theCharValP;
          
  /* count the number of strings */
  for (j = 0; j < theStringLength; j++, theCharWrkP ++)
  {
    /* value separator */
    if (*theCharWrkP == '\\') 
    {
      ioElemP->nb_val++;
      *theCharWrkP = '\0';
    } /* if */
  } /* for ...counting the number of values */
          
  ioElemP->value = (UValue_T *) ecalloc3 ((PapyULong) ioElemP->nb_val, (PapyULong) sizeof (UValue_T));
          	
  /* extraction of the strings */	
  for (j = 0, theCharWrkP = theCharValP; 
       j < (int) ioElemP->nb_val;
       j ++, theCharWrkP += theStringLength + 1)
  {
    theStringLength = strlen (theCharWrkP);
		    
    /* addition to delete the blank if odd string */
    if (ioElemP->vr == UI)
    {
      /* suppress the blank by shifting all the chars to the left */
      /* old was : theCharWrkP [theStringLength - 1] == '0') */
      if (theCharWrkP [theStringLength - 1] == 0x00) 
	      theCharWrkP [theStringLength - 1] = '\0';
    } /* then ...VR = UI */
    else
    {
      if (theCharWrkP [theStringLength - 1] == ' ')
	      theCharWrkP [theStringLength - 1] = '\0';
    } /* else ...VR <> UI */
		    
    ioElemP->value [j].a = theCharWrkP;

  } /* for ...extraction of the strings */
          
} /* endof ExtractString */


										
/********************************************************************************/
/*									 	*/
/*	PutBufferInElement3 : fill_in an element structure (one element) 	*/
/* 	from a buffer made of unsigned chars					*/
/* 	return : standard error message						*/
/*									  	*/
/********************************************************************************/

PapyShort
PutBufferInElement3 (PapyShort inFileNb, unsigned char *ioBuffP, PapyULong inElemLength,
		    SElement *ioElemP, PapyULong *ioBufPosP, PapyLong inInitFilePos)
{
  Item          *theSeqItemP, *theDSitemP;
  Object        *theObjectP, *theObjectP2;  /* MAL */
  SElement      *theSeqGroupP;
  UValue_T      *theValueTP;
  unsigned char *theTmp0P, theTmp1, *theCharP;
  unsigned char theDoubleArr [8], theIncr;
  PapyLong      theCurrFilePos, theInitialFilePos = inInitFilePos;
  PapyULong     ii, i, j, thePosInSeq, thePosInItem, theSeqSize, theSeqGrSize, theImLength;
  PapyULong     theTmpULong, theULong = 0L;
  PapyUShort    theSeqGrNb, theElemNb, *theTmpUsP;   /* *imOW */
  char          *theCharValP, *theCharWrkP; 
  char          *theStringP, *theP;
  int           theEnumSeqNb, theStringLength, theFirstTime, theIsUndefItemLen;
  PapyShort     theErr;


  /* extract the element depending on the value representation */
  switch (ioElemP->vr)
  {
    case RET :
      *ioBufPosP += ioElemP->length;
      break;

    case SS :				/* 16 bits binary signed */
      ioElemP->nb_val = (PapyULong) (inElemLength / 2);
      ioElemP->value = (UValue_T *) ecalloc3 ((PapyULong) ioElemP->nb_val,
					      (PapyULong) sizeof (UValue_T));
      theValueTP = ioElemP->value;
      for (j = 0; j < ioElemP->nb_val; j++, theValueTP++)
      {
        /* points to the right place in the buffer */
        theTmp0P    = ioBuffP;
        theTmp0P   += *ioBufPosP;
        /* updates the current position in the read buffer */
        *ioBufPosP += 2L;  
        /* extract the element according to the little-endian syntax */
        theValueTP->ss  = (PapyUShort) (*(theTmp0P + 1));
        theValueTP->ss  = theValueTP->ss << 8;
        theValueTP->ss |= (PapyUShort) *theTmp0P;
      } /* for */
	    
      break; /* SS */

	  
    case AT :
    case USS :				/* 16 bits binary unsigned */
      ioElemP->nb_val = (PapyULong) (inElemLength / 2);
      ioElemP->value  = (UValue_T *) ecalloc3 ((PapyULong) ioElemP->nb_val,
					       (PapyULong) sizeof (UValue_T));
      theValueTP = ioElemP->value;
      for (j = 0; j < ioElemP->nb_val; j++, theValueTP++)
      {
        /* points to the right place in the buffer */
        theTmp0P    = ioBuffP;
        theTmp0P   += *ioBufPosP;
        /* updates the current position in the read buffer */
        *ioBufPosP += 2L;  
        /* extract the element according to the little-endian syntax */
        theValueTP->us  = (PapyUShort) (*(theTmp0P + 1));
        theValueTP->us  = theValueTP->us << 8;
        theValueTP->us |= (PapyUShort) *theTmp0P;
      } /* for */

      break; /* USS */
	    
        
    case SL :				/* 32 bits binary signed */
      ioElemP->nb_val = (PapyULong) (inElemLength / 4);
      ioElemP->value  = (UValue_T *) ecalloc3 ((PapyULong) ioElemP->nb_val,
					       (PapyULong) sizeof (UValue_T));
      theValueTP = ioElemP->value;
      for (j = 0; j < ioElemP->nb_val; j++, theValueTP++)
      {
        /* points to the right place in the buffer */
        theTmp0P  = ioBuffP;
        theTmp0P += *ioBufPosP;
        /* updates the current position in the read buffer */
        *ioBufPosP += 4L;
        /* extract the element according to the little-endian syntax */
        theTmpULong      = (PapyULong) (*(theTmp0P + 3));
        theTmpULong      = theTmpULong << 24;
        theULong	 = theTmpULong;
        theTmpULong      = (PapyULong) (*(theTmp0P + 2));
        theTmpULong      = theTmpULong << 16;
        theULong	|= theTmpULong;
        theTmpULong      = (PapyULong) (*(theTmp0P + 1));
        theTmpULong      = theTmpULong << 8;
        theULong	|= theTmpULong;
        theTmpULong      = (PapyULong) *theTmp0P;
        theULong        |= theTmpULong;
        theValueTP->sl   = theULong;
      } /* for */

      break; /* SL */
	  
	  
    case UL :				/* 32 bits binary unsigned */
      ioElemP->nb_val = (PapyULong) (inElemLength / 4);
      ioElemP->value  = (UValue_T *) ecalloc3 ((PapyULong) ioElemP->nb_val,
					       (PapyULong) sizeof (UValue_T));
      theValueTP = ioElemP->value;
      for (j = 0; j < ioElemP->nb_val; j++, theValueTP++)
      {
        /* points to the right place in the buffer */
        theTmp0P    = ioBuffP;
        theTmp0P   += *ioBufPosP;
        /* updates the current position in the read buffer */
        *ioBufPosP += 4L;
        /* extract the element according to the little-endian syntax */
        theTmpULong      = (PapyULong) (*(theTmp0P + 3));
        theTmpULong      = theTmpULong << 24;
        theULong	 = theTmpULong;
        theTmpULong      = (PapyULong) (*(theTmp0P + 2));
        theTmpULong      = theTmpULong << 16;
        theULong	|= theTmpULong;
        theTmpULong      = (PapyULong) (*(theTmp0P + 1));
        theTmpULong      = theTmpULong << 8;
        theULong        |= theTmpULong;
        theTmpULong      = (PapyULong) *theTmp0P;
        theULong        |= theTmpULong;
        theValueTP->ul   = theULong;
      } /* for */

      break; /* UL */
	  
	  
    case FL :				/* 32 bits binary floating */
      ioElemP->nb_val = (PapyULong) (inElemLength / 4);
      ioElemP->value  = (UValue_T *) ecalloc3 ((PapyULong) ioElemP->nb_val,
					       (PapyULong) sizeof (UValue_T));
      theValueTP = ioElemP->value;
      for (j = 0; j < ioElemP->nb_val; j++, theValueTP++)
      {
        /* points to the right place in the buffer */
        theTmp0P    = ioBuffP;
        theTmp0P   += *ioBufPosP;
        /* updates the current position in the read buffer */
        *ioBufPosP += 4L;
        /* extract the element according to the little-endian syntax */
        theTmpULong      = (PapyULong) (*(theTmp0P + 3));
        theTmpULong      = theTmpULong << 24;
        theULong	 = theTmpULong;
        theTmpULong      = (PapyULong) (*(theTmp0P + 2));
        theTmpULong      = theTmpULong << 16;
        theULong	|= theTmpULong;
        theTmpULong      = (PapyULong) (*(theTmp0P + 1));
        theTmpULong      = theTmpULong << 8;
        theULong	|= theTmpULong;
        theTmpULong      = (PapyULong) *theTmp0P;
        theULong        |= theTmpULong;
        theValueTP->fl   = (float)theULong;
      } /* for */

      break; /* FL */
	  
	  
    case FD :				/* 64 bits binary floating */
      ioElemP->nb_val = (PapyULong) (inElemLength / 8);
      ioElemP->value  = (UValue_T *) ecalloc3 ((PapyULong) ioElemP->nb_val,
					       (PapyULong) sizeof (UValue_T));
      theValueTP = ioElemP->value;
      for (j = 0; j < ioElemP->nb_val; j++, theValueTP++)
      {
        /* points to the right place in the buffer */
        theTmp0P    = ioBuffP;
        theTmp0P   += *ioBufPosP;
        /* updates the current position in the read buffer */
        *ioBufPosP += 8L;
    
        /* extract the element according to the little-endian syntax */
        for (theIncr = 0; theIncr < 4; theIncr++)
        {
          theDoubleArr [2 * theIncr]       = *theTmp0P;
          theDoubleArr [(2 * theIncr) + 1] = *(theTmp0P + 1);
          theTmp0P += 2;
        } /* for ...extraction of the value */
    
        theValueTP->fd = *((PapyFloatDouble *) &theDoubleArr);
        
      } /* for */

      break; /* FD */
    
    case OB :				/* 1 byte image  */
      ioElemP->nb_val = (PapyULong) 1L;
      ioElemP->value  = (UValue_T *) emalloc3 ((PapyULong) sizeof (UValue_T));
      
      /* allocate room for the element */
      theCharP = (unsigned char *) emalloc3 ((PapyULong) inElemLength);
      
      /* copy the bits of the image to the value */
      theTmp0P = theCharP;
      ioBuffP += *ioBufPosP;
      for (i = 0L; i < inElemLength; theTmp0P++, ioBuffP++, i++)
      {
        *theTmp0P = *ioBuffP;
      } /* for */
      
      ioElemP->value->a = (char *) theCharP;
      *ioBufPosP += inElemLength;
      break; /* OB */
    
    case OW :				/* 2 Bytes image */
      theValueTP = ioElemP->value;

      ioElemP->nb_val = (PapyULong) 1L;
      /*ioElemP->value= (UValue_T *) emalloc3 ((PapyULong) sizeof (UValue_T));*/
      ioElemP->value  = (UValue_T *) ecalloc3 ((PapyULong) ioElemP->nb_val,
					      (PapyULong) sizeof (UValue_T));
      theImLength     = inElemLength / 2;
      
      /* pixel data */
      if (ioElemP->group == 0x7FE0 && ioElemP->element == 0x0010)
      {
//#ifndef __alpha
         /* swap the bytes (little endian) */
        for (i = 0L, theCharP = ioBuffP + (*ioBufPosP); i < theImLength; i++, theCharP += 2)
        {
          theTmp1 	  = *(theCharP + 1);
          *(theCharP + 1) = *theCharP;
          *theCharP       =  theTmp1;
        } /* for */
//#endif
        
        ioElemP->value->ow = (PapyUShort *) (ioBuffP + (*ioBufPosP));
      } /* if ...pixel data */
      else /* not pixel data */
      {
        ioElemP->value->ow = (PapyUShort *) ecalloc3 ((PapyULong) theImLength, 
        			        	      (PapyULong) sizeof (PapyUShort));
        /*ioElemP->value->ow = (PapyUShort *) emalloc3 ((PapyULong)theImLength * sizeof (PapyUShort) + 1L);*/

        for (i = 0L, theTmpUsP = ioElemP->value->ow, ioBuffP += *ioBufPosP; i < theImLength; i++, theTmpUsP++, ioBuffP += 2)
        {
          *theTmpUsP  = (PapyUShort) (*(ioBuffP + 1));
    	    *theTmpUsP  = *theTmpUsP << 8;
    	    *theTmpUsP |= (PapyUShort) *ioBuffP;
        } /* for */
       
        /*ioElemP->value->ow = imOW;*/
      } /* else ...not pixel data */
      
      *ioBufPosP += inElemLength;
      break; /* OW */
	    
	  
    case SQ :				/* sequence */
      /* if not the pointer sequence or the image sequence extract the seq */
      if (!(ioElemP->group == 0x0041 && 
            (ioElemP->element == Papy3EnumToElemNb (ioElemP, papPointerSequenceGr) ||
             ioElemP->element == Papy3EnumToElemNb (ioElemP, papImageSequenceGr)))  &&
          !(ioElemP->group == 0x0088 && 
            ioElemP->element == Papy3EnumToElemNb (ioElemP, papIconImageSequenceGr)))
      {
        ioElemP->nb_val    = 1L;
        ioElemP->value     = (UValue_T *) emalloc3 ((PapyULong) sizeof (UValue_T));
        ioElemP->value->sq = NULL;
	      
	      
        /* loop on the items of the sequence */
        thePosInSeq       = 0L;
        theFirstTime      = TRUE;
        /*theIsUndefItemLen = FALSE;*/
        while (thePosInSeq < inElemLength)
        {
          theIsUndefItemLen = FALSE;

          /* read the basic info on the item */
          if ((theSeqGrNb = Extract2Bytes (ioBuffP, ioBufPosP)) != 0xFFFE)
            RETURN (papGroupNumber);
          if ((theElemNb = Extract2Bytes (ioBuffP, ioBufPosP)) != 0xE000)
            RETURN (papElemNumber);
          theSeqSize = Extract4Bytes (ioBuffP, ioBufPosP);
  	
          thePosInSeq += 8L;	/* size of the item delimiter */
          thePosInItem = 0L;	/* the position in this item of the sequence */
	        
          /* if undefined item length, compute it */
          if (theSeqSize == 0xFFFFFFFF)
          {
            /* set a boolean for futur computing of the seq length */
            theIsUndefItemLen = TRUE;
	          
            /* get the current position of the file pointer */
            theErr = Papy3FTell (gPapyFile [inFileNb], &theCurrFilePos);
            /* position the file pointer at the begining of the item */
            theErr = Papy3FSeek (gPapyFile [inFileNb], SEEK_SET, (theInitialFilePos + (PapyLong) (*ioBufPosP)));
	          
            /* computes the item length from the file */
            theSeqSize = 0L;
            theErr     = ComputeUndefinedItemLength3 (inFileNb, &theSeqSize);
	          
            /* reset the file pointer to its previous position */
            theErr = Papy3FSeek (gPapyFile [inFileNb], SEEK_SET, theCurrFilePos);
          } /* if */

          /* creates an empty object that will point to the list of groups */
          theObjectP = (Object *) emalloc3 ((PapyULong) sizeof (Object));
          theObjectP->whoAmI        = papItem;
          theObjectP->item          = NULL;
          theObjectP->module        = NULL;
          theObjectP->group         = NULL;
          theObjectP->record        = NULL;
          theObjectP->tmpFileLength = 0L;
  
          theSeqItemP = InsertLastInList (&(ioElemP->value->sq), theObjectP);
  	        
          /* keep track of the place where to insert a new object (group) */
  	      /*theSeqItemP = theObjectP->item;  /* problem de conservation du pointeur */
 		      

	        /* loop on the groups of the item */
	        /*while (theSeqSize > (thePosInSeq - 8L))*/
	        /*while (inElemLength > (thePosInSeq - 8L))*/
          /*while (inElemLength > thePosInSeq) last = CHG 5.11.99 */
          /*while (theSeqSize > thePosInSeq) last = CHG 8.11.99 */
          /*while (theSeqSize > (thePosInSeq - 8L)) last = CHG 8.11.99 */
          /* or see PapOldPatch.tar.gzip */
          if (theSeqSize > 8L)
            while (theSeqSize > thePosInItem)
            {
              /* read the basic info on the new group */
              theSeqGrNb  = Extract2Bytes (ioBuffP, ioBufPosP);
              theElemNb   = Extract2Bytes (ioBuffP, ioBufPosP);


              /* test if it is the group length element */
              if (theElemNb == 0x0000)
              {
                /* jump over : implicit : the length of the element (1 * 4 bytes) */
                /*	     explicit : the VR and the length of the element (2 * 2 bytes) */
                *ioBufPosP  += 4L;
                theSeqGrSize = Extract4Bytes (ioBuffP, ioBufPosP);
	            
                /* the theFirstTime ioElemP must be taken into account ... */
                theSeqGrSize += 12L;
	            
                /* reset the ioBuffP pos to begining of the group */
                *ioBufPosP   -= 12L;
              } /* if ...ioElemP = group length */
              /* else, we have to compute the group length */
              else
              {
                /* reset the ioBuffP pos to the begining of the group */
                *ioBufPosP   -= 4L;
	      
                /* get the current position of the file pointer */
                theErr = Papy3FTell (gPapyFile [inFileNb], &theCurrFilePos);
	            
                /* position the file pointer at the begining of the item */
                theErr = Papy3FSeek (gPapyFile [inFileNb], SEEK_SET, (theInitialFilePos + (PapyLong) (*ioBufPosP)));
	          
                /* computes the group length */
                theSeqGrSize = ComputeUndefinedGroupLength3 (inFileNb, (PapyLong) theSeqSize);
	            
                /* then reset the file pointer to its previous position */
                theErr = Papy3FSeek (gPapyFile [inFileNb], SEEK_SET, theCurrFilePos);
	          
              } /* else ...compute the group length */
	          
	          
              /* search the enum group number */
              theEnumSeqNb = Papy3ToEnumGroup (theSeqGrNb);
	          
              /* it could be a private group that has an unknown definition */
              if (theEnumSeqNb < 0)
              {
                /* add the group size plus grNb and elemNb */
                thePosInSeq += theSeqGrSize;
                thePosInItem+= theSeqGrSize;
                *ioBufPosP  += theSeqGrSize;
              } /* if ...private group with unknown definition */
              /* known group => extract it from the buffer */
              else
              {	    
                /* create the group */
                theSeqGroupP = Papy3GroupCreate (theEnumSeqNb);
	            
                /* fill the group struct from the content of the buffer */
                theErr = PutBufferInGroup3 (inFileNb, ioBuffP, theSeqGroupP, theSeqGrNb,
		  			    theSeqGrSize, ioBufPosP, theInitialFilePos);
                if (theErr < 0)
                {
                  efree3 ((void **) &ioBuffP);
                  RETURN (theErr);
                } /* if ...theErr */
		    
                thePosInSeq += theSeqGrSize;  /* add the grNb and elemNb */
                thePosInItem+= theSeqGrSize;	
		    
                /* creation of the object that will encapsulate the group */
                theObjectP2 = (Object *) emalloc3 ((PapyULong) sizeof (Object));
                theObjectP2->whoAmI        = papGroup;
                theObjectP2->objID         = theEnumSeqNb;
                theObjectP2->group         = theSeqGroupP;
                theObjectP2->item          = NULL;
                theObjectP2->module 	     = NULL;
                theObjectP2->tmpFileLength = 0L;
	    	
                /* add the object to the list of this element */
                theDSitemP = InsertLastInList ((Item **) &(theObjectP->item), theObjectP2);
              /*theObjectP->item = theDSitemP;  /* MAL added */
	          
	            /* MAL 
              theDSitemP = InsertLastInList ((Item **) &theSeqItemP, theObjectP);
              
              if (theFirstTime)
	            {
	              theFirstTime = FALSE;
	              ioElemP->value->sq->object->item = theDSitemP;
	            } /* if ...theFirstTime time we are in the loop */

	          } /* else ..known group = extracted group */
	        
	        } /* while ...loop on the groups of the item */
	  
	        /* if it was an item with undefined length move the buffer further the delimiter */
	        if (theIsUndefItemLen)
	        {
	          thePosInSeq += 8L;
	          *ioBufPosP  += 8L;
	        } /* if */
	        
	      } /* while ...loop on the items of the sequence */
	      
      } /* if ...not pointer or image sequence */
	    
      /* pointer or image sequence group 41 */
      else 
      {
	      /* there is a value, but set to NULL */
	      ioElemP->nb_val = 1L;		/* CHG */
	      ioElemP->value  = NULL;
      } /* else ...pointer or image sequence */
      break;
          
    case AE :
    case AS :
    case CS :
    case DA :
    case DS :
    case DT :
    case IS :
    case LO :
    case LT :
    case PN :
    case SH :
    case ST :
    case TM :
    case UI :
    case UN :
    case UT :				/* all kind of strings */
      /*theCharValP = ExtractString (ioBuffP, ioBufPosP, inElemLength);*/
		  				   /* 1 for the string terminator */
      theStringP = (char *) emalloc3 ((PapyULong) (inElemLength + 1L));
      theP       = theStringP;
      theTmp0P   = ioBuffP;
      /* extract the element from the buffer */
      for (ii = 0L; ii < inElemLength; ii++, (*ioBufPosP)++) 
      {
        *(theP++) = theTmp0P [(int) *ioBufPosP];
      }
    
      theStringP [ii] = '\0';
    
      theCharValP = theStringP;
 
      theStringLength = strlen (theCharValP); 
          
      ioElemP->nb_val = 1L;     /* number of strings */
      theCharWrkP = theCharValP;
          
      /* count the number of strings */
      for (j = 0; j < (PapyULong)theStringLength; j++, theCharWrkP ++)
      {
        /* value separator */
        if (*theCharWrkP == '\\') 
	      {
	        ioElemP->nb_val++;
	        *theCharWrkP = '\0';
	      } /* if */
      } /* for ...counting the number of values */
          
      ioElemP->value = (UValue_T *) ecalloc3 ((PapyULong) ioElemP->nb_val,
          			              (PapyULong) sizeof (UValue_T));
          	
      /* extraction of the strings */	
      for (j = 0, theCharWrkP = theCharValP; 
           j < ioElemP->nb_val;
           j ++, theCharWrkP += theStringLength + 1)
      {
	      theStringLength = strlen (theCharWrkP);
		          
	      /* addition to delete the blank if odd string */
        if (ioElemP->vr == UI)
        {
          /* suppress the blank by shifting all the chars to the left */
          /* old was : theCharWrkP [theStringLength - 1] == '0') */
          if (theCharWrkP [theStringLength - 1] == 0x00) 
	          theCharWrkP [theStringLength - 1] = '\0';
        } /* then ...VR = UI */
	      else
	        if (theCharWrkP [theStringLength - 1] == ' ')
	          theCharWrkP [theStringLength - 1] = '\0';
		    
	      ioElemP->value[j].a = PapyStrDup (theCharWrkP);

      } /* for ...extraction of the strings */
          
      efree3 ((void **) &theStringP);
          
      break; /* strings */
          
  } /* switch ...value representation */
  
  return 0;
	
} /* endof PutBufferInElement3 */


										
/********************************************************************************/
/*									 	*/
/*	PutBufferInGroup3 : fill_in a group structure (all the elements) 	*/
/* 	from a buffer made of unsigned chars					*/
/* 	return : the enum group number if successfull				*/
/*		 standard error message otherwise 				*/
/*									  	*/
/********************************************************************************/

PapyShort
PutBufferInGroup3 (PapyShort inFileNb, unsigned char *ioBuffP, SElement *ioGroupP,
		   PapyUShort inPapyGrNb, PapyULong inBytesToRead, PapyULong *ioBufPosP,
		   PapyLong inInitFilePos)
{
  SElement		*theArrElemP;
  PapyULong	 	theElemLength;
  PapyULong		j, theInitialBufPos;
  PapyULong		theTmpULong, theULong = 0L;
  PapyLong		theInitialFilePos, theCurrFilePos;
  PapyUShort	 	theGrNb;
  PapyUShort	 	theElemNb, theElemLengthGr2;
  char			theFoo [3], *theFooP;
  unsigned char		*theCharP; 
  int 			theStructPos, theEnumGrNb, i, theIsOld, theIsUndefSeqLen = FALSE;
  int	 		theShadow, theEnabledShadow [0x00FF], theMaxElem;
  PapyShort		theErr, theCreator;
  
  
  theInitialBufPos  = *ioBufPosP;
  theInitialFilePos = inInitFilePos;
  theIsOld          = TRUE; 
  theEnumGrNb       = Papy3ToEnumGroup (inPapyGrNb);   /* gr_nb papyrus -> enum */
  if (theEnumGrNb < 0)				 /* unknown group number */
  {
    efree3 ((void **) &ioBuffP);
    RETURN (papGroupNumber)
  } /* if */
  
  /* the number of elements of this group */
  theMaxElem = gArrGroup [theEnumGrNb].size;
  


  if (inPapyGrNb >= 0x6000 && inPapyGrNb <= 0x6FFF)	/* overlay or UIN overlay */
  {      
    for (j = 0, theArrElemP = ioGroupP; j < (PapyULong)theMaxElem; j++, theArrElemP++)
      theArrElemP->group = inPapyGrNb;
    
  } /* if ...overlay or UINOverlay group */
  
  if (inPapyGrNb % 2 != 0) 			/* is it a shadow group ? */
  {
    theShadow = TRUE; 
    
    /* disables all elements (initialisation) */
    for (i = 0; i < 0x00FF; i++) theEnabledShadow [i] = FALSE;
  } /* then */
  else theShadow = FALSE;
    
  
  theArrElemP = ioGroupP;
  
  while ((*ioBufPosP - theInitialBufPos) < inBytesToRead)	/* loop on the elements */
  {
    theIsUndefSeqLen = FALSE;
    theStructPos = 0;			  /* pos in the array of elements */
    
    /* points to the right place in the buffer */
    theCharP    = ioBuffP;
    theCharP   += *ioBufPosP;
    /* extract the group number according to the little-endian syntax */
    theGrNb     = (PapyUShort) (*(theCharP + 1));
    theGrNb     = theGrNb << 8;
    theGrNb    |= (PapyUShort) *theCharP;
    /* updates the current position in the read buffer */
    *ioBufPosP += 2L;
    /* points to the right place in the buffer */
    theCharP   += 2;
    
    /* extract the element according to the little-endian syntax */
    theElemNb   = (PapyUShort) (*(theCharP + 1));
    theElemNb   = theElemNb << 8;
    theElemNb  |= (PapyUShort) *theCharP;
    /* updates the current position in the read buffer */
    *ioBufPosP += 2L;
    /* points to the right place in the buffer */
    theCharP   += 2;
    
    /* some special test for the group 2 are necessary */
    if (theGrNb == 0x0002)
    {
      /* test to discover which transfert syntax was used to create the file (implicit or explicit VR) */
      theFooP     = (char *) &theFoo [0];
      theFooP [0] = (char)   *theCharP;
      theFooP [1] = (char) (*(theCharP + 1));
      theFooP [2] = '\0';
      
      /* if the VR is unknown assume the group 2 is using implicit VR */
      if (!(strcmp (theFooP, "AE") == 0 || strcmp (theFooP, "AS") == 0 || strcmp (theFooP, "AT") == 0 ||
            strcmp (theFooP, "CS") == 0 || strcmp (theFooP, "DA") == 0 || strcmp (theFooP, "DS") == 0 ||
            strcmp (theFooP, "DT") == 0 || strcmp (theFooP, "FL") == 0 || strcmp (theFooP, "FD") == 0 ||
            strcmp (theFooP, "IS") == 0 || strcmp (theFooP, "LO") == 0 || strcmp (theFooP, "LT") == 0 ||
            strcmp (theFooP, "OW") == 0 || strcmp (theFooP, "PN") == 0 || strcmp (theFooP, "SH") == 0 ||
            strcmp (theFooP, "SL") == 0 || strcmp (theFooP, "SQ") == 0 || strcmp (theFooP, "SS") == 0 ||
            strcmp (theFooP, "ST") == 0 || strcmp (theFooP, "TM") == 0 || strcmp (theFooP, "UI") == 0 || 
            strcmp (theFooP, "UL") == 0 || strcmp (theFooP, "UN") == 0 || strcmp (theFooP, "US") == 0 ||
            strcmp (theFooP, "UT") == 0 || strcmp (theFooP, "OB") == 0))
        gArrTransfSyntax [inFileNb] = LITTLE_ENDIAN_IMPL;
      
      /* if there are OB values in group 2 it is a recent version of the toolkit ( >  3.3) */
      /* the theIsOld variable will be used later in the code of this routine */
      if (strcmp (theFooP, "OB") == 0) theIsOld = FALSE;
    } /* if ...group 2 */
    
    
    
    /* test wether the transfert syntax is the little-endian explicit VR one */
    if (gArrTransfSyntax [inFileNb] == LITTLE_ENDIAN_EXPL)
    {
      /* extract the VR */
      theFooP     = (char *) &theFoo [0];
      theFooP [0] = (char)   *theCharP;
      theFooP [1] = (char) (*(theCharP + 1));
      theFooP [2] = '\0';
      /* updates the current position in the read buffer */
      *ioBufPosP += 2L;
      /* points to the right place in the buffer */
      theCharP   += 2;
      
      /* extract the element length depending on the extracted VR */
      if (strcmp (theFooP, "OB") == 0 || 
      	  strcmp (theFooP, "OW") == 0 ||
      	  strcmp (theFooP, "SQ") == 0 ||
      	  strcmp (theFooP, "UN") == 0 ||
      	  strcmp (theFooP, "UT") == 0)
      {
        /* updates the current position in the read buffer by jumping over the 2 bytes set to 0 */
        *ioBufPosP += 2L;
        /* points to the right place in the buffer */
        theCharP   += 2;
        
        /* extract the element length according to the little-endian explicit VR syntax */
        theTmpULong      = (PapyULong) (*(theCharP + 3));
        theTmpULong      = theTmpULong << 24;
        theULong	 = theTmpULong;
        theTmpULong      = (PapyULong) (*(theCharP + 2));
        theTmpULong      = theTmpULong << 16;
        theULong	|= theTmpULong;
        theTmpULong      = (PapyULong) (*(theCharP + 1));
        theTmpULong      = theTmpULong << 8;
        theULong	|= theTmpULong;
        theTmpULong      = (PapyULong) *theCharP;
        theULong        |= theTmpULong;
        
        theElemLength    = theULong;
        
        /* updates the current position in the read buffer */
        *ioBufPosP += 4L;
      } /* if ...VR = OB, OW or SQ */
      else
      {
        /* extract the element length according to the little-endian explicit VR syntax */
        theElemLengthGr2  = (PapyUShort) (*(theCharP + 1));
        theElemLengthGr2  = theElemLengthGr2 << 8;
        theElemLengthGr2 |= (PapyUShort) *theCharP;
      
        theElemLength     = (PapyULong) theElemLengthGr2;
        
        /* updates the current position in the read buffer */
        *ioBufPosP += 2L;
      } /* else ...other VRs */
            
    } /* if ...transfert syntax is little_endian explicit VR */
    /* little_endian implicit VR */
    else
    {
      /* extract the element length according to the little-endian implicit VR syntax */
      theTmpULong      = (PapyULong) (*(theCharP + 3));
      theTmpULong      = theTmpULong << 24;
      theULong	       = theTmpULong;
      theTmpULong      = (PapyULong) (*(theCharP + 2));
      theTmpULong      = theTmpULong << 16;
      theULong	      |= theTmpULong;
      theTmpULong      = (PapyULong) (*(theCharP + 1));
      theTmpULong      = theTmpULong << 8;
      theULong	      |= theTmpULong;
      theTmpULong      = (PapyULong) *theCharP;
      theULong 	      |= theTmpULong;
      theElemLength    = theULong;
    
      /* updates the current position in the read buffer */
      *ioBufPosP += 4L;
    } /* else ...little_endian implicit VR */
    

    /* it could be an undefined length, i.e. VR = SQ or VR = UN */
    if (theElemLength == 0xFFFFFFFF)
    {
      /* for futur move of the buffer pointer */
      /* a Sequence Delimiter Item (8 bytes: (FFFE,E0DD) et 00000000) will be added */
      /* at the end of the data element value ...*/
      theIsUndefSeqLen = TRUE;
      
      theElemLength = 0L;
      if (!(theGrNb == 0x7FE0 && theElemNb == 0x0010))
      {
        /* get the current file position */
        theErr = Papy3FTell (gPapyFile [inFileNb], &theCurrFilePos);
        /* position the file pointer to point to the item */
        theErr = Papy3FSeek (gPapyFile [inFileNb], SEEK_SET, theInitialFilePos + (PapyLong) (*ioBufPosP));
        
        if ((theErr = ComputeUndefinedSequenceLength3 (inFileNb, &theElemLength)) < 0)
          RETURN (theErr);
      
        /* theElemLength should not contain the Sequence Delimitation Item
           which exist only if we have a Data Element of undefined length */
        if (theIsUndefSeqLen) theElemLength -= 8L;

        /* reset the file pointer to its previous position */
        theErr = Papy3FSeek (gPapyFile [inFileNb], SEEK_SET, theCurrFilePos);
      
      } /* if ...not image pixel */
      else /* just decide it is the last readable thing */
        theElemLength = inBytesToRead - (*ioBufPosP - theInitialBufPos);
    
    } /* if ...undefined length */
    
    /* odd element length are forbidden */
    if (theElemLength % 2 != 0)
      RETURN (papLengthIsNotEven);

    /* it is a shadow group, so we are looking dynamically for our element range */
    if (theShadow && theElemNb >= 0x0010 && theElemNb <= 0x00FF)
    {
      theCreator = Papy3CheckValidOwnerId (inFileNb, ioBuffP, ioBufPosP, theElemNb, theElemLength, ioGroupP);
      
      /* look for the position in the enum of the group */
      if (theCreator)
      {
        while (theStructPos <= theMaxElem && theElemNb != theArrElemP [theStructPos].element)
          theStructPos++;
        theEnabledShadow [theArrElemP [theStructPos].element] = TRUE;
      } /* if */
    	    
    } /* if */

    else
    {
      if (!theShadow ||
          (theShadow && theElemNb <  0x0010) ||
    	  (theShadow && theElemNb >= 0x1000 && theEnabledShadow [theElemNb >> 8]))
      {    
        /* search the element in the array */
        while (theStructPos <= theMaxElem && theElemNb != theArrElemP [theStructPos].element)
          theStructPos++;
      
        /* element number out of range */
        if (theStructPos >= theMaxElem)
        {
          /* it could be an unknown element (who knows with DICOM ?) */
          /* so we just skip the element */
	        *ioBufPosP += theElemLength;
        } /* if */
      
        else 
        {
          theArrElemP [theStructPos].length = theElemLength;
        
          /* there has been a change in the dictionary. */
          /* This helps the new version to read the old files ( < 3.3) */
          if (theIsOld && theGrNb == 0x0002 && (theElemNb == 0x0001 || theElemNb == 0x0102))
            theArrElemP [theStructPos].vr = USS;
 
      
          /* extract the value of the element from the buffer */      
          if (theElemLength > 0 && 
      	      !(theArrElemP [theStructPos].group   == 0x7FE0 && 
      	        theArrElemP [theStructPos].element == 0x0010))
          {
            /* extract the element depending on the value representation */
	          if ((theErr = PutBufferInElement3 (inFileNb, ioBuffP, theElemLength, 
				               &theArrElemP [theStructPos], ioBufPosP, theInitialFilePos)) < 0)
              RETURN (theErr);  
            
            /* if it was a sequence with an undefined length, move the buffer accordingly */
            if (theIsUndefSeqLen) 
              *ioBufPosP += 8L;
            
          } /* if ...theElemLength > 0 */
        } /* else ...element found */
	
      } /* if 				...we can read this element */
	
      else 		        /* we dont know how to read this element... */
        *ioBufPosP += theElemLength;	 			/* ...so we skip it */
    } /* else ...not creator of a private data element */
    
  } /* while ... loop on the elements */
  
  return (PapyShort) theEnumGrNb;

} /* endof PutBufferInGroup3 */


										
/********************************************************************************/
/*									 	*/
/*	Papy3GroupRead : read the group from the file in a buffer then fill_in 	*/
/*	the group structure (all the elements) from the buffer.			*/
/* 	return : the group number (in the enum_type) if able to fill it		*/
/*		 standard error message otherwise 				*/
/*									  	*/
/********************************************************************************/

PapyShort CALLINGCONV
Papy3GroupRead (PapyShort inFileNb, SElement **ioGroupP)
{

  unsigned char 	*theBuffP;
  PapyULong		theBytesToRead, theGrLength, theBufPos;
  PapyLong		theInitFilePos;
  PapyUShort	 	thePapyGrNb;
  PapyShort		theErr;
  int 			theEnumGrNb;
  enum ETransf_Syntax	thePrevSyntax; 
  
  
  /* get the position in the file for any eventual unknown seq length */
  theErr = Papy3FTell (gPapyFile [inFileNb], &theInitFilePos);
  
  theBufPos = 0L;
  /* read the buffer from the file */
  if (ReadGroup3 (inFileNb, &thePapyGrNb, &theBuffP, &theBytesToRead, &theGrLength) < 0)
  {
    efree3 ((void **) &theBuffP);
    RETURN (papReadGroup)
  } /* if */
  
  /* makes sure we keep the right syntax and set the default syntax instead */
  if (thePapyGrNb == 0x0002)
  {
    thePrevSyntax = gArrTransfSyntax [inFileNb];
    gArrTransfSyntax [inFileNb] = LITTLE_ENDIAN_EXPL;
  } /* if */
    
  theEnumGrNb = Papy3ToEnumGroup (thePapyGrNb);     /* gr_nb papyrus -> enum */
  if (theEnumGrNb < 0)				 /* unknown group number */
  {
    efree3 ((void **) &theBuffP);
    RETURN (papGroupNumber)
  } /* if */
 
  /* allocates the structure of the given group */

  *ioGroupP = Papy3GroupCreate (theEnumGrNb);
  
  /* if the group do not have the group length element, fill it ... */
  if (theGrLength != 0)
  {
    (*ioGroupP)->nb_val    = 1L;
    (*ioGroupP)->value     = (UValue_T *) emalloc3 ((PapyULong) sizeof (UValue_T));
    (*ioGroupP)->value->ul = theGrLength;
  } /* if ...undefined group length */
  
  /* extract the elements of the buffer and put them in the group structure */
  theBufPos = 0L;
  theErr = PutBufferInGroup3 (inFileNb, theBuffP, *ioGroupP, thePapyGrNb, theBytesToRead, &theBufPos, theInitFilePos);
  if (theErr < 0)
  {
    efree3 ((void **) &theBuffP);
    RETURN (theErr);
  } /* if */

  /* frees the read buffer */
  efree3 ((void **) &theBuffP);
  
  /* restore any previous transfert syntax */
  if (thePapyGrNb == 0x0002)
    gArrTransfSyntax [inFileNb] = thePrevSyntax;
  
  RETURN ((PapyShort) thePapyGrNb);
  
} /* endof Papy3GroupRead */
