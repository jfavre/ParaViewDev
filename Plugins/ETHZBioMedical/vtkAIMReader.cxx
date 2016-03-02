/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkAIMReader.cxx
  Language:  C++
  Date:      Tue May 22 09:53:57 MET DST 2001
  Version:   Rev 3, vtk4.0

  Steven K. Boyd
  boyd@biomed.ee.ethz.ch

  Extensive changes by Jean M Favre at CSCS since the earlier version
  Last updated: October 31, 2013 to work under ParaView 4.0
=========================================================================*/

#include "vtkAIMReader.h"

#include "vtkImageData.h"
#include "vtkObjectFactory.h"
#include "vtkShortArray.h"
#include "vtkCharArray.h"
#include "vtkFloatArray.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkStreamingDemandDrivenPipeline.h"

vtkStandardNewMacro(vtkAIMReader);

//-----------------------------------------------------------------------
// Constructor for a vtkAIMReader.
vtkAIMReader::vtkAIMReader()
{
  this->FileName = NULL;
  this->Error = 0;
  this->El_size_mm[0] = 0.01;
  this->El_size_mm[1] = 0.01;
  this->El_size_mm[2] = 0.01;
  this->Offset[0] = 0;
  this->Offset[1] = 0;
  this->Offset[2] = 0;
  this->Position[0] = 0;
  this->Position[1] = 0;
  this->Position[2] = 0;

  D3InitAny(this->image,D1Tundef);
  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(1);
}

//-----------------------------------------------------------------------
vtkAIMReader::~vtkAIMReader()
{
  if (this->FileName)
    {
    delete [] this->FileName;
    }
  D3FreeAnyImage(this->image);
}

//-----------------------------------------------------------------------
// Retrieve the pointer to the processing log.
char *vtkAIMReader::GetProcessingLog()
{
  if ( this->image.proc_log != 0 )
    {
    return this->image.proc_log;
    }
  else
    {
    vtkErrorMacro(<< "No processing log available.\n"
                  << "Make sure vtkAIMReader is Updated() before calling this option.");
    return NULL;
    }
}

//-----------------------------------------------------------------------
// Print AIM information.
void vtkAIMReader::PrintLogInfo()
{
  if ( this->image.proc_log != 0 )
    {
    D3PrintLogInfoAny(this->image, FileName);
    }
  else
    {
    vtkErrorMacro(<< "No image data available.\n"
                  << "Make sure vtkAIMReader is Updated() before calling this option.");
    }

  return;
}

//-----------------------------------------------------------------------
// Print AIM log information.
void vtkAIMReader::PrintProcessingLog()
{
  if ( this->image.proc_log != 0 )
    {
    DnProcLogPrint(this->image);
    }
  else
    {
    vtkErrorMacro(<<"No processing log information available.");
    }

  return;
}

//-----------------------------------------------------------------------
// Determines important information about the input AIM file

int vtkAIMReader::RequestInformation(vtkInformation* request,
    vtkInformationVector** vtkNotUsed(inputVector),
    vtkInformationVector* outputVector)
{
  vtkDebugMacro(<<"\n  Running ExecuteInformation.");
  vtkInformation* outInfo = outputVector->GetInformationObject(0);
  int       status;
  D3int     dim, off, pos;
  D3float   el_size_mm;
  double    spacing[3];

  vtkImageData *output = vtkImageData::SafeDownCast( outInfo->Get(vtkDataObject::DATA_OBJECT()));
      
  if (!this->FileName)
    {
      vtkErrorMacro(<<"A FileName must be specified.");
      return 0;
    }

  AimpackReadImageInfo(&this->image, this->FileName, &status);

  if (!status)
  {
    vtkErrorMacro(<< "File " << this->FileName << " not found.");
    return 0;
  }

  dim = this->image.dim;
  off = this->image.off;
  pos = this->image.pos;
  el_size_mm = this->image.el_size_mm;

  this->El_size_mm[0] = (double)el_size_mm.x;
  this->El_size_mm[1] = (double)el_size_mm.y;
  this->El_size_mm[2] = (double)el_size_mm.z;
  this->Offset[0] = off.x;
  this->Offset[1] = off.y;
  this->Offset[2] = off.z;
  this->Dimension[0] = dim.x;
  this->Dimension[1] = dim.y;
  this->Dimension[2] = dim.z;
  this->Position[0] = pos.x;
  this->Position[1] = pos.y;
  this->Position[2] = pos.z;

  // el_size_mm is originally a float value.  We need to round-off this number.
  // It is assumed that we will never have more than eight significant figures
  // past the decimal point.
  if ( (log10(el_size_mm.x) < -8) ||
       (log10(el_size_mm.y) < -8) ||
       (log10(el_size_mm.z) < -8) )
    vtkWarningMacro(<<"\n  Element precision is less than 1e-8."
                    <<"\n  Errors in el_size_mm conversion may result.");
  spacing[0] = floor(el_size_mm.x*1e8)/1e8;
  spacing[1] = floor(el_size_mm.y*1e8)/1e8;
  spacing[2] = floor(el_size_mm.z*1e8)/1e8;

  // AIM 'pos' data is in unit of voxels, and converted to vtkImageData->Origin
  // which are in real world coordinates.  Therefore, voxels are converted to
  // world coordinates.

  double origin[3]={0,0,0};
  outInfo->Set(vtkDataObject::ORIGIN(), origin, 3);
  outInfo->Set(vtkDataObject::SPACING(), spacing, 3);

  int dEx[6];
  dEx[0] = 0; dEx[1] = dim.x; dEx[2] = 0; dEx[3] = dim.y; dEx[4] = 0; dEx[5] = dim.z;
  outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), dEx, 6);

  return 1;
}

//-----------------------------------------------------------------------
// Reads an AIM file and creates a vtkStructuredPoints dataset.
int vtkAIMReader::RequestData(
  vtkInformation* request,
  vtkInformationVector** vtkNotUsed(inputVector),
  vtkInformationVector* outputVector)
{
  vtkDebugMacro(<<"\n  Running Execute.");
  vtkInformation* outInfo = outputVector->GetInformationObject(0);
  D3int dim;
  int	status, count = 0;

  vtkCharArray *carray;
  vtkShortArray *sarray;
  vtkFloatArray *farray;

  this->Error = 1;
  this->UpdateProgress(.10);
  vtkImageData *output = vtkImageData::GetData(outputVector);

  output->SetExtent(
    outInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT()));

  AimpackReadImageUncompress(&this->image, FileName, &status);
  this->UpdateProgress(.70);

  if (!status) {
    vtkErrorMacro(<< "File could not be read\n  " << FileName);
  }
  switch( this->image.type )
  {
    case D1Tchar:
    case D1TbinCmp:
    case D3Tbit8:
    case D1TcharCmp:
      output->SetScalarType(VTK_CHAR, outInfo);
      output->AllocateScalars(VTK_CHAR, 1);
      break;
    case D1Tshort:
      output->SetScalarType(VTK_SHORT, outInfo);
      output->AllocateScalars(VTK_SHORT, 1);
      cerr << "output->SetScalarType(VTK_SHORT)\n";
      break;

    case D1Tfloat:
    case D1Ti3efloat:
      output->SetScalarType(VTK_FLOAT, outInfo);
      output->AllocateScalars(VTK_FLOAT, 1);
      cerr << "output->SetScalarType(VTK_FLOAT)\n";
      break;
    default:
    vtkErrorMacro(<< "Unknown AIM data type: " << (int)(this->image.type) );
    break;
  }

  dim.x = output->GetDimensions()[0] - 1;
  dim.y = output->GetDimensions()[1] - 1;
  dim.z = output->GetDimensions()[2] - 1;

  // Reads all of the data in the AIM (even any data in the "off" region)

  switch( output->GetScalarType() )
  {
    case VTK_CHAR:
      cerr << "creating new char array\n";
      carray = vtkCharArray::New();
      carray->SetNumberOfComponents(1);
      carray->SetName("aim_data");
      carray->SetArray( (char *) (this->image.dat), dim.x * dim.y * dim.z, 1);
      output->GetCellData()->SetScalars( carray );
      carray->Delete();
      break;

    case VTK_SHORT:
      sarray = vtkShortArray::New();
      sarray->SetNumberOfComponents(1);
      sarray->SetName("aim_data");
      sarray->SetArray( (short *) (this->image.dat), dim.x * dim.y * dim.z, 1);
      output->GetCellData()->SetScalars( sarray );
      sarray->Delete();
      break;

    case VTK_FLOAT:
      farray = vtkFloatArray::New();
      farray->SetNumberOfComponents(1);
      farray->SetName("aim_data");
      farray->SetArray( (float *) (this->image.dat), dim.x * dim.y * dim.z, 1);
      output->GetCellData()->SetScalars( farray );
      farray->Delete();

      break;

    default:
    vtkErrorMacro(<< "Unknown AIM data type: " << (int)(this->image.type) );
    break;    
  }
  this->UpdateProgress(.90);
  //this is extraneous. Probably due to the default allocation scheme. Remove it now.
  output->GetPointData()->RemoveArray("ImageScalars");
  vtkDebugMacro(<<"\n  " << (dim.x * dim.y * dim.z) << " cells read.");

  this->Error = 0;

  return 1;
}

//-----------------------------------------------------------------------
void vtkAIMReader::PrintSelf(ostream& os, vtkIndent indent)
{
  vtkImageAlgorithm::PrintSelf(os,indent);

  os << indent << "Error: " << this->Error << "\n";
  os << indent << "File Name: " 
     << (this->FileName ? this->FileName : "(none)") << "\n";
}



// THIS IS THE OLD WAY THAT AIMS WERE READ

#if 0
  data->GetContinuousIncrements(outExt, outIncX, outIncY, outIncZ);

    
  case VTK_CHAR:

      outPtrC = (char *) data->GetScalarPointer(outExt[0], outExt[2], outExt[4]);

      for (idxZ = 0; idxZ <= maxZ; idxZ++)
      {
        for (idxY = 0; !this->AbortExecute && idxY <= maxY; idxY++)
        {
          for (idxX = 0; idxX <= maxX; idxX++) 
          {
            *outPtrC = ((D1char *)this->image.dat)[count];
            outPtrC++;
            count++;
          }
          outPtrC += outIncY;
        }
        outPtrC += outIncZ;
      }
      break;

    case VTK_SHORT:

      outPtrS = (short *) data->GetScalarPointer(outExt[0], outExt[2], outExt[4]);

      for (idxZ = 0; idxZ <= maxZ; idxZ++)
      {
        for (idxY = 0; !this->AbortExecute && idxY <= maxY; idxY++)
        {
          for (idxX = 0; idxX <= maxX; idxX++) 
          {
            *outPtrS = ((D1short *)this->image.dat)[count];
            outPtrS++;
            count++;
          }
          outPtrS += outIncY;
        }
        outPtrS += outIncZ;
      }
      break;
#endif
