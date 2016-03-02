#include "vtkETHZBioMedicalReader.h"

#include "vtkDataArraySelection.h"
#include "vtkObjectFactory.h"

const char * vtkETHZBioMedicalReader::strain_var_names[] = {
                   "e11",
                   "e22",
                   "e33",
                   "e12", 
                   "e23", 
                   "e31", 
                   "SED", 
                   "EFF"
                   };
const char * vtkETHZBioMedicalReader::stress_var_names[] = { 
                   "s11", 
                   "s22", 
                   "s33", 
                   "s12", 
                   "s23", 
                   "s31",
                   "SVM"
                   };

vtkStandardNewMacro(vtkETHZBioMedicalReader);

vtkETHZBioMedicalReader::vtkETHZBioMedicalReader()
{
  this->FileName = NULL; 
  this->DebugOff();
  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(1);
  this->PointDataArraySelection = vtkDataArraySelection::New();
  this->CellDataArraySelection = vtkDataArraySelection::New();
}

vtkETHZBioMedicalReader::~vtkETHZBioMedicalReader()
{
  vtkDebugMacro(<< "cleaning up inside destructor");
  if (this->FileName)
    delete [] this->FileName;
  this->PointDataArraySelection->Delete();
  this->CellDataArraySelection->Delete();
}

int vtkETHZBioMedicalReader::RequestInformation(vtkInformation* request,
                         vtkInformationVector** inputVector,
                         vtkInformationVector* outputVector)
{
  return 1;
}

int vtkETHZBioMedicalReader::RequestData(
                vtkInformation* vtkNotUsed(request),
                vtkInformationVector** vtkNotUsed(inputVector),
                vtkInformationVector* outputVector)
{
  return 1;
}

void vtkETHZBioMedicalReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "Reader:\n";
  os << indent << indent << "Filename      : " << this->FileName << endl;
}

void vtkETHZBioMedicalReader::EnablePointArray(const char* name)
{
  this->SetPointArrayStatus(name, 1);
}

void vtkETHZBioMedicalReader::DisablePointArray(const char* name)
{
  this->SetPointArrayStatus(name, 0);
}

void vtkETHZBioMedicalReader::EnableAllPointArrays()
{
  this->PointDataArraySelection->EnableAllArrays();
}

void vtkETHZBioMedicalReader::DisableAllPointArrays()
{
  this->PointDataArraySelection->DisableAllArrays();
}

const char* vtkETHZBioMedicalReader::GetPointArrayName(int index)
{
  return this->PointDataArraySelection->GetArrayName(index);
}

int vtkETHZBioMedicalReader::GetPointArrayStatus(const char* name)
{
  return this->PointDataArraySelection->ArrayIsEnabled(name);
}

void vtkETHZBioMedicalReader::SetPointArrayStatus(const char* name, int status)
{
  if(status)
    {
    this->PointDataArraySelection->EnableArray(name);
    }
  else
    {
    this->PointDataArraySelection->DisableArray(name);
    }
}

int vtkETHZBioMedicalReader::GetNumberOfPointArrays()
{
  return this->PointDataArraySelection->GetNumberOfArrays();
}

void vtkETHZBioMedicalReader::EnableCellArray(const char* name)
{
  this->SetCellArrayStatus(name, 1);
}

void vtkETHZBioMedicalReader::DisableCellArray(const char* name)
{
  this->SetCellArrayStatus(name, 0);
}

void vtkETHZBioMedicalReader::EnableAllCellArrays()
{
  this->CellDataArraySelection->EnableAllArrays();
}

void vtkETHZBioMedicalReader::DisableAllCellArrays()
{
  this->CellDataArraySelection->DisableAllArrays();
}

const char* vtkETHZBioMedicalReader::GetCellArrayName(int index)
{
  return this->CellDataArraySelection->GetArrayName(index);
}

int vtkETHZBioMedicalReader::GetCellArrayStatus(const char* name)
{
  return this->CellDataArraySelection->ArrayIsEnabled(name);
}

void vtkETHZBioMedicalReader::SetCellArrayStatus(const char* name, int status)
{
  if(status)
    {
    this->CellDataArraySelection->EnableArray(name);
    }
  else
    {
    this->CellDataArraySelection->DisableArray(name);
    }
}

int vtkETHZBioMedicalReader::GetNumberOfCellArrays()
{
  return this->CellDataArraySelection->GetNumberOfArrays();
}
