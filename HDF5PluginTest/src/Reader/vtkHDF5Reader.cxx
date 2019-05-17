
#include "vtkHDF5Reader.h"

#include "vtkErrorCode.h"
#include "vtkFieldData.h"
#include "vtkFloatArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkUnstructuredGrid.h"

#include <sys/time.h>
#include <algorithm>
using namespace std;

#include <vtk_hdf5.h>

vtkStandardNewMacro(vtkHDF5Reader);

vtkHDF5Reader::vtkHDF5Reader()
{
  this->FileName = nullptr; 
  this->DebugOff();
  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(1);
}
 
vtkHDF5Reader::~vtkHDF5Reader()
{
  vtkDebugMacro(<< "cleaning up inside destructor");
  if (this->FileName)
    delete [] this->FileName;
}

int vtkHDF5Reader::RequestInformation(
                         vtkInformation *vtkNotUsed(request),
                         vtkInformationVector **vtkNotUsed(inputVector),
                         vtkInformationVector* outputVector)
{
  hid_t file_id = H5Fopen(this->FileName, H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t g_id, root_id, mesh_id, dataset_id;
  hid_t filespace;
  hsize_t dimsf[2];

  vtkInformation* outInfo = outputVector->GetInformationObject(0);

  root_id = H5Gopen(file_id, "/", H5P_DEFAULT);

  mesh_id = H5Gopen(root_id, "Unnamed", H5P_DEFAULT);
  dataset_id = H5Dopen(mesh_id, "XYZ", H5P_DEFAULT);
  filespace = H5Dget_space(dataset_id);
  H5Sget_simple_extent_dims(filespace, dimsf, NULL);
  this->NbNodes = dimsf[0];
  cerr << __LINE__ << ": RequestInformation(found " << this->NbNodes << " with " << dimsf[1] << " coordinates values\n";
  H5Sclose(filespace);
  H5Dclose(dataset_id);

  return 1;
}

int vtkHDF5Reader::RequestData(
                vtkInformation* vtkNotUsed(request),
                vtkInformationVector** vtkNotUsed(inputVector),
                vtkInformationVector* outputVector)
{
  hid_t root_id, mesh_id, dataset_id, coords_id, solutions_id=0, dataset_id2, parameters_id, ids_id;
  hid_t filespace;
  hsize_t dimsf[2];
  herr_t   status;
  bool Materials_Found = false;
  struct timeval tv0, tv1, res;

  vtkDebugMacro( << "RequestData(BEGIN)");
  vtkInformation* outInfo = outputVector->GetInformationObject(0);
  vtkDataObject* doOutput = outInfo->Get(vtkDataObject::DATA_OBJECT());
  vtkUnstructuredGrid* mb = vtkUnstructuredGrid::SafeDownCast(doOutput);
  if (!mb)
  {
    return 0;
  }

  this->UpdateProgress(1.0);

  return 1;
}

