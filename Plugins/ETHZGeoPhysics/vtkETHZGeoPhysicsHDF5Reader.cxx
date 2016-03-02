#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkDataArraySelection.h"
#include "vtkDataSetAttributes.h"
#include "vtkDoubleArray.h"
#include "vtkErrorCode.h"
#include "vtkExtentTranslator.h"
#include "vtkFieldData.h"
#include "vtkFloatArray.h"
#include "vtkIdTypeArray.h"
#include "vtkImageData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkIntArray.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkStreamingDemandDrivenPipeline.h"


#include "vtkETHZGeoPhysicsHDF5Reader.h"

#define DATA_LENGTH H5T_NATIVE_FLOAT

//vtkCxxRevisionMacro(vtkETHZGeoPhysicsHDF5Reader, "$Revision: 1.1 $");
vtkStandardNewMacro(vtkETHZGeoPhysicsHDF5Reader);


void vtkETHZGeoPhysicsHDF5Reader::SetWholeExtent(int xMin, int xMax,
                                         int yMin, int yMax,
                                         int zMin, int zMax)
{
  int modified = 0;
cerr << "SetWholeExtent: " <<
xMin << ", " << xMax << ", " <<
yMin << ", " << yMax << ", " <<
zMin << ", " <<  zMax << endl;
  if (this->WholeExtent[0] != xMin)
    {
    modified = 1;
    this->WholeExtent[0] = xMin ;
    }
  if (this->WholeExtent[1] != xMax)
    {
    modified = 1;
    this->WholeExtent[1] = xMax ;
    }
  if (this->WholeExtent[2] != yMin)
    {
    modified = 1;
    this->WholeExtent[2] = yMin ;
    }
  if (this->WholeExtent[3] != yMax)
    {
    modified = 1;
    this->WholeExtent[3] = yMax ;
    }
  if (this->WholeExtent[4] != zMin)
    {
    modified = 1;
    this->WholeExtent[4] = zMin ;
    }
  if (this->WholeExtent[5] != zMax)
    {
    modified = 1;
    this->WholeExtent[5] = zMax ;
    }
  if (modified)
    {
    this->Modified();
    }
}


vtkETHZGeoPhysicsHDF5Reader::vtkETHZGeoPhysicsHDF5Reader()
{
  this->FileName            = NULL; 
  this->DebugOff();
  this->SetNumberOfInputPorts(0);
  //this->SetNumberOfOutputPorts(1);
  this->PointDataArraySelection = vtkDataArraySelection::New();
  //this->CellDataArraySelection = vtkDataArraySelection::New();
  this->ShuffleXZOn();
}

vtkETHZGeoPhysicsHDF5Reader::~vtkETHZGeoPhysicsHDF5Reader()
{
  vtkDebugMacro(<< "cleaning up inside destructor");
  if (this->FileName)
    delete [] this->FileName;
  this->PointDataArraySelection->Delete();
  //this->CellDataArraySelection->Delete();
}

int vtkETHZGeoPhysicsHDF5Reader::RequestInformation(vtkInformation* request,
                         vtkInformationVector** inputVector,
                         vtkInformationVector* outputVector)
{
  hid_t err_id, filespace, dataset_id, file_id, root_id, geo_id;
  int my_extents[6];
  hsize_t dimsf[3];
  vtkInformation* outInfo = outputVector->GetInformationObject(0);

  //outInfo->Set(vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(), -1);

  file_id = H5Fopen(this->FileName, H5F_ACC_RDONLY, H5P_DEFAULT);
  root_id = H5Gopen(file_id, "/", H5P_DEFAULT);
  // before opening, we need to inquire if object really exists and turn off errors.

  herr_t error;
  H5E_auto_t func;
  void *client_data;
  H5Eget_auto2(H5E_DEFAULT, &func, &client_data);
  H5Eset_auto2(H5E_DEFAULT, NULL, NULL);
  error = H5Gget_objinfo(root_id, "Entropy values - UNITS UNKNOWN", 0, NULL);
  H5Eset_auto2(H5E_DEFAULT, func, client_data);

  if(error < 0)
    {
    //cerr << "new HDF5 format\n";
    geo_id = H5Gopen(root_id, "/GeoPhysics", H5P_DEFAULT);
    dataset_id = H5Dopen(geo_id, "Entropy", H5P_DEFAULT);
    this->ShuffleXZOff();
    }
  else
    {
    dataset_id = H5Dopen(root_id, "Entropy values - UNITS UNKNOWN", H5P_DEFAULT);
    this->ShuffleXZOn();
    //cerr << "old HDF5 format\n";
    }
  filespace = H5Dget_space(dataset_id);
  int rank = H5Sget_simple_extent_ndims(filespace);
  H5Sget_simple_extent_dims(filespace, dimsf, NULL);
  H5Sclose(filespace);
  H5Dclose(dataset_id);

  this->PointDataArraySelection->AddArray("Entropy");
  this->PointDataArraySelection->AddArray("Pressure");
  this->PointDataArraySelection->AddArray("Velocity");

//cerr << dimsf[0] << " " << dimsf[1] << " " << dimsf[2] << endl;
  my_extents[0] = 0;
  my_extents[1] = dimsf[0]-1;
  my_extents[2] = 0;
  my_extents[3] = dimsf[1]-1;
  my_extents[4] = 0;
  my_extents[5] = dimsf[2]-1;
  outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), my_extents, 6);
/*
  cerr << "Set WHOLE_EXTENT  extents: " << my_extents[0] << ", " <<
                                                    my_extents[1] << ", " <<
                                                    my_extents[2] << ", " <<
                                                    my_extents[3] << ", " <<
                                                    my_extents[4] << ", " <<
                                                    my_extents[5] << endl;
*/
  double spacing[3];
  double origin[3];

  // set the spacing
  spacing[0] = spacing[1] = spacing[2] = 1;
  outInfo->Set(vtkDataObject::SPACING(), spacing, 3);

  // set the origin.
  origin[0] = -my_extents[1]*0.5;
  origin[1] = -my_extents[3]*0.5;
  origin[2] = -my_extents[5]*0.5;
  outInfo->Set(vtkDataObject::ORIGIN(),  origin, 3);

  if(!this->GetShuffleXZ())
    H5Gclose(geo_id);
  H5Gclose(root_id);
  H5Fclose(file_id);
  vtkDataObject::SetPointDataActiveScalarInfo(outInfo, VTK_FLOAT, 1);
// every pvserver process executes requestInformation but there is no notion of total # of processes
  return 1;
}

int numberOfTuples ;
hid_t    file_id2, root_id2, group_id2;

#define PIECE 0
static void H5LTmake_dataset_float(hid_t group_id, const char *name, int rank, hsize_t dims[], const float *data)
{
   hid_t dataspace_id = H5Screate_simple(rank, dims, NULL);
   hid_t dataset_id = H5Dcreate(group_id, name, H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
   H5Dclose(dataset_id);
   H5Sclose(dataspace_id);
}

void vtkETHZGeoPhysicsHDF5Reader::ExecuteDataWithInformation(vtkDataObject *output, vtkInformation* outInfo)
{
  hid_t    mesh_id, dataset_id, dataset_id2, parameters_id;
  hid_t    aid1, attr1, file_id, root_id, dataspace_id, geo_id;
  herr_t   status;

  vtkDebugMacro( << "RequestData(BEGIN)");
  this->UpdateProgress(0);
  int done=0;
  int tasks = 0;
  if(this->PointDataArraySelection->ArrayIsEnabled("Entropy"))
    tasks++;
  if(this->PointDataArraySelection->ArrayIsEnabled("Pressure"))
    tasks++;
  if(this->PointDataArraySelection->ArrayIsEnabled("Velocity"))
    tasks+= 3;

  int updatePiece, updateNumPieces;

  vtkImageData *data = this->AllocateOutputData(output, outInfo);

  int *outExt = data->GetExtent();
  cerr << "ExecuteData:data GET_EXTENT: " <<
  outExt[0] << ", " << outExt[1] << ", " <<
  outExt[2] << ", " << outExt[3] << ", " <<
  outExt[4] << ", " <<  outExt[5] << endl;
  int extent[6];

  if (!this->FileName)
    {
    vtkErrorMacro(<< "error reading header specified!");
    return;
    }

  numberOfTuples = (outExt[1] - outExt[0] + 1) *
                   (outExt[3] - outExt[2] + 1) *
                   (outExt[5] - outExt[4] + 1);

  hsize_t     dimsm[4];             /* dataset dimensions */
  int dims[3];

  hsize_t stride[3]; /* Stride of hyperslab */
  hsize_t block[3];  /* Block sizes */

  dimsm[0] = outExt[1] - outExt[0] + 1;
  dimsm[1] = outExt[3] - outExt[2] + 1;
  dimsm[2] = outExt[5] - outExt[4] + 1;

  file_id = H5Fopen(this->FileName, H5F_ACC_RDONLY, H5P_DEFAULT);
  root_id = H5Gopen(file_id, "/", H5P_DEFAULT);
  if(!this->GetShuffleXZ())
    {
    geo_id = H5Gopen(root_id, "/GeoPhysics", H5P_DEFAULT);
    }
  if(this->WriteHDF5FileToDisk)
    {
    file_id2 = H5Fcreate("/tmp/geophysics.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    group_id2 = H5Gcreate(file_id2, "/GeoPhysics", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }

  if(this->PointDataArraySelection->ArrayIsEnabled("Entropy"))
    {
    if(this->GetShuffleXZ())
      dataset_id = H5Dopen(root_id, "Entropy values - UNITS UNKNOWN", H5P_DEFAULT);
    else
      {
      geo_id = H5Gopen(root_id, "/GeoPhysics", H5P_DEFAULT);
      dataset_id = H5Dopen(geo_id, "Entropy", H5P_DEFAULT);
      H5Gclose(geo_id);
      }
    this->ReadScalar(data, dataset_id, outExt, "Entropy");
    data->GetPointData()->SetActiveAttribute("Entropy", vtkDataSetAttributes::SCALARS);
    done++;
    this->UpdateProgress((double)done/tasks);
    }
  if(this->PointDataArraySelection->ArrayIsEnabled("Pressure"))
    {
    if(this->GetShuffleXZ())
      dataset_id = H5Dopen(root_id, "Pressure       - UNITS UNKNOWN", H5P_DEFAULT);
    else
      {
      geo_id = H5Gopen(root_id, "/GeoPhysics", H5P_DEFAULT);
      dataset_id = H5Dopen(geo_id, "Pressure", H5P_DEFAULT);
      H5Gclose(geo_id);
      }
    this->ReadScalar(data, dataset_id, outExt, "Pressure");
    done++;
    this->UpdateProgress((double)done/tasks);
    }
  if(this->PointDataArraySelection->ArrayIsEnabled("Velocity"))
    {
  hsize_t     dimsf[3];
#if H5_VERS_RELEASE == 2
  hssize_t     start[3];         /* hyperslab offset in the file */
#else
  hsize_t     start[3];         /* hyperslab offset in the file */
#endif
  hsize_t count[3];
  hid_t       filespace, memspace;
    vtkFloatArray *velocity = vtkFloatArray::New();
    vtkFloatArray *tdata = vtkFloatArray::New();

    if(!this->GetShuffleXZ())
      {
      geo_id = H5Gopen(root_id, "/GeoPhysics", H5P_DEFAULT);
      dataset_id = H5Dopen(geo_id, "Velocity", H5P_DEFAULT);
      this->ReadVector(data, dataset_id, outExt, "Velocity");
      //H5Dread(dataset_id, DATA_LENGTH, H5S_ALL, H5S_ALL, H5P_DEFAULT, velocity->GetPointer(0));
    done+=3;
    this->UpdateProgress((double)done/tasks);
      }
    else
      {
      velocity->SetNumberOfComponents(3);
      velocity->SetNumberOfTuples(numberOfTuples);
      velocity->SetName("Velocity");
      float *vptr = velocity->GetPointer(0);
      cerr << "allocating " << numberOfTuples * 3 * sizeof(float) << " bytes for Velocity\n";
    if(tdata->GetNumberOfTuples() == 0)
      {
      tdata->SetNumberOfComponents(1);
      tdata->SetNumberOfTuples(numberOfTuples);
      cerr << "allocating " << numberOfTuples * 1 * sizeof(float) << " bytes for tdata:Velocity\n";
      }
  dataset_id = H5Dopen(root_id, "X Velocity     - UNITS UNKNOWN", H5P_DEFAULT);
	  
  filespace = H5Dget_space(dataset_id);
  int rank = H5Sget_simple_extent_ndims(filespace);
  H5Sget_simple_extent_dims(filespace, dimsf, NULL);
      //Look at reverse order below!
  start[2] = outExt[0];
  start[1] = outExt[2];
  start[0] = outExt[4];

  count[2]  = outExt[1] - outExt[0] + 1;
  count[1]  = outExt[3] - outExt[2] + 1;
  count[0]  = outExt[5] - outExt[4] + 1;
  
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, NULL, count, NULL);
    
  count[0] = outExt[1] - outExt[0] + 1;
  count[1] = outExt[3] - outExt[2] + 1;
  count[2] = outExt[5] - outExt[4] + 1;

  memspace  = H5Screate_simple(3, count, NULL);
  start[0] = start[1] = start[2] = 0;

  H5Sselect_hyperslab(memspace, H5S_SELECT_SET, start, NULL, count, NULL);

    if(1)//this->GetShuffleXZ())
      {
      status = H5Dread(dataset_id, DATA_LENGTH, memspace, filespace,
                             H5P_DEFAULT, tdata->GetPointer(0));
      dims[0] = dimsf[0]; dims[1] = dimsf[1]; dims[2] = dimsf[2];
      this->ShuffleFromFORTRAN2C_Vector(vptr, tdata, count);
      cerr << "read and shuffled X Velocity\n";
      }
    else
      {
      start[0] = 0; block[0] = 1; stride[0] = 3; count[0]  = numberOfTuples;
      H5Sselect_hyperslab(memspace, H5S_SELECT_SET, start, stride, count, block);
      status = H5Dread(dataset_id, DATA_LENGTH, memspace, filespace,
                             H5P_DEFAULT, velocity->GetPointer(0));
      cerr << "\nread X Velocity\n";
      }
    done++;
    this->UpdateProgress((double)done/tasks);
    H5Dclose(dataset_id);

  dataset_id = H5Dopen(root_id, "Y Velocity     - UNITS UNKNOWN", H5P_DEFAULT);
  filespace = H5Dget_space(dataset_id);
  rank = H5Sget_simple_extent_ndims(filespace);
  H5Sget_simple_extent_dims(filespace, dimsf, NULL);
  
  // remember Look at reverse order below!
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, NULL, count, NULL);
  
    if(1)//this->GetShuffleXZ())
      {
      status = H5Dread(dataset_id, DATA_LENGTH, memspace, filespace,
                             H5P_DEFAULT, tdata->GetPointer(0));
      vptr++; // offset by one to stuff in the 1-st component (0-based)
      this->ShuffleFromFORTRAN2C_Vector(vptr, tdata, count);
      cerr << "read and shuffled Y Velocity\n";
      }
    else
      {
      start[0] = 1;
      H5Sselect_hyperslab(memspace, H5S_SELECT_SET, start, stride, count, block);
      status = H5Dread(dataset_id, DATA_LENGTH, memspace, filespace,
                             H5P_DEFAULT, velocity->GetPointer(0));
      cerr << "\nread Y Velocity\n";
      }
    done++;
    this->UpdateProgress((double)done/tasks);
    H5Dclose(dataset_id);
    
  dataset_id = H5Dopen(root_id, "Z Velocity     - UNITS UNKNOWN", H5P_DEFAULT);
  filespace = H5Dget_space(dataset_id);
  rank = H5Sget_simple_extent_ndims(filespace);
  H5Sget_simple_extent_dims(filespace, dimsf, NULL);

  // remember Look at reverse order below!
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, NULL, count, NULL);
    if(1)//this->GetShuffleXZ())
      {
      status = H5Dread(dataset_id, DATA_LENGTH, memspace, filespace,
                             H5P_DEFAULT, tdata->GetPointer(0));
			         
      vptr++; // offset by one to stuff in the 2-nd component (0-based) 
      this->ShuffleFromFORTRAN2C_Vector(vptr, tdata, count);
      cerr << "read and shuffled Z Velocity\n";
      }
    else
      {
      start[0] = 2;
      H5Sselect_hyperslab(memspace, H5S_SELECT_SET, start,
                                                  stride, count, block);
    
      status = H5Dread(dataset_id, DATA_LENGTH, memspace, filespace,
                             H5P_DEFAULT, velocity->GetPointer(0));
      cerr << "\nread Z Velocity\n";
      }
    data->GetPointData()->SetVectors(velocity);
    //velocity->Delete();
    data->GetPointData()->SetActiveAttribute("Velocity", vtkDataSetAttributes::VECTORS);
    //tdata->FastDelete();
    done++;
    this->UpdateProgress((double)done/tasks);
    H5Dclose(dataset_id);
    H5Sclose(memspace);
    }
    if(this->WriteHDF5FileToDisk)
      {
      dimsm[3] = 3; //set the 4th dimension
      H5LTmake_dataset_float(group_id2, "Velocity", 4, dimsm, velocity->GetPointer(0));
      }
    //data->GetPointData()->SetVectors(velocity);
    velocity->Delete();
    //data->GetPointData()->SetActiveAttribute("Velocity", vtkDataSetAttributes::VECTORS);
     tdata->Delete();
    }

  if(this->WriteHDF5FileToDisk)
    {
    H5Gclose(group_id2);
    H5Fclose(file_id2);
    }
  //if(!this->GetShuffleXZ())
    //H5Gclose(geo_id);

  vtkDebugMacro( << "RequestData(END)");
  //cerr << "vtkETHZGeoPhysicsHDF5Reader::RequestData(END)\n";
  this->UpdateProgress(1);
  return ;
}


void
vtkETHZGeoPhysicsHDF5Reader::ShuffleFromFORTRAN2C(vtkFloatArray *dest, vtkFloatArray *src,   hsize_t *dimsm)
{
  int t = dimsm[2] * dimsm[1];
  for(int z=0; z < dimsm[0]; z++)
    for(int y=0; y < dimsm[1]; y++)
      for(int x=0; x < dimsm[2]; x++)
        dest->SetValue(z*t+y*dimsm[2]+x, src->GetValue(x*t+y*dimsm[2]+z));
}

void
vtkETHZGeoPhysicsHDF5Reader::ReadScalar(vtkImageData *data, hid_t dataset_id, int *outExt, const char *name)
{
  hsize_t     dimsf[4];
#if H5_VERS_RELEASE == 2
  hssize_t     start[4];         /* hyperslab offset in the file */
#else
  hsize_t     start[4];         /* hyperslab offset in the file */
#endif
  hsize_t count[4];
  hid_t       filespace, memspace;

  vtkFloatArray *tdata = vtkFloatArray::New();
  vtkFloatArray *fa = vtkFloatArray::New();
  fa->SetNumberOfComponents(1);
  fa->SetNumberOfTuples(numberOfTuples);

  if(this->GetShuffleXZ())
    {
    tdata->SetNumberOfComponents(1);
    tdata->SetNumberOfTuples(numberOfTuples);
    }

  fa->SetName(name);
  //cerr << "Proc[" << "?" << "/" << "?" << "] allocating " << numberOfTuples  << " sizeof(float) * bytes for " << name << endl;

  filespace = H5Dget_space(dataset_id);
  int rank = H5Sget_simple_extent_ndims(filespace);
  H5Sget_simple_extent_dims(filespace, dimsf, NULL);

  //cerr << "filespace[] = "; for(int u=0; u < rank; u++) cerr << dimsf[u] << " "; cerr << endl;
//Look at reverse order below!
  start[2] = outExt[0];
  start[1] = outExt[2];
  start[0] = outExt[4];

  count[2]  = outExt[1] - outExt[0] + 1;
  count[1]  = outExt[3] - outExt[2] + 1;
  count[0]  = outExt[5] - outExt[4] + 1;

  //cerr << "File Start[] "; for(int u=0; u < rank; u++) cerr << start[u] << " "; cerr << endl;
  //cerr << "File Count[] "; for(int u=0; u < rank; u++) cerr << count[u] << " "; cerr <<endl;

  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, NULL, count, NULL);

  count[0] = outExt[1] - outExt[0] + 1;
  count[1] = outExt[3] - outExt[2] + 1;
  count[2] = outExt[5] - outExt[4] + 1;

  memspace  = H5Screate_simple(3, count, NULL);
  start[0] = start[1] = start[2] = 0;

  //cerr << "Mem Start[] "; for(int u=0; u < rank; u++) cerr << start[u] << " "; cerr << endl;
  //cerr << "Mem Count[] "; for(int u=0; u < rank; u++) cerr << count[u] << " "; cerr << endl;

  H5Sselect_hyperslab(memspace, H5S_SELECT_SET, start, NULL, count, NULL);

  if(this->GetShuffleXZ())
    {
    H5Dread(dataset_id, DATA_LENGTH, memspace, filespace, H5P_DEFAULT, tdata->GetPointer(0));

    this->ShuffleFromFORTRAN2C(fa, tdata, count);

    // tdata->Delete(); do not delete. We re-use the array for each var to load
    cerr << "read and shuffled " << name << endl;
    }
  else
    {
    H5Dread(dataset_id, DATA_LENGTH, memspace, filespace, H5P_DEFAULT, fa->GetPointer(0));
    cerr << "read " << name << endl;
    }

  H5Dclose(dataset_id);
  if(this->WriteHDF5FileToDisk)
    {
    H5LTmake_dataset_float(group_id2, name, 3, count, fa->GetPointer(0));
    }
  data->GetPointData()->AddArray(fa);
  fa->Delete();
  //data->GetPointData()->SetActiveAttribute(name, vtkDataSetAttributes::SCALARS);
  tdata->FastDelete();
}

void
vtkETHZGeoPhysicsHDF5Reader::ReadVector(vtkImageData *data, hid_t dataset_id, int *outExt, const char *name)
{
  hsize_t     dimsf[4];
#if H5_VERS_RELEASE == 2
  hssize_t     start[4];         /* hyperslab offset in the file */
#else
  hsize_t     start[4];         /* hyperslab offset in the file */
#endif
  hsize_t count[4];
  hid_t       filespace, memspace;

  vtkFloatArray *fa = vtkFloatArray::New();
  fa->SetNumberOfComponents(3);
  fa->SetNumberOfTuples(numberOfTuples);

  fa->SetName(name);
  //cerr << "Proc[" << "?" << "/" << "?" << "] allocating " << numberOfTuples  << " sizeof(float) * bytes for " << name << endl;

  filespace = H5Dget_space(dataset_id);
  int rank = H5Sget_simple_extent_ndims(filespace);
  H5Sget_simple_extent_dims(filespace, dimsf, NULL);

  //cerr << "filespace[] = "; for(int u=0; u < rank; u++) cerr << dimsf[u] << " "; cerr << endl;
//Look at reverse order below!
  start[3] = 0;
  start[2] = outExt[0];
  start[1] = outExt[2];
  start[0] = outExt[4];

  count[3] = 3;
  count[2]  = outExt[1] - outExt[0] + 1;
  count[1]  = outExt[3] - outExt[2] + 1;
  count[0]  = outExt[5] - outExt[4] + 1;

  //cerr << "File Start[] "; for(int u=0; u < rank; u++) cerr << start[u] << " "; cerr << endl;
  //cerr << "File Count[] "; for(int u=0; u < rank; u++) cerr << count[u] << " "; cerr <<endl;

  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, NULL, count, NULL);

  count[0] = outExt[1] - outExt[0] + 1;
  count[1] = outExt[3] - outExt[2] + 1;
  count[2] = outExt[5] - outExt[4] + 1;
  count[3] = 3;
  memspace  = H5Screate_simple(4, count, NULL);
  start[0] = start[1] = start[2] = 0;
  start[3] = 0;
  //cerr << "Mem Start[] "; for(int u=0; u < rank; u++) cerr << start[u] << " "; cerr << endl;
  //cerr << "Mem Count[] "; for(int u=0; u < rank; u++) cerr << count[u] << " "; cerr << endl;

  H5Sselect_hyperslab(memspace, H5S_SELECT_SET, start, NULL, count, NULL);
  
  H5Dread(dataset_id, DATA_LENGTH, memspace, filespace, H5P_DEFAULT, fa->GetPointer(0));
  cerr << "read " << name << endl;

  H5Dclose(dataset_id);
  if(this->WriteHDF5FileToDisk)
    {
    H5LTmake_dataset_float(group_id2, name, 3, count, fa->GetPointer(0));
    }
  //data->GetPointData()->AddArray(fa);
  data->GetPointData()->SetVectors(fa);
  fa->Delete();
  data->GetPointData()->SetActiveAttribute(name, vtkDataSetAttributes::VECTORS);
  //tdata->FastDelete();
}

void
vtkETHZGeoPhysicsHDF5Reader::ShuffleFromFORTRAN2C_Vector(float *dest, vtkFloatArray *src, hsize_t *dimsm)
{
  int t = dimsm[2] * dimsm[1];
  int Zoffset, offset0, offset1;
  for(int z=0; z < dimsm[0]; z++)
  {
    Zoffset = z * t;
    for(int y=0; y < dimsm[1]; y++)
      {
      offset0 = y * dimsm[2];
      offset1 = offset0 + z;
      offset0 += Zoffset;
      for(int x=0; x < dimsm[2]; x++)
        {
        *(dest + 3*(offset0 + x)) = src->GetValue(x*dimsm[1]*dimsm[2]+offset1);
	}
      }
  }
}
void vtkETHZGeoPhysicsHDF5Reader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "Reader:\n";
  os << indent << indent << "Geom file      : " << this->FileName << endl;
}

void vtkETHZGeoPhysicsHDF5Reader::EnablePointArray(const char* name)
{
  this->SetPointArrayStatus(name, 1);
}

void vtkETHZGeoPhysicsHDF5Reader::DisablePointArray(const char* name)
{
  this->SetPointArrayStatus(name, 0);
}

void vtkETHZGeoPhysicsHDF5Reader::EnableAllPointArrays()
{
  this->PointDataArraySelection->EnableAllArrays();
}

void vtkETHZGeoPhysicsHDF5Reader::DisableAllPointArrays()
{
  this->PointDataArraySelection->DisableAllArrays();
}

const char* vtkETHZGeoPhysicsHDF5Reader::GetPointArrayName(int index)
{
  return this->PointDataArraySelection->GetArrayName(index);
}

int vtkETHZGeoPhysicsHDF5Reader::GetPointArrayStatus(const char* name)
{
  return this->PointDataArraySelection->ArrayIsEnabled(name);
}

void vtkETHZGeoPhysicsHDF5Reader::SetPointArrayStatus(const char* name, int status)
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

int vtkETHZGeoPhysicsHDF5Reader::GetNumberOfPointArrays()
{
  return this->PointDataArraySelection->GetNumberOfArrays();
}

#ifdef CELL
void vtkETHZGeoPhysicsHDF5Reader::EnableCellArray(const char* name)
{
  this->SetCellArrayStatus(name, 1);
}

void vtkETHZGeoPhysicsHDF5Reader::DisableCellArray(const char* name)
{
  this->SetCellArrayStatus(name, 0);
}

void vtkETHZGeoPhysicsHDF5Reader::EnableAllCellArrays()
{
  this->CellDataArraySelection->EnableAllArrays();
}

void vtkETHZGeoPhysicsHDF5Reader::DisableAllCellArrays()
{
  this->CellDataArraySelection->DisableAllArrays();
}

const char* vtkETHZGeoPhysicsHDF5Reader::GetCellArrayName(int index)
{
  return this->CellDataArraySelection->GetArrayName(index);
}

int vtkETHZGeoPhysicsHDF5Reader::GetCellArrayStatus(const char* name)
{
  return this->CellDataArraySelection->ArrayIsEnabled(name);
}

void vtkETHZGeoPhysicsHDF5Reader::SetCellArrayStatus(const char* name, int status)
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

int vtkETHZGeoPhysicsHDF5Reader::GetNumberOfCellArrays()
{
  return this->CellDataArraySelection->GetNumberOfArrays();
}
#endif
