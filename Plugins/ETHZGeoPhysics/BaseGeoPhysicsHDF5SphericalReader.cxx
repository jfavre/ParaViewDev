#include "vtkDataArraySelection.h"
#include "vtkDataSetAttributes.h"
#include "vtkDoubleArray.h"
#include "vtkErrorCode.h"
#include "vtkExtentTranslator.h"
#include "vtkFieldData.h"
#include "vtkMath.h"
#include "vtkFloatArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkIntArray.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include <math.h>


#define DATA_LENGTH H5T_NATIVE_FLOAT

#include "BaseGeoPhysicsHDF5SphericalReader.h"

vtkPoints*
ConstructPoints(int *ext,
                double *SinTheta, double *SinPhi, double *CosTheta, double *CosPhi,
                float *Radius, float *Theta, float *Phi,
                bool xform)
{
  vtkPoints *points;
  double x, y, z, depth, radius;
  double theta, phi, st, ct, temp;
  int i, j, k;
  int id, num;

  points = vtkPoints::New();
  num = (ext[1]-ext[0]+1)*(ext[3]-ext[2]+1)*(ext[5]-ext[4]+1);
  //cerr << "ConstructPoints() allocates " << num << " points\n";
  points->Allocate(num);
  points->SetNumberOfPoints(num);

  id = 0;
  //cerr << "From Radius[" << ext[4] << "-" << ext[5] << "]\n";
  //cerr << "From Theta[" << ext[2] << "-" << ext[3] << "]\n";
  //cerr << "From Phi[" << ext[0] << "-" << ext[1] << "]\n";
 // radius increasing slowest

  for (i = ext[4]; i <= ext[5]; ++i)
    {
    radius = Radius[i];
    for (j = ext[2]; j <= ext[3]; ++j) // theta index
      {
      theta = Theta[j];
      if(xform)
        {
        st = SinTheta[j];
        ct = CosTheta[j];
        z = ct * radius; // outside of innermost loop for efficiency
        temp = st * radius;
        }
      for (k = ext[0]; k <= ext[1]; ++k, ++id) // phi index
        {
        phi = Phi[k];
        if(xform)
          {
          x = temp * CosPhi[k];
          y = temp * SinPhi[k];
          points->SetPoint(id, x, y, z);
          }
        else
          points->SetPoint(id, radius, theta, phi);
        }
      }
    }
  return points;
}


vtkFloatArray *
ReadScalarVar(hid_t *basegroup_id, int *ext, const char *varname)
{
  hid_t   dataset_id, filespace, memspace;
  hsize_t start[3], lstart[3];
  hsize_t count[3], dimsf[4];
  herr_t err;

  //Look at reverse order below!
  start[2] = ext[0];
  start[1] = ext[2];
  start[0] = ext[4];
  count[2]  = ext[1] - ext[0] + 1;
  count[1]  = ext[3] - ext[2] + 1;
  count[0]  = ext[5] - ext[4] + 1;

  if((dataset_id = H5Dopen(*basegroup_id, varname, H5P_DEFAULT)) < 0)
    {
    cerr << "error opening dataset " << varname << endl;
    return NULL;
    }
  else
    {
    int numberOfTuples = (ext[1]-ext[0]+1) * (ext[3]-ext[2]+1) * (ext[5]-ext[4]+1);
    vtkFloatArray *tdata = vtkFloatArray::New();
    tdata->SetNumberOfComponents(1);
    tdata->SetName(varname);
    tdata->SetNumberOfTuples(numberOfTuples);
    filespace = H5Dget_space(dataset_id);
    int rank = H5Sget_simple_extent_ndims(filespace);
    H5Sget_simple_extent_dims(filespace, dimsf, NULL);

    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, NULL, count, NULL);

    memspace  = H5Screate_simple(3, count, NULL);
    lstart[0] = lstart[1] = lstart[2] = 0;
    H5Sselect_hyperslab(memspace, H5S_SELECT_SET, lstart, NULL, count, NULL);
	
    err = H5Dread(dataset_id, DATA_LENGTH, memspace, filespace, H5P_DEFAULT, tdata->GetPointer(0));
    if(err < 0)
      {
      cerr << "H5Dread returned negative value for variable " << varname << endl;
      }
    H5Dclose(dataset_id);
    H5Sclose(memspace);
    H5Sclose(filespace);
    return tdata;
    }
}

vtkFloatArray *
ReadVector(hid_t *root_id, int *ext, const char *varname,
           double *SinTheta, double *SinPhi, double *CosTheta, double *CosPhi,
           bool xform)
{
  hsize_t     dimsf[4];
  hsize_t     start[4];         /* hyperslab offset in the file */
  hsize_t     count[4], block[4];
  hid_t  dataset_id, filespace, memspace, group_id2;
  int rank;
  std::vector<std::string> comp;
  comp.clear();
  comp.push_back("p");
  comp.push_back("t");
  comp.push_back("r");
  std::vector<std::string>::iterator varit;

  //cerr << "allocating " << numberOfTuples  << " sizeof(float) * bytes for " << name << endl;
  int offset = 0;
  int numberOfTuples = (ext[1] - ext[0] + 1) * (ext[3] - ext[2] + 1) * (ext[5] - ext[4] + 1);
  vtkFloatArray *fa = vtkFloatArray::New();
  fa->SetNumberOfComponents(3);
  fa->SetNumberOfTuples(numberOfTuples);
  fa->SetName(varname);
  for(varit = comp.begin(); varit != comp.end(); varit++)
    {
    dataset_id = H5Dopen(*root_id, (varname  + ("/" + (varname + *varit))).c_str(), H5P_DEFAULT);
    //cerr << "vector look for = " << name  + ("/" + (name + *varit)) << endl;
    filespace = H5Dget_space(dataset_id);
    rank = H5Sget_simple_extent_ndims(filespace);
    H5Sget_simple_extent_dims(filespace, dimsf, NULL);
      //Look at reverse order below!
    start[0] = ext[4];
    start[1] = ext[2];
    start[2] = ext[0];

    count[0]  = ext[5] - ext[4] + 1;
    count[1]  = ext[3] - ext[2] + 1;
    count[2]  = ext[1] - ext[0] + 1;
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, NULL, count, NULL);

    count[0] = numberOfTuples;
    count[1] = 3;
    memspace  = H5Screate_simple(2, count, NULL);
  
    start[0] = 0;
    start[1] = offset++;
    count[0] = numberOfTuples;
    count[1] = 1;
    H5Sselect_hyperslab(memspace, H5S_SELECT_SET, start, NULL, count, NULL);
  
    H5Dread(dataset_id, DATA_LENGTH, memspace, filespace, H5P_DEFAULT, fa->GetPointer(0));
    H5Dclose(dataset_id);
    H5Sclose(memspace);
    H5Sclose(filespace);
    }
  if(xform)
    {
    double *d3, d33[3];
    int id = 0;
    for (int i = ext[4]; i <= ext[5]; ++i)
      {
       //cerr << sin(this->Phi->GetValue(k)) << endl;
      //radius = this->Radius->GetValue(i);
      for (int j = ext[2]; j <= ext[3]; ++j) // theta index
        {
        double st = SinTheta[j];
        double ct = CosTheta[j];
        // phi index
        for (int k = ext[0]; k <= ext[1]; ++k)
          { 
	  d3 = fa->GetTuple3(id);
	  d33[0] = st * CosPhi[k] * d3[2] +
	                  ct * CosPhi[k] * d3[1] -
			  SinPhi[k] * d3[0];
	  d33[1] = st * SinPhi[k] * d3[2] +
			 ct * SinPhi[k] * d3[1] +
			 CosPhi[k] * d3[0];
	  d33[2] = ct * d3[2] - st * d3[1];
	  fa->SetTuple3(id, d33[0], d33[1], d33[2]);
          ++id;
          }
        }
      }
    }
  return fa;
}
