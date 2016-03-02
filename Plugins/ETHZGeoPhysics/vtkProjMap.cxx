//#define PARALLEL_DEBUG 1
#include "vtkDataSetAttributes.h"
#include "vtkDoubleArray.h"
#include "vtkErrorCode.h"
#include "vtkFieldData.h"
#include "vtkFloatArray.h"
#include "vtkStructuredGrid.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkPlaneSource.h"
#include "vtkGeoProjectionSource.h"
#include "vtkGeoProjection.h"
#include "vtkSphericalTransform.h"

#include <math.h>
#include "vtkProjMap.h"

//vtkCxxRevisionMacro(vtkProjMap, "$Revision: 1.1 $");
vtkStandardNewMacro(vtkProjMap);

vtkProjMap::vtkProjMap()
{
  this->DebugOff();
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(2);
  this->Projection = 70;
}

int vtkProjMap::FillInputPortInformation(
  int port, vtkInformation* info)
{
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkStructuredGrid");
  return 1;
}

int vtkProjMap::FillOutputPortInformation(
  int port, vtkInformation* info)
{
  if(port == 0)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkStructuredGrid");
  else if (port == 1)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
  return 1;
}

vtkProjMap::~vtkProjMap()
{
}

int vtkProjMap::RequestInformation(vtkInformation* request,
                         vtkInformationVector** inputVector,
                         vtkInformationVector* outputVector)
{
  vtkInformation* outInfo = outputVector->GetInformationObject(0);
  //outInfo->Set(vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(), -1);
  int wholeExt[6];
  wholeExt[0] = 0;
  wholeExt[1] = -1;
  wholeExt[2] = 0;
  wholeExt[3] = -1;
  wholeExt[4] = 0;
  wholeExt[5] = -1;
  outInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), wholeExt, 6);

  outInfo = outputVector->GetInformationObject(1);
  outInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), wholeExt, 6);
  return 1;
}

int vtkProjMap::RequestData(
                vtkInformation* vtkNotUsed(request),
                vtkInformationVector** inputVector,
                vtkInformationVector* outputVector)
{
  vtkDebugMacro( << "vtkProjMap::RequestData(BEGIN)");
  this->UpdateProgress(0);
  vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation* outInfo = outputVector->GetInformationObject(0);

  int updatePiece     = inInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
  int updateNumPieces = inInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());

#ifdef PARALLEL_DEBUG
  char fname[23];
  sprintf(fname,"/scratch/daint/jfavre/out.%02d.txt", updatePiece);
  FILE *errs = fopen(fname, "w");
  fprintf(errs, "updatePiece = %d, updateNumPieces = %d\n\n", updatePiece, updateNumPieces);
#endif

  vtkStructuredGrid *inData = static_cast<vtkStructuredGrid*>(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  int *inExt = inData->GetExtent();

  int ext[6], wext[6];

  for (int i=0; i < 6; i++)
    ext[i] = inExt[i];

  //outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), ext);
  outInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), wext);
  int resolution = (ext[5] - ext[4] + 1) * (ext[3] - ext[2] + 1) * (ext[1] - ext[0] + 1);
#ifdef PARALLEL_DEBUG
  fprintf(errs,"updateExt[6] = %d,%d,%d,%d,%d,%d\n",ext[0] ,ext[1] , ext[2] ,  ext[3] , ext[4] , ext[5] );
  fprintf(errs,"wholeExt[6] = %d,%d,%d,%d,%d,%d\n",wext[0] ,wext[1] , wext[2] ,  wext[3] , wext[4] , wext[5] );
  fprintf(errs, "resolution = %d\n", resolution);
#endif
  
  vtkPointData *inPD = inData->GetPointData();
  //cerr << "input() has " << resolution << " points " << endl;

  vtkStructuredGrid *outData = vtkStructuredGrid::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));
  outData->SetExtent(ext);

  outData->ShallowCopy((vtkDataObject*)inData);

  //vtkGeoProjection *proj1 = vtkGeoProjection::New();
  //for (int i=0; i < 180; i++)
    //cerr << i << " = " << proj1->GetProjectionName(i) << endl;
  //proj1->Delete();

  vtkGeoProjectionSource *projsrc = vtkGeoProjectionSource::New();
  projsrc->SetProjection(this->Projection);
  //projsrc->Update();

  // do not use the input's Points since they have been transformed already to spherical space.
  // rebuild a 2D map of points
  vtkPoints *inPts = vtkPoints::New();
  //inPts->DeepCopy(inData->GetPoints());
  vtkDataArray *Rarray, *Tarray, *Parray;
  inPts->Allocate(resolution);
  inPts->SetNumberOfPoints(resolution);
  int id=0;
  double theta, theta0, phi, phi0, previous_phi, radius, x[3], bounds[6];
// the angle phi calculated is non-continuous, ranging from 0 to PI/2, then from -PI/2 to PI/2, then again from -PI/2 to 0. I need to shift it to range from 0 to 2*PI
  for (int i = ext[4]; i <= ext[5]; ++i)
    {
    Rarray = inData->GetFieldData()->GetArray("R");
    if(Rarray)
      radius = Rarray->GetTuple1(i);
    else
      radius = 1.0 + i;
    //cerr << "ConstructPoints() allocates " << resolution << " points for radius " << radius << endl;
    Tarray = inData->GetFieldData()->GetArray("theta");
    for (int j = ext[2]; j <= ext[3]; ++j) // theta index
      {
      theta0 = j * -180.0 /(ext[3]-ext[2])+90.0;
      if(Tarray)
        theta = Tarray->GetTuple1(j);
      else
        theta = (j-wext[2])* M_PI/(wext[3] - wext[2]);
//can optimize later to avoid multiply and then divide by ppi
      theta = (theta/M_PI)*180.0;

      if(theta > 180.0) theta=180;
      theta -= 90.0;
      theta *= -1.0;
      previous_phi = -M_PI;
      Parray = inData->GetFieldData()->GetArray("phi");
      for (int k = ext[0]; k <= ext[1]; ++k) // phi index
        {
        phi0 = k * 360.0 /(ext[1]-ext[0])-180.0;

        inData->GetPoints()->GetPoint(id, x); // expressed in world coords, must be converted to degrees
 
        //theta = 90.0  - 180.0*acos(x[2]/radius)/M_PI;
        //if(theta < 0) { cerr <<"negative theta = " << theta << endl; theta=0;}
        //if(theta > 180) { cerr <<"super-positive theta = " << theta << endl; theta=180;}
        phi = atan(x[1]/x[0]);
//if(j == (ext[3]-ext[2])/2)cerr << phi0 << ", " << phi << " ";
        if(phi < previous_phi){ phi += M_PI;} // to push over the 2nd ramp to join with 1st ramp
        if(phi < previous_phi){ phi += M_PI;} // to push over the 3rd ramp to join with 2nd ramp
//if(j == (ext[3]-ext[2])/2)cerr << phi << " ";
        previous_phi = phi; 
        // phi now ranges from 0 to 2*PI
        //phi = (phi/M_PI - 1.0)*180.0;
        if(Parray)
          phi = Parray->GetTuple1(k);
        else
          phi = (k-wext[0])* 2.0 * M_PI/(wext[1] - wext[0]);
        phi = (phi/M_PI)*180.0;
        if(phi < 0) phi=0;
        if(phi > 360.0) phi=360.0;
//if(j == (ext[3]-ext[2])/2)cerr << phi << endl;
//cerr << ph2i0 << ", " << phi << endl;
        phi -= 180.0;
        inPts->SetPoint(id, phi, theta, radius);
        ++id;
        }
      }
    }
  vtkPoints *newPoints = vtkPoints::New();
  projsrc->GetTransform()->TransformPoints(inPts, newPoints);
  outData->SetPoints(newPoints);
  newPoints->Modified();
  newPoints->GetBounds(bounds);
  std::cerr << "Output Bounds[6] = " << bounds[0] << " "  << bounds[1] << " " << bounds[2] << " " << bounds[3] << " " << bounds[4] << " " << bounds[5] <<  endl;
  newPoints->Delete();
  inPts->Delete();
/////////////////////////////////////////////////////////
// do now the second output for the Land texture
/////////////////////////////////////////////////////////

  outInfo = outputVector->GetInformationObject(1);

  vtkPolyData *land = vtkPolyData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkPlaneSource *plane = vtkPlaneSource::New();
  plane->SetOrigin(-180, -90, 0);
  plane->SetPoint1(180, -90, 0);
  plane->SetPoint2(-180, 90, 0.0);
  plane->SetXResolution(121);
  plane->SetYResolution(121);
  plane->Update();
  land->ShallowCopy(plane->GetOutput());
  vtkPoints *landPoints = vtkPoints::New();
  projsrc->GetTransform()->TransformPoints(plane->GetOutput()->GetPoints(), landPoints);
  land->SetPoints(landPoints);
  landPoints->Delete();
  plane->Delete();
  projsrc->Delete();
  vtkDebugMacro( << "vtkProjMap::RequestData(END)");
  this->UpdateProgress(1);
#ifdef PARALLEL_DEBUG
  fclose(errs);
#endif
  return 1;
}

void vtkProjMap::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}
