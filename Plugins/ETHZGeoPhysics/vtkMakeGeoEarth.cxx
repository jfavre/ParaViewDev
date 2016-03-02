#include "vtkPlaneSource.h"
#include "vtkXMLPolyDataReader.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h" 
#include "vtkObjectFactory.h"
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkPolyDataNormals.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkGlobeSource.h"
#include "vtkMakeGeoEarth.h"

#include "vtkPVConfig.h" // For PARAVIEW_DATA_ROOT
//#include "vtkPQConfig.h" // For PARAVIEW_DATA_ROOT

//vtkCxxRevisionMacro(vtkMakeGeoEarth, "$Revision: 1.1 $");
vtkStandardNewMacro(vtkMakeGeoEarth);

vtkMakeGeoEarth::vtkMakeGeoEarth()
{
  this->DebugOff();
  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(2);
  this->Radius = 1.0;
  this->Xres = 90;
  this->Yres = 45;
  vtkPolyData *pd = vtkPolyData::New();
  pd->ReleaseData();
  this->GetExecutive()->SetOutputData(1, pd);
  pd->Delete();
}

vtkMakeGeoEarth::~vtkMakeGeoEarth()
{

}

int vtkMakeGeoEarth::RequestInformation(vtkInformation* request,
                         vtkInformationVector** inputVector,
                         vtkInformationVector* outputVector)
{
  vtkInformation* info;
  info = outputVector->GetInformationObject(0);
  //info->Set(vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(), -1);

  info = outputVector->GetInformationObject(1);
  //info->Set(vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(), -1);
  return 1;
}

int vtkMakeGeoEarth::FillOutputPortInformation(int port, vtkInformation* info)
{
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
  return 1;
}

int vtkMakeGeoEarth::RequestData(vtkInformation* request,
                         vtkInformationVector** inputVector,
                         vtkInformationVector* outputVector)
{
  int i, j;
  vtkInformation* info;
  info = outputVector->GetInformationObject(0);
  vtkPolyData* output0 = vtkPolyData::SafeDownCast(info->Get(vtkDataObject::DATA_OBJECT()));
  info = outputVector->GetInformationObject(1);
  vtkPolyData* output1 = vtkPolyData::SafeDownCast(info->Get(vtkDataObject::DATA_OBJECT()));
  int numPoints;
  vtkPointSet *psInput;
  double step, scaled_tex_x, endTheta = 90.; // default = 180
  step = (endTheta +180.0)/(this->Xres);
  vtkPlaneSource *plane = vtkPlaneSource::New();
  plane->SetOrigin(-180, -90, 0);
  plane->SetPoint1(endTheta, -90, 0);
  plane->SetPoint2(-180, 90, 0.0);
  plane->SetXResolution(this->Xres);
  plane->SetYResolution(this->Yres);
  plane->Update();
  // modify the X texture coordinates to go from 0.0 to 0.75 (if we use 90 degrees in SetPoint1(90, -90, 0))
  vtkDataArray *tcoords = plane->GetOutput()->GetPointData()->GetTCoords();
  double w0, w1;
  for (j = 0; j <= this->Yres; j++)
    {
    for (i = 0; i <= this->Xres; i++)
      {
      // Y texture coord is unchanged
      scaled_tex_x = (i*step)/360.0;
      tcoords->SetComponent(j*(this->Xres+1) + i, 0, scaled_tex_x);
      }
    }
  output0->ShallowCopy(plane->GetOutput());
  vtkPoints *newPoints0 = vtkPoints::New();
  vtkPoints *oldPoints0 = 0;
  psInput = vtkPointSet::SafeDownCast(plane->GetOutput());
  oldPoints0 = psInput->GetPoints();
  newPoints0->DeepCopy(oldPoints0);
  output0->SetPoints(newPoints0);
  numPoints = psInput->GetNumberOfPoints();
  newPoints0->Delete();
  vtkDataArray *normals = (output0->GetPointData()->GetNormals());

  // Convert the points to global coordinates
  for (i = 0; i < numPoints; i++)
    {
    double theta, phi;
    double a[3];
    oldPoints0->GetPoint(i, a);
    theta = a[0];
    phi = a[1];
    double x[3];
    vtkGlobeSource::ComputeGlobePoint(theta, phi, 1.0, x); // use radius = 1, and post multiply to avoid
                                                           // having to divide in the line below for Normals()

    newPoints0->SetPoint(i, x[0]*this->Radius, x[1]*this->Radius, x[2]*this->Radius);
// at this point, normals have not been transformed as JB pointed out.
    normals->SetTuple3(i, x[0], x[1], x[2]);
    }


  vtkXMLPolyDataReader *reader = vtkXMLPolyDataReader::New();
  reader->SetFileName(PARAVIEW_DATA_ROOT "/Data/political.vtp");
  reader->Update();

  output1->ShallowCopy(reader->GetOutput());
  vtkPoints *newPoints = vtkPoints::New();
  vtkPoints *oldPoints = 0;
  psInput = vtkPointSet::SafeDownCast(reader->GetOutput());
  oldPoints = psInput->GetPoints();
  newPoints->DeepCopy(oldPoints);
  output1->SetPoints(newPoints);
  numPoints = psInput->GetNumberOfPoints();
  newPoints->Delete();
  vtkDoubleArray *thetas = vtkDoubleArray::New(); vtkDoubleArray *phis = vtkDoubleArray::New();
  thetas->SetNumberOfTuples(numPoints); phis->SetNumberOfTuples(numPoints);
  thetas->SetName("theta"); phis->SetName("phi");
  oldPoints->GetData()->GetData(0, numPoints-1, 0, 0, thetas);
  oldPoints->GetData()->GetData(0, numPoints-1, 1, 1, phis);
  output1->GetPointData()->AddArray(thetas); output1->GetPointData()->AddArray(phis); thetas->Delete(); phis->Delete();
  // Convert the points to global coordinates
  for (i = 0; i < numPoints; i++)
    {
    double theta, phi;
    double a[3];
    oldPoints->GetPoint(i, a);
    theta = a[0];
    phi = a[1];
    double x[3];
    vtkGlobeSource::ComputeGlobePoint(theta, phi, this->Radius, x);
    newPoints->SetPoint(i, x[0], x[1], x[2]);
    }
  reader->Delete();
  plane->Delete();
  return 1;
}

void vtkMakeGeoEarth::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "Radius:\n" << this->Radius << endl;
}
