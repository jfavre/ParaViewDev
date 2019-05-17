#ifndef __vtkHDF5Reader_h
#define __vtkHDF5Reader_h

#include "HDF5Module.h" // for export macro
#include "vtkUnstructuredGridAlgorithm.h"

class HDF5_EXPORT vtkHDF5Reader : public vtkUnstructuredGridAlgorithm
{
public:
  static vtkHDF5Reader *New();
  vtkTypeMacro(vtkHDF5Reader, vtkUnstructuredGridAlgorithm);

  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

protected:
  vtkHDF5Reader();
  ~vtkHDF5Reader();

  int RequestData(vtkInformation*,
                  vtkInformationVector**,
                  vtkInformationVector*);
  int RequestInformation(vtkInformation*,
                         vtkInformationVector**,
                         vtkInformationVector*);
  
  char *FileName;
  int NbNodes;
 private:
  vtkHDF5Reader(const vtkHDF5Reader&) = delete;
  void operator=(const vtkHDF5Reader&) = delete;
};

#endif
