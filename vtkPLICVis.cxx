#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkSmartPointer.h"
#include "vtkRectilinearGrid.h"
#include "vtkDataArray.h"
#include "vtkFloatArray.h"
#include "vtkCellData.h"
#include "vtkCellArray.h"
#include "vtkPointData.h"
#include "vtkRectilinearGrid.h"

#include "vtkPLICVis.h"
#include "plicvis_impl.h"

#include <vector>

vtkStandardNewMacro(vtkPLICVis);

//----------------------------------------------------------------------------
vtkPLICVis::vtkPLICVis()
{
  this->SetNumberOfInputPorts(1);
}
//----------------------------------------------------------------------------
int vtkPLICVis::FillInputPortInformation(int port, vtkInformation* info)
{
  if (!this->Superclass::FillInputPortInformation(port, info)) {
    return 0;
  }
  if (port == 0) {
    info->Set( vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkRectilinearGrid");
    return 1;
  }
  return 0;
}
//----------------------------------------------------------------------------
int vtkPLICVis::RequestUpdateExtent(vtkInformation *vtkNotUsed(request),
				    vtkInformationVector **inputVector,
				    vtkInformationVector *outputVector)
{
  // set one ghost level -----------------------------------------------------
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS(), 1);
  
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  outInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS(), 1);

  return 1;
}
//----------------------------------------------------------------------------
int vtkPLICVis::RequestData(vtkInformation *vtkNotUsed(request),
			    vtkInformationVector **inputVector,
			    vtkInformationVector *outputVector)
{
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkSmartPointer<vtkRectilinearGrid> input = vtkRectilinearGrid::
    SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkPolyData *output = vtkPolyData::
      SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkDataArray *coords[3] = {input->GetXCoordinates(),
			     input->GetYCoordinates(),
			     input->GetZCoordinates()};
  std::vector<float> dx;
  std::vector<float> dy;
  std::vector<float> dz;

  nodeCoordsToEdgeLengths(coords[0], dx);
  nodeCoordsToEdgeLengths(coords[1], dy);
  nodeCoordsToEdgeLengths(coords[2], dz);

  std::vector<float3> vertices(0);
  std::vector<int> indices(0);

  int nodeRes[3];
  input->GetDimensions(nodeRes);
  int cellRes[3] = {nodeRes[0]-1,nodeRes[1]-1,nodeRes[2]-1};

  int extent[6];
  input->GetExtent(extent);

  int imin = 1;
  int imax = cellRes[0]-1;
  int jmin = 1;
  int jmax = cellRes[1]-1;
  int kmin = 1;
  int kmax = cellRes[2]-1;

  vtkDataArray *data = input->GetCellData()->GetArray("Data");
  int vertexID = 0;
  for (int k = kmin; k < kmax; ++k) {
    float oz = coords[2]->GetComponent(k,0);
    for (int j = jmin; j < jmax; ++j) {
      float oy = coords[1]->GetComponent(j,0);
      for (int i = imin; i < imax; ++i) {
	float ox = coords[0]->GetComponent(i,0);

	int idx = i + j*cellRes[0] + k*cellRes[0]*cellRes[1];
	float f = data->GetComponent(idx, 0);
	if (f <= EMF0) {
	  continue;
	}

	if (f >= EMF1 && !interfaceCell(data, i, j, k, cellRes)) {
	  continue;
	}

	float3 grad = computeGradient(data, i, j, k, cellRes, dx, dy, dz);

	generatePLIC(f, grad, dx[i], dy[j], dz[k], 
	 	     ox, oy, oz, vertices, indices, vertexID);
      }
    }
  }

  std::vector<std::vector<int>> borders(0);
  
  //  extractPLICBorders(vertices, indices, borders);

  const int numPoints = vertices.size();
  const int numTriangles = indices.size()/3;

  vtkPoints *points = vtkPoints::New();
  points->SetNumberOfPoints(numPoints);
  for (int i = 0; i < numPoints; ++i) {
    float p[3] = {vertices[i].x, vertices[i].y, vertices[i].z};
    points->SetPoint(i, p);
  }

  vtkIdTypeArray *cellIndices = vtkIdTypeArray::New();
  cellIndices->SetNumberOfComponents(1);
  cellIndices->SetNumberOfTuples(numTriangles*4);
  for (int i = 0; i < numTriangles; ++i) {
    cellIndices->SetValue(i*4, 3);
    cellIndices->SetValue(i*4+1, indices[i*3+0]);
    cellIndices->SetValue(i*4+2, indices[i*3+1]);
    cellIndices->SetValue(i*4+3, indices[i*3+2]);
  }
  vtkCellArray *cells = vtkCellArray::New();
  cells->SetCells(numTriangles, cellIndices);

  output->SetPoints(points);
  output->SetPolys(cells);

  return 1;
}

////////// External Operators /////////////
void vtkPLICVis::PrintSelf(ostream &os, vtkIndent indent)
{
}
