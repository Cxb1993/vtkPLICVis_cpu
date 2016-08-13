#ifndef PLICVIS_IMPL_H
#define PLICVIS_IMPL_H

#include "vtkDataArray.h"
#include <vector>

#define EMF0 0.000001f
#define EMF1 0.999999f


float3 computeGradient(vtkDataArray *data, int cell_i, int cell_j, int cell_k, int cellRes[3],
		       std::vector<float> &dx,std::vector<float> &dy,std::vector<float> &dz);

void generatePLIC(float f, float3 grad, 
		  float dx, float dy, float dz,
		  float ox, float oy, float oz, 
		  std::vector<float3> &vertices,
		  std::vector<int> &indices,
		  int &vertexID);

void nodeCoordsToEdgeLengths(vtkDataArray *coords, std::vector<float> &dx);

int interfaceCell(vtkDataArray *data, int cell_i, int cell_j, int cell_k, int cellRes[3]);

void extractPLICBorders(std::vector<float3> &vertices,
			std::vector<int> &indices,
			std::vector<std::vector<int>> &borders);

#endif//PLICVIS_IMPL_H
