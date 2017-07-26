#ifndef OPTIMIZATION_H
#define OPTIMIZATION_H

#include "centerline.h"
#include <Eigen/Sparse>
#include <unsupported/Eigen/IterativeSolvers>
#include <vtkExtractEdges.h>
#include <vector>
#include <omp.h>

#define REGWEIGHT 0
#define BENDWEIGHT 1

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;

void constructAandb(std::map<vtkIdType, double> &coefficientMap, double *b,
                vtkSmartPointer<vtkPolyData> t_colon, vtkSmartPointer<vtkPolyData> SurfaceLineUp,
                bool *Is_Fixed, vtkSmartPointer<vtkIdList> Ids,int* InvertIds, RenderManager* t_rendermanager);

vtkSmartPointer<vtkPolyData> Optimize(vtkSmartPointer<vtkPolyData> t_colon,
                                      vtkSmartPointer<vtkPolyData> SurfaceLineUp, bool *Is_Fixed, RenderManager* t_rendermanager);

void updateA(int m, int i1, int i2, double weight, std::map<vtkIdType, double> &coefficientMap);

vtkSmartPointer<vtkIdList> GetFacetsOfEdge(vtkSmartPointer<vtkPolyData> mesh, int idx1, int idx2);

double ComputeStretchWeight(vtkSmartPointer<vtkPolyData> t_colon, int idx1, int idx2);

vtkSmartPointer<vtkIdList> GetConnectedVertices(vtkSmartPointer<vtkPolyData> mesh, vtkIdType id);

double ComputeFacetArea(vtkSmartPointer<vtkPolyData> t_colon,  vtkIdType cellid);
#endif // OPTIMIZATION_H
