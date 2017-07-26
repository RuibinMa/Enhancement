#ifndef CENTERLINE_H
#define CENTERLINE_H
#include "vtkpolydatagroup.h"
#include "object.h"
#include "rendermanager.h"
#include "filemanager.h"
#include "optimization.h"
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <vtkPlane.h>
#include <vtkPlanes.h>
#include <vtkPlaneSource.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkSplineFilter.h>
#include <vtkCardinalSpline.h>
#include <vtkGlyph3D.h>
#include <vtkParametricSpline.h>
#include <vtkParametricFunctionSource.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkCellArray.h>
#include <vtkLine.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkThinPlateSplineTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkPointLocator.h>
#include <vtkCutter.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkUnsignedCharArray.h>
#include <vtkPointData.h>
#include <vtkAppendPolyData.h>
#include <vtkCleanPolyData.h>
#include <vtkMatrix3x3.h>
#include <vtkTriangleFilter.h>
#include <vtkFeatureEdges.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkTriangle.h>
#include <vtkClipPolyData.h>
#include <vtkDecimatePro.h>
#include <vtkDijkstraGraphGeodesicPath.h>
#include <vtkIntersectionPolyDataFilter.h>
#include <vtkLandmarkTransform.h>

class Centerline : public Object
{
private:
    double* Gaussian1DKernel(int len, double std);
    double* GaussianDerivative1DKernel(int len, double std);
    double* conv(double *A, double *B, int lenA, int lenB, int *lenC);
public:
    Centerline();
    ~Centerline(){}
    int GetNumberOfPoints();
    vtkSmartPointer<vtkPlane> GetVerticalPlane(vtkIdType i, double* normalDirection);
    void SmoothCenterline(double std, vtkSmartPointer<vtkIntArray> ViolationNums = NULL);
    void Smooth(vtkSmartPointer<vtkDoubleArray> Candidate, double std);
    void GaussianTangents(vtkSmartPointer<vtkDoubleArray> Tangents, double std);
    void GaussianNormals(vtkSmartPointer<vtkDoubleArray> Normals, vtkSmartPointer<vtkDoubleArray> Curvatures, vtkSmartPointer<vtkDoubleArray> Tangents, double std);
    void UniformSample(int resolution, vtkSmartPointer<vtkPolyData> line = NULL);
    void splineTangent(double *tangent, vtkSmartPointer<vtkParametricSpline> spline, double t_u, double stepsize);
    void splineNormal(double *normal, vtkSmartPointer<vtkParametricSpline> spline, double t_u, double stepsize, double &curvature);
    void PutNormalsOnSameSide(vtkSmartPointer<vtkDoubleArray> Normals, vtkSmartPointer<vtkDoubleArray> Curvatures);
    double splineTorsion(vtkSmartPointer<vtkParametricSpline> spline, double t_u, double stepsize);

    vtkSmartPointer<vtkPolyData> EliminateTorsion(RenderManager* t_rendermanager,RenderManager* t_rendermanager_right, vtkSmartPointer<vtkPolyData> t_colon, FileManager* t_filemanager);
    vtkSmartPointer<vtkPolyData> Deformation(vtkSmartPointer<vtkDoubleArray> S, vtkSmartPointer<vtkDoubleArray> Curvatures,
                                             vtkSmartPointer<vtkDoubleArray> Tangents, vtkSmartPointer<vtkDoubleArray> Normals,
                                             vtkSmartPointer<vtkPolyData> t_colon, RenderManager *t_rendermanager,
                                             vtkSmartPointer<vtkDoubleArray> PlaneOriginals, vtkSmartPointer<vtkDoubleArray> PlaneNormals,
                                             vtkSmartPointer<vtkDoubleArray> RefDirections, FileManager* t_filemanager);
    vtkSmartPointer<vtkPolyData> Deformation_v2(vtkSmartPointer<vtkDoubleArray> S, vtkSmartPointer<vtkDoubleArray> Curvatures,
                                             vtkSmartPointer<vtkDoubleArray> Tangents, vtkSmartPointer<vtkDoubleArray> Normals,
                                             vtkSmartPointer<vtkPolyData> t_colon, RenderManager *t_rendermanager,
                                             vtkSmartPointer<vtkDoubleArray> PlaneOriginals, vtkSmartPointer<vtkDoubleArray> PlaneNormals,
                                             vtkSmartPointer<vtkDoubleArray> RefDirections, FileManager* t_filemanager);
    vtkSmartPointer<vtkPolyData> Deformation_v3(vtkSmartPointer<vtkDoubleArray> S, vtkSmartPointer<vtkDoubleArray> Curvatures, vtkSmartPointer<vtkIdList> CurvaturePointIds,
                                             vtkSmartPointer<vtkDoubleArray> Tangents, vtkSmartPointer<vtkDoubleArray> Normals,
                                             vtkSmartPointer<vtkPolyData> t_colon, RenderManager *t_rendermanager,
                                             vtkSmartPointer<vtkDoubleArray> PlaneOriginals, vtkSmartPointer<vtkDoubleArray> PlaneNormals,
                                             vtkSmartPointer<vtkDoubleArray> RefDirections, FileManager* t_filemanager);
    vtkSmartPointer<vtkPolyData> Deformation_v3_1(vtkSmartPointer<vtkDoubleArray> S, vtkSmartPointer<vtkDoubleArray> Curvatures, vtkSmartPointer<vtkIdList> CurvaturePointIds,
                                             vtkSmartPointer<vtkDoubleArray> Tangents, vtkSmartPointer<vtkDoubleArray> Normals,
                                             vtkSmartPointer<vtkPolyData> t_colon, RenderManager *t_rendermanager, RenderManager *t_rendermanager_right,
                                             vtkSmartPointer<vtkDoubleArray> PlaneOriginals, vtkSmartPointer<vtkDoubleArray> PlaneNormals,
                                             vtkSmartPointer<vtkDoubleArray> RefDirections, FileManager* t_filemanager);
    vtkSmartPointer<vtkPolyData> Deformation_v3_2(vtkSmartPointer<vtkDoubleArray> S, vtkSmartPointer<vtkDoubleArray> Curvatures, vtkSmartPointer<vtkIdList> CurvaturePointIds,
                                             vtkSmartPointer<vtkDoubleArray> Tangents, vtkSmartPointer<vtkDoubleArray> Normals,
                                             vtkSmartPointer<vtkPolyData> t_colon, RenderManager *t_rendermanager,
                                             vtkSmartPointer<vtkDoubleArray> PlaneOriginals, vtkSmartPointer<vtkDoubleArray> PlaneNormals,
                                             vtkSmartPointer<vtkDoubleArray> RefDirections, FileManager* t_filemanager);
    vtkSmartPointer<vtkPolyData> Deformation_v4(vtkSmartPointer<vtkDoubleArray> S, vtkSmartPointer<vtkDoubleArray> Curvatures,
                                             vtkSmartPointer<vtkDoubleArray> Tangents, vtkSmartPointer<vtkDoubleArray> Normals,
                                             vtkSmartPointer<vtkPolyData> t_colon, RenderManager *t_rendermanager,
                                             vtkSmartPointer<vtkDoubleArray> PlaneOriginals, vtkSmartPointer<vtkDoubleArray> PlaneNormals,
                                             vtkSmartPointer<vtkDoubleArray> RefDirections, FileManager* t_filemanager);
    void VisualizeTNB(vtkSmartPointer<vtkDoubleArray> S, vtkSmartPointer<vtkDoubleArray> Curvatures,
                      vtkSmartPointer<vtkDoubleArray> Tangents, vtkSmartPointer<vtkDoubleArray> Normals, vtkSmartPointer<vtkDoubleArray> Binormals,
                      RenderManager* t_rendermanager);
    void VisualizeSpoke(vtkSmartPointer<vtkPoints> CurvaturePoints, vtkSmartPointer<vtkIntArray> ViolationNums, RenderManager *t_rendermanager);
    void VisualizeOriginalCurve(RenderManager* t_rendermanager);


    vtkSmartPointer<vtkPolyData> ReorderContour(vtkSmartPointer<vtkPolyData> cutCircle);
    vtkSmartPointer<vtkPolyData> FormPlate(vtkSmartPointer<vtkPolyData> cutCircle);
    double SinglePath(double **costVrt, double **costHrz,int sx, int sy, int tx, int ty, int *steps);
    void ContoursToSurface(RenderManager* t_rendermanager, FileManager* t_filemanager);
    vtkSmartPointer<vtkPolyData> ConnectTwoContours(vtkSmartPointer<vtkPolyData> circle1, vtkSmartPointer<vtkPolyData> circle2);
    void ConnectTwoContoursTest(RenderManager* t_rendermanager, FileManager* t_filemanager);
    vtkSmartPointer<vtkPolyData> PieceBetweenPlanes(vtkSmartPointer<vtkPolyData> t_colon,
                                                    vtkSmartPointer<vtkDoubleArray> PlaneOriginals, vtkSmartPointer<vtkDoubleArray> PlaneNormals,
                                                    vtkIdType left, vtkIdType right, vtkSmartPointer<vtkPolyData> lastpiece = NULL);
    vtkIdType GetTheSectionIdOfAPoint(vtkSmartPointer<vtkPolyData> t_colon, vtkIdType pointid,
                                      vtkSmartPointer<vtkDoubleArray> PlaneOriginals, vtkSmartPointer<vtkDoubleArray> PlaneNormals, vtkSmartPointer<vtkIdList> CurvaturePointIds);
    vtkIdType GetTheSectionIdOfAPoint_v2(vtkSmartPointer<vtkPolyData> t_colon, vtkIdType pointid,
                                         vtkSmartPointer<vtkDoubleArray> PlaneOriginals, vtkSmartPointer<vtkDoubleArray> PlaneNormals, vtkSmartPointer<vtkIdList> CurvaturePointIds);
    void GetSectionIds(vtkPolyData* t_colon, vtkIdType pointid, vtkIdList* SectionIds,
                       vtkDoubleArray* PlaneOriginals, vtkDoubleArray* PlaneNormals);
    void GetSectionIds_loop(vtkPolyData* t_colon, vtkIdType seed, vtkIdList* SectionIds,
                       vtkDoubleArray* PlaneOriginals, vtkDoubleArray* PlaneNormals);
    void GetSectionIds_loop_combinehighcurvatures(vtkPolyData* t_colon, vtkIdType seed, vtkIdList* SectionIds,
                       vtkDoubleArray* PlaneOriginals, vtkDoubleArray* PlaneNormals, vtkDoubleArray* Curvatures);
    void GetSectionIds_loop_v2(vtkPolyData* t_colon, vtkIdType seed, vtkIdList* SectionIds,
                       vtkDoubleArray* PlaneOriginals, vtkDoubleArray* PlaneNormals, vtkPolyDataGroup* CircleGroup);
    void GetRegionIds(vtkPolyData* t_colon, bool* Is_Fixed, int* RegionIds, int *NumberOfRegions);
    vtkSmartPointer<vtkIntArray> GetSectionHeadAndTail(vtkSmartPointer<vtkDoubleArray> Curvatures, double thres);
    double LorentzianInterpolationFactor(double l, double angle);

};

void VisualizePoints(vtkSmartPointer<vtkPoints> points, double r, double g, double b, float size, RenderManager* t_rendermanager);
#endif // CENTERLINE_H
