#ifndef SWEEPINGPLANE_H
#define SWEEPINGPLANE_H


#include "object.h"
#include <vtkSmoothPolyDataFilter.h>
#include <vtkDecimatePro.h>
#include <vtkFillHolesFilter.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkExtractSelection.h>
#include <vtkSelection.h>
#include <vtkSelectionNode.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkExtractEdges.h>

class  SweepingPlane : public Object
{
private:
    double color[3];
    double linewidth;
    bool boundary;
public:
    SweepingPlane();
    ~SweepingPlane(){}
    void InputData(vtkSmartPointer<vtkPolyData> poly);
};

#endif // SWEEPINGPLANE_H
