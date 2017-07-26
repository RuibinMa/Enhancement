#ifndef TRACERMARKER_H
#define TRACERMARKER_H


#include "object.h"
#include <vtkSmoothPolyDataFilter.h>
#include <vtkDecimatePro.h>
#include <vtkFillHolesFilter.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkExtractSelection.h>
#include <vtkSelection.h>
#include <vtkSelectionNode.h>
#include <vtkVertexGlyphFilter.h>


class  TracerMarker : public Object
{
private:
    double color[3];
    double size;
public:
    TracerMarker();
    ~TracerMarker(){}
    void InputPoints(vtkSmartPointer<vtkPoints> points);
    void SetColor(double r, double g, double b);
    void SetSize(double s);
};


#endif // TRACERMARKER_H
