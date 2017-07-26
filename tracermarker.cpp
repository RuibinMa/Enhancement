#include "tracermarker.h"

TracerMarker::TracerMarker()
{
    color[0] = 0;
    color[1] = 1;
    color[2] = 0;
    size = 5;
    actor->GetProperty()->SetColor(color);
    actor->GetProperty()->SetPointSize(size);
}
void TracerMarker::InputPoints(vtkSmartPointer<vtkPoints> points)
{
    vtkSmartPointer<vtkPolyData> poly = vtkSmartPointer<vtkPolyData>::New();
    poly->SetPoints(points);
    vtkSmartPointer<vtkVertexGlyphFilter> vertexFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
    vertexFilter->SetInputData(poly);
    vertexFilter->Update();
    model->Reset();
    model->DeepCopy(vertexFilter->GetOutput());
    mapper->SetInputData(model);
    mapper->Update();
    actor->SetMapper(mapper);
}
void TracerMarker::SetColor(double r, double g, double b)
{
    color[0] = r;
    color[1] = g;
    color[2] = b;
    actor->GetProperty()->SetColor(color);
}
void TracerMarker::SetSize(double s)
{
    size = s;
    actor->GetProperty()->SetPointSize(size);
}
