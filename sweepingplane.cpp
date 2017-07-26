#include "sweepingplane.h"

SweepingPlane::SweepingPlane()
{
    color[0] = 0;
    color[1] = 1;
    color[2] = 1;
    linewidth = 3;
    actor->GetProperty()->SetColor(color);
    actor->GetProperty()->SetLineWidth(linewidth);
    actor->GetProperty()->SetOpacity(0.5);
    boundary = false;
}
void SweepingPlane::InputData(vtkSmartPointer<vtkPolyData> poly)
{
    if(boundary)
    {
        vtkSmartPointer<vtkExtractEdges> extract = vtkSmartPointer<vtkExtractEdges>::New();
        extract->SetInputData(poly);
        extract->Update();
        SetInput(extract->GetOutput());
    }
    else
    {
        SetInput(poly);
    }
}
