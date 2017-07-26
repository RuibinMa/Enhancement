#ifndef VTKPOLYDATAGROUP_H
#define VTKPOLYDATAGROUP_H
#include <vtkPolyData.h>
#include <vtkIdList.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkVertexGlyphFilter.h>

struct PolydataMember{ // currently, can only store the point locations. can be extended to cell information
    double ** points; // should be a Nx3 array
    int numofpoints;
    struct PolydataMember* next;
};

class vtkPolyDataGroup
{
public:
    vtkPolyDataGroup();
    ~vtkPolyDataGroup();
    void AddMember(vtkPolyData* object);
    vtkSmartPointer<vtkPolyData> GetMember(int id);
    int GetNumOfMembers(){return NumOfMembers;}

private:
    int maxPolyDataNumber = 2000;
    int NumOfMembers;
    struct PolydataMember* head;

};

#endif // VTKPOLYDATAGROUP_H
