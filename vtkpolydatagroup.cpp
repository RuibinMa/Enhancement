#include "vtkpolydatagroup.h"

vtkPolyDataGroup::vtkPolyDataGroup()
{
    NumOfMembers = 0;
    head = NULL;
}
vtkPolyDataGroup::~vtkPolyDataGroup()
{
    if(head)
    {
        struct PolydataMember* p = head;
        while(head)
        {
            head = head->next;
            if(p->numofpoints!=0)
            {
                for(int i = 0; i < p->numofpoints; i++)
                {
                    free(p->points[i]);
                }
                free(p->points);
            }
            delete p;
            p = head;
        }
    }
}

void vtkPolyDataGroup::AddMember(vtkPolyData *object)
{
    if(NumOfMembers >= maxPolyDataNumber)
    {
        std::cerr<<"PolyData Group has reached the max member number"<<endl;
        return;
    }

    // fill the newmember
    struct PolydataMember* newmember = new PolydataMember;
    newmember->next = NULL;
    int num = object->GetNumberOfPoints();
    newmember->numofpoints = num;
    newmember->points = (double **) malloc(num * sizeof(double*));
    for(int i = 0; i < num; i++)
    {
        newmember->points[i] = (double*) malloc(3 * sizeof(double));
        double p[3];
        object->GetPoint(i, p);
        newmember->points[i][0] = p[0];
        newmember->points[i][1] = p[1];
        newmember->points[i][2] = p[2];
    }

    // link the newmember
    struct PolydataMember* p = head;
    if(!head)
        head = newmember;
    else
    {
        while(p->next!=NULL)
        {
            p = p->next;
        }
        p->next = newmember;
    }
    NumOfMembers++;
}
vtkSmartPointer<vtkPolyData> vtkPolyDataGroup::GetMember(int id)
{
    assert(id < NumOfMembers);
    assert(head);
    struct PolydataMember* p = head;
    for(int i = 0; i < id; i++)
    {
        p = p->next;
    }

    int num = p->numofpoints;
    vtkSmartPointer<vtkPoints> Points = vtkSmartPointer<vtkPoints>::New();
    for(int i = 0; i< num; i++)
    {
        Points->InsertNextPoint(p->points[i][0], p->points[i][1], p->points[i][2]);
    }
    vtkSmartPointer<vtkPolyData> Poly = vtkSmartPointer<vtkPolyData>::New();
    Poly->SetPoints(Points);
    vtkSmartPointer<vtkVertexGlyphFilter> vertexFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
    vertexFilter->SetInputData(Poly);
    vertexFilter->Update();
    Poly->DeepCopy(vertexFilter->GetOutput());
    return Poly;
}
