#include "centerline.h"
#include <iostream>
#include <fstream>
#include <vtkMath.h>
#include <unistd.h>

Centerline::Centerline()
{
    actor->GetProperty()->SetColor(1,1,0);
    actor->GetProperty()->SetLineWidth(5);
}
int Centerline::GetNumberOfPoints()
{
    return model->GetNumberOfPoints();
}
vtkSmartPointer<vtkPlane> Centerline::GetVerticalPlane(vtkIdType i, double *normalDirection)
{
    if(i>GetNumberOfPoints())
    {
        std::cerr<<"index exceed the number of points("<<GetNumberOfPoints()<<")"<<endl;
        exit(1);
    }
    vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();
    double p[3], plast[3], pnext[3], direction[3];
    model->GetPoint(i,p);
    if(i == 0)
    {
        model->GetPoint(i+1 , pnext);
        vtkMath::Subtract(pnext, p, direction);
        vtkMath::Normalize(direction);
    }
    else if(i == GetNumberOfPoints()-1)
    {
        model->GetPoint(i-1 , plast);
        vtkMath::Subtract(p, plast, direction);
        vtkMath::Normalize(direction);
    }
    else
    {
        model->GetPoint(i+1, pnext);
        model->GetPoint(i-1, plast);
        double direction_next[3];
        double direction_last[3];
        vtkMath::Subtract(pnext, p, direction_next);
        vtkMath::Normalize(direction_next);
        vtkMath::Subtract(p, plast, direction_last);
        vtkMath::Normalize(direction_last);
        vtkMath::Add(direction_next, direction_last, direction);
        vtkMath::Normalize(direction);
    }
    plane->SetNormal(direction);
    plane->SetOrigin(p);
    normalDirection[0] = direction[0];
    normalDirection[1] = direction[1];
    normalDirection[2] = direction[2];
    return plane;
}
double* Centerline::Gaussian1DKernel(int len, double std)
{
    int i, ctr;
    double *G;
    G = (double *)calloc(len, sizeof(double));
    ctr = (len - 1) / 2;
    double sum;
    for(i=0; i<len; i++)
    {
        G[i] = 1/(sqrt(2 * 3.1415927) * std) * exp(- (i-ctr)*(i-ctr) / (2 * std * std));
        sum+=G[i];
    }
    for(i=0; i<len; i++)
    {
        G[i] /= sum;
    }
    //std::cout<<"len "<<len<<" ctr "<<ctr<<std::endl;
    return G;
}

double* Centerline::GaussianDerivative1DKernel(int len, double std)
{
    int i, ctr;
    double *DG;
    DG = (double *)calloc(len, sizeof(double));
    ctr = (len - 1) / 2;
    double absMax = 0;
    double temp;
    for(i=0; i<len; i++)
    {
        DG[i] = -(i-ctr) * 1/(sqrt(2 * 3.1415927) * std * std * std) * exp(- (i-ctr)*(i-ctr) / (2 * std * std));
        temp = DG[i];
        if(temp < 0) temp = -temp;
        if(absMax < temp)
            absMax = temp;
    }
    for(i=0; i<len; i++)
    {
        DG[i] /= absMax;
    }
    std::cout<<"len "<<len<<" ctr "<<ctr<<std::endl;
    return DG;
}

double* Centerline::conv(double *A, double *B, int lenA, int lenB, int *lenC)
{
    int nconv;
    int i, j, i1;
    double tmp;
    double *C;

    //allocated convolution array
    nconv = lenA+lenB-1;
    C = (double*) calloc(nconv, sizeof(double));

    //convolution process
    for (i=0; i<nconv; i++)
    {
        i1 = i;
        tmp = 0.0;
        for (j=0; j<lenB; j++)
        {
            if(i1>=0 && i1<lenA)
                tmp = tmp + (A[i1]*B[j]);

            i1 = i1-1;
            C[i] = tmp;
        }
    }
    //get length of convolution array
    (*lenC) = nconv;
    //return convolution array
    return(C);
}
void Centerline::SmoothCenterline(double std, vtkSmartPointer<vtkIntArray> ViolationNums)
{
    int halfLength = (int)round(std*4.5);
    int GaussianKernelLen = 2*halfLength + 1;
    double GaussianStd = std;
    double *G = Gaussian1DKernel(GaussianKernelLen, GaussianStd);
    //for(int i = 0; i < GaussianKernelLen; i++)
    //{
    //    std::cout<<G[i]<<std::endl;
    //}
    int NumPoints, NumSmoothedPoints;
    NumPoints = model->GetNumberOfPoints();
    double *x, *y, *z, *gx, *gy, *gz;
    x = (double *) calloc(NumPoints + GaussianKernelLen-1, sizeof(double));
    y = (double *) calloc(NumPoints + GaussianKernelLen-1, sizeof(double));
    z = (double *) calloc(NumPoints + GaussianKernelLen-1, sizeof(double));
    double p[3];
    for(vtkIdType i = 0; i< NumPoints + GaussianKernelLen - 1; i++)
    {
        if(i<(GaussianKernelLen-1)/2)
        {
            model->GetPoint(0, p);
            double pright[3];
            model->GetPoint((GaussianKernelLen - 1)/2 - i, pright);
            vtkMath::MultiplyScalar(p, 2);
            vtkMath::Subtract(p, pright, p);
        }
        else if(i >= (GaussianKernelLen-1)/2 + NumPoints)
        {
            model->GetPoint((vtkIdType)NumPoints -1, p);
            double pleft[3];
            model->GetPoint((GaussianKernelLen-1)/2 + 2*NumPoints - 2 -i, pleft);
            vtkMath::MultiplyScalar(p, 2);
            vtkMath::Subtract(p, pleft, p);
        }
        else
            model->GetPoint(i - (GaussianKernelLen-1)/2, p);
        x[i] = p[0];
        y[i] = p[1];
        z[i] = p[2];
    }
    gx = conv(x, G, NumPoints + GaussianKernelLen -1, GaussianKernelLen, &NumSmoothedPoints);
    gy = conv(y, G, NumPoints + GaussianKernelLen -1, GaussianKernelLen, &NumSmoothedPoints);
    gz = conv(z, G, NumPoints + GaussianKernelLen -1, GaussianKernelLen, &NumSmoothedPoints);
    double gp[3];
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
    for(vtkIdType i = 0; i <NumPoints; i++)
    {
        if(ViolationNums == NULL || ViolationNums->GetValue(i) > 0) // selective smoothing
        {
            gp[0] = gx[i + GaussianKernelLen -1];
            gp[1] = gy[i + GaussianKernelLen -1];
            gp[2] = gz[i + GaussianKernelLen -1];
            points->InsertNextPoint(gp);
        }
        else
        {
            double p[3];
            model->GetPoint(i, p);
            points->InsertNextPoint(p);
        }
    }
    for(int i=0; i<NumPoints-1; i++)
    {
        vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
        line->GetPointIds()->SetId(0, i);
        line->GetPointIds()->SetId(1, i+1);
        lines->InsertNextCell(line);
    }
    vtkSmartPointer<vtkPolyData> polyline = vtkSmartPointer<vtkPolyData>::New();
    polyline->SetPoints(points);
    polyline->SetLines(lines);

    model->DeepCopy(polyline);
    mapper->Update();
    actor->GetProperty()->SetPointSize(5);
    //std::cout<<"centerline smoothed:"<<model->GetNumberOfPoints()<<std::endl;
    free(x);
    free(y);
    free(z);
    free(gx);
    free(gy);
    free(gz);
    free(G);
}
void Centerline::Smooth(vtkSmartPointer<vtkDoubleArray> Candidate, double std)
{
    int halfLength = (int)round(std*4.5);
    int GaussianKernelLen = 2*halfLength + 1;
    double GaussianStd = std;
    double *G = Gaussian1DKernel(GaussianKernelLen, GaussianStd);
    //for(int i = 0; i < GaussianKernelLen; i++)
    //{
    //    std::cout<<G[i]<<std::endl;
    //}
    int NumPoints, NumSmoothedPoints;
    NumPoints = Candidate->GetNumberOfTuples();
    double *x, *y, *z, *gx, *gy, *gz;
    x = (double *) calloc(NumPoints + GaussianKernelLen-1, sizeof(double));
    y = (double *) calloc(NumPoints + GaussianKernelLen-1, sizeof(double));
    z = (double *) calloc(NumPoints + GaussianKernelLen-1, sizeof(double));
    double p[3];
    for(vtkIdType i = 0; i< NumPoints + GaussianKernelLen - 1; i++)
    {
        if(i<(GaussianKernelLen-1)/2)
        {
            Candidate->GetTuple(0, p);
            double pright[3];
            Candidate->GetTuple((GaussianKernelLen - 1)/2 - i, pright);
            vtkMath::MultiplyScalar(p, 2);
            vtkMath::Subtract(p, pright, p);
        }
        else if(i >= (GaussianKernelLen-1)/2 + NumPoints)
        {
            Candidate->GetTuple((vtkIdType)NumPoints -1, p);
            double pleft[3];
            Candidate->GetTuple((GaussianKernelLen-1)/2 + 2*NumPoints - 2 -i, pleft);
            vtkMath::MultiplyScalar(p, 2);
            vtkMath::Subtract(p, pleft, p);
        }
        else
            Candidate->GetTuple(i - (GaussianKernelLen-1)/2, p);
        x[i] = p[0];
        y[i] = p[1];
        z[i] = p[2];
    }
    gx = conv(x, G, NumPoints + GaussianKernelLen -1, GaussianKernelLen, &NumSmoothedPoints);
    gy = conv(y, G, NumPoints + GaussianKernelLen -1, GaussianKernelLen, &NumSmoothedPoints);
    gz = conv(z, G, NumPoints + GaussianKernelLen -1, GaussianKernelLen, &NumSmoothedPoints);
    //std::cout<<NumSmoothedPoints<<endl;
    double gp[3];
    for(vtkIdType i = 0; i <NumPoints; i++)
    {
        gp[0] = gx[i + GaussianKernelLen -1];
        gp[1] = gy[i + GaussianKernelLen -1];
        gp[2] = gz[i + GaussianKernelLen -1];
        Candidate->SetTuple(i, gp);
    }
    free(x);
    free(y);
    free(z);
    free(gx);
    free(gy);
    free(gz);
    free(G);
}

void Centerline::UniformSample(int resolution, vtkSmartPointer<vtkPolyData> line)
{
    vtkSmartPointer<vtkParametricSpline> spline = vtkSmartPointer<vtkParametricSpline>::New();
    if(line == NULL)
    {
        spline->SetPoints(model->GetPoints());
    }
    else
    {
        spline->SetPoints(line->GetPoints());
    }
    vtkSmartPointer<vtkParametricFunctionSource> functionSource = vtkSmartPointer<vtkParametricFunctionSource>::New();
    functionSource->SetParametricFunction(spline);
    functionSource->Update();
    functionSource->SetUResolution(resolution);
    functionSource->Update();
    if(line == NULL)
    {
        model->DeepCopy(functionSource->GetOutput());
    }
    else
    {
        line->DeepCopy(functionSource->GetOutput());
    }
}

void Centerline::splineTangent(double *tangent, vtkSmartPointer<vtkParametricSpline> spline, double t_u, double stepsize)
{
    double u[3], Du[9];
    u[0] = t_u;
    double p[3];
    spline->Evaluate(u, p, Du);
    if(u[0] < stepsize)
    {
        double pnext[3], unext[3];
        unext[0] = u[0] + stepsize;
        spline->Evaluate(unext, pnext, Du);
        vtkMath::Subtract(pnext, p, tangent);
    }
    else if(u[0] > 1 - stepsize)
    {
        double plast[3], ulast[3];
        ulast[0] = u[0] - stepsize;
        spline->Evaluate(ulast, plast, Du);
        vtkMath::Subtract(p, plast, tangent);
    }
    else
    {
        double pnext[3], plast[3], unext[3], ulast[3];
        unext[0] = u[0] + stepsize;
        spline->Evaluate(unext, pnext, Du);
        ulast[0] = u[0] - stepsize;
        spline->Evaluate(ulast, plast, Du);
        vtkMath::Subtract(pnext, plast, tangent);
    }
    vtkMath::Normalize(tangent);
}
void Centerline::splineNormal(double *normal, vtkSmartPointer<vtkParametricSpline> spline, double t_u, double stepsize, double &curvature)
{
    double u[3];
    u[0] = t_u;
    double tangent[3];
    splineTangent(tangent, spline, t_u, stepsize);
    if(u[0] < stepsize)
    {
        double tangentnext[3];
        splineTangent(tangentnext, spline, t_u + stepsize, stepsize);
        vtkMath::Subtract(tangentnext, tangent, normal);
        curvature = vtkMath::Normalize(normal) / stepsize;
    }
    else if(u[0] > 1 - stepsize)
    {
        double tangentlast[3];
        splineTangent(tangentlast, spline, t_u - stepsize, stepsize);
        vtkMath::Subtract(tangent, tangentlast, normal);
        curvature = vtkMath::Normalize(normal) / stepsize;
    }
    else
    {
        double tangentnext[3], tangentlast[3];
        splineTangent(tangentnext, spline, t_u + stepsize, stepsize);
        splineTangent(tangentlast, spline, t_u - stepsize, stepsize);
        vtkMath::Subtract(tangentnext, tangentlast, normal);
        curvature = vtkMath::Normalize(normal) / (2*stepsize);
    }
}
double Centerline::splineTorsion(vtkSmartPointer<vtkParametricSpline> spline, double t_u, double stepsize)
{
    double u[3];
    u[0] = t_u;
    double tangent[3], normal[3], binormal[3], dbinormal[3];
    double curvature, torsion;
    splineTangent(tangent, spline, t_u, stepsize);
    splineNormal(normal, spline, t_u, stepsize, curvature);
    vtkMath::Cross(tangent, normal, binormal);
    //std::cout<<vtkMath::Norm(binormal)<<endl;
    if(u[0] < stepsize)
    {
        double tangentnext[3];
        splineTangent(tangentnext, spline, t_u + stepsize, stepsize);
        double normalnext[3];
        splineNormal(normalnext, spline, t_u + stepsize, stepsize, curvature);
        double binormalnext[3];
        vtkMath::Cross(tangentnext, normalnext, binormalnext);
        vtkMath::Subtract(binormalnext, binormal, dbinormal);
        torsion = - vtkMath::Dot(normal, dbinormal) / stepsize;
    }
    else if(u[0] > 1 - stepsize)
    {
        double tangentlast[3];
        splineTangent(tangentlast, spline, t_u - stepsize, stepsize);
        double normallast[3];
        splineNormal(normallast, spline, t_u - stepsize, stepsize, curvature);
        double binormallast[3];
        vtkMath::Cross(tangentlast, normallast, binormallast);
        vtkMath::Subtract(binormal, binormallast, dbinormal);
        torsion = - vtkMath::Dot(normal, dbinormal) / stepsize;
    }
    else
    {
        double tangentlast[3], tangentnext[3], normallast[3], normalnext[3], binormallast[3], binormalnext[3];
        splineTangent(tangentnext, spline, t_u + stepsize, stepsize);
        splineNormal(normalnext, spline, t_u + stepsize, stepsize, curvature);
        vtkMath::Cross(tangentnext, normalnext, binormalnext);
        splineTangent(tangentlast, spline, t_u - stepsize, stepsize);
        splineNormal(normallast, spline, t_u - stepsize, stepsize, curvature);
        vtkMath::Cross(tangentlast, normallast, binormallast);
        vtkMath::Subtract(binormalnext, binormallast, dbinormal);
        torsion = - vtkMath::Dot(normal, dbinormal) / (2*stepsize);
    }
    return torsion;
}
void Centerline::PutNormalsOnSameSide(vtkSmartPointer<vtkDoubleArray> Normals, vtkSmartPointer<vtkDoubleArray> Curvatures)
{
    for(vtkIdType i = 1; i< Normals->GetNumberOfTuples(); i++)
    {
        double normal[3], lastnormal[3];
        Normals->GetTuple(i, normal);
        Normals->GetTuple(i-1, lastnormal);
        int lastagreement;
        lastagreement = (vtkMath::Dot(normal, lastnormal)>0)?1:0;

        if(lastagreement == 0)
        {
            vtkMath::MultiplyScalar(normal, -1.0);
            Normals->SetTuple(i, normal);
            double curvature;
            curvature = -1 * Curvatures->GetValue(i);
            Curvatures->SetValue(i, curvature);
        }
    }
}

void Centerline::GaussianTangents(vtkSmartPointer<vtkDoubleArray> Tangents, double std)
{
    int halfLength = (int)round(std*4.5);
    int GaussianKernelLen = 2*halfLength + 1;
    double GaussianStd = std;
    double *DG = GaussianDerivative1DKernel(GaussianKernelLen, GaussianStd);

    for(int i = 0; i < GaussianKernelLen; i++)
    {
        std::cout<<DG[i]<<std::endl;
    }
    int NumPoints, NumSmoothedPoints;
    NumPoints = model->GetNumberOfPoints();
    double *x, *y, *z, *gtx, *gty, *gtz;
    x = (double *) calloc(NumPoints + GaussianKernelLen-1, sizeof(double));
    y = (double *) calloc(NumPoints + GaussianKernelLen-1, sizeof(double));
    z = (double *) calloc(NumPoints + GaussianKernelLen-1, sizeof(double));
    double p[3];
    for(vtkIdType i = 0; i< NumPoints + GaussianKernelLen - 1; i++)
    {
        if(i<(GaussianKernelLen-1)/2)
            model->GetPoint(0, p);
        else if(i >= (GaussianKernelLen-1)/2 + NumPoints)
            model->GetPoint((vtkIdType)NumPoints -1, p);
        else
            model->GetPoint(i - (GaussianKernelLen-1)/2, p);
        x[i] = p[0];
        y[i] = p[1];
        z[i] = p[2];
    }
    gtx = conv(x, DG, NumPoints + GaussianKernelLen -1, GaussianKernelLen, &NumSmoothedPoints);
    gty = conv(y, DG, NumPoints + GaussianKernelLen -1, GaussianKernelLen, &NumSmoothedPoints);
    gtz = conv(z, DG, NumPoints + GaussianKernelLen -1, GaussianKernelLen, &NumSmoothedPoints);
    std::cout<<NumSmoothedPoints<<endl;
    double gt[3];
    for(vtkIdType i = 0; i <NumPoints; i++)
    {

        gt[0] = gtx[i + GaussianKernelLen -1];
        gt[1] = gty[i + GaussianKernelLen -1];
        gt[2] = gtz[i + GaussianKernelLen -1];
        vtkMath::Normalize(gt);
        std::cout<<i<<" "<<gt[0]<<" "<<gt[1]<<" "<<gt[2]<<" "<<std::endl;
        //std::cout<<i<<" "<<x[i + (GaussianKernelLen-1)/2]<<" "<<y[i + (GaussianKernelLen-1)/2]<<" "<<z[i + (GaussianKernelLen-1)/2]<<" "<<std::endl;
        //std::cout<<i<<" "<<model->GetPoint(i)[0]<<" "<<model->GetPoint(i)[1]<<" "<<model->GetPoint(i)[2]<<" "<<std::endl;

        Tangents->InsertNextTuple(gt);
    }
    free(x);
    free(y);
    free(z);
    free(gtx);
    free(gty);
    free(gtz);
    free(DG);
}

vtkSmartPointer<vtkPolyData> Centerline::Deformation(vtkSmartPointer<vtkDoubleArray> S, vtkSmartPointer<vtkDoubleArray> Curvatures,
                                                     vtkSmartPointer<vtkDoubleArray> Tangents, vtkSmartPointer<vtkDoubleArray> Normals,
                                                     vtkSmartPointer<vtkPolyData> t_colon, RenderManager *t_rendermanager,
                                                     vtkSmartPointer<vtkDoubleArray> PlaneOriginals, vtkSmartPointer<vtkDoubleArray> PlaneNormals,
                                                     vtkSmartPointer<vtkDoubleArray> RefDirections, FileManager *t_filemanager)
{
    //PutNormalsOnSameSide(Normals, Curvatures);
    std::cout<<"Deformation"<<endl;
    bool straight = true;
    // Eliminate the torsion by growing the curve on a plane, according to: -dNnew/dSnew = -k*Tnew
    double point[3], nextpoint[3];
    double tangent[3], nexttangent[3];
    double normal[3], nextnormal[3];
    double binormal[3];
    double ds, curvature;
    vtkSmartPointer<vtkDoubleArray> NewTangents = vtkSmartPointer<vtkDoubleArray>::New();
    NewTangents->SetNumberOfComponents(3);
    vtkSmartPointer<vtkDoubleArray> NewNormals = vtkSmartPointer<vtkDoubleArray>::New();
    NewNormals->SetNumberOfComponents(3);
    vtkSmartPointer<vtkDoubleArray> NewBinormals = vtkSmartPointer<vtkDoubleArray>::New();
    NewBinormals->SetNumberOfComponents(3);
    vtkSmartPointer<vtkPoints> newpoints = vtkSmartPointer<vtkPoints>::New();
    for(vtkIdType i=0; i<model->GetNumberOfPoints(); i++)
    {
        if(i == 0)
        {
            model->GetPoint(0, point);
            point[0] = point[0] + 60;
            binormal[0] = 0;
            binormal[1] = 0;
            binormal[2] = 1;
            tangent[0] = 1;
            tangent[1] = 0;
            tangent[2] = 0;
            normal[0] = 0;
            normal[1] = 1;
            normal[2] = 0;
        }
        ds = (i != model->GetNumberOfPoints()-1)?(S->GetValue(i + 1) - S->GetValue(i)):0;
        curvature = Curvatures->GetValue(i);

        // Record the New Axis System
        NewTangents->InsertNextTuple(tangent);
        NewNormals->InsertNextTuple(normal);
        NewBinormals->InsertNextTuple(binormal);

        double dnormal[3];
        if(straight)
        {
            dnormal[0] = 0; dnormal[1] = 0; dnormal[2] = 0;
        }
        else
        {
            dnormal[0] = tangent[0]; dnormal[1] = tangent[1]; dnormal[2] = tangent[2];
            double RelaxationFactor = 2;
            vtkMath::MultiplyScalar(dnormal, -curvature / RelaxationFactor);
        }
        vtkMath::Add(normal, dnormal, nextnormal);
        vtkMath::Normalize(nextnormal);
        //std::cout<<i<<"\ts="<<S->GetValue(i)<<"("<<ds<<")"<<"\tk="<<curvature<<endl;
        vtkMath::Cross(nextnormal, binormal, nexttangent);
        vtkMath::Normalize(nexttangent);
        vtkMath::MultiplyScalar(tangent, ds);
        vtkMath::Add(point, tangent, nextpoint);
        newpoints->InsertNextPoint(point);

        point[0] = nextpoint[0]; point[1] = nextpoint[1]; point[2] = nextpoint[2];
        normal[0] = nextnormal[0]; normal[1] = nextnormal[1]; normal[2] = nextnormal[2];
        tangent[0] = nexttangent[0]; tangent[1] = nexttangent[1]; tangent[2] = nexttangent[2];
    }

    vtkSmartPointer<vtkPolyData> newcenterline = vtkSmartPointer<vtkPolyData>::New();
    newcenterline->SetPoints(newpoints);
    vtkSmartPointer<vtkVertexGlyphFilter> newVertexFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
    newVertexFilter->SetInputData(newcenterline);
    newVertexFilter->Update();
    vtkSmartPointer<vtkPolyDataMapper> newMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    newMapper->SetInputConnection(newVertexFilter->GetOutputPort());
    vtkSmartPointer<vtkActor> newActor = vtkSmartPointer<vtkActor>::New();
    newActor->SetMapper(newMapper);
    t_rendermanager->renderModel(newActor);

    // Simple Deformation
    // Get the old coordinates
    /*
    vtkSmartPointer<vtkPoints> newColonPoints = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkDoubleArray> Coordinates = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkPointLocator> pointLocator = vtkSmartPointer<vtkPointLocator>::New();
    pointLocator->SetDataSet(model);
    pointLocator->BuildLocator();
    for(vtkIdType i = 0; i< t_colon->GetNumberOfPoints(); i++)
    {
        double p[3], centerp[3];
        t_colon->GetPoint(i, p);
        vtkIdType id;
        id = pointLocator->FindClosestPoint(p);

        model->GetPoint(id, centerp);
        double vector[3];
        vtkMath::Subtract(p, centerp, vector);
        double coordinate[3];
        coordinate[0] = vtkMath::Dot(vector, Tangents->GetTuple(id));
        coordinate[1] = vtkMath::Dot(vector, Normals->GetTuple(id));
        coordinate[2] = vtkMath::Dot(vector, Binormals->GetTuple(id));
        Coordinates->InsertNextTuple(coordinate);

        double newp[3], newcenterp[3];
        newcenterline->GetPoint(id, newcenterp);
        double x[3], y[3], z[3], tmp1[3], tmp2[3];
        NewTangents->GetTuple(id, x);
        NewNormals->GetTuple(id, y);
        NewBinormals->GetTuple(id, z);
        vtkMath::MultiplyScalar(x, coordinate[0]);
        vtkMath::MultiplyScalar(y, coordinate[1]);
        vtkMath::MultiplyScalar(z, coordinate[2]);
        vtkMath::Add(newcenterp, x, tmp1);
        vtkMath::Add(tmp1, y, tmp2);
        vtkMath::Add(tmp2, z, newp);

        newColonPoints->InsertNextPoint(newp);
    }

    vtkSmartPointer<vtkPolyData> newColonPoly = vtkSmartPointer<vtkPolyData>::New();
    newColonPoly->SetPoints(newColonPoints);
    vtkSmartPointer<vtkVertexGlyphFilter> newColonVertexFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
    newColonVertexFilter->SetInputData(newColonPoly);

    vtkSmartPointer<vtkPolyDataMapper> newColonMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    newColonMapper->SetInputConnection(newColonVertexFilter->GetOutputPort());
    vtkSmartPointer<vtkActor> newColonActor = vtkSmartPointer<vtkActor>::New();
    newColonActor->SetMapper(newColonMapper);
    //t_rendermanager->renderModel(newColonActor);
    return newColonPoly;
    */
    // Line Up the Cross Sections

    vtkSmartPointer<vtkCutter> cutter = vtkSmartPointer<vtkCutter>::New();
    cutter->SetInputData(t_colon);
    vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();

    vtkSmartPointer<vtkPolyDataConnectivityFilter> connectivityFilter = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
    connectivityFilter->SetInputConnection(cutter->GetOutputPort());
    connectivityFilter->SetExtractionModeToClosestPointRegion();

    //vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();

    vtkSmartPointer<vtkAppendPolyData> appendFilter = vtkSmartPointer<vtkAppendPolyData>::New();
    vtkSmartPointer<vtkCleanPolyData> cleanFilter = vtkSmartPointer<vtkCleanPolyData>::New();

    vtkSmartPointer<vtkPolyData> CutCircleLineUp = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPolyData> lastCircle = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPolyData> lastNewCircle = vtkSmartPointer<vtkPolyData>::New();

    vtkSmartPointer<vtkAppendPolyData> appendSurfaceFilter = vtkSmartPointer<vtkAppendPolyData>::New();
    vtkSmartPointer<vtkCleanPolyData> cleanSurfaceFilter = vtkSmartPointer<vtkCleanPolyData>::New();

    vtkSmartPointer<vtkPolyData> SurfaceLineUp = vtkSmartPointer<vtkPolyData>::New();

    for(vtkIdType i = 0; i<model->GetNumberOfPoints(); i++)
    {
        //vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
        double newp[3], oldp[3];
        newcenterline->GetPoint(i, newp);
        model->GetPoint(i, oldp);
        double told[3], nold[3], bold[3];
        Tangents->GetTuple(i, told);

        //Normals->GetTuple(i, nold);
        RefDirections->GetTuple(i, nold);

        vtkMath::Cross(told, nold, bold);

        plane->SetOrigin(PlaneOriginals->GetTuple(i));
        plane->SetNormal(PlaneNormals->GetTuple(i));
        cutter->SetCutFunction(plane);
        cutter->Update();

        connectivityFilter->SetClosestPoint(oldp);
        connectivityFilter->SetInputData(cutter->GetOutput());
        connectivityFilter->Update();

        vtkSmartPointer<vtkPolyData> cutCircle = vtkSmartPointer<vtkPolyData>::New();
        cutCircle = connectivityFilter->GetOutput();
        //std::cout<<i<<" "<<"cutline points:"<<cutCircle->GetNumberOfPoints()<<endl;

        if(i != 0)
        {
            for(vtkIdType j = 0; j<lastCircle->GetNumberOfPoints(); j++)
            {
                if(cutCircle->GetNumberOfPoints() >= 100)
                    break;
                connectivityFilter->SetClosestPoint(lastCircle->GetPoint(j));
                connectivityFilter->Update();
                cutCircle = connectivityFilter->GetOutput();
            }
        }
        lastCircle->DeepCopy(cutCircle);
        vtkSmartPointer<vtkPoints> newCutCircle = vtkSmartPointer<vtkPoints>::New();

        for(vtkIdType j=0; j<cutCircle->GetNumberOfPoints(); j++)
        {
            double p[3];
            cutCircle->GetPoint(j, p);
            double vector[3];
            vtkMath::Subtract(p, oldp, vector);
            double coordinate[3];
            coordinate[0] = vtkMath::Dot(vector, told);
            coordinate[1] = vtkMath::Dot(vector, nold);
            coordinate[2] = vtkMath::Dot(vector, bold);
            double pp[3];
            double vx[3], vy[3], vz[3], tmp1[3], tmp2[3];
            NewTangents->GetTuple(i, vx);
            NewNormals->GetTuple(i, vy);
            NewBinormals->GetTuple(i, vz);
            vtkMath::MultiplyScalar(vx, coordinate[0]);
            vtkMath::MultiplyScalar(vy, coordinate[1]);
            vtkMath::MultiplyScalar(vz, coordinate[2]);
            vtkMath::Add(newp, vx, tmp1);
            vtkMath::Add(tmp1, vy, tmp2);
            vtkMath::Add(tmp2, vz, pp);
            newCutCircle->InsertNextPoint(pp);

        }

        // check the new CutCircle's order
        /*
        vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
        for(vtkIdType y = 0; y < newCutCircle->GetNumberOfPoints()-1; y++)
        {
            vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
            line->GetPointIds()->SetId(0, y);
            line->GetPointIds()->SetId(1, y+1);
            lines->InsertNextCell(line);
            points->InsertNextPoint(newCutCircle->GetPoint(y));
        }
        vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
        line->GetPointIds()->SetId(0, newCutCircle->GetNumberOfPoints()-1);
        line->GetPointIds()->SetId(1, 0);
        lines->InsertNextCell(line);
        points->InsertNextPoint(newCutCircle->GetPoint(newCutCircle->GetNumberOfPoints()-1));
        vtkSmartPointer<vtkPolyData> poly = vtkSmartPointer<vtkPolyData>::New();
        poly->SetPoints(points);
        poly->SetLines(lines);
        */
        cutCircle->SetPoints(newCutCircle);

        if(i > 0)
        {
            std::cout<<"Connect Contours "<<i-1<<" and "<<i<<" "<<lastNewCircle->GetNumberOfPoints()<<"->"<<cutCircle->GetNumberOfPoints()<<endl;;
            vtkSmartPointer<vtkPolyData> surface = vtkSmartPointer<vtkPolyData>::New();
            surface->DeepCopy(ConnectTwoContours(lastNewCircle, cutCircle));
            appendSurfaceFilter->RemoveAllInputs();
            appendSurfaceFilter->AddInputData(SurfaceLineUp);
            appendSurfaceFilter->AddInputData(surface);
            appendFilter->Update();
            cleanSurfaceFilter->SetInputConnection(appendSurfaceFilter->GetOutputPort());
            cleanSurfaceFilter->Update();
            SurfaceLineUp->DeepCopy(cleanSurfaceFilter->GetOutput());
        }

        lastNewCircle->DeepCopy(cutCircle);

        appendFilter->RemoveAllInputs();
        appendFilter->AddInputData(CutCircleLineUp);
        appendFilter->AddInputData(cutCircle);
        appendFilter->Update();
        cleanFilter->SetInputConnection(appendFilter->GetOutputPort());
        cleanFilter->Update();
        CutCircleLineUp->DeepCopy(cleanFilter->GetOutput());
    }

    t_filemanager->SaveFile(CutCircleLineUp, "CutCircleLineUp.vtp");
    t_filemanager->SaveFile(SurfaceLineUp, "SurfaceLineUp.vtp");
    vtkSmartPointer<vtkPolyDataMapper> CutCircleLineUpMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    CutCircleLineUpMapper->SetInputData(CutCircleLineUp);
    CutCircleLineUpMapper->Update();
    vtkSmartPointer<vtkActor> CutCircleLineUpActor = vtkSmartPointer<vtkActor>::New();
    CutCircleLineUpActor->SetMapper(CutCircleLineUpMapper);
    CutCircleLineUpActor->GetProperty()->SetColor(1, 0, 0);
    t_rendermanager->renderModel(CutCircleLineUpActor);

    vtkSmartPointer<vtkPolyDataMapper> SurfaceLineUpMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    SurfaceLineUpMapper->SetInputData(SurfaceLineUp);
    SurfaceLineUpMapper->Update();
    vtkSmartPointer<vtkActor> SurfaceLineUpActor = vtkSmartPointer<vtkActor>::New();
    SurfaceLineUpActor->SetMapper(SurfaceLineUpMapper);
    t_rendermanager->renderModel(SurfaceLineUpActor);

    std::cout<<"Deformation End"<<endl;
    return SurfaceLineUp;
}
void Centerline::VisualizeTNB(vtkSmartPointer<vtkDoubleArray> S, vtkSmartPointer<vtkDoubleArray> Curvatures,
                              vtkSmartPointer<vtkDoubleArray> Tangents, vtkSmartPointer<vtkDoubleArray> Normals, vtkSmartPointer<vtkDoubleArray> Binormals,
                              RenderManager* t_rendermanager)
{
    // visualize the tangent, normal, binormal. The length of normal plotted is propotional to the curvature;
    vtkSmartPointer<vtkPoints> PointsT = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkPoints> PointsN = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkPoints> PointsB = vtkSmartPointer<vtkPoints>::New();
    for(vtkIdType i = 0; i< model->GetNumberOfPoints(); i++)
    {
        double tangent[3], normal[3], binormal[3];
        Tangents->GetTuple(i, tangent);
        Normals->GetTuple(i, normal);
        Binormals->GetTuple(i, binormal);
        //std::cout<<tangent[0]<<" "<<tangent[1]<<" "<<tangent[2]<<" "
        //         <<normal[0]<<" "<<normal[1]<<" "<<normal[2]<<" "
        //         <<binormal[0]<<" "<<binormal[1]<<" "<<binormal[2]<<std::endl;
        double p[3], pt[3], pn[3], pb[3];
        double ds = S->GetValue(i+1) - S->GetValue(i);
        if(i == model->GetNumberOfPoints()-1)
            ds = S->GetValue(i) - S->GetValue(i-1);
        model->GetPoint(i, p);
        vtkMath::MultiplyScalar(tangent, ds);
        vtkMath::Add(p, tangent, pt);
        PointsT->InsertNextPoint(pt);
        vtkMath::MultiplyScalar(normal, 5*Curvatures->GetValue(i) /* 5 Torsions->GetValue(i)*10*/);
        vtkMath::Add(p, normal, pn);
        PointsN->InsertNextPoint(pn);
        vtkMath::MultiplyScalar(binormal, 3);
        vtkMath::Add(p, binormal, pb);
        PointsB->InsertNextPoint(pb);
    }
    vtkSmartPointer<vtkPolyData> PolyT = vtkSmartPointer<vtkPolyData>::New();
    PolyT->SetPoints(PointsT);
    vtkSmartPointer<vtkVertexGlyphFilter> vertexFilterT = vtkSmartPointer<vtkVertexGlyphFilter>::New();
    vertexFilterT->SetInputData(PolyT);
    vertexFilterT->Update();
    vtkSmartPointer<vtkPolyDataMapper> PolyTMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    PolyTMapper->SetInputConnection(vertexFilterT->GetOutputPort());
    vtkSmartPointer<vtkActor> PolyTActor = vtkSmartPointer<vtkActor>::New();
    PolyTActor->SetMapper(PolyTMapper);
    PolyTActor->GetProperty()->SetColor(0,1,0);
    PolyTActor->GetProperty()->SetPointSize(3);

    vtkSmartPointer<vtkPolyData> PolyN = vtkSmartPointer<vtkPolyData>::New();
    PolyN->SetPoints(PointsN);
    vtkSmartPointer<vtkVertexGlyphFilter> vertexFilterN = vtkSmartPointer<vtkVertexGlyphFilter>::New();
    vertexFilterN->SetInputData(PolyN);
    vertexFilterN->Update();
    vtkSmartPointer<vtkPolyDataMapper> PolyNMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    PolyNMapper->SetInputConnection(vertexFilterN->GetOutputPort());
    vtkSmartPointer<vtkActor> PolyNActor = vtkSmartPointer<vtkActor>::New();
    PolyNActor->SetMapper(PolyNMapper);
    PolyNActor->GetProperty()->SetColor(1,0,0);
    PolyNActor->GetProperty()->SetPointSize(3);

    vtkSmartPointer<vtkPolyData> PolyB = vtkSmartPointer<vtkPolyData>::New();
    PolyB->SetPoints(PointsB);
    vtkSmartPointer<vtkVertexGlyphFilter> vertexFilterB = vtkSmartPointer<vtkVertexGlyphFilter>::New();
    vertexFilterB->SetInputData(PolyB);
    vertexFilterB->Update();
    vtkSmartPointer<vtkPolyDataMapper> PolyBMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    PolyBMapper->SetInputConnection(vertexFilterB->GetOutputPort());
    vtkSmartPointer<vtkActor> PolyBActor = vtkSmartPointer<vtkActor>::New();
    PolyBActor->SetMapper(PolyBMapper);
    PolyBActor->GetProperty()->SetColor(0,0,1);
    PolyBActor->GetProperty()->SetPointSize(3);

    t_rendermanager->renderModel(PolyTActor);
    t_rendermanager->renderModel(PolyNActor);
    t_rendermanager->renderModel(PolyBActor);
}
void Centerline::VisualizeOriginalCurve(RenderManager *t_rendermanager)
{
    vtkSmartPointer<vtkPolyData> OriginalCurve = vtkSmartPointer<vtkPolyData>::New();
    OriginalCurve->DeepCopy(model);
    vtkSmartPointer<vtkPolyDataMapper> OriginalCurveMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    OriginalCurveMapper->SetInputData(OriginalCurve);
    OriginalCurveMapper->Update();
    vtkSmartPointer<vtkActor> OriginalCurveActor = vtkSmartPointer<vtkActor>::New();
    OriginalCurveActor->SetMapper(OriginalCurveMapper);
    OriginalCurveActor->GetProperty()->SetColor(0,0,1);
    t_rendermanager->renderModel(OriginalCurveActor);
}

void Centerline::VisualizeSpoke(vtkSmartPointer<vtkPoints> CurvaturePoints, vtkSmartPointer<vtkIntArray> ViolationNums, RenderManager *t_rendermanager)
{
    vtkSmartPointer<vtkPoints> spokePoints = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> spokeLines = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
    colors->SetNumberOfComponents(3);
    unsigned char red[3] = {255, 0, 0};
    unsigned char blue[3] = {0, 0, 255};
    for(vtkIdType i=0; i<model->GetNumberOfPoints(); i++)
    {
        double centerPoint[3];
        model->GetPoint(i, centerPoint);
        double curvaturePoint[3];
        CurvaturePoints->GetPoint(i, curvaturePoint);

        spokePoints->InsertNextPoint(centerPoint);
        spokePoints->InsertNextPoint(curvaturePoint);
        vtkSmartPointer<vtkLine> spokeLine = vtkSmartPointer<vtkLine>::New();
        spokeLine->GetPointIds()->SetId(0, 2*i);
        spokeLine->GetPointIds()->SetId(1, 2*i+1);
        spokeLines->InsertNextCell(spokeLine);
        if(ViolationNums->GetValue(i) > 0)
        {
            colors->InsertNextTypedTuple(red);
        }
        else
        {
            colors->InsertNextTypedTuple(blue);
        }
    }
    vtkSmartPointer<vtkPolyData> PolyCurvaturePoints = vtkSmartPointer<vtkPolyData>::New();
    PolyCurvaturePoints->SetPoints(CurvaturePoints);
    vtkSmartPointer<vtkVertexGlyphFilter> CurvaturePointsVertexFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
    CurvaturePointsVertexFilter->SetInputData(PolyCurvaturePoints);
    CurvaturePointsVertexFilter->Update();
    vtkSmartPointer<vtkPolyData> CurvaturePointsPoly = vtkSmartPointer<vtkPolyData>::New();
    CurvaturePointsPoly->ShallowCopy(CurvaturePointsVertexFilter->GetOutput());
    CurvaturePointsPoly->GetPointData()->SetScalars(colors);

    vtkSmartPointer<vtkPolyDataMapper> CurvaturePointsMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    CurvaturePointsMapper->SetInputData(CurvaturePointsPoly);
    CurvaturePointsMapper->Update();

    vtkSmartPointer<vtkActor> CurvaturePointsActor = vtkSmartPointer<vtkActor>::New();
    CurvaturePointsActor->SetMapper(CurvaturePointsMapper);
    CurvaturePointsActor->GetProperty()->SetPointSize(3);
    t_rendermanager->renderModel(CurvaturePointsActor);

    vtkSmartPointer<vtkPolyData> spokePoly = vtkSmartPointer<vtkPolyData>::New();
    spokePoly->SetPoints(spokePoints);
    spokePoly->SetLines(spokeLines);

    vtkSmartPointer<vtkPolyDataMapper> spokeMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    spokeMapper->SetInputData(spokePoly);
    spokeMapper->Update();

    vtkSmartPointer<vtkActor> spokeActor = vtkSmartPointer<vtkActor>::New();
    spokeActor->SetMapper(spokeMapper);

    t_rendermanager->renderModel(spokeActor);
}
void VisualizePoints(vtkSmartPointer<vtkPoints> points, double r, double g, double b, float size, RenderManager *t_rendermanager)
{
    vtkSmartPointer<vtkPolyData> poly = vtkSmartPointer<vtkPolyData>::New();
    poly->SetPoints(points);
    vtkSmartPointer<vtkVertexGlyphFilter> vertexFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
    vertexFilter->SetInputData(poly);
    vertexFilter->Update();
    vtkSmartPointer<vtkPolyDataMapper> pointsMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    pointsMapper->SetInputConnection(vertexFilter->GetOutputPort());
    pointsMapper->Update();
    vtkSmartPointer<vtkActor> pointsActor = vtkSmartPointer<vtkActor>::New();
    pointsActor->SetMapper(pointsMapper);
    pointsActor->GetProperty()->SetColor(r, g, b);
    pointsActor->GetProperty()->SetPointSize(size);

    t_rendermanager->renderModel(pointsActor);
}

vtkSmartPointer<vtkPolyData> Centerline::EliminateTorsion(RenderManager* t_rendermanager, RenderManager *t_rendermanager_right, vtkSmartPointer<vtkPolyData> t_colon, FileManager *t_filemanager)
{
    bool use_spline = true;  // whether use spline
    double stepSize = 0.0005;

    VisualizeOriginalCurve(t_rendermanager);
    vtkSmartPointer<vtkPointLocator> CurvaturePointsLocator = vtkSmartPointer<vtkPointLocator>::New();
    CurvaturePointsLocator->SetDataSet(t_colon);
    CurvaturePointsLocator->BuildLocator();

    vtkSmartPointer<vtkParametricSpline> spline = vtkSmartPointer<vtkParametricSpline>::New();

    // Arrays for tangents, normals, binormals, curvatures, torsions, length parameters, unified parameters
    vtkSmartPointer<vtkDoubleArray> Tangents = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> Normals = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> Binormals = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> Curvatures = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> Torsions = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> S = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> U = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkPoints> CurvaturePoints = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkIdList> CurvaturePointIds = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkPoints> ViolationPoints = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkDoubleArray> Directions = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkIntArray> ViolationNums = vtkSmartPointer<vtkIntArray>::New();
    vtkSmartPointer<vtkAppendPolyData> appendFilter = vtkSmartPointer<vtkAppendPolyData>::New();
    vtkSmartPointer<vtkCleanPolyData> cleanFilter = vtkSmartPointer<vtkCleanPolyData>::New();
    vtkSmartPointer<vtkPolyData> IllCutCircles = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPolyData> NormalCutCircles = vtkSmartPointer<vtkPolyData>::New();

    vtkSmartPointer<vtkDoubleArray> PlaneOriginals = vtkSmartPointer<vtkDoubleArray>::New(); PlaneOriginals->SetNumberOfComponents(3); PlaneOriginals->SetNumberOfTuples(model->GetNumberOfPoints());
    vtkSmartPointer<vtkDoubleArray> PlaneNormals = vtkSmartPointer<vtkDoubleArray>::New(); PlaneNormals->SetNumberOfComponents(3); PlaneNormals->SetNumberOfTuples(model->GetNumberOfPoints());

    // Smooth Centerline
    /*
    for(int i = 0; i < 20; i++)
    {
        SmoothCenterline(3, NULL);
    }
    */


    int MaxIter = 1; int modify = 0; // if modify==1 do one more loop to visualize the effect of the very last modification
    for(int iter = 0; iter < MaxIter + modify; iter++)
    {
        int N = model->GetNumberOfPoints();
        if(iter < MaxIter)
            std::cout<<"Iteration: "<<iter<<endl;
        else
            std::cout<<"Iteration: Visualize the effect of the last modification"<<endl;
        Tangents->Reset();          Tangents->SetNumberOfComponents(3);        Tangents->SetNumberOfTuples(N);
        Normals->Reset();           Normals->SetNumberOfComponents(3);         Normals->SetNumberOfTuples(N);
        Binormals->Reset();         Binormals->SetNumberOfComponents(3);       Binormals->SetNumberOfTuples(N);
        Curvatures->Reset();                                                   Curvatures->SetNumberOfValues(N);
        Torsions->Reset();                                                     Torsions->SetNumberOfValues(N);
        S->Reset();                                                            S->SetNumberOfValues(N);
        U->Reset();                                                            U->SetNumberOfValues(N);
        CurvaturePoints->Reset();                                              CurvaturePoints->SetNumberOfPoints(N);
        CurvaturePointIds->Reset();                                            CurvaturePointIds->SetNumberOfIds(N);
        Directions->Reset();        Directions->SetNumberOfComponents(3);      Directions->SetNumberOfTuples(N);
        ViolationNums->Reset();                                                ViolationNums->SetNumberOfValues(N);
        spline->RemoveAllObservers();
        spline->SetPoints(model->GetPoints());

        // Get the length as parameter on each point
        double cumS = 0;
        S->InsertValue(0, cumS);
        for(vtkIdType i = 1; i< model->GetNumberOfPoints(); i++)
        {
            double p[3], lastp[3], distance;
            model->GetPoint(i-1, lastp);
            model->GetPoint(i, p);
            distance = sqrt(vtkMath::Distance2BetweenPoints(p, lastp));
            cumS += distance;
            S->InsertValue(i, cumS);
            //std::cout<<i<<" distance "<<distance<<endl;
        }
        std::cout<<"Curve Length:"<<S->GetValue(model->GetNumberOfPoints()-1)<<endl;
        // Get the normalized parameter, used in spline
        for(vtkIdType i = 0; i< model->GetNumberOfPoints(); i++)
        {
            double u;
            u = S->GetValue(i) / cumS;
            U->InsertValue(i, u);
        }
        // Get the tangent vector on each point
        for(vtkIdType i = 0; i< model->GetNumberOfPoints(); i++)
        {
            double tangent[3];
            if(use_spline)
            {
                splineTangent(tangent, spline, U->GetValue(i), stepSize);
            }
            else
            {
                GetVerticalPlane(i, tangent);
            }
            Tangents->InsertTuple(i, tangent);
        }
        // Get the normal vector and curvature on each point
        for(vtkIdType i = 0; i< model->GetNumberOfPoints(); i++)
        {
            double normal[3];
            double curvature;
            if(use_spline)
            {
                splineNormal(normal, spline, U->GetValue(i), stepSize, curvature);
                curvature /= cumS;
            }
            else
            {
                double tangent[3];
                Tangents->GetTuple(i, tangent);
                if(i == 0)
                {
                    double nexttangent[3];
                    Tangents->GetTuple(i + 1, nexttangent);
                    double snext, s, delta_s;
                    s = S->GetValue(i);
                    snext = S->GetValue(i+1);
                    delta_s = snext - s;
                    vtkMath::Subtract(nexttangent, tangent, normal);
                    curvature = vtkMath::Norm(normal) / delta_s;
                }
                else if(i == model->GetNumberOfPoints()-1)
                {
                    double lasttangent[3];
                    Tangents->GetTuple(i - 1, lasttangent);
                    double s, slast, delta_s;
                    s = S->GetValue(i);
                    slast = S->GetValue(i - 1);
                    delta_s = s - slast;
                    vtkMath::Subtract(tangent, lasttangent, normal);
                    curvature = vtkMath::Norm(normal) / delta_s;
                }
                else
                {
                    double lasttangent[3], nexttangent[3], normal1[3], normal2[3], tmp[3];
                    Tangents->GetTuple(i + 1, nexttangent);
                    Tangents->GetTuple(i - 1, lasttangent);
                    double snext, slast, delta_s;
                    snext = S->GetValue(i + 1);
                    slast = S->GetValue(i - 1);
                    delta_s = snext - slast;

                    vtkMath::Subtract(nexttangent, tangent, normal1);
                    //vtkMath::Normalize(normal1);
                    vtkMath::Subtract(tangent, lasttangent, normal2);
                    //vtkMath::Normalize(normal2);
                    vtkMath::Add(normal1, normal2, normal);
                    vtkMath::Subtract(nexttangent, lasttangent, tmp);
                    curvature = vtkMath::Norm(tmp) / delta_s;
                }
                vtkMath::Normalize(normal);
            }
            Normals->InsertTuple(i, normal);
            Curvatures->InsertValue(i, curvature);
        }

        // Get the binormal on each point

        for(vtkIdType i = 0; i<model->GetNumberOfPoints(); i++)
        {
            double tangent[3], normal[3], binormal[3];
            Tangents->GetTuple(i, tangent);
            Normals->GetTuple(i, normal);
            vtkMath::Cross(tangent, normal, binormal);
            vtkMath::Normalize(binormal);
            Binormals->InsertTuple(i, binormal);
        }
        // Get the torsion on each point
        for(vtkIdType i = 0; i<model->GetNumberOfPoints(); i++)
        {
            double torsion;
            if(use_spline)
            {
                torsion = splineTorsion(spline, U->GetValue(i), stepSize) / cumS;
            }
            else
            {
                double binormal[3], normal[3];
                Binormals->GetTuple(i, binormal);
                Normals->GetTuple(i, normal);
                double s;
                s = S->GetValue(i);
                if(i == 0)
                {
                    double nextbinormal[3], dbinormal[3];
                    Binormals->GetTuple(i + 1, nextbinormal);
                    double snext = S->GetValue(i + 1);
                    vtkMath::Subtract(nextbinormal, binormal, dbinormal);
                    torsion = - 1/(snext - s) * vtkMath::Dot(dbinormal, normal);
                }
                else if(i == model->GetNumberOfPoints()-1)
                {
                    double lastbinormal[3], dbinormal[3];
                    Binormals->GetTuple(i - 1, lastbinormal);
                    double slast = S->GetValue(i - 1);
                    vtkMath::Subtract(binormal, lastbinormal, dbinormal);
                    torsion = - 1/(s - slast) * vtkMath::Dot(dbinormal, normal);
                }
                else
                {
                    double nextbinormal[3], lastbinormal[3], dbinormal[3];
                    Binormals->GetTuple(i + 1, nextbinormal);
                    Binormals->GetTuple(i - 1, lastbinormal);
                    double snext, slast;
                    snext = S->GetValue(i + 1);
                    slast = S->GetValue(i - 1);
                    vtkMath::Subtract(nextbinormal, lastbinormal, dbinormal);
                    torsion = -1/(snext - slast) * vtkMath::Dot(dbinormal, normal);
                }
            }
            Torsions->InsertValue(i, torsion);
            //std::cout<<"torsion: "<<i<<" "<<torsion<<endl;
        }


        // cut circle
        vtkSmartPointer<vtkCutter> cutter = vtkSmartPointer<vtkCutter>::New();
        cutter->SetInputData(t_colon);
        vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();
        vtkSmartPointer<vtkPolyData> lastCircle = vtkSmartPointer<vtkPolyData>::New();
        vtkSmartPointer<vtkPolyData> cutline = vtkSmartPointer<vtkPolyData>::New();

        vtkSmartPointer<vtkPolyDataConnectivityFilter> connectivityFilter = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
        connectivityFilter->SetInputConnection(cutter->GetOutputPort());
        connectivityFilter->SetExtractionModeToClosestPointRegion();

        vtkSmartPointer<vtkDoubleArray> Radius = vtkSmartPointer<vtkDoubleArray>::New();

        // calculate the cut circles and curvature points

        for(vtkIdType i=0; i<model->GetNumberOfPoints(); i++)
        {
            double centerPoint[3];
            model->GetPoint(i, centerPoint);
            connectivityFilter->SetClosestPoint(centerPoint);

            plane->SetOrigin(model->GetPoint(i));
            plane->SetNormal(Tangents->GetTuple(i));

            PlaneOriginals->SetTuple(i, model->GetPoint(i)); // record
            PlaneNormals->SetTuple(i, Tangents->GetTuple(i));

            cutter->SetCutFunction(plane);
            cutter->Update();
            connectivityFilter->Update();

            cutline = connectivityFilter->GetOutput();

            std::cout<<i<<" "<<"cutline points:"<<cutline->GetNumberOfPoints()<<endl;

            if(i == 0)
            {
                cutline->DeepCopy(ReorderContour(cutline));
                UniformSample(cutline->GetNumberOfPoints()*2, cutline);
                lastCircle->DeepCopy(cutline);
            }
            else
            {
                for(vtkIdType j = 0; j < lastCircle->GetNumberOfPoints(); j++)
                {
                    if(cutline->GetNumberOfPoints() >= 80)
                        break;
                    connectivityFilter->SetClosestPoint(lastCircle->GetPoint(j));
                    connectivityFilter->Update();
                    cutline = connectivityFilter->GetOutput();
                }
                cutline->DeepCopy(ReorderContour(cutline));
                UniformSample(cutline->GetNumberOfPoints()*(int)(249866 / t_colon->GetNumberOfPoints()), cutline);
                lastCircle->DeepCopy(cutline);
            }

            /*
                if(i == 1)
                {
                vtkSmartPointer<vtkXMLPolyDataWriter> writer =
                  vtkSmartPointer<vtkXMLPolyDataWriter>::New();
                writer->SetFileName("cutcircle1.vtp");
                writer->SetInputData(cutline);
                writer->Write();
                for(vtkIdType y = 0; y < cutline->GetNumberOfPoints()-1; y++)
                {
                    double pt[3], npt[3], v[3];
                    cutline->GetPoint(y, pt);
                    cutline->GetPoint(y+1, npt);
                    //std::cout<<pt[0]<<" "<<pt[1]<<" "<<pt[2]<<endl;
                    vtkMath::Subtract(npt, pt, v);
                    std::cout<<vtkMath::Norm(v)<<endl;
                }
                exit(0);
                }

                */
            double curvaturePoint[3], p[3], v[3], normal[3], angleCos, maxScore = -INFINITY;
            int violationNum = 0, innerNum = 0;
            double direction[3] = {0,0,0};

            Normals->GetTuple(i, normal);
            for(vtkIdType j=0; j<cutline->GetNumberOfPoints(); j++)
            {
                cutline->GetPoint(j, p);
                vtkMath::Subtract(p, centerPoint, v);
                double r = vtkMath::Normalize(v);
                angleCos = vtkMath::Dot(v, normal);
                double sign = 0;
                if(angleCos > 0)
                {
                    innerNum++;
                    double Krel = angleCos * Curvatures->GetValue(j);
                    if(r*Krel >= 1)
                    {
                        violationNum++;
                        sign = 1;
                        if(iter == MaxIter - 1)
                            ViolationPoints->InsertNextPoint(p);
                    }
                    //double rr = (r < 15)? r : 15;
                    //double step = 0.1 * pow(rr, 3) * (sign * pow((Krel - 1/r),2) + Krel * Krel);
                    //vtkMath::MultiplyScalar(v, step);
                    //vtkMath::Add(direction, v, direction);
                }
                if(angleCos > maxScore)
                {
                    maxScore = angleCos;
                    curvaturePoint[0] = p[0];
                    curvaturePoint[1] = p[1];
                    curvaturePoint[2] = p[2];
                }
            }

            Normals->GetTuple(i, direction);
            double step = 0;
            if(violationNum!=0)
                step = 0.5 * 15 * 15 * Curvatures->GetValue(i);
            //step = 0.5 * pow(15, 3) * pow(Curvatures->GetValue(i), 2);
            vtkMath::MultiplyScalar(direction, step);

            /*
                if(innerNum > 0)
                {
                    vtkMath::MultiplyScalar(direction, 1/(double)innerNum);
                    //std::cout<<i<<" "<<vtkMath::Norm(direction)<<endl;
                }
                */

            Directions->InsertTuple(i, direction);
            ViolationNums->InsertValue(i, violationNum);
            CurvaturePoints->InsertPoint(i, curvaturePoint);
            vtkIdType curvaturePointid = CurvaturePointsLocator->FindClosestPoint(curvaturePoint);
            CurvaturePointIds->InsertId(i, curvaturePointid);

            double r = (sqrt(vtkMath::Distance2BetweenPoints(curvaturePoint, centerPoint)));
            Radius->InsertValue(i, r);

            /*
            if(violationNum > 0 && iter == MaxIter + modify - 1)
            {
                appendFilter->RemoveAllInputs();
                appendFilter->AddInputData(cutline);
                appendFilter->AddInputData(IllCutCircles);
                appendFilter->Update();
                cleanFilter->SetInputConnection(appendFilter->GetOutputPort());
                cleanFilter->Update();
                IllCutCircles->DeepCopy(cleanFilter->GetOutput());
            }
            else if(iter == MaxIter + modify - 1)
            {
                appendFilter->RemoveAllInputs();
                appendFilter->AddInputData(cutline);
                appendFilter->AddInputData(NormalCutCircles);
                appendFilter->Update();
                cleanFilter->SetInputConnection(appendFilter->GetOutputPort());
                cleanFilter->Update();
                NormalCutCircles->DeepCopy(cleanFilter->GetOutput());
            }
            */

        }

        if(iter == MaxIter - 1)
        {
            //VisualizeTNB(S, Curvatures, Tangents, Normals, Binormals, t_rendermanager);
            //VisualizeSpoke(CurvaturePoints, ViolationNums, t_rendermanager);
        }
        if(modify && iter < MaxIter)
        {
            for(int n = 0; n < 20; n++)
            {
                Smooth(Directions, 3);
            }
            vtkSmartPointer<vtkPoints> NewPoints = vtkSmartPointer<vtkPoints>::New();
            for(vtkIdType i = 0; i < model->GetNumberOfPoints(); i++)
            {
                double direction[3];
                Directions->GetTuple(i, direction);
                double p[3], newp[3];
                model->GetPoint(i, p);
                vtkMath::Add(p, direction, newp);
                NewPoints->InsertPoint(i, newp);
                //std::cout<<vtkMath::Norm(direction)<<endl;
            }
            model->SetPoints(NewPoints);
            //SmoothCenterline(3, NULL);
        }
        std::cout<<"Iteration Ends"<<endl;
    }

    // visualize the ill cut circles
    /*
    vtkSmartPointer<vtkPolyDataMapper> IllCutCirclesMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    IllCutCirclesMapper->SetInputData(IllCutCircles);
    IllCutCirclesMapper->Update();
    vtkSmartPointer<vtkActor> IllCutCirclesActor = vtkSmartPointer<vtkActor>::New();
    IllCutCirclesActor->SetMapper(IllCutCirclesMapper);
    IllCutCirclesActor->GetProperty()->SetColor(1, 0, 0);
    t_rendermanager->renderModel(IllCutCirclesActor);
    */
    // visualize the normal cut circles
    /*
    vtkSmartPointer<vtkPolyDataMapper> NormalCutCirclesMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    NormalCutCirclesMapper->SetInputData(NormalCutCircles);
    NormalCutCirclesMapper->Update();
    vtkSmartPointer<vtkActor> NormalCutCirclesActor = vtkSmartPointer<vtkActor>::New();
    NormalCutCirclesActor->SetMapper(NormalCutCirclesMapper);
    NormalCutCirclesActor->GetProperty()->SetColor(0, 1, 1);
    t_rendermanager->renderModel(NormalCutCirclesActor);
    */

    // visualize the violation points
    /*
    vtkSmartPointer<vtkPolyData> PolyViolationPoints = vtkSmartPointer<vtkPolyData>::New();
    PolyViolationPoints->SetPoints(ViolationPoints);
    vtkSmartPointer<vtkVertexGlyphFilter> ViolationPointsVertexFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
    ViolationPointsVertexFilter->SetInputData(PolyViolationPoints);
    ViolationPointsVertexFilter->Update();
    vtkSmartPointer<vtkPolyDataMapper> ViolationPointsMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    ViolationPointsMapper->SetInputConnection(ViolationPointsVertexFilter->GetOutputPort());
    ViolationPointsMapper->Update();
    vtkSmartPointer<vtkActor> ViolationPointsActor = vtkSmartPointer<vtkActor>::New();
    ViolationPointsActor->SetMapper(ViolationPointsMapper);
    ViolationPointsActor->GetProperty()->SetPointSize(3);
    ViolationPointsActor->GetProperty()->SetColor(1,0,0);
    t_rendermanager->renderModel(ViolationPointsActor);
    */


    // relax the Orthogonality
    if(0)
    {
    vtkSmartPointer<vtkIdList> tempViolationTail = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkIdList> tempViolationHead = vtkSmartPointer<vtkIdList>::New();
    for(vtkIdType i = 0; i < model->GetNumberOfPoints()-1; i++)
    {
        if(ViolationNums->GetValue(i)==0 && ViolationNums->GetValue(i+1)!=0)
            tempViolationHead->InsertNextId(i+1);
        else if(ViolationNums->GetValue(i)!=0 && ViolationNums->GetValue(i+1)==0)
            tempViolationTail->InsertNextId(i);
    }

    vtkSmartPointer<vtkIdList> ViolationTail = vtkSmartPointer<vtkIdList>::New();
    std::cout<<tempViolationTail->GetNumberOfIds()<<endl;
    if(tempViolationHead->GetId(0) > tempViolationTail->GetId(0))
    {
        for(vtkIdType i = 1; i<tempViolationTail->GetNumberOfIds(); i++)
        {
            ViolationTail->InsertNextId(tempViolationTail->GetId(i));
        }
    }
    else
        ViolationTail->DeepCopy(tempViolationTail);
    std::cout<<ViolationTail->GetNumberOfIds()<<endl;

    vtkSmartPointer<vtkIdList> ViolationHead = vtkSmartPointer<vtkIdList>::New();
    std::cout<<tempViolationHead->GetNumberOfIds()<<endl;
    if(tempViolationHead->GetId(tempViolationHead->GetNumberOfIds()-1) > ViolationTail->GetId(ViolationTail->GetNumberOfIds()-1))
    {
        for(vtkIdType i = 0; i<tempViolationHead->GetNumberOfIds()-1; i++)
        {
            ViolationHead->InsertNextId(tempViolationHead->GetId(i));
        }
    }
    else
        ViolationHead->DeepCopy(tempViolationHead);
    std::cout<<ViolationHead->GetNumberOfIds()<<endl;

    assert(ViolationHead->GetNumberOfIds() == ViolationTail->GetNumberOfIds());

    vtkSmartPointer<vtkCutter> newcutter = vtkSmartPointer<vtkCutter>::New();
    newcutter->SetInputData(t_colon);
    vtkSmartPointer<vtkPolyDataConnectivityFilter> newconnectivityFilter = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
    newconnectivityFilter->SetInputConnection(newcutter->GetOutputPort());
    newconnectivityFilter->SetExtractionModeToClosestPointRegion();

    //vtkSmartPointer<vtkPolyData> NewCutlines = vtkSmartPointer<vtkPolyData>::New();
    //vtkSmartPointer<vtkAppendPolyData> newappendFilter = vtkSmartPointer<vtkAppendPolyData>::New();
    //vtkSmartPointer<vtkCleanPolyData> newcleanFilter = vtkSmartPointer<vtkCleanPolyData>::New();

    for(vtkIdType i = 0; i < ViolationHead->GetNumberOfIds(); i++)
    {
        std::cout<<ViolationHead->GetId(i)<<"->"<<ViolationTail->GetId(i)<<endl;
        double t1[3], t2[3], p1[3], p2[3];
        Tangents->GetTuple(ViolationHead->GetId(i)-1, t1);
        Tangents->GetTuple(ViolationTail->GetId(i)+1, t2);
        model->GetPoint(ViolationHead->GetId(i)-1, p1);
        model->GetPoint(ViolationTail->GetId(i)+1, p2);
        double **A;
        A = (double **)malloc(2*sizeof(double *));
        A[0] = (double *)malloc(2*sizeof(double));
        A[1] = (double *)malloc(2*sizeof(double));
        A[0][0] = t1[0]; A[0][1] = t1[1];
        A[1][0] = t2[0]; A[1][1] = t2[1];
        double *x;
        x = (double *)malloc(2*sizeof(double));
        x[0] = vtkMath::Dot(t1, p1);
        x[1] = vtkMath::Dot(t2, p2);

        vtkMath::SolveLinearSystem(A, x, 2);

        double P[3];
        P[0] = x[0]; P[1] = x[1]; P[2] = 0;
        free(A[0]);
        free(A[1]);
        free(A);
        free(x);

        double sections = ViolationTail->GetId(i) - ViolationHead->GetId(i) + 2;
        for(vtkIdType j = 1; j < sections; j++)
        {
            double t[3], l, lambda1, lambda2, angle;
            l = j/sections;
            angle = vtkMath::AngleBetweenVectors(t1, t2);
            lambda1 = LorentzianInterpolationFactor(1-l, angle);
            lambda2 = LorentzianInterpolationFactor(l, angle);
            //t[0] = (sections - j)/sections * t1[0] + j/sections * t2[0];
            //t[1] = (sections - j)/sections * t1[1] + j/sections * t2[1];
            //t[2] = (sections - j)/sections * t1[2] + j/sections * t2[2];
            t[0] = lambda1 * t1[0] + lambda2 * t2[0];
            t[1] = lambda1 * t1[1] + lambda2 * t2[1];
            t[2] = lambda1 * t1[2] + lambda2 * t2[2];
            vtkMath::Normalize(t);
            vtkSmartPointer<vtkPlane> newplane = vtkSmartPointer<vtkPlane>::New();
            newplane->SetNormal(t);
            newplane->SetOrigin(P);

            PlaneOriginals->SetTuple(ViolationHead->GetId(i) + j - 1, P); // record
            PlaneNormals->SetTuple(ViolationHead->GetId(i) + j - 1, t);
            //std::cout<<"set "<<ViolationHead->GetId(i)+j-1<<endl;

            /* //visualize the modified cutcircle in yellow
            newcutter->SetCutFunction(newplane);
            newcutter->Update();
            double centerPoint[3];
            model->GetPoint(ViolationHead->GetId(i) + j - 1, centerPoint);
            newconnectivityFilter->SetClosestPoint(centerPoint);
            newconnectivityFilter->Update();
            vtkSmartPointer<vtkPolyData> newcutline = vtkSmartPointer<vtkPolyData>::New();
            newcutline = newconnectivityFilter->GetOutput();

            newappendFilter->RemoveAllInputs();
            newappendFilter->AddInputData(newcutline);
            newappendFilter->AddInputData(NewCutlines);
            newappendFilter->Update();
            newcleanFilter->SetInputConnection(newappendFilter->GetOutputPort());
            newcleanFilter->Update();
            NewCutlines->DeepCopy(newcleanFilter->GetOutput());
            */

        }
    }
}
    /* // visualize the new cutcircles in yellow
    vtkSmartPointer<vtkPolyDataMapper> NewCutlinesMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    NewCutlinesMapper->SetInputData(NewCutlines);
    NewCutlinesMapper->Update();
    vtkSmartPointer<vtkActor> NewCutlinesActor = vtkSmartPointer<vtkActor>::New();
    NewCutlinesActor->SetMapper(NewCutlinesMapper);
    NewCutlinesActor->GetProperty()->SetColor(2, 1, 0);
    t_rendermanager->renderModel(NewCutlinesActor);
    */
    // Rotation minimization frames

    vtkSmartPointer<vtkDoubleArray> RefDirections = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkPoints> EndPoints = vtkSmartPointer<vtkPoints>::New();
    EndPoints->SetNumberOfPoints(model->GetNumberOfPoints());
    EndPoints->SetPoint(0, model->GetPoint(0));
    RefDirections->SetNumberOfComponents(3);
    RefDirections->InsertNextTuple(Normals->GetTuple(0));
    for(vtkIdType i = 1; i < model->GetNumberOfPoints(); i++)
    {
        double lastx[3], x[3], v1[3], v1r[3], v1t[3];
        model->GetPoint(i-1, lastx);
        model->GetPoint(i, x);
        double lastr[3], lastt[3];
        RefDirections->GetTuple(i-1, lastr);
        PlaneNormals->GetTuple(i-1, lastt);

        vtkMath::Subtract(x, lastx, v1);
        double c1 = vtkMath::Dot(v1, v1);
        double c1r = 2/c1 * vtkMath::Dot(v1, lastr);
        double c1t = 2/c1 * vtkMath::Dot(v1, lastt);
        v1r[0] = c1r * v1[0]; v1r[1] = c1r * v1[1]; v1r[2] = c1r * v1[2];
        v1t[0] = c1t * v1[0]; v1t[1] = c1t * v1[1]; v1t[2] = c1t * v1[2];

        double rL[3], tL[3];
        vtkMath::Subtract(lastr, v1r, rL);
        vtkMath::Subtract(lastt, v1t, tL);

        double t[3], v2[3], v2r[3];
        PlaneNormals->GetTuple(i, t);
        vtkMath::Subtract(t, tL, v2);
        double c2 = vtkMath::Dot(v2, v2);
        double c2r = 2/c2 * vtkMath::Dot(v2, rL);
        v2r[0] = c2r * v2[0]; v2r[1] = c2r * v2[1]; v2r[2] = c2r * v2[2];

        double r[3];
        vtkMath::Subtract(rL, v2r, r);
        RefDirections->InsertNextTuple(r);

        // visualize the end points
        double end[3], cp[3], et[3];
        model->GetPoint(i, cp);
        RefDirections->GetTuple(i, et);
        vtkMath::MultiplyScalar(et, 2);
        vtkMath::Add(cp, et, end);
        EndPoints->SetPoint(i, end);
    }
    //VisualizeSpoke(EndPoints, ViolationNums, t_rendermanager);

    // entrance to deformation versions
    //return Deformation_v2(S, Curvatures,Tangents, Normals, t_colon, t_rendermanager, PlaneOriginals, PlaneNormals, RefDirections, t_filemanager);
    return Deformation_v3_1(S, Curvatures, CurvaturePointIds,Tangents, Normals, t_colon, t_rendermanager, t_rendermanager_right, PlaneOriginals, PlaneNormals, RefDirections, t_filemanager);
    //return Deformation_v3(S, Curvatures, CurvaturePointIds,Tangents, Normals, t_colon, t_rendermanager, PlaneOriginals, PlaneNormals, InterpolatedRefDirections, t_filemanager);
    //return t_colon;
}

void CreateCircle( const double& z, const double& radius, const int& resolution, vtkPolyData* polyData )
{
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

  points->SetNumberOfPoints( resolution );

  for( int i = 0 ; i < resolution; ++i )
    {
    double theta = vtkMath::RadiansFromDegrees(360.*i/double(resolution));
    double x = radius*cos(theta);
    double y = radius*sin(theta);
    points->SetPoint( i, x, y, z );
    }

  polyData->Initialize();
  polyData->SetPoints( points );

  vtkSmartPointer<vtkVertexGlyphFilter> vertexFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
  vertexFilter->SetInputData(polyData);
  vertexFilter->Update();
  polyData->DeepCopy(vertexFilter->GetOutput());
}
vtkSmartPointer<vtkPolyData> Centerline::ReorderContour(vtkSmartPointer<vtkPolyData> cutCircle)
{
    int * Ids = (int *)malloc(cutCircle->GetNumberOfPoints()*sizeof(int));

    //std::cout<<"c,p: "<<cutCircle->GetNumberOfCells()<<" "<<cutCircle->GetNumberOfPoints()<<endl;
    for(vtkIdType i=0; i<cutCircle->GetNumberOfPoints(); i++)
    {
        Ids[i] = -1;
    }

    for(vtkIdType i=0; i<cutCircle->GetNumberOfCells(); i++)
    {
        vtkSmartPointer<vtkIdList> ids = vtkSmartPointer<vtkIdList>::New();
        cutCircle->GetCellPoints(i, ids);

        Ids[ids->GetId(0)]= ids->GetId(1); // assume the id pairs are connected head to tail, in order

    }

    if(cutCircle->GetNumberOfCells() < cutCircle->GetNumberOfPoints())
    {
        int * connected = (int*)malloc(cutCircle->GetNumberOfPoints() * sizeof(int));
        vtkIdType head, tail;
        for(vtkIdType i=0; i < cutCircle->GetNumberOfPoints(); i++)
        {
            connected[i] = 0;
            if(Ids[i] < 0)
            {
                head = i;
            }
        }
        for(vtkIdType i=0; i < cutCircle->GetNumberOfPoints(); i++)
        {
            if(Ids[i] >= 0)
            {
                connected[Ids[i]] = 1;
            }
        }
        for(vtkIdType i=0; i < cutCircle->GetNumberOfPoints(); i++)
        {
            if(connected[i] == 0)
            {
                tail = i;
                break;
            }
        }
        Ids[head] = tail;
        //std::cout<<"connected "<<head<<" to "<<tail<<endl;
        free(connected);
    }

    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();

    vtkIdType lastId = 0;
    for(vtkIdType i=0; i<cutCircle->GetNumberOfPoints()-1; i++)
    {
        points->InsertNextPoint(cutCircle->GetPoint(lastId));

        //std::cout<<i<<"   "<<lastId<<"->";
        lastId = Ids[lastId];
        //std::cout<<lastId<<endl;

        vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
        line->GetPointIds()->SetId(0, i);
        line->GetPointIds()->SetId(1, i+1);
        lines->InsertNextCell(line);
    }
    points->InsertNextPoint(cutCircle->GetPoint(lastId));
    vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
    line->GetPointIds()->SetId(0, cutCircle->GetNumberOfPoints()-1);
    line->GetPointIds()->SetId(1, 0);
    lines->InsertNextCell(line);

    vtkSmartPointer<vtkPolyData> ReorderedCircle = vtkSmartPointer<vtkPolyData>::New();

    ReorderedCircle->SetPoints(points);
    ReorderedCircle->SetLines(lines);

    //std::cout<<cutCircle->GetNumberOfPoints()<<"->"<<ReorderedCircle->GetNumberOfPoints()<<endl;
    free(Ids);
    return ReorderedCircle;
}
vtkSmartPointer<vtkPolyData> Centerline::FormPlate(vtkSmartPointer<vtkPolyData> cutCircle)
{
    // assume the cutcircle is processed by ReorderContour
    int N = cutCircle->GetNumberOfPoints();
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();
    double cp[3] = {0};
    for(int i = 0;  i<N; i++)
    {
        double p[3];
        cutCircle->GetPoint(i, p);
        points->InsertNextPoint(p);
        cp[0] += p[0]; cp[1] += p[1]; cp[2] += p[2];
        vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
        triangle->GetPointIds()->SetId(0, i);
        triangle->GetPointIds()->SetId(1, (i+1)%N);
        triangle->GetPointIds()->SetId(2, N);
        triangles->InsertNextCell(triangle);
    }
    cp[0] /= N; cp[1] /= N; cp[2] /= N;
    points->InsertNextPoint(cp);
    vtkSmartPointer<vtkPolyData> plate = vtkSmartPointer<vtkPolyData>::New();
    plate->SetPoints(points);
    plate->SetPolys(triangles);




    vtkSmartPointer<vtkPoints> testedgepoints = vtkSmartPointer<vtkPoints>::New();
    double p1[3] = {0,0,-1}, p2[3] = {0,0,1};
    testedgepoints->InsertNextPoint(p1);
    testedgepoints->InsertNextPoint(p2);
    vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
    line->GetPointIds()->SetId(0, 0);
    line->GetPointIds()->SetId(1,1);
    vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
    lines->InsertNextCell(line);
    vtkSmartPointer<vtkPolyData> testedge = vtkSmartPointer<vtkPolyData>::New();
    testedge->SetPoints(testedgepoints);
    testedge->SetLines(lines);
    vtkSmartPointer<vtkIntersectionPolyDataFilter> interfilter = vtkSmartPointer<vtkIntersectionPolyDataFilter>::New();
    interfilter->AddInputData(testedge);
    interfilter->AddInputData(plate);
    interfilter->GetOutput();
    vtkSmartPointer<vtkPolyData> inter = vtkSmartPointer<vtkPolyData>::New();


    return plate;
}

double Centerline::SinglePath(double **costVrt, double **costHrz, int sx, int sy, int tx, int ty, int *steps)
{
    int m = tx - sx;
    int n = ty - sy;
    double** cost = (double**)malloc((m+1)*sizeof(double *));
    for(int i=0; i< m+1; i++)
    {
        cost[i] = (double*)malloc((n+1)*sizeof(double));
    }
    int** step = (int**)malloc((m+1)*sizeof(int*));
    for(int i=0; i<m+1; i++)
    {
        step[i] = (int*)malloc((n+1)*sizeof(int));
    }
    cost[0][0] = 0;
    step[0][0] = -1; // from nowhere
    for(int i=1; i<m+1; i++)
    {
        cost[i][0] = costVrt[sx + i-1][sy] + cost[i-1][0];
        step[i][0] = 0; // from top
    }
    for(int i=1; i<n+1; i++)
    {
        cost[0][i] = costHrz[sx][sy + i-1] + cost[0][i-1];
        step[0][i] = 1; // from left
    }
    for(int i=1; i<m+1; i++)
    {
        for(int j=1; j<n+1; j++)
        {
            double c1 = costVrt[sx + i-1][sy + j] + cost[i-1][j];
            double c2 = costHrz[sx + i][sy + j-1] + cost[i][j-1];
            if(c1<c2)
            {
                cost[i][j] = c1;
                step[i][j] = 0; // from top
            }
            else
            {
                cost[i][j] = c2;
                step[i][j] = 1; // from left
            }
        }
    }

    int x = m, y = n;
    for(int i=0; i < m+n; i++)
    {
        //std::cout<<"x:"<<x<<"  y:"<<y<<" "<<step[x][y]<<endl;
        steps[m+n-1-i] = step[x][y];
        if(step[x][y] == 0)
            x--;
        else if(step[x][y] == 1)
            y--;
    }
    double result = cost[m][n];
    for(int i=0; i< m+1; i++)
    {
        free(cost[i]);
    }
    free(cost);
    for(int i=0; i<m+1; i++)
    {
        free(step[i]);
    }
    free(step);
    return result;
}

void Centerline::ContoursToSurface(RenderManager* t_rendermanager, FileManager *t_filemanager)
{
    std::cout<<"Contours To Surface"<<endl;
    vtkSmartPointer<vtkPolyData> circle1 = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPolyData> circle2 = vtkSmartPointer<vtkPolyData>::New();

    t_filemanager->LoadNewFile("cutcircle0.vtp");
    circle1->DeepCopy( t_filemanager->getfile());
    t_filemanager->LoadNewFile("cutcircle1.vtp");
    circle2->DeepCopy( t_filemanager->getfile());
    std::cout<<circle1->GetNumberOfPoints()<<" "<<circle2->GetNumberOfPoints()<<endl;
    vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
    transform->Translate(0, 0, 3);
    transform->Update();
    vtkSmartPointer<vtkTransformFilter> transformFilter = vtkSmartPointer<vtkTransformFilter>::New();
    transformFilter->SetInputData(circle2);
    transformFilter->SetTransform(transform);
    transformFilter->Update();
    circle2->DeepCopy(transformFilter->GetOutput());
    circle1->DeepCopy(ReorderContour(circle1));
    circle2->DeepCopy(ReorderContour(circle2));

    vtkSmartPointer<vtkPolyDataMapper> circleMapper1 = vtkSmartPointer<vtkPolyDataMapper>::New();
    circleMapper1->SetInputData(circle1);
    circleMapper1->Update();
    vtkSmartPointer<vtkActor> circleActor1 = vtkSmartPointer<vtkActor>::New();
    circleActor1->SetMapper(circleMapper1);
    circleActor1->GetProperty()->SetColor(0, 1, 0);
    circleActor1->GetProperty()->SetRepresentationToPoints();
    t_rendermanager->renderModel(circleActor1);

    vtkSmartPointer<vtkPolyDataMapper> circleMapper2 = vtkSmartPointer<vtkPolyDataMapper>::New();
    circleMapper2->SetInputData(circle2);
    circleMapper2->Update();
    vtkSmartPointer<vtkActor> circleActor2 = vtkSmartPointer<vtkActor>::New();
    circleActor2->SetMapper(circleMapper2);
    circleActor2->GetProperty()->SetColor(1, 0, 0);
    circleActor2->GetProperty()->SetRepresentationToPoints();
    t_rendermanager->renderModel(circleActor2);

    // connect two neighboring contours
    // calculate the cost array
    int m = circle1->GetNumberOfPoints();
    int n = circle2->GetNumberOfPoints();
    std::cout<<m<<"           "<<n<<endl;
    double** costVrt = (double**)malloc((2*m)*sizeof(double *));
    for(int i=0; i<2*m; i++)
    {
        costVrt[i] = (double*)malloc((n+1)*sizeof(double));
    }
    double** costHrz = (double**)malloc((2*m+1)*sizeof(double *));
    for(int i=0; i<2*m+1; i++)
    {
        costHrz[i] = (double*)malloc(n*sizeof(double));
    }

    for(int i=0; i<2*m; i++)
    {
        for(int j=0; j<=n; j++)
        {
            vtkIdType l = i%m, r = (i+1)%m, down = j%n;
            double a, b, c;
            double pl[3], pr[3], pdown[3];
            circle1->GetPoint(l, pl);
            circle1->GetPoint(r, pr);
            circle2->GetPoint(down, pdown);
            a = sqrt(vtkMath::Distance2BetweenPoints(pl, pdown));
            b = sqrt(vtkMath::Distance2BetweenPoints(pr, pdown));
            c = sqrt(vtkMath::Distance2BetweenPoints(pl, pr));
            double s = (a+b+c)/2;
            costVrt[i][j] = sqrt(s*(s-a)*(s-b)*(s-c));
            //std::cout<<i<<" "<<j<<" "<<costVrt[i][j]<<endl;
        }
    }
    for(int i=0; i<=2*m; i++)
    {
        for(int j=0; j<n; j++)
        {
            vtkIdType up = i%m, l = j%n, r = (j+1)%n;
            double a, b, c;
            double pup[3], pl[3], pr[3];
            circle1->GetPoint(up, pup);
            circle2->GetPoint(l, pl);
            circle2->GetPoint(r, pr);
            a = sqrt(vtkMath::Distance2BetweenPoints(pl, pup));
            b = sqrt(vtkMath::Distance2BetweenPoints(pr, pup));
            c = sqrt(vtkMath::Distance2BetweenPoints(pl, pr));
            double s = (a+b+c)/2;
            costHrz[i][j] = sqrt(s*(s-a)*(s-b)*(s-c));
            //std::cout<<i<<" "<<j<<" "<<s*(s-a)*(s-b)*(s-c)<<endl;
        }
    }
    std::cout<<"Cost Array Calculated"<<endl;

    int** steps = (int**)malloc(m*sizeof(int*));
    for(int i = 0; i < m; i++)
    {
        steps[i] = (int *)malloc((m+n)*sizeof(int));
    }
    double * costs = (double*)malloc(m*sizeof(double));

    double minCost = INFINITY;
    int minx = 0;
    for(int x = 0; x < m; x++)
    {
        costs[x] = SinglePath(costVrt, costHrz, x, 0, x+m, n, steps[x]);
        //std::cout<<"x = "<<x<<" cost = "<<costs[x]<<endl;
        if(costs[x] < minCost)
        {
            minCost = costs[x];
            minx = x;
        }
    }
    std::cout<<"minCost = "<<minCost<<" at x="<<minx<<endl;

    vtkSmartPointer<vtkPoints> surfacePoints = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> surfaceTriangles = vtkSmartPointer<vtkCellArray>::New();

    for(vtkIdType i = 0; i < m; i++)
    {
        double p[3];
        circle1->GetPoint(i, p);
        surfacePoints->InsertNextPoint(p);
    }
    for(vtkIdType i = 0; i < n; i++)
    {
        double p[3];
        circle2->GetPoint(i, p);
        surfacePoints->InsertNextPoint(p);
    }


    vtkIdType count = 0;
    vtkIdType x, y;
    x = minx; y = 0;

    for(int i=0; i < m+n; i++)
    {
        //std::cout<<i<<"step  "<<steps[i]<<"    "<<x<<","<<y<<endl;
        vtkIdType id1, id2, id3;
        id1 = x%m + count;
        id2 = y%n + count + m;
        if(steps[minx][i] == 0)
        {
            id3 = (x + 1)%m + count;
            //std::cout<<m<<" "<<n<<" "<<x%m<<" "<<y%n<<" "<<(x+1)%m<<endl;
            x++;
        }
        else
        {
            id3 = (y + 1)%n + count + m;
            //std::cout<<m<<" "<<n<<" "<<x%m<<" "<<y%n<<" "<<(y+1)%m<<endl;
            y++;
        }
        vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
        triangle->GetPointIds()->SetId(0, id1);
        triangle->GetPointIds()->SetId(1, id2);
        triangle->GetPointIds()->SetId(2, id3);
        surfaceTriangles->InsertNextCell(triangle);
    }
    vtkSmartPointer<vtkPolyData> surface = vtkSmartPointer<vtkPolyData>::New();
    surface->SetPoints(surfacePoints);
    surface->SetPolys(surfaceTriangles);
    vtkSmartPointer<vtkPolyDataMapper> surfaceMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    surfaceMapper->SetInputData(surface);
    surfaceMapper->Update();
    vtkSmartPointer<vtkActor> surfaceActor = vtkSmartPointer<vtkActor>::New();
    surfaceActor->SetMapper(surfaceMapper);
    //surfaceActor->GetProperty()->SetRepresentationToWireframe();
    t_rendermanager->renderModel(surfaceActor);

    // free the cost array
    for(int i=0; i<2*m; i++)
    {
        free(costVrt[i]);
    }
    free(costVrt);
    for(int i=0; i<2*m+1; i++)
    {
        free(costHrz[i]);
    }
    free(costHrz);
    for(int i = 0; i < m; i++)
    {
        free(steps[i]);
    }
    free(steps);
    free(costs);
}
vtkSmartPointer<vtkPolyData> Centerline::ConnectTwoContours(vtkSmartPointer<vtkPolyData> circle1, vtkSmartPointer<vtkPolyData> circle2)
{

    circle1->DeepCopy(ReorderContour(circle1));
    circle2->DeepCopy(ReorderContour(circle2));

    // connect two neighboring contours
    // calculate the cost array
    int m = circle1->GetNumberOfPoints();
    int n = circle2->GetNumberOfPoints();
    double** costVrt = (double**)malloc((2*m)*sizeof(double *));
    for(int i=0; i<2*m; i++)
    {
        costVrt[i] = (double*)malloc((n+1)*sizeof(double));
    }
    double** costHrz = (double**)malloc((2*m+1)*sizeof(double *));
    for(int i=0; i<2*m+1; i++)
    {
        costHrz[i] = (double*)malloc(n*sizeof(double));
    }

    for(int i=0; i<2*m; i++)
    {
        for(int j=0; j<=n; j++)
        {
            vtkIdType l = i%m, r = (i+1)%m, down = j%n;
            double a, b, c;
            double pl[3], pr[3], pdown[3];
            circle1->GetPoint(l, pl);
            circle1->GetPoint(r, pr);
            circle2->GetPoint(down, pdown);
            a = sqrt(vtkMath::Distance2BetweenPoints(pl, pdown));
            b = sqrt(vtkMath::Distance2BetweenPoints(pr, pdown));
            c = sqrt(vtkMath::Distance2BetweenPoints(pl, pr));
            double s = (a+b+c)/2;
            costVrt[i][j] = sqrt(s*(s-a)*(s-b)*(s-c));
            //std::cout<<i<<" "<<j<<" "<<costVrt[i][j]<<endl;
        }
    }
    for(int i=0; i<=2*m; i++)
    {
        for(int j=0; j<n; j++)
        {
            vtkIdType up = i%m, l = j%n, r = (j+1)%n;
            double a, b, c;
            double pup[3], pl[3], pr[3];
            circle1->GetPoint(up, pup);
            circle2->GetPoint(l, pl);
            circle2->GetPoint(r, pr);
            a = sqrt(vtkMath::Distance2BetweenPoints(pl, pup));
            b = sqrt(vtkMath::Distance2BetweenPoints(pr, pup));
            c = sqrt(vtkMath::Distance2BetweenPoints(pl, pr));
            double s = (a+b+c)/2;
            costHrz[i][j] = sqrt(s*(s-a)*(s-b)*(s-c));
            //std::cout<<i<<" "<<j<<" "<<s*(s-a)*(s-b)*(s-c)<<endl;
        }
    }

    int** steps = (int**)malloc(m*sizeof(int*));
    for(int i = 0; i < m; i++)
    {
        steps[i] = (int *)malloc((m+n)*sizeof(int));
    }
    double * costs = (double*)malloc(m*sizeof(double));

    double minCost = INFINITY;
    int minx = 0;
    for(int x = 0; x < m; x++)
    {
        costs[x] = SinglePath(costVrt, costHrz, x, 0, x+m, n, steps[x]);
        //std::cout<<"x = "<<x<<" cost = "<<costs[x]<<endl;
        if(costs[x] < minCost)
        {
            minCost = costs[x];
            minx = x;
        }
    }
    std::cout<<"minCost = "<<minCost<<" at x="<<minx<<endl;

    vtkSmartPointer<vtkPoints> surfacePoints = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> surfaceTriangles = vtkSmartPointer<vtkCellArray>::New();

    for(vtkIdType i = 0; i < m; i++)
    {
        double p[3];
        circle1->GetPoint(i, p);
        surfacePoints->InsertNextPoint(p);
    }
    for(vtkIdType i = 0; i < n; i++)
    {
        double p[3];
        circle2->GetPoint(i, p);
        surfacePoints->InsertNextPoint(p);
    }

    vtkIdType count = 0;
    vtkIdType x, y;
    x = minx; y = 0;

    for(int i=0; i < m+n; i++)
    {
        //std::cout<<i<<"step  "<<steps[i]<<"    "<<x<<","<<y<<endl;
        vtkIdType id1, id2, id3;
        id1 = x%m + count;
        id2 = y%n + count + m;
        if(steps[minx][i] == 0)
        {
            id3 = (x + 1)%m + count;
            //std::cout<<m<<" "<<n<<" "<<x%m<<" "<<y%n<<" "<<(x+1)%m<<endl;
            x++;
        }
        else
        {
            id3 = (y + 1)%n + count + m;
            //std::cout<<m<<" "<<n<<" "<<x%m<<" "<<y%n<<" "<<(y+1)%m<<endl;
            y++;
        }
        vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
        triangle->GetPointIds()->SetId(0, id1);
        triangle->GetPointIds()->SetId(1, id2);
        triangle->GetPointIds()->SetId(2, id3);
        surfaceTriangles->InsertNextCell(triangle);
    }
    vtkSmartPointer<vtkPolyData> surface = vtkSmartPointer<vtkPolyData>::New();
    surface->SetPoints(surfacePoints);
    surface->SetPolys(surfaceTriangles);

    // free the cost array
    for(int i=0; i<2*m; i++)
    {
        free(costVrt[i]);
    }
    free(costVrt);
    for(int i=0; i<2*m+1; i++)
    {
        free(costHrz[i]);
    }
    free(costHrz);
    for(int i = 0; i < m; i++)
    {
        free(steps[i]);
    }
    free(steps);
    free(costs);

    return surface;
}
void Centerline::ConnectTwoContoursTest(RenderManager *t_rendermanager, FileManager *t_filemanager)
{
    std::cout<<"Connect Two Contours Test"<<endl;
    vtkSmartPointer<vtkPolyData> circle1 = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPolyData> circle2 = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPolyData> circle3 = vtkSmartPointer<vtkPolyData>::New();

    t_filemanager->LoadNewFile("cutcircle0.vtp");
    circle1->DeepCopy( t_filemanager->getfile());
    t_filemanager->LoadNewFile("cutcircle1.vtp");
    circle2->DeepCopy( t_filemanager->getfile());
    t_filemanager->LoadNewFile("cutcircle0.vtp");
    circle3->DeepCopy( t_filemanager->getfile());

    vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
    transform->Translate(0, 0, 3);
    transform->Update();
    vtkSmartPointer<vtkTransformFilter> transformFilter = vtkSmartPointer<vtkTransformFilter>::New();
    transformFilter->SetInputData(circle2);
    transformFilter->SetTransform(transform);
    transformFilter->Update();
    circle2->DeepCopy(transformFilter->GetOutput());

    transform->Translate(0, 0, 6);
    transform->Update();
    transformFilter->RemoveAllInputs();
    transformFilter->SetInputData(circle3);
    transformFilter->Update();
    circle3->DeepCopy(transformFilter->GetOutput());

    circle1->DeepCopy(ReorderContour(circle1));
    circle2->DeepCopy(ReorderContour(circle2));
    circle3->DeepCopy(ReorderContour(circle3));

    vtkSmartPointer<vtkPolyDataMapper> circleMapper1 = vtkSmartPointer<vtkPolyDataMapper>::New();
    circleMapper1->SetInputData(circle1);
    circleMapper1->Update();
    vtkSmartPointer<vtkActor> circleActor1 = vtkSmartPointer<vtkActor>::New();
    circleActor1->SetMapper(circleMapper1);
    circleActor1->GetProperty()->SetColor(0, 1, 0);
    circleActor1->GetProperty()->SetRepresentationToPoints();
    //t_rendermanager->renderModel(circleActor1);

    vtkSmartPointer<vtkPolyDataMapper> circleMapper2 = vtkSmartPointer<vtkPolyDataMapper>::New();
    circleMapper2->SetInputData(circle2);
    circleMapper2->Update();
    vtkSmartPointer<vtkActor> circleActor2 = vtkSmartPointer<vtkActor>::New();
    circleActor2->SetMapper(circleMapper2);
    circleActor2->GetProperty()->SetColor(1, 0, 0);
    circleActor2->GetProperty()->SetRepresentationToPoints();
    //t_rendermanager->renderModel(circleActor2);

    vtkSmartPointer<vtkPolyDataMapper> circleMapper3 = vtkSmartPointer<vtkPolyDataMapper>::New();
    circleMapper3->SetInputData(circle3);
    circleMapper3->Update();
    vtkSmartPointer<vtkActor> circleActor3 = vtkSmartPointer<vtkActor>::New();
    circleActor3->SetMapper(circleMapper3);
    circleActor3->GetProperty()->SetColor(0, 0, 1);
    circleActor3->GetProperty()->SetRepresentationToPoints();
    //t_rendermanager->renderModel(circleActor3);

    vtkSmartPointer<vtkPolyData> surface1 = vtkSmartPointer<vtkPolyData>::New();
    surface1->DeepCopy(ConnectTwoContours(circle1, circle2));

    vtkSmartPointer<vtkPolyData> surface2 = vtkSmartPointer<vtkPolyData>::New();
    surface2->DeepCopy(ConnectTwoContours(circle2, circle3));

    std::cout<<circle1->GetNumberOfPoints()<<" "<<circle2->GetNumberOfPoints()<<" "<<circle3->GetNumberOfPoints()<<endl;
    std::cout<<surface1->GetNumberOfPoints()<<" "<<surface2->GetNumberOfPoints()<<endl;
    vtkSmartPointer<vtkAppendPolyData> appendFilter = vtkSmartPointer<vtkAppendPolyData>::New();
    appendFilter->AddInputData(surface1);
    appendFilter->AddInputData(surface2);
    appendFilter->Update();
    vtkSmartPointer<vtkCleanPolyData> cleanFilter = vtkSmartPointer<vtkCleanPolyData>::New();
    cleanFilter->SetInputConnection(appendFilter->GetOutputPort());
    cleanFilter->Update();
    std::cout<<cleanFilter->GetOutput()->GetNumberOfPoints()<<endl;
    vtkSmartPointer<vtkPolyDataMapper> surfaceMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    surfaceMapper->SetInputConnection(cleanFilter->GetOutputPort());
    surfaceMapper->Update();
    vtkSmartPointer<vtkActor> surfaceActor = vtkSmartPointer<vtkActor>::New();
    surfaceActor->SetMapper(surfaceMapper);
    t_rendermanager->renderModel(surfaceActor);
}

vtkSmartPointer<vtkPolyData> Centerline::Deformation_v2(vtkSmartPointer<vtkDoubleArray> S, vtkSmartPointer<vtkDoubleArray> Curvatures,
                                                     vtkSmartPointer<vtkDoubleArray> Tangents, vtkSmartPointer<vtkDoubleArray> Normals,
                                                     vtkSmartPointer<vtkPolyData> t_colon, RenderManager *t_rendermanager,
                                                     vtkSmartPointer<vtkDoubleArray> PlaneOriginals, vtkSmartPointer<vtkDoubleArray> PlaneNormals,
                                                     vtkSmartPointer<vtkDoubleArray> RefDirections, FileManager *t_filemanager)
{
    //PutNormalsOnSameSide(Normals, Curvatures);
    std::cout<<"Deformation"<<endl;
    bool straight = true;
    // Eliminate the torsion by growing the curve on a plane, according to: -dNnew/dSnew = -k*Tnew
    double point[3], nextpoint[3];
    double tangent[3], nexttangent[3];
    double normal[3], nextnormal[3];
    double binormal[3];
    double ds, curvature;
    vtkSmartPointer<vtkDoubleArray> NewTangents = vtkSmartPointer<vtkDoubleArray>::New();
    NewTangents->SetNumberOfComponents(3);
    vtkSmartPointer<vtkDoubleArray> NewNormals = vtkSmartPointer<vtkDoubleArray>::New();
    NewNormals->SetNumberOfComponents(3);
    vtkSmartPointer<vtkDoubleArray> NewBinormals = vtkSmartPointer<vtkDoubleArray>::New();
    NewBinormals->SetNumberOfComponents(3);
    vtkSmartPointer<vtkPoints> newpoints = vtkSmartPointer<vtkPoints>::New();
    for(vtkIdType i=0; i<model->GetNumberOfPoints(); i++)
    {
        if(i == 0)
        {
            model->GetPoint(0, point);
            point[0] = point[0] + 60;
            binormal[0] = 0;
            binormal[1] = 0;
            binormal[2] = 1;
            tangent[0] = 1;
            tangent[1] = 0;
            tangent[2] = 0;
            normal[0] = 0;
            normal[1] = 1;
            normal[2] = 0;
        }
        ds = (i != model->GetNumberOfPoints()-1)?(S->GetValue(i + 1) - S->GetValue(i)):0;
        curvature = Curvatures->GetValue(i);

        // Record the New Axis System
        NewTangents->InsertNextTuple(tangent);
        NewNormals->InsertNextTuple(normal);
        NewBinormals->InsertNextTuple(binormal);

        double dnormal[3];
        if(straight)
        {
            dnormal[0] = 0; dnormal[1] = 0; dnormal[2] = 0;
        }
        else
        {
            dnormal[0] = tangent[0]; dnormal[1] = tangent[1]; dnormal[2] = tangent[2];
            double RelaxationFactor = 2;
            vtkMath::MultiplyScalar(dnormal, -curvature / RelaxationFactor);
        }
        vtkMath::Add(normal, dnormal, nextnormal);
        vtkMath::Normalize(nextnormal);
        //std::cout<<i<<"\ts="<<S->GetValue(i)<<"("<<ds<<")"<<"\tk="<<curvature<<endl;
        vtkMath::Cross(nextnormal, binormal, nexttangent);
        vtkMath::Normalize(nexttangent);
        vtkMath::MultiplyScalar(tangent, ds);
        vtkMath::Add(point, tangent, nextpoint);
        newpoints->InsertNextPoint(point);

        point[0] = nextpoint[0]; point[1] = nextpoint[1]; point[2] = nextpoint[2];
        normal[0] = nextnormal[0]; normal[1] = nextnormal[1]; normal[2] = nextnormal[2];
        tangent[0] = nexttangent[0]; tangent[1] = nexttangent[1]; tangent[2] = nexttangent[2];
    }

    vtkSmartPointer<vtkPolyData> newcenterline = vtkSmartPointer<vtkPolyData>::New();
    newcenterline->SetPoints(newpoints);
    vtkSmartPointer<vtkVertexGlyphFilter> newVertexFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
    newVertexFilter->SetInputData(newcenterline);
    newVertexFilter->Update();
    vtkSmartPointer<vtkPolyDataMapper> newMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    newMapper->SetInputConnection(newVertexFilter->GetOutputPort());
    vtkSmartPointer<vtkActor> newActor = vtkSmartPointer<vtkActor>::New();
    newActor->SetMapper(newMapper);
    t_rendermanager->renderModel(newActor);

    // Line Up the Cross Sections
    vtkSmartPointer<vtkCutter> cutter = vtkSmartPointer<vtkCutter>::New();
    cutter->SetInputData(t_colon);
    vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();

    vtkSmartPointer<vtkPolyDataConnectivityFilter> connectivityFilter = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
    connectivityFilter->SetInputConnection(cutter->GetOutputPort());
    connectivityFilter->SetExtractionModeToClosestPointRegion();

    //vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();

    vtkSmartPointer<vtkAppendPolyData> appendFilter = vtkSmartPointer<vtkAppendPolyData>::New();
    vtkSmartPointer<vtkCleanPolyData> cleanFilter = vtkSmartPointer<vtkCleanPolyData>::New();

    vtkSmartPointer<vtkPolyData> CutCircleLineUp = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPolyData> lastCircle = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPolyData> lastNewCircle = vtkSmartPointer<vtkPolyData>::New();

    vtkSmartPointer<vtkPolyData> SurfaceLineUp = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPointLocator> pointLocator = vtkSmartPointer<vtkPointLocator>::New();
    pointLocator->SetDataSet(t_colon);
    pointLocator->BuildLocator();

    bool* Is_Fixed = (bool*)malloc(sizeof(bool) * t_colon->GetNumberOfPoints());
    memset(Is_Fixed, 0, sizeof(bool) * t_colon->GetNumberOfPoints());
    vtkSmartPointer<vtkDoubleArray> DeformationField = vtkSmartPointer<vtkDoubleArray>::New();
    DeformationField->SetNumberOfComponents(3);
    DeformationField->SetNumberOfTuples(t_colon->GetNumberOfPoints());

    for(vtkIdType i = 0; i < model->GetNumberOfPoints(); i++)
    {
        //vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
        double newp[3], oldp[3];
        newcenterline->GetPoint(i, newp);
        model->GetPoint(i, oldp);
        double told[3], nold[3], bold[3];
        Tangents->GetTuple(i, told);

        //Normals->GetTuple(i, nold);
        RefDirections->GetTuple(i, nold);

        vtkMath::Cross(told, nold, bold);

        plane->SetOrigin(PlaneOriginals->GetTuple(i));
        plane->SetNormal(PlaneNormals->GetTuple(i));
        cutter->SetCutFunction(plane);
        cutter->Update();

        connectivityFilter->SetClosestPoint(oldp);
        connectivityFilter->SetInputData(cutter->GetOutput());
        connectivityFilter->Update();

        vtkSmartPointer<vtkPolyData> cutCircle = vtkSmartPointer<vtkPolyData>::New();
        cutCircle = connectivityFilter->GetOutput();
        //std::cout<<i<<" "<<"cutline points:"<<cutCircle->GetNumberOfPoints()<<endl;

        if(i != 0)
        {
            for(vtkIdType j = 0; j<lastCircle->GetNumberOfPoints(); j++)
            {
                if(cutCircle->GetNumberOfPoints() >= 100)
                    break;
                connectivityFilter->SetClosestPoint(lastCircle->GetPoint(j));
                connectivityFilter->Update();
                cutCircle = connectivityFilter->GetOutput();
            }
        }

        // down sample the cutCircle
        cutCircle->DeepCopy(ReorderContour(cutCircle));
        UniformSample(30, cutCircle);

        lastCircle->DeepCopy(cutCircle);
        vtkSmartPointer<vtkPoints> newCutCircle = vtkSmartPointer<vtkPoints>::New();

        for(vtkIdType j=0; j<cutCircle->GetNumberOfPoints(); j++)
        {
            double p[3];
            cutCircle->GetPoint(j, p);
            double vector[3];
            vtkMath::Subtract(p, oldp, vector);
            double coordinate[3];
            coordinate[0] = vtkMath::Dot(vector, told);
            coordinate[1] = vtkMath::Dot(vector, nold);
            coordinate[2] = vtkMath::Dot(vector, bold);
            double pp[3];
            double vx[3], vy[3], vz[3], tmp1[3], tmp2[3];
            NewTangents->GetTuple(i, vx);
            NewNormals->GetTuple(i, vy);
            NewBinormals->GetTuple(i, vz);
            vtkMath::MultiplyScalar(vx, coordinate[0]);
            vtkMath::MultiplyScalar(vy, coordinate[1]);
            vtkMath::MultiplyScalar(vz, coordinate[2]);
            vtkMath::Add(newp, vx, tmp1);
            vtkMath::Add(tmp1, vy, tmp2);
            vtkMath::Add(tmp2, vz, pp);
            newCutCircle->InsertNextPoint(pp);

            // find the deformation of fixed points
            vtkIdType id = pointLocator->FindClosestPoint(p);
            Is_Fixed[id] = true;
            double deformationVector[3];
            vtkMath::Subtract(pp, p, deformationVector);
            DeformationField->SetTuple(id, deformationVector);
        }

        cutCircle->SetPoints(newCutCircle);
        lastNewCircle->DeepCopy(cutCircle);

        appendFilter->RemoveAllInputs();
        appendFilter->AddInputData(CutCircleLineUp);
        appendFilter->AddInputData(cutCircle);
        appendFilter->Update();
        cleanFilter->SetInputConnection(appendFilter->GetOutputPort());
        cleanFilter->Update();
        CutCircleLineUp->DeepCopy(cleanFilter->GetOutput());
    }
    // set the non-fixed points' Is_Fixed and DeformationField
    vtkSmartPointer<vtkPointLocator> centerlinePointLocator = vtkSmartPointer<vtkPointLocator>::New();
    centerlinePointLocator->SetDataSet(model);
    centerlinePointLocator->BuildLocator();

    ofstream file;
    file.open("/home/ruibinma/Desktop/deformation_fixed.txt");
    for(vtkIdType i = 0; i < t_colon->GetNumberOfPoints(); i++)
    {
        if(!Is_Fixed[i])
        {
            double p[3], cp[3], ncp[3];
            t_colon->GetPoint(i, p);
            vtkIdType centerid = centerlinePointLocator->FindClosestPoint(p);
            model->GetPoint(centerid, cp);
            newcenterline->GetPoint(centerid, ncp);

            double told[3], nold[3], bold[3];
            Tangents->GetTuple(centerid, told);
            //Normals->GetTuple(i, nold);
            RefDirections->GetTuple(centerid, nold);
            vtkMath::Cross(told, nold, bold);

            double vector[3];
            vtkMath::Subtract(p, cp, vector);
            double coordinate[3];
            coordinate[0] = vtkMath::Dot(vector, told);
            coordinate[1] = vtkMath::Dot(vector, nold);
            coordinate[2] = vtkMath::Dot(vector, bold);
            double pp[3];
            double vx[3], vy[3], vz[3], tmp1[3], tmp2[3];
            NewTangents->GetTuple(centerid, vx);
            NewNormals->GetTuple(centerid, vy);
            NewBinormals->GetTuple(centerid, vz);
            vtkMath::MultiplyScalar(vx, coordinate[0]);
            vtkMath::MultiplyScalar(vy, coordinate[1]);
            vtkMath::MultiplyScalar(vz, coordinate[2]);
            vtkMath::Add(ncp, vx, tmp1);
            vtkMath::Add(tmp1, vy, tmp2);
            vtkMath::Add(tmp2, vz, pp);

            double deformationVector[3];
            vtkMath::Subtract(pp, p, deformationVector);
            DeformationField->SetTuple(i, deformationVector);
        }

        double v[3];
        DeformationField->GetTuple(i, v);
        double p[3];
        t_colon->GetPoint(i, p);
        std::cout<<i<<" -> "<<p[0]<<" "<<p[1]<<" "<<p[2]<<" + "<<v[0]<<" "<<v[1]<<" "<<v[2]<<" "<<Is_Fixed[i]<<std::endl;
        file<<v[0]<<" "<<v[1]<<" "<<v[2]<<" "<<Is_Fixed[i]<<std::endl;
    }
    file.close();

    vtkSmartPointer<vtkPolyDataMapper> CutCircleLineUpMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    CutCircleLineUpMapper->SetInputData(CutCircleLineUp);
    CutCircleLineUpMapper->Update();
    vtkSmartPointer<vtkActor> CutCircleLineUpActor = vtkSmartPointer<vtkActor>::New();
    CutCircleLineUpActor->SetMapper(CutCircleLineUpMapper);
    CutCircleLineUpActor->GetProperty()->SetColor(1, 0, 0);
    t_rendermanager->renderModel(CutCircleLineUpActor);
    std::cout<<"cutcircles have points"<<CutCircleLineUp->GetNumberOfPoints()<<endl;


    vtkSmartPointer<vtkPoints> lineuppoints = vtkSmartPointer<vtkPoints>::New();
    for(vtkIdType i = 0; i < t_colon->GetNumberOfPoints(); i++)
    {
        double v[3];
        DeformationField->GetTuple(i, v);
        double p[3], pp[3];
        t_colon->GetPoint(i, p);
        vtkMath::Add(p, v, pp);
        lineuppoints->InsertNextPoint(pp);
    }
    SurfaceLineUp->DeepCopy(t_colon);
    SurfaceLineUp->SetPoints(lineuppoints);
    vtkSmartPointer<vtkPolyDataMapper> SurfaceLineUpMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    SurfaceLineUpMapper->SetInputData(SurfaceLineUp);
    SurfaceLineUpMapper->Update();
    vtkSmartPointer<vtkActor> SurfaceLineUpActor = vtkSmartPointer<vtkActor>::New();
    SurfaceLineUpActor->SetMapper(SurfaceLineUpMapper);
    t_rendermanager->renderModel(SurfaceLineUpActor);


    free(Is_Fixed);
    std::cout<<"Deformation End"<<endl;
    return SurfaceLineUp;
}
vtkSmartPointer<vtkPolyData> Centerline::Deformation_v3(vtkSmartPointer<vtkDoubleArray> S, vtkSmartPointer<vtkDoubleArray> Curvatures, vtkSmartPointer<vtkIdList> CurvaturePointIds,
                                                     vtkSmartPointer<vtkDoubleArray> Tangents, vtkSmartPointer<vtkDoubleArray> Normals,
                                                     vtkSmartPointer<vtkPolyData> t_colon, RenderManager *t_rendermanager,
                                                     vtkSmartPointer<vtkDoubleArray> PlaneOriginals, vtkSmartPointer<vtkDoubleArray> PlaneNormals,
                                                     vtkSmartPointer<vtkDoubleArray> RefDirections, FileManager *t_filemanager)
{
    //PutNormalsOnSameSide(Normals, Curvatures);
    std::cout<<"Deformation"<<endl;
    bool straight = true;
    // Eliminate the torsion by growing the curve on a plane, according to: -dNnew/dSnew = -k*Tnew
    double point[3], nextpoint[3];
    double tangent[3], nexttangent[3];
    double normal[3], nextnormal[3];
    double binormal[3];
    double ds, curvature;
    vtkSmartPointer<vtkDoubleArray> NewTangents = vtkSmartPointer<vtkDoubleArray>::New();
    NewTangents->SetNumberOfComponents(3);
    vtkSmartPointer<vtkDoubleArray> NewNormals = vtkSmartPointer<vtkDoubleArray>::New();
    NewNormals->SetNumberOfComponents(3);
    vtkSmartPointer<vtkDoubleArray> NewBinormals = vtkSmartPointer<vtkDoubleArray>::New();
    NewBinormals->SetNumberOfComponents(3);
    vtkSmartPointer<vtkPoints> newpoints = vtkSmartPointer<vtkPoints>::New();
    for(vtkIdType i=0; i<model->GetNumberOfPoints(); i++)
    {
        if(i == 0)
        {
            model->GetPoint(0, point);
            point[0] = point[0] + 60;
            binormal[0] = 0;
            binormal[1] = 0;
            binormal[2] = 1;
            tangent[0] = 1;
            tangent[1] = 0;
            tangent[2] = 0;
            normal[0] = 0;
            normal[1] = 1;
            normal[2] = 0;
        }
        ds = (i != model->GetNumberOfPoints()-1)?(S->GetValue(i + 1) - S->GetValue(i)):0;
        curvature = Curvatures->GetValue(i);

        // Record the New Axis System
        NewTangents->InsertNextTuple(tangent);
        NewNormals->InsertNextTuple(normal);
        NewBinormals->InsertNextTuple(binormal);

        double dnormal[3];
        if(straight)
        {
            dnormal[0] = 0; dnormal[1] = 0; dnormal[2] = 0;
        }
        else
        {
            dnormal[0] = tangent[0]; dnormal[1] = tangent[1]; dnormal[2] = tangent[2];
            double RelaxationFactor = 2;
            vtkMath::MultiplyScalar(dnormal, -curvature / RelaxationFactor);
        }
        vtkMath::Add(normal, dnormal, nextnormal);
        vtkMath::Normalize(nextnormal);
        //std::cout<<i<<"\ts="<<S->GetValue(i)<<"("<<ds<<")"<<"\tk="<<curvature<<endl;
        vtkMath::Cross(nextnormal, binormal, nexttangent);
        vtkMath::Normalize(nexttangent);
        vtkMath::MultiplyScalar(tangent, ds);
        vtkMath::Add(point, tangent, nextpoint);
        newpoints->InsertNextPoint(point);

        point[0] = nextpoint[0]; point[1] = nextpoint[1]; point[2] = nextpoint[2];
        normal[0] = nextnormal[0]; normal[1] = nextnormal[1]; normal[2] = nextnormal[2];
        tangent[0] = nexttangent[0]; tangent[1] = nexttangent[1]; tangent[2] = nexttangent[2];
    }

    vtkSmartPointer<vtkPolyData> newcenterline = vtkSmartPointer<vtkPolyData>::New();
    newcenterline->SetPoints(newpoints);
    vtkSmartPointer<vtkVertexGlyphFilter> newVertexFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
    newVertexFilter->SetInputData(newcenterline);
    newVertexFilter->Update();
    vtkSmartPointer<vtkPolyDataMapper> newMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    newMapper->SetInputConnection(newVertexFilter->GetOutputPort());
    vtkSmartPointer<vtkActor> newActor = vtkSmartPointer<vtkActor>::New();
    newActor->SetMapper(newMapper);
    t_rendermanager->renderModel(newActor);

    // Line Up the Cross Sections
    vtkSmartPointer<vtkCutter> cutter = vtkSmartPointer<vtkCutter>::New();
    cutter->SetInputData(t_colon);

    vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();

    vtkSmartPointer<vtkPolyDataConnectivityFilter> connectivityFilter = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
    connectivityFilter->SetInputConnection(cutter->GetOutputPort());
    connectivityFilter->SetExtractionModeToClosestPointRegion();

    vtkSmartPointer<vtkAppendPolyData> appendFilter = vtkSmartPointer<vtkAppendPolyData>::New();

    vtkSmartPointer<vtkPolyData> OriginCutCircle = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPolyData> CutCircleLineUp = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPolyData> lastCircle = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPolyData> reflastCircle = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPolyData> lastNewCircle = vtkSmartPointer<vtkPolyData>::New();

    vtkSmartPointer<vtkPolyData> SurfaceLineUp = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPointLocator> pointLocator = vtkSmartPointer<vtkPointLocator>::New();
    pointLocator->SetDataSet(t_colon);
    pointLocator->BuildLocator();

    vtkSmartPointer<vtkDoubleArray> DeformationField = vtkSmartPointer<vtkDoubleArray>::New();
    DeformationField->SetNumberOfComponents(3);
    DeformationField->SetNumberOfTuples(t_colon->GetNumberOfPoints());
    vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    vtkSmartPointer<vtkPolyData> source = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPolyData> target = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(t_colon->GetNumberOfPoints());
    int count = 0;

    /*
    vtkSmartPointer<vtkIdList> Sections = vtkSmartPointer<vtkIdList>::New();
    for(vtkIdType i = 0; i < t_colon->GetNumberOfPoints(); i++)
    {
        vtkIdType section = GetTheSectionIdOfAPoint_v2(t_colon, i, PlaneOriginals, PlaneNormals, CurvaturePointIds);
        Sections->InsertNextId(section);
    }
    */

    double p[3], t[3], pp[3], pseed[3];
    model->GetPoint(0, p);
    Tangents->GetTuple(0, t);
    vtkMath::MultiplyScalar(t, -50);
    vtkMath::Add(p, t, pp);
    vtkSmartPointer<vtkPointLocator> seedLocator = vtkSmartPointer<vtkPointLocator>::New();
    seedLocator->SetDataSet(t_colon);
    seedLocator->BuildLocator();
    vtkIdType seed = seedLocator->FindClosestPoint(pp);
    t_colon->GetPoint(seed, pseed);
    vtkSmartPointer<vtkPoints> seedpoints = vtkSmartPointer<vtkPoints>::New();
    seedpoints->InsertNextPoint(pseed);
    vtkSmartPointer<vtkPolyData> seedpoly = vtkSmartPointer<vtkPolyData>::New();
    seedpoly->SetPoints(seedpoints);
    vtkSmartPointer<vtkVertexGlyphFilter> seedVertex = vtkSmartPointer<vtkVertexGlyphFilter>::New();
    seedVertex->SetInputData(seedpoly);
    seedVertex->Update();
    vtkSmartPointer<vtkPolyDataMapper> seedMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    seedMapper->SetInputConnection(seedVertex->GetOutputPort());
    seedMapper->Update();
    vtkSmartPointer<vtkActor> seedActor = vtkSmartPointer<vtkActor>::New();
    seedActor->SetMapper(seedMapper);
    seedActor->GetProperty()->SetPointSize(10);
    seedActor->GetProperty()->SetColor(0,1,0);
    t_rendermanager->renderModel(seedActor);

    vtkSmartPointer<vtkIdList> Sections = vtkSmartPointer<vtkIdList>::New();
    for(vtkIdType i = 0; i < t_colon->GetNumberOfPoints(); i++)
    {
        Sections->InsertNextId(-1);
    }
    Sections->SetId(seed, 0);
    GetSectionIds_loop(t_colon, seed, Sections, PlaneOriginals, PlaneNormals);

    for(vtkIdType i = 0; i<model->GetNumberOfPoints(); i++)
    {
        //vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
        double newp[3], oldp[3];
        newcenterline->GetPoint(i, newp);
        model->GetPoint(i, oldp);
        double told[3], nold[3], bold[3];
        Tangents->GetTuple(i, told);

        //Normals->GetTuple(i, nold);
        RefDirections->GetTuple(i, nold);

        vtkMath::Cross(told, nold, bold);

        plane->SetOrigin(PlaneOriginals->GetTuple(i));
        plane->SetNormal(PlaneNormals->GetTuple(i));
        cutter->SetCutFunction(plane);
        cutter->Update();

        connectivityFilter->SetClosestPoint(oldp);
        connectivityFilter->SetInputData(cutter->GetOutput());
        connectivityFilter->Update();

        vtkSmartPointer<vtkPolyData> cutCircle = vtkSmartPointer<vtkPolyData>::New();
        cutCircle = connectivityFilter->GetOutput();
        //std::cout<<i<<" "<<"cutline points:"<<cutCircle->GetNumberOfPoints()<<endl;
        if(i != 0)
        {
            for(vtkIdType j = 0; j<reflastCircle->GetNumberOfPoints(); j++)
            {

                if(cutCircle->GetNumberOfPoints() >= 100)
                    break;
                connectivityFilter->SetClosestPoint(lastCircle->GetPoint(j));
                connectivityFilter->Update();
                cutCircle = connectivityFilter->GetOutput();
            }
        }
        reflastCircle->DeepCopy(cutCircle);
        appendFilter->RemoveAllInputs();
        appendFilter->AddInputData(cutCircle);
        appendFilter->AddInputData(OriginCutCircle);
        appendFilter->Update();
        OriginCutCircle->DeepCopy(appendFilter->GetOutput());

        // build source landmarks
        cutCircle->DeepCopy(ReorderContour(cutCircle));
        UniformSample(40, cutCircle);
        appendFilter->RemoveAllInputs();
        appendFilter->AddInputData(cutCircle);
        appendFilter->AddInputData(lastCircle);
        appendFilter->Update();
        source->ShallowCopy(appendFilter->GetOutput());
        lastCircle->DeepCopy(cutCircle);

        vtkSmartPointer<vtkPoints> newCutCircle = vtkSmartPointer<vtkPoints>::New();

        for(vtkIdType j=0; j<cutCircle->GetNumberOfPoints(); j++)
        {
            double p[3];
            cutCircle->GetPoint(j, p);
            double vector[3];
            vtkMath::Subtract(p, oldp, vector);
            double coordinate[3];
            coordinate[0] = vtkMath::Dot(vector, told);
            coordinate[1] = vtkMath::Dot(vector, nold);
            coordinate[2] = vtkMath::Dot(vector, bold);
            double pp[3];
            double vx[3], vy[3], vz[3], tmp1[3], tmp2[3];
            NewTangents->GetTuple(i, vx);
            NewNormals->GetTuple(i, vy);
            NewBinormals->GetTuple(i, vz);
            vtkMath::MultiplyScalar(vx, coordinate[0]);
            vtkMath::MultiplyScalar(vy, coordinate[1]);
            vtkMath::MultiplyScalar(vz, coordinate[2]);
            vtkMath::Add(newp, vx, tmp1);
            vtkMath::Add(tmp1, vy, tmp2);
            vtkMath::Add(tmp2, vz, pp);
            newCutCircle->InsertNextPoint(pp);
        }
        cutCircle->SetPoints(newCutCircle);
        // build target landmarks
        appendFilter->RemoveAllInputs();
        appendFilter->AddInputData(cutCircle);
        appendFilter->AddInputData(lastNewCircle);
        appendFilter->Update();
        target->ShallowCopy(appendFilter->GetOutput());

        lastNewCircle->DeepCopy(cutCircle);

        appendFilter->RemoveAllInputs();
        appendFilter->AddInputData(CutCircleLineUp);
        appendFilter->AddInputData(cutCircle);
        appendFilter->Update();
        //cleanFilter->SetInputConnection(appendFilter->GetOutputPort());
        //cleanFilter->Update();
        CutCircleLineUp->DeepCopy(appendFilter->GetOutput());

        // process the points in this section
        vtkSmartPointer<vtkPoints> sectionpoints = vtkSmartPointer<vtkPoints>::New();
        vtkSmartPointer<vtkPolyData> sectionpoly = vtkSmartPointer<vtkPolyData>::New();
        vtkSmartPointer<vtkIdList> sectionids = vtkSmartPointer<vtkIdList>::New();
        for(vtkIdType j =0; j < t_colon->GetNumberOfPoints(); j++)
        {
            vtkIdType section = Sections->GetId(j);
            if(section == i)
            {
                double p[3];
                t_colon->GetPoint(j, p);
                sectionids->InsertNextId(j);
                sectionpoints->InsertNextPoint(p);
            }
        }
        sectionpoly->SetPoints(sectionpoints);
        count += sectionpoly->GetNumberOfPoints();
        std::cout<<count<<endl;
        vtkSmartPointer<vtkThinPlateSplineTransform> transform = vtkSmartPointer<vtkThinPlateSplineTransform>::New();
        transform->SetSourceLandmarks(source->GetPoints());
        transform->SetTargetLandmarks(target->GetPoints());
        transform->SetBasisToR();
        transformFilter->SetTransform(transform);
        transformFilter->SetInputData(sectionpoly);
        transformFilter->Update();
        sectionpoly->DeepCopy(transformFilter->GetOutput());
        assert(sectionpoly->GetNumberOfPoints() == sectionids->GetNumberOfIds());
        for(vtkIdType j = 0; j < sectionpoly->GetNumberOfPoints(); j++)
        {
            double p[3];
            sectionpoly->GetPoint(j, p);
            points->SetPoint(sectionids->GetId(j), p);
        }
    }
    // very last piece
    source->ShallowCopy(lastCircle);
    target->ShallowCopy(lastNewCircle);
    vtkSmartPointer<vtkPoints> sectionpoints = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkPolyData> sectionpoly = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkIdList> sectionids = vtkSmartPointer<vtkIdList>::New();
    for(vtkIdType j =0; j < t_colon->GetNumberOfPoints(); j++)
    {
        vtkIdType section = Sections->GetId(j);
        if(section == model->GetNumberOfPoints())
        {
            double p[3];
            t_colon->GetPoint(j, p);
            sectionids->InsertNextId(j);
            sectionpoints->InsertNextPoint(p);
        }
    }
    sectionpoly->SetPoints(sectionpoints);
    count += sectionpoly->GetNumberOfPoints();
    std::cout<<count<<endl;
    vtkSmartPointer<vtkThinPlateSplineTransform> transform = vtkSmartPointer<vtkThinPlateSplineTransform>::New();
    transform->SetSourceLandmarks(source->GetPoints());
    transform->SetTargetLandmarks(target->GetPoints());
    transform->SetBasisToR();
    transformFilter->RemoveAllInputs();
    transformFilter->SetInputData(sectionpoly);
    transformFilter->SetTransform(transform);
    transformFilter->Update();
    sectionpoly->DeepCopy(transformFilter->GetOutput());
    assert(sectionpoly->GetNumberOfPoints() == sectionids->GetNumberOfIds());
    for(vtkIdType j = 0; j < sectionpoly->GetNumberOfPoints(); j++)
    {
        double p[3];
        sectionpoly->GetPoint(j, p);
        points->SetPoint(sectionids->GetId(j), p);
    }

    SurfaceLineUp->DeepCopy(t_colon);
    SurfaceLineUp->SetPoints(points);

    vtkSmartPointer<vtkPolyDataMapper> OriginCutCircleMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    OriginCutCircleMapper->SetInputData(OriginCutCircle);
    OriginCutCircleMapper->Update();
    vtkSmartPointer<vtkActor> OriginCutCircleActor = vtkSmartPointer<vtkActor>::New();
    OriginCutCircleActor->SetMapper(OriginCutCircleMapper);
    OriginCutCircleActor->GetProperty()->SetColor(0, 0, 1);
    t_rendermanager->renderModel(OriginCutCircleActor);

    vtkSmartPointer<vtkPolyDataMapper> CutCircleLineUpMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    CutCircleLineUpMapper->SetInputData(CutCircleLineUp);
    CutCircleLineUpMapper->Update();
    vtkSmartPointer<vtkActor> CutCircleLineUpActor = vtkSmartPointer<vtkActor>::New();
    CutCircleLineUpActor->SetMapper(CutCircleLineUpMapper);
    CutCircleLineUpActor->GetProperty()->SetColor(1, 0, 0);
    t_rendermanager->renderModel(CutCircleLineUpActor);

    vtkSmartPointer<vtkPolyDataMapper> SurfaceLineUpMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    SurfaceLineUpMapper->SetInputData(SurfaceLineUp);
    SurfaceLineUpMapper->Update();
    vtkSmartPointer<vtkActor> SurfaceLineUpActor = vtkSmartPointer<vtkActor>::New();
    SurfaceLineUpActor->SetMapper(SurfaceLineUpMapper);
    t_rendermanager->renderModel(SurfaceLineUpActor);
    t_filemanager->SaveFile(SurfaceLineUp, "SurfaceLineUp_v3.stl");

    std::cout<<"Deformation End"<<endl;
    return SurfaceLineUp;
}
vtkSmartPointer<vtkPolyData> Centerline::Deformation_v3_1(vtkSmartPointer<vtkDoubleArray> S, vtkSmartPointer<vtkDoubleArray> Curvatures, vtkSmartPointer<vtkIdList> CurvaturePointIds,
                                                     vtkSmartPointer<vtkDoubleArray> Tangents, vtkSmartPointer<vtkDoubleArray> Normals,
                                                     vtkSmartPointer<vtkPolyData> t_colon, RenderManager *t_rendermanager, RenderManager *t_rendermanager_right,
                                                     vtkSmartPointer<vtkDoubleArray> PlaneOriginals, vtkSmartPointer<vtkDoubleArray> PlaneNormals,
                                                     vtkSmartPointer<vtkDoubleArray> RefDirections, FileManager *t_filemanager)
{

    //PutNormalsOnSameSide(Normals, Curvatures);
    std::cout<<"Deformation"<<endl;
    int choice = 2; // 0-straight(stretch or press); 1-sin; 2-circle; 3-helix; 4-twist; 5-L-shape; 6-waterpipe
    double translate = 0;
    // Eliminate the torsion by growing the curve on a plane, according to: -dNnew/dSnew = -k*Tnew
    double point[3], nextpoint[3];
    double tangent[3], nexttangent[3];
    double normal[3], nextnormal[3];
    double binormal[3];
    double ds, curvature;
    vtkSmartPointer<vtkDoubleArray> NewTangents = vtkSmartPointer<vtkDoubleArray>::New();
    NewTangents->SetNumberOfComponents(3);
    vtkSmartPointer<vtkDoubleArray> NewNormals = vtkSmartPointer<vtkDoubleArray>::New();
    NewNormals->SetNumberOfComponents(3);
    vtkSmartPointer<vtkDoubleArray> NewBinormals = vtkSmartPointer<vtkDoubleArray>::New();
    NewBinormals->SetNumberOfComponents(3);
    vtkSmartPointer<vtkPoints> newpoints = vtkSmartPointer<vtkPoints>::New();
    for(int i=0; i<model->GetNumberOfPoints(); i++)
    {
        double stretch = 1;
        if(choice == 0)
        {
            if(i == 0)
            {
                model->GetPoint(0, point);
                point[0] = point[0] + translate;
                binormal[0] = 0;
                binormal[1] = 0;
                binormal[2] = 1;
                normal[0] = 0;
                normal[1] = 1;
                normal[2] = 0;
            }
            ds = (i != model->GetNumberOfPoints()-1)?(S->GetValue(i + 1) - S->GetValue(i)):0;
            ds = ds*stretch;

            tangent[0] = 1;
            tangent[1] = 0;
            tangent[2] = 0;

            NewTangents->InsertNextTuple(tangent);
            NewNormals->InsertNextTuple(normal);
            NewBinormals->InsertNextTuple(binormal);

            vtkMath::MultiplyScalar(tangent, ds);
            vtkMath::Add(point, tangent, nextpoint);
            newpoints->InsertNextPoint(point);
            point[0] = nextpoint[0]; point[1] = nextpoint[1]; point[2] = nextpoint[2];
        }
        else if(choice == 1)
        {

            if(i == 0)
            {
                model->GetPoint(0, point);
                point[0] = point[0] + translate;
                binormal[0] = 0;
                binormal[1] = 0;
                binormal[2] = 1;
            }
            ds = (i != model->GetNumberOfPoints()-1)?(S->GetValue(i + 1) - S->GetValue(i)):0;

            double p = S->GetValue(i) / S->GetValue(model->GetNumberOfPoints()-1);
            double x = cos(4*3.1415926*p)*sqrt(2)/2;
            double y = sqrt(1-x*x);
            normal[0] = x; normal[1] = y; normal[2] = 0;
            vtkMath::Cross(normal, binormal, tangent);

            NewTangents->InsertNextTuple(tangent);
            NewNormals->InsertNextTuple(normal);
            NewBinormals->InsertNextTuple(binormal);

            vtkMath::MultiplyScalar(tangent, ds);
            newpoints->InsertNextPoint(point);
            vtkMath::Add(point, tangent, nextpoint);
            point[0] = nextpoint[0]; point[1] = nextpoint[1]; point[2] = nextpoint[2];
        }
        else if(choice == 2)
        {

            if(i == 0)
            {
                model->GetPoint(0, point);
                point[0] = point[0] + translate;
                binormal[0] = 0;
                binormal[1] = 0;
                binormal[2] = 1;
            }
            ds = (i != model->GetNumberOfPoints()-1)?(S->GetValue(i + 1) - S->GetValue(i)):0;

            double p = S->GetValue(i) / S->GetValue(model->GetNumberOfPoints()-1) / 2; // divide by 2 to be half a circle
            double x = sin(2*3.1415926*p);
            double y = cos(2*3.1415926*p);
            tangent[0] = x; tangent[1] = y; tangent[2] = 0;
            vtkMath::Cross(binormal, tangent, normal);

            NewTangents->InsertNextTuple(tangent);
            NewNormals->InsertNextTuple(normal);
            NewBinormals->InsertNextTuple(binormal);

            vtkMath::MultiplyScalar(tangent, ds);
            newpoints->InsertNextPoint(point);
            vtkMath::Add(point, tangent, nextpoint);
            point[0] = nextpoint[0]; point[1] = nextpoint[1]; point[2] = nextpoint[2];
        }
        else if(choice == 3)
        {
            double r = 30, c = 10, k=2;
            if(i == 0)
            {
                model->GetPoint(0, point);
                point[0] = point[0] + translate;
            }
            ds = (i != model->GetNumberOfPoints()-1)?(S->GetValue(i + 1) - S->GetValue(i)):0;

            double p = S->GetValue(i) / S->GetValue(model->GetNumberOfPoints()-1);
            tangent[0] = -r*sin(2*k*3.1415926*p); tangent[1] = r*cos(2*k*3.1415926*p); tangent[2] = c;
            vtkMath::Normalize(tangent);
            normal[0] = -cos(2*k*3.1415926*p); normal[1] = -sin(2*k*3.1415926*p); normal[2] = 0;
            vtkMath::Cross(tangent, normal, binormal);

            NewTangents->InsertNextTuple(tangent);
            NewNormals->InsertNextTuple(normal);
            NewBinormals->InsertNextTuple(binormal);

            vtkMath::MultiplyScalar(tangent, ds);
            newpoints->InsertNextPoint(point);
            vtkMath::Add(point, tangent, nextpoint);
            point[0] = nextpoint[0]; point[1] = nextpoint[1]; point[2] = nextpoint[2];
        }
        if(choice == 4)
        {
            double k=8;
            double stretch = 1;
            if(i == 0)
            {
                model->GetPoint(0, point);
                point[0] = point[0] + translate;
            }
            ds = (i != model->GetNumberOfPoints()-1)?(S->GetValue(i + 1) - S->GetValue(i)):0;
            ds = ds*stretch;

            tangent[0] = 1;
            tangent[1] = 0;
            tangent[2] = 0;
            double p = S->GetValue(i) / S->GetValue(model->GetNumberOfPoints()-1);
            normal[0] = 0; normal[1] = cos(2*k*3.1415926*p); normal[2] = sin(2*k*3.1415926*p);
            vtkMath::Cross(tangent, normal, binormal);

            NewTangents->InsertNextTuple(tangent);
            NewNormals->InsertNextTuple(normal);
            NewBinormals->InsertNextTuple(binormal);

            vtkMath::MultiplyScalar(tangent, ds);
            vtkMath::Add(point, tangent, nextpoint);
            newpoints->InsertNextPoint(point);
            point[0] = nextpoint[0]; point[1] = nextpoint[1]; point[2] = nextpoint[2];
        }
        else if(choice == 5)
        {
            if(i == 0)
            {
                model->GetPoint(0, point);
                point[0] = point[0] + translate;
                binormal[0] = 0;
                binormal[1] = 0;
                binormal[2] = 1;
            }
            ds = (i != model->GetNumberOfPoints()-1)?(S->GetValue(i + 1) - S->GetValue(i)):0;

            if(i < 450)
            {
                tangent[0] = 1;
                tangent[1] = 0;
                tangent[2] = 0;
            }
            else if(i >= 550)
            {
                tangent[0] = 0;
                tangent[1] = 1;
                tangent[2] = 0;
            }
            else
            {
                double p = (S->GetValue(i) - S->GetValue(449))/(S->GetValue(550) - S->GetValue(449));
                tangent[0] = cos(3.1415926/2*p);
                tangent[1] = sin(3.1415926/2*p);
                tangent[2] = 0;
            }
            vtkMath::Cross(binormal, tangent, normal);


            NewTangents->InsertNextTuple(tangent);
            NewNormals->InsertNextTuple(normal);
            NewBinormals->InsertNextTuple(binormal);

            vtkMath::MultiplyScalar(tangent, ds);
            vtkMath::Add(point, tangent, nextpoint);
            newpoints->InsertNextPoint(point);
            point[0] = nextpoint[0]; point[1] = nextpoint[1]; point[2] = nextpoint[2];
        }
        else if(choice == 6)
        {
            if(i == 0)
            {
                model->GetPoint(0, point);
                point[0] = point[0] + translate;
                binormal[0] = 0;
                binormal[1] = 0;
                binormal[2] = 1;
            }
            ds = (i != model->GetNumberOfPoints()-1)?(S->GetValue(i + 1) - S->GetValue(i)):0;

            if(i < 500)
            {
                if(i < 293)
                {
                    tangent[0] = 1;
                    tangent[1] = 0;
                    tangent[2] = 0;
                }
                else if(i >= 293 && i < 373)
                {
                    double p = (S->GetValue(i) - S->GetValue(292))/(S->GetValue(373) - S->GetValue(292));
                    tangent[0] = cos(3.1415926/2*p);
                    tangent[1] = sin(3.1415926/2*p);
                    tangent[2] = 0;
                }
                else
                {
                    tangent[0] = 0;
                    tangent[1] = 1;
                    tangent[0] = 0;
                }
                vtkMath::Cross(binormal, tangent, normal);
            }
            else
            {
                normal[0] = -1; normal[1] = 0; normal[2] = 0;
                if(i < 626)
                {
                    tangent[0] = 0;
                    tangent[1] = 1;
                    tangent[2] = 0;
                }
                else if(i >= 706)
                {
                    tangent[0] = 0;
                    tangent[1] = 0;
                    tangent[2] = 1;
                }
                else
                {
                    double p = (S->GetValue(i) - S->GetValue(625))/(S->GetValue(706) - S->GetValue(625));
                    tangent[0] = 0;
                    tangent[1] = cos(3.1415926/2*p);
                    tangent[2] = sin(3.1415926/2*p);
                }
                vtkMath::Cross(tangent, normal, binormal);
            }


            NewTangents->InsertNextTuple(tangent);
            NewNormals->InsertNextTuple(normal);
            NewBinormals->InsertNextTuple(binormal);

            vtkMath::MultiplyScalar(tangent, ds);
            vtkMath::Add(point, tangent, nextpoint);
            newpoints->InsertNextPoint(point);
            point[0] = nextpoint[0]; point[1] = nextpoint[1]; point[2] = nextpoint[2];
        }
    }
    /*
    for(vtkIdType i=0; i<model->GetNumberOfPoints(); i++)
    {
        if(i == 0)
        {
            model->GetPoint(0, point);
            point[0] = point[0] + 60;
            binormal[0] = 0;
            binormal[1] = 0;
            binormal[2] = 1;
            tangent[0] = 1;
            tangent[1] = 0;
            tangent[2] = 0;
            normal[0] = 0;
            normal[1] = 1;
            normal[2] = 0;
        }
        ds = (i != model->GetNumberOfPoints()-1)?(S->GetValue(i + 1) - S->GetValue(i)):0;
        curvature = Curvatures->GetValue(i);

        // Record the New Axis System
        NewTangents->InsertNextTuple(tangent);
        NewNormals->InsertNextTuple(normal);
        NewBinormals->InsertNextTuple(binormal);


        double dnormal[3];
        if(straight)
        {
            dnormal[0] = 0; dnormal[1] = 0; dnormal[2] = 0;
        }
        else
        {
            dnormal[0] = tangent[0]; dnormal[1] = tangent[1]; dnormal[2] = tangent[2];
            double RelaxationFactor = 2;
            vtkMath::MultiplyScalar(dnormal, -curvature / RelaxationFactor);
        }
        vtkMath::Add(normal, dnormal, nextnormal);
        vtkMath::Normalize(nextnormal);

        //std::cout<<i<<"\ts="<<S->GetValue(i)<<"("<<ds<<")"<<"\tk="<<curvature<<endl;
        vtkMath::Cross(nextnormal, binormal, nexttangent);
        vtkMath::Normalize(nexttangent);
        vtkMath::MultiplyScalar(tangent, ds);
        vtkMath::Add(point, tangent, nextpoint);
        newpoints->InsertNextPoint(point);

        point[0] = nextpoint[0]; point[1] = nextpoint[1]; point[2] = nextpoint[2];
        normal[0] = nextnormal[0]; normal[1] = nextnormal[1]; normal[2] = nextnormal[2];
        tangent[0] = nexttangent[0]; tangent[1] = nexttangent[1]; tangent[2] = nexttangent[2];
    }
    */

    vtkSmartPointer<vtkPolyData> newcenterline = vtkSmartPointer<vtkPolyData>::New();
    newcenterline->DeepCopy(model);
    newcenterline->SetPoints(newpoints);
    vtkSmartPointer<vtkPolyDataMapper> newMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    newMapper->SetInputData(newcenterline);
    vtkSmartPointer<vtkActor> newActor = vtkSmartPointer<vtkActor>::New();
    newActor->SetMapper(newMapper);
    newActor->GetProperty()->SetLineWidth(5);
    newActor->GetProperty()->SetColor(1, 1, 0);
    t_rendermanager_right->renderModel(newActor);

    // Line Up the Cross Sections
    vtkSmartPointer<vtkCutter> cutter = vtkSmartPointer<vtkCutter>::New();
    cutter->SetInputData(t_colon);
    vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();
    vtkSmartPointer<vtkPolyDataConnectivityFilter> connectivityFilter = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
    connectivityFilter->SetInputConnection(cutter->GetOutputPort());
    connectivityFilter->SetExtractionModeToClosestPointRegion();
    vtkSmartPointer<vtkAppendPolyData> appendFilter = vtkSmartPointer<vtkAppendPolyData>::New();

    vtkSmartPointer<vtkPolyData> OriginCutCircle = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPolyData> CutCircleLineUp = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPolyData> lastCircle = vtkSmartPointer<vtkPolyData>::New();

    vtkSmartPointer<vtkPolyData> SurfaceLineUp = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPointLocator> pointLocator = vtkSmartPointer<vtkPointLocator>::New();
    pointLocator->SetDataSet(t_colon);
    pointLocator->BuildLocator();

    vtkPolyDataGroup* CircleGroup = new vtkPolyDataGroup;
    vtkPolyDataGroup* ResampledCircleGroup = new vtkPolyDataGroup;
    vtkPolyDataGroup* ResampledLineUpGroup = new vtkPolyDataGroup;

    for(vtkIdType i = 0; i<model->GetNumberOfPoints(); i++)
    {
        double newp[3], oldp[3];
        newcenterline->GetPoint(i, newp);
        model->GetPoint(i, oldp);
        double told[3], nold[3], bold[3];
        Tangents->GetTuple(i, told);

        //Normals->GetTuple(i, nold);
        RefDirections->GetTuple(i, nold);

        vtkMath::Cross(told, nold, bold);

        plane->SetOrigin(PlaneOriginals->GetTuple(i));
        plane->SetNormal(PlaneNormals->GetTuple(i));
        cutter->SetCutFunction(plane);
        cutter->Update();

        connectivityFilter->SetClosestPoint(oldp);
        connectivityFilter->SetInputData(cutter->GetOutput());
        connectivityFilter->Update();

        vtkSmartPointer<vtkPolyData> cutCircle = vtkSmartPointer<vtkPolyData>::New();
        cutCircle = connectivityFilter->GetOutput();
        if(i != 0)
        {
            for(vtkIdType j = 0; j<lastCircle->GetNumberOfPoints(); j++)
            {

                if(cutCircle->GetNumberOfPoints() >= 70)
                    break;
                connectivityFilter->SetClosestPoint(lastCircle->GetPoint(j));
                connectivityFilter->Update();
                cutCircle = connectivityFilter->GetOutput();
            }
        }
        lastCircle->DeepCopy(cutCircle);

        CircleGroup->AddMember(cutCircle);// record the cutCircle into group

        appendFilter->RemoveAllInputs();
        appendFilter->AddInputData(cutCircle);
        appendFilter->AddInputData(OriginCutCircle);
        appendFilter->Update();
        OriginCutCircle->DeepCopy(appendFilter->GetOutput());

        vtkSmartPointer<vtkPoints> newCutCircle = vtkSmartPointer<vtkPoints>::New();
        vtkSmartPointer<vtkPolyData> newCutCirclePoly = vtkSmartPointer<vtkPolyData>::New();
        newCutCirclePoly->DeepCopy(cutCircle);
        for(vtkIdType j=0; j<cutCircle->GetNumberOfPoints(); j++)
        {
            double p[3];
            cutCircle->GetPoint(j, p);
            double vector[3];
            vtkMath::Subtract(p, oldp, vector);
            double coordinate[3];
            coordinate[0] = vtkMath::Dot(vector, told);
            coordinate[1] = vtkMath::Dot(vector, nold);
            coordinate[2] = vtkMath::Dot(vector, bold);
            double pp[3];
            double vx[3], vy[3], vz[3], tmp1[3], tmp2[3];
            NewTangents->GetTuple(i, vx);
            NewNormals->GetTuple(i, vy);
            NewBinormals->GetTuple(i, vz);
            vtkMath::MultiplyScalar(vx, coordinate[0]);
            vtkMath::MultiplyScalar(vy, coordinate[1]);
            vtkMath::MultiplyScalar(vz, coordinate[2]);
            vtkMath::Add(newp, vx, tmp1);
            vtkMath::Add(tmp1, vy, tmp2);
            vtkMath::Add(tmp2, vz, pp);
            newCutCircle->InsertNextPoint(pp);
        }
        newCutCirclePoly->SetPoints(newCutCircle);
        appendFilter->RemoveAllInputs();
        appendFilter->AddInputData(CutCircleLineUp);
        appendFilter->AddInputData(newCutCirclePoly);
        appendFilter->Update();
        CutCircleLineUp->DeepCopy(appendFilter->GetOutput());


        cutCircle->DeepCopy(ReorderContour(cutCircle));
        UniformSample(40, cutCircle);
        ResampledCircleGroup->AddMember(cutCircle);

        vtkSmartPointer<vtkPoints> newResampledCutCircle = vtkSmartPointer<vtkPoints>::New();
        for(vtkIdType j=0; j<cutCircle->GetNumberOfPoints(); j++)
        {
            double p[3];
            cutCircle->GetPoint(j, p);
            double vector[3];
            vtkMath::Subtract(p, oldp, vector);
            double coordinate[3];
            coordinate[0] = vtkMath::Dot(vector, told);
            coordinate[1] = vtkMath::Dot(vector, nold);
            coordinate[2] = vtkMath::Dot(vector, bold);
            double pp[3];
            double vx[3], vy[3], vz[3], tmp1[3], tmp2[3];
            NewTangents->GetTuple(i, vx);
            NewNormals->GetTuple(i, vy);
            NewBinormals->GetTuple(i, vz);
            vtkMath::MultiplyScalar(vx, coordinate[0]);
            vtkMath::MultiplyScalar(vy, coordinate[1]);
            vtkMath::MultiplyScalar(vz, coordinate[2]);
            vtkMath::Add(newp, vx, tmp1);
            vtkMath::Add(tmp1, vy, tmp2);
            vtkMath::Add(tmp2, vz, pp);
            newResampledCutCircle->InsertNextPoint(pp);
        }
        cutCircle->SetPoints(newResampledCutCircle);

        ResampledLineUpGroup->AddMember(cutCircle);
    }

    /*
    vtkSmartPointer<vtkPolyData> test1 = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPolyData> test2 = vtkSmartPointer<vtkPolyData>::New();
    test1->DeepCopy(ResampledCircleGroup->GetMember(0));
    test2->DeepCopy(ResampledCircleGroup->GetMember(1));

    vtkSmartPointer<vtkVertexGlyphFilter> testvertex1 = vtkSmartPointer<vtkVertexGlyphFilter>::New();
    testvertex1->SetInputData(test1);
    testvertex1->Update();
    vtkSmartPointer<vtkVertexGlyphFilter> testvertex2 = vtkSmartPointer<vtkVertexGlyphFilter>::New();
    testvertex2->SetInputData(test2);
    testvertex2->Update();

    appendFilter->RemoveAllInputs();
    appendFilter->AddInputConnection(testvertex1->GetOutputPort());
    appendFilter->AddInputConnection(testvertex2->GetOutputPort());
    appendFilter->Update();

    vtkSmartPointer<vtkPolyDataMapper> testMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    testMapper->SetInputConnection(appendFilter->GetOutputPort());
    testMapper->Update();
    vtkSmartPointer<vtkActor> testActor = vtkSmartPointer<vtkActor>::New();
    testActor->SetMapper(testMapper);
    testActor->GetProperty()->SetColor(1,0,0);
    testActor->GetProperty()->SetPointSize(5);
    //std::cout<<"test has "<<test->GetNumberOfPoints()<<endl;
    t_rendermanager->renderModel(testActor);
    */


    vtkSmartPointer<vtkPolyDataMapper> OriginCutCircleMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    OriginCutCircleMapper->SetInputData(OriginCutCircle);
    OriginCutCircleMapper->Update();
    vtkSmartPointer<vtkActor> OriginCutCircleActor = vtkSmartPointer<vtkActor>::New();
    OriginCutCircleActor->SetMapper(OriginCutCircleMapper);
    OriginCutCircleActor->GetProperty()->SetColor(0, 0.5, 1);
    OriginCutCircleActor->GetProperty()->SetLineWidth(2);
    t_rendermanager->renderModel(OriginCutCircleActor);

    vtkSmartPointer<vtkPolyDataMapper> CutCircleLineUpMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    CutCircleLineUpMapper->SetInputData(CutCircleLineUp);
    CutCircleLineUpMapper->Update();
    vtkSmartPointer<vtkActor> CutCircleLineUpActor = vtkSmartPointer<vtkActor>::New();
    CutCircleLineUpActor->SetMapper(CutCircleLineUpMapper);
    CutCircleLineUpActor->GetProperty()->SetColor(1, 0, 0);
    CutCircleLineUpActor->GetProperty()->SetLineWidth(2);
    t_rendermanager_right->renderModel(CutCircleLineUpActor);


    // deformation of surface begin
    // place the seed
    double p[3], t[3], pp[3], pseed[3];
    model->GetPoint(0, p);
    Tangents->GetTuple(0, t);
    vtkMath::MultiplyScalar(t, -50);
    vtkMath::Add(p, t, pp);
    vtkSmartPointer<vtkPointLocator> seedLocator = vtkSmartPointer<vtkPointLocator>::New();
    seedLocator->SetDataSet(t_colon);
    seedLocator->BuildLocator();
    vtkIdType seed = seedLocator->FindClosestPoint(pp);
    t_colon->GetPoint(seed, pseed);
    std::cout<<"Seed Point ID: "<<seed<<std::endl;
    vtkSmartPointer<vtkPoints> seedpoints = vtkSmartPointer<vtkPoints>::New();
    seedpoints->InsertNextPoint(pseed);
    seedpoints->InsertNextPoint(p);
    VisualizePoints(seedpoints, 0, 1, 0, 10, t_rendermanager);

    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(t_colon->GetNumberOfPoints());

    /*
    for(int i = 0; i < CircleGroup->GetNumOfMembers(); i++)
    {
        std::cout<<"circle "<<i<<" has "<<CircleGroup->GetMember(i)->GetNumberOfPoints()<<
                   "(resampled to "<<ResampledCircleGroup->GetMember(i)->GetNumberOfPoints()<<")"<<endl;
    }
    */

    vtkSmartPointer<vtkIdList> Sections = vtkSmartPointer<vtkIdList>::New();
    for(vtkIdType i = 0; i < t_colon->GetNumberOfPoints(); i++)
    {
        Sections->InsertNextId(-1);
    }
    Sections->SetId(seed, 0);
    GetSectionIds_loop(t_colon, seed, Sections, PlaneOriginals, PlaneNormals);
    //GetSectionIds_loop_v2(t_colon, seed, Sections, PlaneOriginals, PlaneNormals, CircleGroup);
    //GetSectionIds_loop_combinehighcurvatures(t_colon, seed, Sections, PlaneOriginals, PlaneNormals, Curvatures);

    vtkSmartPointer<vtkPolyData> source = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPolyData> target = vtkSmartPointer<vtkPolyData>::New();
    for(int i = 0; i <= model->GetNumberOfPoints(); i++)
    {
        if(i == 0)
        {
            appendFilter->RemoveAllInputs();
            appendFilter->AddInputData(ResampledCircleGroup->GetMember(0));
            appendFilter->AddInputData(ResampledCircleGroup->GetMember(1));
            appendFilter->Update();
            source->DeepCopy(appendFilter->GetOutput());
            appendFilter->RemoveAllInputs();
            appendFilter->AddInputData(ResampledLineUpGroup->GetMember(0));
            appendFilter->AddInputData(ResampledLineUpGroup->GetMember(1));
            appendFilter->Update();
            target->DeepCopy(appendFilter->GetOutput());
            //source->DeepCopy(ResampledCircleGroup->GetMember(i));
            //target->DeepCopy(ResampledLineUpGroup->GetMember(i));

        }
        else if(i == model->GetNumberOfPoints())
        {
            appendFilter->RemoveAllInputs();
            appendFilter->AddInputData(ResampledCircleGroup->GetMember(i-2));
            appendFilter->AddInputData(ResampledCircleGroup->GetMember(i-1));
            appendFilter->Update();
            source->DeepCopy(appendFilter->GetOutput());
            appendFilter->RemoveAllInputs();
            appendFilter->AddInputData(ResampledLineUpGroup->GetMember(i-2));
            appendFilter->AddInputData(ResampledLineUpGroup->GetMember(i-1));
            appendFilter->Update();
            target->DeepCopy(appendFilter->GetOutput());
            //source->DeepCopy(ResampledCircleGroup->GetMember(i-1));
            //target->DeepCopy(ResampledLineUpGroup->GetMember(i-1));
        }
        else
        {
            appendFilter->RemoveAllInputs();
            appendFilter->AddInputData(ResampledCircleGroup->GetMember(i-1));
            appendFilter->AddInputData(ResampledCircleGroup->GetMember(i));
            appendFilter->Update();
            source->DeepCopy(appendFilter->GetOutput());
            appendFilter->RemoveAllInputs();
            appendFilter->AddInputData(ResampledLineUpGroup->GetMember(i-1));
            appendFilter->AddInputData(ResampledLineUpGroup->GetMember(i));
            appendFilter->Update();
            target->DeepCopy(appendFilter->GetOutput());
        }
        vtkSmartPointer<vtkPoints> sectionpoints = vtkSmartPointer<vtkPoints>::New();
        vtkSmartPointer<vtkPolyData> sectionpoly = vtkSmartPointer<vtkPolyData>::New();
        vtkSmartPointer<vtkIdList> sectionids = vtkSmartPointer<vtkIdList>::New();
        for(vtkIdType j =0; j < t_colon->GetNumberOfPoints(); j++)
        {
            vtkIdType section = Sections->GetId(j);
            if(section == i)
            {
                double p[3];
                t_colon->GetPoint(j, p);
                sectionids->InsertNextId(j);
                sectionpoints->InsertNextPoint(p);
            }
        }
        sectionpoly->SetPoints(sectionpoints);
        vtkSmartPointer<vtkThinPlateSplineTransform> transform = vtkSmartPointer<vtkThinPlateSplineTransform>::New();
        transform->SetSourceLandmarks(source->GetPoints());
        transform->SetTargetLandmarks(target->GetPoints());
        transform->SetBasisToR();
        transform->Update();
        vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
        transformFilter->RemoveAllInputs();
        transformFilter->SetInputData(sectionpoly);
        transformFilter->SetTransform(transform);
        transformFilter->Update();
        sectionpoly->DeepCopy(transformFilter->GetOutput());

        for(vtkIdType j = 0; j < sectionpoly->GetNumberOfPoints(); j++)
        {
            double p[3];
            sectionpoly->GetPoint(j, p);
            points->SetPoint(sectionids->GetId(j), p);
        }

        //std::cout<<"section "<<i<<" "<<sectionpoly->GetNumberOfPoints()<<endl;
    }
    SurfaceLineUp->DeepCopy(t_colon);
    SurfaceLineUp->SetPoints(points);

    vtkSmartPointer<vtkPolyDataMapper> SurfaceLineUpMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    SurfaceLineUpMapper->SetInputData(SurfaceLineUp);
    SurfaceLineUpMapper->Update();
    vtkSmartPointer<vtkActor> SurfaceLineUpActor = vtkSmartPointer<vtkActor>::New();
    SurfaceLineUpActor->SetMapper(SurfaceLineUpMapper);
    //t_rendermanager->renderModel(SurfaceLineUpActor);
    t_filemanager->SaveFile(SurfaceLineUp, "SurfaceLineUp_v3_1.off");

    // output a deformation file
    pointLocator->RemoveAllObservers();
    pointLocator->SetDataSet(SurfaceLineUp);
    pointLocator->BuildLocator();
    bool* Is_Fixed = (bool*)malloc(sizeof(bool) * t_colon->GetNumberOfPoints());
    memset(Is_Fixed, 0, sizeof(bool) * t_colon->GetNumberOfPoints());
    //vtkSmartPointer<vtkDoubleArray> DeformationField = vtkSmartPointer<vtkDoubleArray>::New();
    for(vtkIdType i = 0; i < OriginCutCircle->GetNumberOfPoints(); i++)
    {
        double p[3];
        CutCircleLineUp->GetPoint(i, p);
        vtkIdType id = pointLocator->FindClosestPoint(p);
        Is_Fixed[id] = true;
    }

    bool * marked = (bool*)malloc(t_colon->GetNumberOfPoints()*sizeof(bool));
    memset(marked, 0, t_colon->GetNumberOfPoints()*sizeof(bool));
    vtkSmartPointer<vtkIntArray> lines = vtkSmartPointer<vtkIntArray>::New();
    lines->SetNumberOfComponents(2);
    for(int i=0; i<t_colon->GetNumberOfPoints(); i++)
    {
        vtkSmartPointer<vtkIdList> linkedpts = vtkSmartPointer<vtkIdList>::New();
        vtkSmartPointer<vtkIdList> cellids = vtkSmartPointer<vtkIdList>::New();
        t_colon->GetPointCells(i, cellids);
        for(int j=0; j<cellids->GetNumberOfIds(); j++)
        {
            vtkSmartPointer<vtkIdList> ptids = vtkSmartPointer<vtkIdList>::New();
            t_colon->GetCellPoints(cellids->GetId(j), ptids);
            assert(ptids->GetNumberOfIds() <= 3);
            for(int k = 0; k < ptids->GetNumberOfIds(); k++)
            {
                bool linked = false;
                int ptid = ptids->GetId(k);
                if(ptid == i || marked[ptid])
                    continue;
                for(int l = 0; l < linkedpts->GetNumberOfIds(); l++)
                {
                    if(linkedpts->GetId(l) == ptid)
                    {
                        linked = true;
                        break;
                    }
                }
                if(linked)
                    continue;
                int tuple[2];
                tuple[0] = i;
                tuple[1] = ptid;
                //std::cout<<i<<" "<<ptid<<endl;
                linkedpts->InsertNextId(ptid);
                lines->InsertNextTypedTuple(tuple);
            }
        }
        marked[i] = true;
    }
    for(int i=0; i < lines->GetNumberOfTuples(); i++)
    {
        int tuple[2];
        lines->GetTypedTuple(i, tuple);
        int idx1 = tuple[0];
        int idx2 = tuple[1];
        if(GetFacetsOfEdge(t_colon, idx1, idx2)->GetNumberOfIds() <= 1)
        {
            Is_Fixed[idx1] = true;
            Is_Fixed[idx2] = true;
        }
    }
    free(marked);

    ofstream file;
    file.open("./centerline_config.txt");
    assert(model->GetNumberOfPoints() == newcenterline->GetNumberOfPoints());
    assert(PlaneNormals->GetNumberOfTuples() == NewTangents->GetNumberOfTuples());
    file<<model->GetNumberOfPoints()<<endl;
    for(int i=0; i<model->GetNumberOfPoints(); i++)
    {
        double pold[3], nold[3], pnew[3], nnew[3];
        model->GetPoint(i, pold);
        newcenterline->GetPoint(i, pnew);
        PlaneNormals->GetTuple(i, nold);
        NewTangents->GetTuple(i, nnew);
        file<<pold[0]<<" "<<pold[1]<<" "<<pold[2]<<" ";
        file<<nold[0]<<" "<<nold[1]<<" "<<nold[2]<<" ";
        file<<pnew[0]<<" "<<pnew[1]<<" "<<pnew[2]<<" ";
        file<<nnew[0]<<" "<<nnew[1]<<" "<<nnew[2]<<std::endl;
    }
    file.close();

    //
    //visualize the fixed points
    /*
    vtkSmartPointer<vtkPoints> fixedPoints = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkPoints> fixedPoints_mapped = vtkSmartPointer<vtkPoints>::New();
    for(vtkIdType i = 0; i < t_colon->GetNumberOfPoints(); i++)
    {
        if(Is_Fixed[i])
        {
            double p[3], q[3];
            SurfaceLineUp->GetPoint(i, p);
            fixedPoints_mapped->InsertNextPoint(p);
            t_colon->GetPoint(i, q);
            fixedPoints->InsertNextPoint(q);
        }
    }
    VisualizePoints(fixedPoints, 1, 1, 0, 5, t_rendermanager);
    VisualizePoints(fixedPoints_mapped, 1, 1, 0, 5, t_rendermanager_right);
    */
    //
    // optimization

    // test
    /*
    vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
    //transform->RotateWXYZ(double angle, double x, double y, double z);
    //transform->RotateWXYZ(180, 0, 1, 0);
    transform->Translate(40, 0, 0);

    vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter =
            vtkSmartPointer<vtkTransformPolyDataFilter>::New();

    transformFilter->SetTransform(transform);
    transformFilter->SetInputData(t_colon);
    transformFilter->Update();
    SurfaceLineUp->DeepCopy(transformFilter->GetOutput());
    */

    int NumberOfRegions = 0;
    int* RegionIds = (int *)malloc(sizeof(int)*t_colon->GetNumberOfPoints());
    for(int i=0; i < t_colon->GetNumberOfPoints(); i++)
    {
        RegionIds[i] = -1;
    }
    GetRegionIds(t_colon, Is_Fixed, RegionIds, &NumberOfRegions);
    int count = 0;
    for(int i=0; i < t_colon->GetNumberOfPoints(); i++)
    {
        if(RegionIds[i] > -0.5)
        {
            count++;
        }
    }
    std::cout<<"count="<<count<<" Number Of Unfixed Regions="<<NumberOfRegions<<endl;

    // Affine Transform
    vtkSmartPointer<vtkPolyData> OptimizationInitial = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPoints> OptimizationInitialPoints = vtkSmartPointer<vtkPoints>::New();
    OptimizationInitial->DeepCopy(SurfaceLineUp);
    OptimizationInitialPoints->DeepCopy(OptimizationInitial->GetPoints());

    for(int i=0; i<NumberOfRegions; i++)
    {
        //std::cout<<"region "<<i<<endl;
        vtkSmartPointer<vtkIdList> regionIds = vtkSmartPointer<vtkIdList>::New();
        vtkSmartPointer<vtkPoints> regionPoints = vtkSmartPointer<vtkPoints>::New();
        vtkSmartPointer<vtkIdList> sourceIds = vtkSmartPointer<vtkIdList>::New();
        vtkSmartPointer<vtkPoints> sourcePoints = vtkSmartPointer<vtkPoints>::New();
        vtkSmartPointer<vtkPoints> targetPoints = vtkSmartPointer<vtkPoints>::New();
        for(int j=0; j<t_colon->GetNumberOfPoints(); j++)
        {
            if(RegionIds[j] == i)
            {
                double p[3];
                t_colon->GetPoint(j, p);
                regionIds->InsertNextId(j);
                regionPoints->InsertNextPoint(p);
                vtkSmartPointer<vtkIdList> ptids = vtkSmartPointer<vtkIdList>::New();
                ptids = GetConnectedVertices(t_colon, j);
                for(int k=0; k < ptids->GetNumberOfIds(); k++)
                {
                    int id = ptids->GetId(k);
                    if(!Is_Fixed[id])
                        continue;

                    bool exist = false;

                    for(int l=0; l < sourceIds->GetNumberOfIds(); l++)
                    {
                        if(sourceIds->GetId(l) == id)
                        {
                            exist = true;
                            break;
                        }
                    }
                    if(!exist)
                    {
                        double ps[3], pt[3];
                        t_colon->GetPoint(id, ps);
                        SurfaceLineUp->GetPoint(id, pt);

                        sourceIds->InsertNextId(id);
                        sourcePoints->InsertNextPoint(ps);
                        targetPoints->InsertNextPoint(pt);
                        regionIds->InsertNextId(id);
                        regionPoints->InsertNextPoint(ps);
                    }
                }
            }
        }
        vtkSmartPointer<vtkPolyData> regionPoly = vtkSmartPointer<vtkPolyData>::New();
        regionPoly->SetPoints(regionPoints);
        vtkSmartPointer<vtkVertexGlyphFilter> vertexFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
        vertexFilter->SetInputData(regionPoly);
        vertexFilter->Update();

        vtkSmartPointer<vtkLandmarkTransform> landmarkTransform = vtkSmartPointer<vtkLandmarkTransform>::New();
        //std::cout<<"numberoflandmarks="<<sourcePoints->GetNumberOfPoints()<<" "<<targetPoints->GetNumberOfPoints()<<endl;
        landmarkTransform->SetSourceLandmarks(sourcePoints);
        landmarkTransform->SetTargetLandmarks(targetPoints);
        landmarkTransform->SetModeToAffine();
        //landmarkTransform->SetModeToSimilarity();
        //landmarkTransform->SetModeToRigidBody();
        landmarkTransform->Update();
        vtkSmartPointer<vtkTransformFilter> transformFilter = vtkSmartPointer<vtkTransformFilter>::New();
        transformFilter->SetTransform(landmarkTransform);
        transformFilter->SetInputConnection(vertexFilter->GetOutputPort());
        transformFilter->Update();
        vtkSmartPointer<vtkPoints> affinePoints = vtkSmartPointer<vtkPoints>::New();
        affinePoints = transformFilter->GetOutput()->GetPoints();

        /*
        if(i==165)
        {
            VisualizePoints(regionPoints, 0, 1, 0, 5, t_rendermanager);
            VisualizePoints(affinePoints, 0, 1, 0, 5, t_rendermanager_right);
            VisualizePoints(sourcePoints, 1, 0, 0, 5, t_rendermanager);
            VisualizePoints(targetPoints, 1, 0, 0, 5, t_rendermanager_right);
        }
        double z[3];

        targetPoints->GetPoint(0, z);
        for(int well=0; well<affinePoints->GetNumberOfPoints(); well++)
        {
            double p[3];
            affinePoints->GetPoint(well, p);
            if(vtkMath::Distance2BetweenPoints(z,p) > 50*50)
                std::cout<<i<<endl;
        }
*/

        for(int ii=0; ii < affinePoints->GetNumberOfPoints(); ii++)
        {
            double p[3];
            int id = regionIds->GetId(ii);
            affinePoints->GetPoint(ii, p);
            //regionPoly->GetPoint(ii, p);
            OptimizationInitialPoints->SetPoint(id, p);
        }

        //std::cout<<"region: "<<i<<endl;
    }

    OptimizationInitial->SetPoints(OptimizationInitialPoints);
    t_filemanager->SaveFile(OptimizationInitial, "OptimizationInitial.off");
    vtkSmartPointer<vtkPolyData> OptimizedSurface = vtkSmartPointer<vtkPolyData>::New();
    //OptimizedSurface = Optimize(OptimizationInitial, SurfaceLineUp, Is_Fixed, t_rendermanager);
    OptimizedSurface->DeepCopy(OptimizationInitial);

/*
    vtkSmartPointer<vtkSmoothPolyDataFilter> smoothFilter = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
    smoothFilter->SetInputData(OptimizedSurface);
    smoothFilter->SetNumberOfIterations(20);
    smoothFilter->SetRelaxationFactor(0.1);
    smoothFilter->FeatureEdgeSmoothingOff();
    smoothFilter->BoundarySmoothingOn();
    smoothFilter->Update();
    OptimizedSurface->DeepCopy(smoothFilter->GetOutput());
*/

    vtkSmartPointer<vtkPolyDataMapper> optimizedMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    optimizedMapper->SetInputData(OptimizedSurface);//
    optimizedMapper->Update();
    vtkSmartPointer<vtkActor> optimizedActor = vtkSmartPointer<vtkActor>::New();
    optimizedActor->SetMapper(optimizedMapper);
    optimizedActor->GetProperty()->SetOpacity(0.5);
    //t_rendermanager_right->renderModel(optimizedActor);

    t_filemanager->SaveFile(OptimizedSurface, "OptimizedSurface.off");

    free(Is_Fixed);
    free(RegionIds);
    delete CircleGroup;
    delete ResampledCircleGroup;
    delete ResampledLineUpGroup;
    std::cout<<"Deformation End"<<endl;
    return OptimizedSurface;
}
vtkSmartPointer<vtkPolyData> Centerline::Deformation_v3_2(vtkSmartPointer<vtkDoubleArray> S, vtkSmartPointer<vtkDoubleArray> Curvatures, vtkSmartPointer<vtkIdList> CurvaturePointIds,
                                                     vtkSmartPointer<vtkDoubleArray> Tangents, vtkSmartPointer<vtkDoubleArray> Normals,
                                                     vtkSmartPointer<vtkPolyData> t_colon, RenderManager *t_rendermanager,
                                                     vtkSmartPointer<vtkDoubleArray> PlaneOriginals, vtkSmartPointer<vtkDoubleArray> PlaneNormals,
                                                     vtkSmartPointer<vtkDoubleArray> RefDirections, FileManager *t_filemanager)
{
    //PutNormalsOnSameSide(Normals, Curvatures);
    std::cout<<"Deformation"<<endl;
    bool straight = true;
    // Eliminate the torsion by growing the curve on a plane, according to: -dNnew/dSnew = -k*Tnew
    double point[3], nextpoint[3];
    double tangent[3], nexttangent[3];
    double normal[3], nextnormal[3];
    double binormal[3];
    double ds, curvature;
    vtkSmartPointer<vtkDoubleArray> NewTangents = vtkSmartPointer<vtkDoubleArray>::New();
    NewTangents->SetNumberOfComponents(3);
    vtkSmartPointer<vtkDoubleArray> NewNormals = vtkSmartPointer<vtkDoubleArray>::New();
    NewNormals->SetNumberOfComponents(3);
    vtkSmartPointer<vtkDoubleArray> NewBinormals = vtkSmartPointer<vtkDoubleArray>::New();
    NewBinormals->SetNumberOfComponents(3);
    vtkSmartPointer<vtkPoints> newpoints = vtkSmartPointer<vtkPoints>::New();
    for(vtkIdType i=0; i<model->GetNumberOfPoints(); i++)
    {
        if(i == 0)
        {
            model->GetPoint(0, point);
            point[0] = point[0] + 60;
            binormal[0] = 0;
            binormal[1] = 0;
            binormal[2] = 1;
            tangent[0] = 1;
            tangent[1] = 0;
            tangent[2] = 0;
            normal[0] = 0;
            normal[1] = 1;
            normal[2] = 0;
        }
        ds = (i != model->GetNumberOfPoints()-1)?(S->GetValue(i + 1) - S->GetValue(i)):0;
        curvature = Curvatures->GetValue(i);

        // Record the New Axis System
        NewTangents->InsertNextTuple(tangent);
        NewNormals->InsertNextTuple(normal);
        NewBinormals->InsertNextTuple(binormal);

        double dnormal[3];
        if(straight)
        {
            dnormal[0] = 0; dnormal[1] = 0; dnormal[2] = 0;
        }
        else
        {
            dnormal[0] = tangent[0]; dnormal[1] = tangent[1]; dnormal[2] = tangent[2];
            double RelaxationFactor = 2;
            vtkMath::MultiplyScalar(dnormal, -curvature / RelaxationFactor);
        }
        vtkMath::Add(normal, dnormal, nextnormal);
        vtkMath::Normalize(nextnormal);
        //std::cout<<i<<"\ts="<<S->GetValue(i)<<"("<<ds<<")"<<"\tk="<<curvature<<endl;
        vtkMath::Cross(nextnormal, binormal, nexttangent);
        vtkMath::Normalize(nexttangent);
        vtkMath::MultiplyScalar(tangent, ds);
        vtkMath::Add(point, tangent, nextpoint);
        newpoints->InsertNextPoint(point);

        point[0] = nextpoint[0]; point[1] = nextpoint[1]; point[2] = nextpoint[2];
        normal[0] = nextnormal[0]; normal[1] = nextnormal[1]; normal[2] = nextnormal[2];
        tangent[0] = nexttangent[0]; tangent[1] = nexttangent[1]; tangent[2] = nexttangent[2];
    }

    vtkSmartPointer<vtkPolyData> newcenterline = vtkSmartPointer<vtkPolyData>::New();
    newcenterline->SetPoints(newpoints);
    vtkSmartPointer<vtkVertexGlyphFilter> newVertexFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
    newVertexFilter->SetInputData(newcenterline);
    newVertexFilter->Update();
    vtkSmartPointer<vtkPolyDataMapper> newMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    newMapper->SetInputConnection(newVertexFilter->GetOutputPort());
    vtkSmartPointer<vtkActor> newActor = vtkSmartPointer<vtkActor>::New();
    newActor->SetMapper(newMapper);
    t_rendermanager->renderModel(newActor);

    // Line Up the Cross Sections
    vtkSmartPointer<vtkCutter> cutter = vtkSmartPointer<vtkCutter>::New();
    cutter->SetInputData(t_colon);
    vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();
    vtkSmartPointer<vtkPolyDataConnectivityFilter> connectivityFilter = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
    connectivityFilter->SetInputConnection(cutter->GetOutputPort());
    connectivityFilter->SetExtractionModeToClosestPointRegion();
    vtkSmartPointer<vtkAppendPolyData> appendFilter = vtkSmartPointer<vtkAppendPolyData>::New();

    vtkSmartPointer<vtkPolyData> OriginCutCircle = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPolyData> CutCircleLineUp = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPolyData> lastCircle = vtkSmartPointer<vtkPolyData>::New();

    vtkSmartPointer<vtkPolyData> SurfaceLineUp = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPointLocator> pointLocator = vtkSmartPointer<vtkPointLocator>::New();
    pointLocator->SetDataSet(t_colon);
    pointLocator->BuildLocator();

    vtkPolyDataGroup* CircleGroup = new vtkPolyDataGroup;
    vtkPolyDataGroup* ResampledCircleGroup = new vtkPolyDataGroup;
    vtkPolyDataGroup* ResampledLineUpGroup = new vtkPolyDataGroup;

    for(vtkIdType i = 0; i<model->GetNumberOfPoints(); i++)
    {
        double newp[3], oldp[3];
        newcenterline->GetPoint(i, newp);
        model->GetPoint(i, oldp);
        double told[3], nold[3], bold[3];
        Tangents->GetTuple(i, told);

        //Normals->GetTuple(i, nold);
        RefDirections->GetTuple(i, nold);

        vtkMath::Cross(told, nold, bold);

        plane->SetOrigin(PlaneOriginals->GetTuple(i));
        plane->SetNormal(PlaneNormals->GetTuple(i));
        cutter->SetCutFunction(plane);
        cutter->Update();

        connectivityFilter->SetClosestPoint(oldp);
        connectivityFilter->SetInputData(cutter->GetOutput());
        connectivityFilter->Update();

        vtkSmartPointer<vtkPolyData> cutCircle = vtkSmartPointer<vtkPolyData>::New();
        cutCircle = connectivityFilter->GetOutput();
        if(i != 0)
        {
            for(vtkIdType j = 0; j<lastCircle->GetNumberOfPoints(); j++)
            {

                if(cutCircle->GetNumberOfPoints() >= 70)
                    break;
                connectivityFilter->SetClosestPoint(lastCircle->GetPoint(j));
                connectivityFilter->Update();
                cutCircle = connectivityFilter->GetOutput();
            }
        }
        lastCircle->DeepCopy(cutCircle);

        CircleGroup->AddMember(cutCircle);// record the cutCircle into group

        appendFilter->RemoveAllInputs();
        appendFilter->AddInputData(cutCircle);
        appendFilter->AddInputData(OriginCutCircle);
        appendFilter->Update();
        OriginCutCircle->DeepCopy(appendFilter->GetOutput());

        vtkSmartPointer<vtkPoints> newCutCircle = vtkSmartPointer<vtkPoints>::New();
        vtkSmartPointer<vtkPolyData> newCutCirclePoly = vtkSmartPointer<vtkPolyData>::New();
        newCutCirclePoly->DeepCopy(cutCircle);
        for(vtkIdType j=0; j<cutCircle->GetNumberOfPoints(); j++)
        {
            double p[3];
            cutCircle->GetPoint(j, p);
            double vector[3];
            vtkMath::Subtract(p, oldp, vector);
            double coordinate[3];
            coordinate[0] = vtkMath::Dot(vector, told);
            coordinate[1] = vtkMath::Dot(vector, nold);
            coordinate[2] = vtkMath::Dot(vector, bold);
            double pp[3];
            double vx[3], vy[3], vz[3], tmp1[3], tmp2[3];
            NewTangents->GetTuple(i, vx);
            NewNormals->GetTuple(i, vy);
            NewBinormals->GetTuple(i, vz);
            vtkMath::MultiplyScalar(vx, coordinate[0]);
            vtkMath::MultiplyScalar(vy, coordinate[1]);
            vtkMath::MultiplyScalar(vz, coordinate[2]);
            vtkMath::Add(newp, vx, tmp1);
            vtkMath::Add(tmp1, vy, tmp2);
            vtkMath::Add(tmp2, vz, pp);
            newCutCircle->InsertNextPoint(pp);
        }
        newCutCirclePoly->SetPoints(newCutCircle);
        appendFilter->RemoveAllInputs();
        appendFilter->AddInputData(CutCircleLineUp);
        appendFilter->AddInputData(newCutCirclePoly);
        appendFilter->Update();
        CutCircleLineUp->DeepCopy(appendFilter->GetOutput());


        cutCircle->DeepCopy(ReorderContour(cutCircle));
        UniformSample(40, cutCircle);
        ResampledCircleGroup->AddMember(cutCircle);

        vtkSmartPointer<vtkPoints> newResampledCutCircle = vtkSmartPointer<vtkPoints>::New();
        for(vtkIdType j=0; j<cutCircle->GetNumberOfPoints(); j++)
        {
            double p[3];
            cutCircle->GetPoint(j, p);
            double vector[3];
            vtkMath::Subtract(p, oldp, vector);
            double coordinate[3];
            coordinate[0] = vtkMath::Dot(vector, told);
            coordinate[1] = vtkMath::Dot(vector, nold);
            coordinate[2] = vtkMath::Dot(vector, bold);
            double pp[3];
            double vx[3], vy[3], vz[3], tmp1[3], tmp2[3];
            NewTangents->GetTuple(i, vx);
            NewNormals->GetTuple(i, vy);
            NewBinormals->GetTuple(i, vz);
            vtkMath::MultiplyScalar(vx, coordinate[0]);
            vtkMath::MultiplyScalar(vy, coordinate[1]);
            vtkMath::MultiplyScalar(vz, coordinate[2]);
            vtkMath::Add(newp, vx, tmp1);
            vtkMath::Add(tmp1, vy, tmp2);
            vtkMath::Add(tmp2, vz, pp);
            newResampledCutCircle->InsertNextPoint(pp);
        }
        cutCircle->SetPoints(newResampledCutCircle);

        ResampledLineUpGroup->AddMember(cutCircle);
    }

    vtkSmartPointer<vtkPolyDataMapper> OriginCutCircleMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    OriginCutCircleMapper->SetInputData(OriginCutCircle);
    OriginCutCircleMapper->Update();
    vtkSmartPointer<vtkActor> OriginCutCircleActor = vtkSmartPointer<vtkActor>::New();
    OriginCutCircleActor->SetMapper(OriginCutCircleMapper);
    OriginCutCircleActor->GetProperty()->SetColor(0, 0, 1);
    t_rendermanager->renderModel(OriginCutCircleActor);

    vtkSmartPointer<vtkPolyDataMapper> CutCircleLineUpMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    CutCircleLineUpMapper->SetInputData(CutCircleLineUp);
    CutCircleLineUpMapper->Update();
    vtkSmartPointer<vtkActor> CutCircleLineUpActor = vtkSmartPointer<vtkActor>::New();
    CutCircleLineUpActor->SetMapper(CutCircleLineUpMapper);
    CutCircleLineUpActor->GetProperty()->SetColor(1, 0, 0);
    t_rendermanager->renderModel(CutCircleLineUpActor);

    // deformation of surface begin
    // place the seed
    double p[3], t[3], pp[3], pseed[3];
    model->GetPoint(0, p);
    Tangents->GetTuple(0, t);
    vtkMath::MultiplyScalar(t, -50);
    vtkMath::Add(p, t, pp);
    vtkSmartPointer<vtkPointLocator> seedLocator = vtkSmartPointer<vtkPointLocator>::New();
    seedLocator->SetDataSet(t_colon);
    seedLocator->BuildLocator();
    vtkIdType seed = seedLocator->FindClosestPoint(pp);
    t_colon->GetPoint(seed, pseed);
    std::cout<<"Seed Point ID: "<<seed<<std::endl;
    vtkSmartPointer<vtkPoints> seedpoints = vtkSmartPointer<vtkPoints>::New();
    seedpoints->InsertNextPoint(pseed);
    vtkSmartPointer<vtkPolyData> seedpoly = vtkSmartPointer<vtkPolyData>::New();
    seedpoly->SetPoints(seedpoints);
    vtkSmartPointer<vtkVertexGlyphFilter> seedVertex = vtkSmartPointer<vtkVertexGlyphFilter>::New();
    seedVertex->SetInputData(seedpoly);
    seedVertex->Update();
    vtkSmartPointer<vtkPolyDataMapper> seedMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    seedMapper->SetInputConnection(seedVertex->GetOutputPort());
    seedMapper->Update();
    vtkSmartPointer<vtkActor> seedActor = vtkSmartPointer<vtkActor>::New();
    seedActor->SetMapper(seedMapper);
    seedActor->GetProperty()->SetPointSize(10);
    seedActor->GetProperty()->SetColor(0,1,0);
    t_rendermanager->renderModel(seedActor); // place the seed

    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(t_colon->GetNumberOfPoints());

    /*
    for(int i = 0; i < CircleGroup->GetNumOfMembers(); i++)
    {
        std::cout<<"circle "<<i<<" has "<<CircleGroup->GetMember(i)->GetNumberOfPoints()<<
                   "(resampled to "<<ResampledCircleGroup->GetMember(i)->GetNumberOfPoints()<<")"<<endl;
    }
    */

    vtkSmartPointer<vtkIdList> Sections = vtkSmartPointer<vtkIdList>::New();
    for(vtkIdType i = 0; i < t_colon->GetNumberOfPoints(); i++)
    {
        Sections->InsertNextId(-1);
    }
    Sections->SetId(seed, 0);
    GetSectionIds_loop(t_colon, seed, Sections, PlaneOriginals, PlaneNormals);
    //GetSectionIds_loop_v2(t_colon, seed, Sections, PlaneOriginals, PlaneNormals, CircleGroup);
    //GetSectionIds_loop_combinehighcurvatures(t_colon, seed, Sections, PlaneOriginals, PlaneNormals, Curvatures);

    vtkSmartPointer<vtkPolyData> source = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPolyData> target = vtkSmartPointer<vtkPolyData>::New();
    for(int i = 0; i <= model->GetNumberOfPoints(); i++)
    {
        if(i == 0)
        {
            appendFilter->RemoveAllInputs();
            appendFilter->AddInputData(ResampledCircleGroup->GetMember(0));
            appendFilter->AddInputData(ResampledCircleGroup->GetMember(1));
            appendFilter->Update();
            source->DeepCopy(appendFilter->GetOutput());
            appendFilter->RemoveAllInputs();
            appendFilter->AddInputData(ResampledLineUpGroup->GetMember(0));
            appendFilter->AddInputData(ResampledLineUpGroup->GetMember(1));
            appendFilter->Update();
            target->DeepCopy(appendFilter->GetOutput());
            //source->DeepCopy(ResampledCircleGroup->GetMember(i));
            //target->DeepCopy(ResampledLineUpGroup->GetMember(i));

        }
        else if(i == model->GetNumberOfPoints())
        {
            appendFilter->RemoveAllInputs();
            appendFilter->AddInputData(ResampledCircleGroup->GetMember(i-2));
            appendFilter->AddInputData(ResampledCircleGroup->GetMember(i-1));
            appendFilter->Update();
            source->DeepCopy(appendFilter->GetOutput());
            appendFilter->RemoveAllInputs();
            appendFilter->AddInputData(ResampledLineUpGroup->GetMember(i-2));
            appendFilter->AddInputData(ResampledLineUpGroup->GetMember(i-1));
            appendFilter->Update();
            target->DeepCopy(appendFilter->GetOutput());
            //source->DeepCopy(ResampledCircleGroup->GetMember(i-1));
            //target->DeepCopy(ResampledLineUpGroup->GetMember(i-1));
        }
        else
        {
            appendFilter->RemoveAllInputs();
            appendFilter->AddInputData(ResampledCircleGroup->GetMember(i-1));
            appendFilter->AddInputData(ResampledCircleGroup->GetMember(i));
            appendFilter->Update();
            source->DeepCopy(appendFilter->GetOutput());
            appendFilter->RemoveAllInputs();
            appendFilter->AddInputData(ResampledLineUpGroup->GetMember(i-1));
            appendFilter->AddInputData(ResampledLineUpGroup->GetMember(i));
            appendFilter->Update();
            target->DeepCopy(appendFilter->GetOutput());
        }
        vtkSmartPointer<vtkPoints> sectionpoints = vtkSmartPointer<vtkPoints>::New();
        vtkSmartPointer<vtkPolyData> sectionpoly = vtkSmartPointer<vtkPolyData>::New();
        vtkSmartPointer<vtkIdList> sectionids = vtkSmartPointer<vtkIdList>::New();
        for(vtkIdType j =0; j < t_colon->GetNumberOfPoints(); j++)
        {
            vtkIdType section = Sections->GetId(j);
            if(section == i)
            {
                double p[3];
                t_colon->GetPoint(j, p);
                sectionids->InsertNextId(j);
                sectionpoints->InsertNextPoint(p);
            }
        }
        sectionpoly->SetPoints(sectionpoints);
        vtkSmartPointer<vtkThinPlateSplineTransform> transform = vtkSmartPointer<vtkThinPlateSplineTransform>::New();
        transform->SetSourceLandmarks(source->GetPoints());
        transform->SetTargetLandmarks(target->GetPoints());
        transform->SetBasisToR();
        transform->Update();
        vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
        transformFilter->RemoveAllInputs();
        transformFilter->SetInputData(sectionpoly);
        transformFilter->SetTransform(transform);
        transformFilter->Update();
        sectionpoly->DeepCopy(transformFilter->GetOutput());


        for(vtkIdType j = 0; j < sectionpoly->GetNumberOfPoints(); j++)
        {
            double p[3];
            sectionpoly->GetPoint(j, p);
            points->SetPoint(sectionids->GetId(j), p);
        }

        //std::cout<<"section "<<i<<" "<<sectionpoly->GetNumberOfPoints()<<endl;
    }
    SurfaceLineUp->DeepCopy(t_colon);
    SurfaceLineUp->SetPoints(points);

    vtkSmartPointer<vtkPolyDataMapper> SurfaceLineUpMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    SurfaceLineUpMapper->SetInputData(SurfaceLineUp);
    SurfaceLineUpMapper->Update();
    vtkSmartPointer<vtkActor> SurfaceLineUpActor = vtkSmartPointer<vtkActor>::New();
    SurfaceLineUpActor->SetMapper(SurfaceLineUpMapper);
    //t_rendermanager->renderModel(SurfaceLineUpActor);
    //t_filemanager->SaveFile(SurfaceLineUp, "SurfaceLineUp_v3_1.stl");

    pointLocator->RemoveAllObservers();
    pointLocator->SetDataSet(SurfaceLineUp);
    pointLocator->BuildLocator();
    bool* Is_Fixed = (bool*)malloc(sizeof(bool) * t_colon->GetNumberOfPoints());
    memset(Is_Fixed, 0, sizeof(bool) * t_colon->GetNumberOfPoints());
    //vtkSmartPointer<vtkDoubleArray> DeformationField = vtkSmartPointer<vtkDoubleArray>::New();
    for(vtkIdType i = 0; i < OriginCutCircle->GetNumberOfPoints(); i++)
    {
        double p[3], dist2;
        CutCircleLineUp->GetPoint(i, p);
        //vtkIdType id = pointLocator->FindClosestPoint(p);
        vtkIdType id = pointLocator->FindClosestPointWithinRadius(1.5, p, dist2);
        if(id >=0 && id < t_colon->GetNumberOfPoints())
            Is_Fixed[id] = true;
    }

    // test
    /*
    for(int i=0; i< Sections->GetNumberOfIds(); i++)
    {
        if(Sections->GetId(i) != 0)
        {
            Is_Fixed[i] = true;
        }
    }
    */

    bool * marked = (bool*)malloc(t_colon->GetNumberOfPoints()*sizeof(bool));
    memset(marked, 0, t_colon->GetNumberOfPoints()*sizeof(bool));
    vtkSmartPointer<vtkIntArray> lines = vtkSmartPointer<vtkIntArray>::New();
    lines->SetNumberOfComponents(2);
    for(int i=0; i<t_colon->GetNumberOfPoints(); i++)
    {
        vtkSmartPointer<vtkIdList> linkedpts = vtkSmartPointer<vtkIdList>::New();
        vtkSmartPointer<vtkIdList> cellids = vtkSmartPointer<vtkIdList>::New();
        t_colon->GetPointCells(i, cellids);
        for(int j=0; j<cellids->GetNumberOfIds(); j++)
        {
            vtkSmartPointer<vtkIdList> ptids = vtkSmartPointer<vtkIdList>::New();
            t_colon->GetCellPoints(cellids->GetId(j), ptids);
            assert(ptids->GetNumberOfIds() <= 3);
            for(int k = 0; k < ptids->GetNumberOfIds(); k++)
            {
                bool linked = false;
                int ptid = ptids->GetId(k);
                if(ptid == i || marked[ptid])
                    continue;
                for(int l = 0; l < linkedpts->GetNumberOfIds(); l++)
                {
                    if(linkedpts->GetId(l) == ptid)
                    {
                        linked = true;
                        break;
                    }
                }
                if(linked)
                    continue;
                int tuple[2];
                tuple[0] = i;
                tuple[1] = ptid;
                //std::cout<<i<<" "<<ptid<<endl;
                linkedpts->InsertNextId(ptid);
                lines->InsertNextTypedTuple(tuple);
            }
        }
        marked[i] = true;
    }
    for(int i=0; i < lines->GetNumberOfTuples(); i++)
    {
        int tuple[2];
        lines->GetTypedTuple(i, tuple);
        int idx1 = tuple[0];
        int idx2 = tuple[1];
        if(GetFacetsOfEdge(t_colon, idx1, idx2)->GetNumberOfIds() <= 1)
        {
            Is_Fixed[idx1] = true;
            Is_Fixed[idx2] = true;
        }
    }
    free(marked);

    //visualize the fixed points

    vtkSmartPointer<vtkPoints> fixedPoints = vtkSmartPointer<vtkPoints>::New();
    for(vtkIdType i = 1; i < t_colon->GetNumberOfPoints(); i++)
    {
        if(Is_Fixed[i])
        {
            double p[3], q[3];
            SurfaceLineUp->GetPoint(i, p);
            fixedPoints->InsertNextPoint(p);
            t_colon->GetPoint(i, q);
            fixedPoints->InsertNextPoint(q);
        }
    }
    vtkSmartPointer<vtkPolyData> fixedPoly = vtkSmartPointer<vtkPolyData>::New();
    fixedPoly->SetPoints(fixedPoints);
    vtkSmartPointer<vtkVertexGlyphFilter> vertexFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
    vertexFilter->SetInputData(fixedPoly);
    vertexFilter->Update();
    vtkSmartPointer<vtkPolyDataMapper> fixedMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    fixedMapper->SetInputConnection(vertexFilter->GetOutputPort());
    fixedMapper->Update();
    vtkSmartPointer<vtkActor> fixedActor = vtkSmartPointer<vtkActor>::New();
    fixedActor->SetMapper(fixedMapper);
    fixedActor->GetProperty()->SetPointSize(5);
    fixedActor->GetProperty()->SetColor(1,1,0);
    t_rendermanager->renderModel(fixedActor);

    //
    // optimization
    vtkSmartPointer<vtkPolyData> OptimizedSurface = vtkSmartPointer<vtkPolyData>::New();
    OptimizedSurface = Optimize(t_colon, SurfaceLineUp, Is_Fixed, t_rendermanager);

    // Modify the Section based on the Optimized Surface
    vtkSmartPointer<vtkIdList> NewSections = vtkSmartPointer<vtkIdList>::New();
    for(vtkIdType i = 0; i < t_colon->GetNumberOfPoints(); i++)
    {
        NewSections->InsertNextId(-1);
    }
    NewSections->SetId(seed, 0);
    vtkSmartPointer<vtkDoubleArray> NewPlaneOriginals = vtkSmartPointer<vtkDoubleArray>::New();
    for(int i=0; i<newcenterline->GetNumberOfPoints(); i++)
    {
        double p[3];
        newcenterline->GetPoint(i, p);
        NewPlaneOriginals->InsertNextTuple(p);
    }
    GetSectionIds_loop(OptimizedSurface, seed, NewSections, NewPlaneOriginals, NewTangents);

    // New Piece-wise thin plate spline transform
    vtkSmartPointer<vtkPolyData> Result = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPoints> resultpoints = vtkSmartPointer<vtkPoints>::New();
    resultpoints->SetNumberOfPoints(t_colon->GetNumberOfPoints());
    for(int i = 0; i <= model->GetNumberOfPoints(); i++)
    {
        if(i == 0)
        {
            appendFilter->RemoveAllInputs();
            appendFilter->AddInputData(ResampledCircleGroup->GetMember(0));
            appendFilter->AddInputData(ResampledCircleGroup->GetMember(1));
            appendFilter->Update();
            source->DeepCopy(appendFilter->GetOutput());
            appendFilter->RemoveAllInputs();
            appendFilter->AddInputData(ResampledLineUpGroup->GetMember(0));
            appendFilter->AddInputData(ResampledLineUpGroup->GetMember(1));
            appendFilter->Update();
            target->DeepCopy(appendFilter->GetOutput());
            //source->DeepCopy(ResampledCircleGroup->GetMember(i));
            //target->DeepCopy(ResampledLineUpGroup->GetMember(i));
        }
        else if(i == model->GetNumberOfPoints())
        {
            appendFilter->RemoveAllInputs();
            appendFilter->AddInputData(ResampledCircleGroup->GetMember(i-2));
            appendFilter->AddInputData(ResampledCircleGroup->GetMember(i-1));
            appendFilter->Update();
            source->DeepCopy(appendFilter->GetOutput());
            appendFilter->RemoveAllInputs();
            appendFilter->AddInputData(ResampledLineUpGroup->GetMember(i-2));
            appendFilter->AddInputData(ResampledLineUpGroup->GetMember(i-1));
            appendFilter->Update();
            target->DeepCopy(appendFilter->GetOutput());
            //source->DeepCopy(ResampledCircleGroup->GetMember(i-1));
            //target->DeepCopy(ResampledLineUpGroup->GetMember(i-1));
        }
        else
        {
            int range = 2;
            appendFilter->RemoveAllInputs();
            for(int k=-range; k <= range-1; k++)
            {
                if( i + k >= 0 && i + k < model->GetNumberOfPoints())
                {
                    appendFilter->AddInputData(ResampledCircleGroup->GetMember(i + k));
                }
            }
            appendFilter->Update();
            source->DeepCopy(appendFilter->GetOutput());

            appendFilter->RemoveAllInputs();
            for(int k=-range; k <= range-1; k++)
            {
                if( i + k >= 0 && i + k < model->GetNumberOfPoints())
                {
                    appendFilter->AddInputData(ResampledLineUpGroup->GetMember(i + k));
                }
            }
            appendFilter->Update();
            target->DeepCopy(appendFilter->GetOutput());
        }
        vtkSmartPointer<vtkPoints> sectionpoints = vtkSmartPointer<vtkPoints>::New();
        vtkSmartPointer<vtkPolyData> sectionpoly = vtkSmartPointer<vtkPolyData>::New();
        vtkSmartPointer<vtkIdList> sectionids = vtkSmartPointer<vtkIdList>::New();
        for(vtkIdType j =0; j < t_colon->GetNumberOfPoints(); j++)
        {
            vtkIdType section = NewSections->GetId(j);
            if(section == i)
            {
                double p[3];
                t_colon->GetPoint(j, p);
                sectionids->InsertNextId(j);
                sectionpoints->InsertNextPoint(p);
            }
        }
        sectionpoly->SetPoints(sectionpoints);
        vtkSmartPointer<vtkThinPlateSplineTransform> transform = vtkSmartPointer<vtkThinPlateSplineTransform>::New();
        transform->SetSourceLandmarks(source->GetPoints());
        transform->SetTargetLandmarks(target->GetPoints());
        transform->SetBasisToR();
        transform->Update();
        vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
        transformFilter->RemoveAllInputs();
        transformFilter->SetInputData(sectionpoly);
        transformFilter->SetTransform(transform);
        transformFilter->Update();
        sectionpoly->DeepCopy(transformFilter->GetOutput());
        for(vtkIdType j = 0; j < sectionpoly->GetNumberOfPoints(); j++)
        {
            double p[3];
            sectionpoly->GetPoint(j, p);
            resultpoints->SetPoint(sectionids->GetId(j), p);
        }
        std::cout<<"section "<<i<<" "<<sectionpoly->GetNumberOfPoints()<<endl;
    }
    Result->DeepCopy(t_colon);
    Result->SetPoints(resultpoints);

    vtkSmartPointer<vtkSmoothPolyDataFilter> smoothFilter = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
    smoothFilter->SetInputData(Result);
    smoothFilter->SetNumberOfIterations(15);
    smoothFilter->SetRelaxationFactor(0.1);
    smoothFilter->FeatureEdgeSmoothingOff();
    smoothFilter->BoundarySmoothingOn();
    smoothFilter->Update();
    Result = smoothFilter->GetOutput();

    vtkSmartPointer<vtkPolyDataMapper> ResultMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    ResultMapper->SetInputData(Result);
    ResultMapper->Update();
    vtkSmartPointer<vtkActor> ResultActor = vtkSmartPointer<vtkActor>::New();
    ResultActor->SetMapper(ResultMapper);
    t_rendermanager->renderModel(ResultActor);

    t_filemanager->SaveFile(Result, "Result_v3_2.off");


    free(Is_Fixed);
    delete CircleGroup;
    delete ResampledCircleGroup;
    delete ResampledLineUpGroup;
    std::cout<<"Deformation End"<<endl;
    return SurfaceLineUp;
}
vtkSmartPointer<vtkPolyData> Centerline::Deformation_v4(vtkSmartPointer<vtkDoubleArray> S, vtkSmartPointer<vtkDoubleArray> Curvatures,
                                                     vtkSmartPointer<vtkDoubleArray> Tangents, vtkSmartPointer<vtkDoubleArray> Normals,
                                                     vtkSmartPointer<vtkPolyData> t_colon, RenderManager *t_rendermanager,
                                                     vtkSmartPointer<vtkDoubleArray> PlaneOriginals, vtkSmartPointer<vtkDoubleArray> PlaneNormals,
                                                     vtkSmartPointer<vtkDoubleArray> RefDirections, FileManager *t_filemanager)
{
    //PutNormalsOnSameSide(Normals, Curvatures);
    std::cout<<"Deformation"<<endl;
    bool straight = true;
    // Eliminate the torsion by growing the curve on a plane, according to: -dNnew/dSnew = -k*Tnew
    double point[3], nextpoint[3];
    double tangent[3], nexttangent[3];
    double normal[3], nextnormal[3];
    double binormal[3];
    double ds, curvature;
    vtkSmartPointer<vtkDoubleArray> NewTangents = vtkSmartPointer<vtkDoubleArray>::New();
    NewTangents->SetNumberOfComponents(3);
    vtkSmartPointer<vtkDoubleArray> NewNormals = vtkSmartPointer<vtkDoubleArray>::New();
    NewNormals->SetNumberOfComponents(3);
    vtkSmartPointer<vtkDoubleArray> NewBinormals = vtkSmartPointer<vtkDoubleArray>::New();
    NewBinormals->SetNumberOfComponents(3);
    vtkSmartPointer<vtkPoints> newpoints = vtkSmartPointer<vtkPoints>::New();
    for(vtkIdType i=0; i<model->GetNumberOfPoints(); i++)
    {
        if(i == 0)
        {
            model->GetPoint(0, point);
            point[0] = point[0] + 100;
            binormal[0] = 0;
            binormal[1] = 0;
            binormal[2] = 1;
            tangent[0] = 1;
            tangent[1] = 0;
            tangent[2] = 0;
            normal[0] = 0;
            normal[1] = 1;
            normal[2] = 0;
        }
        ds = (i != model->GetNumberOfPoints()-1)?(S->GetValue(i + 1) - S->GetValue(i)):0;
        curvature = Curvatures->GetValue(i);

        // Record the New Axis System
        NewTangents->InsertNextTuple(tangent);
        NewNormals->InsertNextTuple(normal);
        NewBinormals->InsertNextTuple(binormal);

        double dnormal[3];
        if(straight)
        {
            dnormal[0] = 0; dnormal[1] = 0; dnormal[2] = 0;
        }
        else
        {
            dnormal[0] = tangent[0]; dnormal[1] = tangent[1]; dnormal[2] = tangent[2];
            double RelaxationFactor = 2;
            vtkMath::MultiplyScalar(dnormal, -curvature / RelaxationFactor);
        }
        vtkMath::Add(normal, dnormal, nextnormal);
        vtkMath::Normalize(nextnormal);
        //std::cout<<i<<"\ts="<<S->GetValue(i)<<"("<<ds<<")"<<"\tk="<<curvature<<endl;
        vtkMath::Cross(nextnormal, binormal, nexttangent);
        vtkMath::Normalize(nexttangent);
        vtkMath::MultiplyScalar(tangent, ds);
        vtkMath::Add(point, tangent, nextpoint);
        newpoints->InsertNextPoint(point);

        point[0] = nextpoint[0]; point[1] = nextpoint[1]; point[2] = nextpoint[2];
        normal[0] = nextnormal[0]; normal[1] = nextnormal[1]; normal[2] = nextnormal[2];
        tangent[0] = nexttangent[0]; tangent[1] = nexttangent[1]; tangent[2] = nexttangent[2];
    }

    vtkSmartPointer<vtkPolyData> newcenterline = vtkSmartPointer<vtkPolyData>::New();
    newcenterline->SetPoints(newpoints);
    vtkSmartPointer<vtkVertexGlyphFilter> newVertexFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
    newVertexFilter->SetInputData(newcenterline);
    newVertexFilter->Update();
    vtkSmartPointer<vtkPolyDataMapper> newMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    newMapper->SetInputConnection(newVertexFilter->GetOutputPort());
    vtkSmartPointer<vtkActor> newActor = vtkSmartPointer<vtkActor>::New();
    newActor->SetMapper(newMapper);
    t_rendermanager->renderModel(newActor);

    // Line Up the Cross Sections
    vtkSmartPointer<vtkCutter> cutter = vtkSmartPointer<vtkCutter>::New();
    cutter->SetInputData(t_colon);

    vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();

    vtkSmartPointer<vtkPolyDataConnectivityFilter> connectivityFilter = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
    connectivityFilter->SetInputConnection(cutter->GetOutputPort());
    connectivityFilter->SetExtractionModeToClosestPointRegion();

    //vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();

    vtkSmartPointer<vtkAppendPolyData> appendFilter = vtkSmartPointer<vtkAppendPolyData>::New();
    vtkSmartPointer<vtkCleanPolyData> cleanFilter = vtkSmartPointer<vtkCleanPolyData>::New();

    vtkSmartPointer<vtkPolyData> CutCircleLineUp = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPolyData> lastCircle = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPolyData> lastNewCircle = vtkSmartPointer<vtkPolyData>::New();

    vtkSmartPointer<vtkPolyData> CutCircleOrigin = vtkSmartPointer<vtkPolyData>::New();

    vtkSmartPointer<vtkPolyData> SurfaceLineUp = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPointLocator> pointLocator = vtkSmartPointer<vtkPointLocator>::New();
    pointLocator->SetDataSet(t_colon);
    pointLocator->BuildLocator();

    bool* Is_Fixed = (bool*)malloc(sizeof(bool) * t_colon->GetNumberOfPoints());
    memset(Is_Fixed, 0, sizeof(bool) * t_colon->GetNumberOfPoints());
    vtkSmartPointer<vtkDoubleArray> DeformationField = vtkSmartPointer<vtkDoubleArray>::New();
    DeformationField->SetNumberOfComponents(3);
    DeformationField->SetNumberOfTuples(t_colon->GetNumberOfPoints());
    vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    vtkSmartPointer<vtkPolyData> piece = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPolyData> transformedpiece = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPolyData> lastpiece = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPolyData> source = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPolyData> target = vtkSmartPointer<vtkPolyData>::New();

    for(vtkIdType i = 0; i<model->GetNumberOfPoints(); i++)
    {
        //vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
        double newp[3], oldp[3];
        newcenterline->GetPoint(i, newp);
        model->GetPoint(i, oldp);
        double told[3], nold[3], bold[3];
        Tangents->GetTuple(i, told);

        //Normals->GetTuple(i, nold);
        RefDirections->GetTuple(i, nold);

        vtkMath::Cross(told, nold, bold);

        plane->SetOrigin(PlaneOriginals->GetTuple(i));
        plane->SetNormal(PlaneNormals->GetTuple(i));
        cutter->SetCutFunction(plane);
        cutter->Update();

        connectivityFilter->SetClosestPoint(oldp);
        connectivityFilter->SetInputData(cutter->GetOutput());
        connectivityFilter->Update();

        vtkSmartPointer<vtkPolyData> cutCircle = vtkSmartPointer<vtkPolyData>::New();
        cutCircle = connectivityFilter->GetOutput();
        //std::cout<<i<<" "<<"cutline points:"<<cutCircle->GetNumberOfPoints()<<endl;

        if(i != 0)
        {
            for(vtkIdType j = 0; j<lastCircle->GetNumberOfPoints(); j++)
            {
                if(cutCircle->GetNumberOfPoints() >= 100)
                    break;
                connectivityFilter->SetClosestPoint(lastCircle->GetPoint(j));
                connectivityFilter->Update();
                cutCircle = connectivityFilter->GetOutput();
            }
        }

        // build source landmarks
        //cutCircle->DeepCopy(ReorderContour(cutCircle));
        //UniformSample(80, cutCircle);

        appendFilter->RemoveAllInputs();
        appendFilter->AddInputData(cutCircle);
        appendFilter->AddInputData(lastCircle);
        appendFilter->Update();
        source->ShallowCopy(appendFilter->GetOutput());


        lastCircle->DeepCopy(cutCircle);

        vtkSmartPointer<vtkPoints> newCutCircle = vtkSmartPointer<vtkPoints>::New();

        for(vtkIdType j=0; j<cutCircle->GetNumberOfPoints(); j++)
        {
            double p[3];
            cutCircle->GetPoint(j, p);
            double vector[3];
            vtkMath::Subtract(p, oldp, vector);
            double coordinate[3];
            coordinate[0] = vtkMath::Dot(vector, told);
            coordinate[1] = vtkMath::Dot(vector, nold);
            coordinate[2] = vtkMath::Dot(vector, bold);
            double pp[3];
            double vx[3], vy[3], vz[3], tmp1[3], tmp2[3];
            NewTangents->GetTuple(i, vx);
            NewNormals->GetTuple(i, vy);
            NewBinormals->GetTuple(i, vz);
            vtkMath::MultiplyScalar(vx, coordinate[0]);
            vtkMath::MultiplyScalar(vy, coordinate[1]);
            vtkMath::MultiplyScalar(vz, coordinate[2]);
            vtkMath::Add(newp, vx, tmp1);
            vtkMath::Add(tmp1, vy, tmp2);
            vtkMath::Add(tmp2, vz, pp);
            newCutCircle->InsertNextPoint(pp);
        }
        cutCircle->SetPoints(newCutCircle);
        // build target landmarks

        appendFilter->RemoveAllInputs();
        appendFilter->AddInputData(cutCircle);
        appendFilter->AddInputData(lastNewCircle);
        appendFilter->Update();
        target->ShallowCopy(appendFilter->GetOutput());

        lastNewCircle->DeepCopy(cutCircle);

        appendFilter->RemoveAllInputs();
        appendFilter->AddInputData(CutCircleLineUp);
        appendFilter->AddInputData(cutCircle);
        appendFilter->Update();
        //cleanFilter->SetInputConnection(appendFilter->GetOutputPort());
        //cleanFilter->Update();
        CutCircleLineUp->DeepCopy(appendFilter->GetOutput());

        if(1)
        {
        piece = PieceBetweenPlanes(t_colon, PlaneOriginals, PlaneNormals, i-1, i, lastpiece);
        lastpiece->DeepCopy(piece);
        vtkSmartPointer<vtkThinPlateSplineTransform> transform = vtkSmartPointer<vtkThinPlateSplineTransform>::New();
        transform->SetSourceLandmarks(source->GetPoints());
        transform->SetTargetLandmarks(target->GetPoints());
        transform->SetBasisToR();
        transformFilter->RemoveAllInputs();
        transformFilter->SetInputData(piece);
        transformFilter->SetTransform(transform);
        transformFilter->Update();
        transformedpiece->ShallowCopy(transformFilter->GetOutput());

        appendFilter->RemoveAllInputs();
        appendFilter->AddInputData(SurfaceLineUp);
        appendFilter->AddInputData(transformedpiece);
        appendFilter->Update();
        cleanFilter->SetInputConnection(appendFilter->GetOutputPort());
        cleanFilter->Update();
        SurfaceLineUp->DeepCopy(cleanFilter->GetOutput());
        std::cout<<i<<" "<<SurfaceLineUp->GetNumberOfPoints()<<" "<<piece->GetNumberOfPoints()<<endl;
        }
    }

    // the very last piece
    if(1)
    {
    source->ShallowCopy(lastCircle);
    target->ShallowCopy(lastNewCircle);
    piece = PieceBetweenPlanes(t_colon, PlaneOriginals, PlaneNormals, model->GetNumberOfPoints()-1, model->GetNumberOfPoints(), lastpiece);
    lastpiece->DeepCopy(piece);
    vtkSmartPointer<vtkThinPlateSplineTransform> transform = vtkSmartPointer<vtkThinPlateSplineTransform>::New();
    transform->SetSourceLandmarks(source->GetPoints());
    transform->SetTargetLandmarks(target->GetPoints());
    transform->SetBasisToR();
    transformFilter->RemoveAllInputs();
    transformFilter->SetInputData(piece);
    transformFilter->SetTransform(transform);
    transformFilter->Update();
    transformedpiece->ShallowCopy(transformFilter->GetOutput());
    appendFilter->RemoveAllInputs();
    appendFilter->AddInputData(SurfaceLineUp);
    appendFilter->AddInputData(transformedpiece);
    appendFilter->Update();
    cleanFilter->SetInputConnection(appendFilter->GetOutputPort());
    cleanFilter->Update();
    SurfaceLineUp->DeepCopy(cleanFilter->GetOutput());
    std::cout<<model->GetNumberOfPoints()<<" "<<SurfaceLineUp->GetNumberOfPoints()<<" "<<piece->GetNumberOfPoints()<<endl;
    }
    //

    std::cout<<"cutcircles have points"<<CutCircleOrigin->GetNumberOfPoints()<<endl;
    std::cout<<"cutcircles have points"<<CutCircleLineUp->GetNumberOfPoints()<<endl;

    vtkSmartPointer<vtkPolyDataMapper> CutCircleLineUpMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    CutCircleLineUpMapper->SetInputData(CutCircleLineUp);
    CutCircleLineUpMapper->Update();
    vtkSmartPointer<vtkActor> CutCircleLineUpActor = vtkSmartPointer<vtkActor>::New();
    CutCircleLineUpActor->SetMapper(CutCircleLineUpMapper);
    CutCircleLineUpActor->GetProperty()->SetColor(1, 0, 0);
    //t_rendermanager->renderModel(CutCircleLineUpActor);

    vtkSmartPointer<vtkPolyDataMapper> SurfaceLineUpMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    SurfaceLineUpMapper->SetInputData(SurfaceLineUp);
    SurfaceLineUpMapper->Update();
    vtkSmartPointer<vtkActor> SurfaceLineUpActor = vtkSmartPointer<vtkActor>::New();
    SurfaceLineUpActor->SetMapper(SurfaceLineUpMapper);
    t_rendermanager->renderModel(SurfaceLineUpActor);

    t_filemanager->SaveFile(SurfaceLineUp, "SurfaceLineUp_v31.stl");

    std::cout<<"Deformation End"<<endl;
    return SurfaceLineUp;
}

vtkSmartPointer<vtkPolyData> Centerline::PieceBetweenPlanes(vtkSmartPointer<vtkPolyData> t_colon,
                                                            vtkSmartPointer<vtkDoubleArray> PlaneOriginals, vtkSmartPointer<vtkDoubleArray> PlaneNormals,
                                                            vtkIdType left, vtkIdType right, vtkSmartPointer<vtkPolyData> lastpiece)
{
    std::cout<<left<<" to "<<right<<endl;
    assert(left == right - 1);
    int threshold = 300;
    vtkSmartPointer<vtkClipPolyData> clipper = vtkSmartPointer<vtkClipPolyData>::New();
    clipper->SetInputData(t_colon);
    vtkSmartPointer<vtkPolyDataConnectivityFilter> connectivityFilter = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
    connectivityFilter->SetInputConnection(clipper->GetOutputPort());
    connectivityFilter->SetExtractionModeToClosestPointRegion();
    vtkSmartPointer<vtkPolyData> output = vtkSmartPointer<vtkPolyData>::New();
    if(left < 0) // the section before the first plane
    {
        vtkSmartPointer<vtkPlane> planeright = vtkSmartPointer<vtkPlane>::New();
        double originright[3], normalright[3];
        PlaneOriginals->GetTuple(right, originright);
        PlaneNormals->GetTuple(right, normalright);
        vtkMath::MultiplyScalar(normalright, -1);
        planeright->SetOrigin(originright);
        planeright->SetNormal(normalright);
        clipper->SetClipFunction(planeright);
        clipper->Update();
        double centerpoint[3];
        model->GetPoint(right, centerpoint);
        connectivityFilter->SetClosestPoint(centerpoint);
        connectivityFilter->Update();
        output->DeepCopy(connectivityFilter->GetOutput());
    }
    else if(right >= model->GetNumberOfPoints())
    {
        vtkSmartPointer<vtkPlane> planeleft = vtkSmartPointer<vtkPlane>::New();
        double originleft[3], normalleft[3];
        PlaneOriginals->GetTuple(left, originleft);
        PlaneNormals->GetTuple(left, normalleft);
        planeleft->SetOrigin(originleft);
        planeleft->SetNormal(normalleft);
        clipper->SetClipFunction(planeleft);
        clipper->Update();
        double centerpoint[3];
        model->GetPoint(left, centerpoint);
        connectivityFilter->SetClosestPoint(centerpoint);
        connectivityFilter->Update();
        output->DeepCopy(connectivityFilter->GetOutput());

        if(output->GetNumberOfPoints() < threshold)
        {
            double p[3];
            for(vtkIdType j = 0; j < lastpiece->GetNumberOfPoints(); j++)
            {
                lastpiece->GetPoint(j, p);
                connectivityFilter->SetClosestPoint(p);
                connectivityFilter->Update();
                output->DeepCopy(connectivityFilter->GetOutput());
                if(output->GetNumberOfPoints() >= threshold)
                    break;
            }
        }
    }
    else
    {
        vtkSmartPointer<vtkPlane> planeright = vtkSmartPointer<vtkPlane>::New();
        double originright[3], normalright[3];
        PlaneOriginals->GetTuple(right, originright);
        PlaneNormals->GetTuple(right, normalright);
        vtkMath::MultiplyScalar(normalright, -1);
        planeright->SetOrigin(originright);
        planeright->SetNormal(normalright);
        clipper->SetClipFunction(planeright);
        clipper->Update();

        vtkSmartPointer<vtkPolyData> temp = vtkSmartPointer<vtkPolyData>::New();
        temp->DeepCopy(clipper->GetOutput());
        clipper->RemoveAllInputs();
        clipper->SetInputData(temp);

        vtkSmartPointer<vtkPlane> planeleft = vtkSmartPointer<vtkPlane>::New();
        double originleft[3], normalleft[3];
        PlaneOriginals->GetTuple(left, originleft);
        PlaneNormals->GetTuple(left, normalleft);
        planeleft->SetOrigin(originleft);
        planeleft->SetNormal(normalleft);
        clipper->SetClipFunction(planeleft);
        clipper->Update();

        double centerpointleft[3], centerpointright[3];
        model->GetPoint(left, centerpointleft);
        model->GetPoint(right, centerpointright);
        double centerpoint[3];
        vtkMath::Add(centerpointleft, centerpointright, centerpoint);
        vtkMath::MultiplyScalar(centerpoint, 0.5);
        connectivityFilter->SetClosestPoint(centerpoint);
        connectivityFilter->Update();

        output->DeepCopy(connectivityFilter->GetOutput());

        if(output->GetNumberOfPoints() < threshold)
        {
            double p[3];
            for(vtkIdType j = 0; j < lastpiece->GetNumberOfPoints(); j++)
            {
                lastpiece->GetPoint(j, p);
                connectivityFilter->SetClosestPoint(p);
                connectivityFilter->Update();
                output->DeepCopy(connectivityFilter->GetOutput());
                if(output->GetNumberOfPoints() >= threshold)
                    break;
            }
        }
    }
    return output;
}
vtkIdType Centerline::GetTheSectionIdOfAPoint(vtkSmartPointer<vtkPolyData> t_colon, vtkIdType pointid,
                                              vtkSmartPointer<vtkDoubleArray> PlaneOriginals, vtkSmartPointer<vtkDoubleArray> PlaneNormals, vtkSmartPointer<vtkIdList> CurvaturePointIds)
{
    std::cout<<"Get the section id of "<<pointid<<endl;
    double p[3], vl[3], vr[3], cpl[3], cpr[3], cp[3];
    t_colon->GetPoint(pointid, p);
    double originright[3], normalright[3], originleft[3], normalleft[3];
    double minDis2 = INFINITY, dis2;
    vtkIdType output;
    vtkIdType left, right;

    for(vtkIdType i = 0; i <= model->GetNumberOfPoints(); i++)
    {
        left = i-1; right = i;
        if(left < 0)
        {
            PlaneOriginals->GetTuple(right, originright);
            PlaneNormals->GetTuple(right, normalright);
            vtkMath::Subtract(p, originright, vr);
            //std::cout<<i<<" "<<vtkMath::Dot(vr, normalright)<<endl;
            if(vtkMath::Dot(vr, normalright) < 0)
            {
                model->GetPoint(right, cpr);
                dis2 = vtkMath::Distance2BetweenPoints(p, cpr);
                if(dis2 < minDis2)
                {
                    minDis2 = dis2;
                    output = i;
                }
            }
        }
        else if(right >= model->GetNumberOfPoints())
        {
            PlaneOriginals->GetTuple(left, originleft);
            PlaneNormals->GetTuple(left, normalleft);
            vtkMath::Subtract(p, originleft, vl);
            //std::cout<<i<<" "<<vtkMath::Dot(vl, normalleft)<<endl;
            if(vtkMath::Dot(vl, normalleft) >= 0)
            {
                model->GetPoint(left, cpl);
                dis2 = vtkMath::Distance2BetweenPoints(p, cpl);
                if(dis2 < minDis2)
                {
                    minDis2 = dis2;
                    output = i;
                }
            }
        }
        else
        {
            PlaneOriginals->GetTuple(right, originright);
            PlaneNormals->GetTuple(right, normalright);
            PlaneOriginals->GetTuple(left, originleft);
            PlaneNormals->GetTuple(left, normalleft);
            vtkMath::Subtract(p, originright, vr);
            vtkMath::Subtract(p, originleft, vl);

            //std::cout<<i<<" "<<vtkMath::Dot(vl, normalleft)<<" "<<vtkMath::Dot(vr, normalright)<<endl;
            if(vtkMath::Dot(vl, normalleft) >= 0 && vtkMath::Dot(vr, normalright) < 0)
            {
                model->GetPoint(left, cpl);
                model->GetPoint(right, cpr);
                vtkMath::Add(cpl, cpr, cp);
                vtkMath::MultiplyScalar(cp, 0.5);
                dis2 = vtkMath::Distance2BetweenPoints(p, cp);
                if(dis2 < minDis2)
                {
                    minDis2 = dis2;
                    output = i;
                }
            }
        }
    }
    return output;
}
vtkIdType Centerline::GetTheSectionIdOfAPoint_v2(vtkSmartPointer<vtkPolyData> t_colon, vtkIdType pointid,
                                              vtkSmartPointer<vtkDoubleArray> PlaneOriginals, vtkSmartPointer<vtkDoubleArray> PlaneNormals, vtkSmartPointer<vtkIdList> CurvaturePointIds)
{
    std::cout<<"Get the section id of "<<pointid<<endl;
    double p[3], vl[3], vr[3], cpl[3], cpr[3], cp[3];
    //double curvaturePointl[3], curvaturePointr[3];

    t_colon->GetPoint(pointid, p);
    double originright[3], normalright[3], originleft[3], normalleft[3];
    double dis2, dis2Threshold = 10000;
    vtkIdType output;
    vtkIdType left, right;
    double len, minLen = INFINITY;

    vtkSmartPointer<vtkDijkstraGraphGeodesicPath> geodesicPath = vtkSmartPointer<vtkDijkstraGraphGeodesicPath>::New();
    geodesicPath->SetInputData(t_colon);
    geodesicPath->SetStartVertex(pointid);

    for(vtkIdType i = 0; i <= model->GetNumberOfPoints(); i++)
    {
        left = i-1; right = i;
        if(left < 0)
        {
            PlaneOriginals->GetTuple(right, originright);
            PlaneNormals->GetTuple(right, normalright);
            vtkMath::Subtract(p, originright, vr);
            //std::cout<<i<<" "<<vtkMath::Dot(vr, normalright)<<endl;
            if(vtkMath::Dot(vr, normalright) < 0)
            {
                model->GetPoint(right, cpr);
                dis2 = vtkMath::Distance2BetweenPoints(p, cpr);
                if(dis2 < dis2Threshold)
                {

                    geodesicPath->SetEndVertex(CurvaturePointIds->GetId(right));
                    geodesicPath->Update();
                    len = geodesicPath->GetOutput()->GetNumberOfPoints();
                    if(len < minLen)
                    {
                        minLen = len;
                        output = i;
                    }
                }
            }
        }
        else if(right >= model->GetNumberOfPoints())
        {
            PlaneOriginals->GetTuple(left, originleft);
            PlaneNormals->GetTuple(left, normalleft);
            vtkMath::Subtract(p, originleft, vl);
            //std::cout<<i<<" "<<vtkMath::Dot(vl, normalleft)<<endl;
            if(vtkMath::Dot(vl, normalleft) >= 0)
            {
                model->GetPoint(left, cpl);
                dis2 = vtkMath::Distance2BetweenPoints(p, cpl);
                if(dis2 < dis2Threshold)
                {
                    geodesicPath->SetEndVertex(CurvaturePointIds->GetId(left));
                    geodesicPath->Update();
                    len = geodesicPath->GetOutput()->GetNumberOfPoints();
                    if(len < minLen)
                    {
                        minLen = len;
                        output = i;
                    }
                }
            }
        }
        else
        {
            PlaneOriginals->GetTuple(right, originright);
            PlaneNormals->GetTuple(right, normalright);
            PlaneOriginals->GetTuple(left, originleft);
            PlaneNormals->GetTuple(left, normalleft);
            vtkMath::Subtract(p, originright, vr);
            vtkMath::Subtract(p, originleft, vl);

            //std::cout<<i<<" "<<vtkMath::Dot(vl, normalleft)<<" "<<vtkMath::Dot(vr, normalright)<<endl;
            if(vtkMath::Dot(vl, normalleft) >= 0 && vtkMath::Dot(vr, normalright) < 0)
            {
                model->GetPoint(left, cpl);
                model->GetPoint(right, cpr);
                vtkMath::Add(cpl, cpr, cp);
                vtkMath::MultiplyScalar(cp, 0.5);
                dis2 = vtkMath::Distance2BetweenPoints(p, cp);
                if(dis2 < dis2Threshold)
                {
                    geodesicPath->SetEndVertex(CurvaturePointIds->GetId(left));
                    geodesicPath->Update();
                    len = geodesicPath->GetOutput()->GetNumberOfPoints();
                    if(len < minLen)
                    {
                        minLen = len;
                        output = i;
                    }
                }
            }
        }
    }
    return output;
}
void Centerline::GetSectionIds(vtkPolyData *t_colon, vtkIdType pointid, vtkIdList *SectionIds,
                               vtkDoubleArray *PlaneOriginals, vtkDoubleArray *PlaneNormals)
{


    vtkSmartPointer<vtkIdList> connectedVertices = vtkSmartPointer<vtkIdList>::New();
    connectedVertices = GetConnectedVertices(t_colon, pointid);
    vtkIdType currentSectionId = SectionIds->GetId(pointid);

    std::cout<<"seed : "<<pointid<<"(current section "<<currentSectionId<<") "<<connectedVertices->GetNumberOfIds()<<": ";
    for(vtkIdType i = 0; i < connectedVertices->GetNumberOfIds(); i++)
    {
        std::cout<<connectedVertices->GetId(i)<<" ";
    }
    std::cout<<endl;

    for(vtkIdType i = 0; i < connectedVertices->GetNumberOfIds(); i++)
    {
        double p[3];
        double originright[3], normalright[3], vr[3];
        double originleft[3], normalleft[3], vl[3];
        double dotl, dotr;
        vtkIdType id = connectedVertices->GetId(i);
        t_colon->GetPoint(id, p);
        if(SectionIds->GetId(id) >= 0) // this point has already been assigned
            continue; // do nothing for this point
        vtkIdType left = currentSectionId-1, right = currentSectionId;
        if(left < 0) // id is at the left end
        {
            PlaneOriginals->GetTuple(right, originright);
            PlaneNormals->GetTuple(right, normalright);
            vtkMath::Subtract(p, originright, vr);
            if(vtkMath::Dot(vr, normalright) <= 0) // id satisfy the section criteron
            {
                SectionIds->SetId(id, right);
            }
            // else go right
            else
            {
                while(1){
                    left++; right ++;
                    if(right > model->GetNumberOfPoints())
                    {
                        std::cerr<<"failed to find section id (go right to the end)"<<endl;
                        exit(1);
                    }
                    PlaneOriginals->GetTuple(left, originleft);
                    PlaneNormals->GetTuple(left, normalleft);
                    vtkMath::Subtract(p, originleft, vl);
                    dotl = vtkMath::Dot(vl, normalleft);
                    if(right == model->GetNumberOfPoints())
                    {
                        if(dotl > 0)
                        {
                            SectionIds->SetId(id, right);
                            break;
                        }
                    }
                    else
                    {
                        PlaneOriginals->GetTuple(right, originright);
                        PlaneNormals->GetTuple(right, normalright);
                        vtkMath::Subtract(p, originright, vr);
                        dotr = vtkMath::Dot(vr, normalright);
                        {
                            if(dotl > 0 && dotr <=0)
                            {
                                SectionIds->SetId(id, right);
                                break;
                            }
                        }
                    }
                } // end of go right
            }
        }
        else if(right >= model->GetNumberOfPoints()) // id is at the right end
        {
            PlaneOriginals->GetTuple(left, originleft);
            PlaneNormals->GetTuple(left, normalleft);
            vtkMath::Subtract(p, originleft, vl);
            if(vtkMath::Dot(vl, normalleft) > 0) // id satisfy the section criteron
            {
                SectionIds->SetId(id, right);
            }
            // else go left
            else
            {
                while(1){
                    left--; right--;
                    if(left < -1)
                    {
                        std::cerr<<"failed to find section id (go left to the end)"<<endl;
                        exit(1);
                    }
                    PlaneOriginals->GetTuple(right, originright);
                    PlaneNormals->GetTuple(right, normalright);
                    vtkMath::Subtract(p, originright, vr);
                    dotr = vtkMath::Dot(vr, normalright);
                    if(left == -1)
                    {
                        if(dotr <= 0)
                        {
                            SectionIds->SetId(id, right);
                            break;
                        }
                    }
                    else
                    {
                        PlaneOriginals->GetTuple(left, originleft);
                        PlaneNormals->GetTuple(left, normalleft);
                        vtkMath::Subtract(p, originleft, vl);
                        dotl = vtkMath::Dot(vl, normalleft);
                        {
                            if(dotl > 0 && dotr <=0)
                            {
                                SectionIds->SetId(id, right);
                                break;
                            }
                        }
                    }
                } // end of go left
            }
        }
        else // id is in the middle
        {
            PlaneOriginals->GetTuple(right, originright);
            PlaneNormals->GetTuple(right, normalright);
            vtkMath::Subtract(p, originright, vr);
            PlaneOriginals->GetTuple(left, originleft);
            PlaneNormals->GetTuple(left, normalleft);
            vtkMath::Subtract(p, originleft, vl);
            dotr = vtkMath::Dot(vr, normalright);
            dotl = vtkMath::Dot(vl, normalleft);
            //std::cout<<"dotr: "<<dotr<<" dotl: "<<dotl<<std::endl;
            if(dotr <= 0 && dotl > 0) // satisfy
            {
                SectionIds->SetId(id, right);
            }
            else if(dotr <= 0 && dotl <= 0) // should go left(the point is on the left side of the current section)
            {
                while(1){
                    left--; right--;
                    if(left < -1)
                    {
                        std::cerr<<"failed to find section id (go left to the end)"<<endl;
                        exit(1);
                    }
                    PlaneOriginals->GetTuple(right, originright);
                    PlaneNormals->GetTuple(right, normalright);
                    vtkMath::Subtract(p, originright, vr);
                    dotr = vtkMath::Dot(vr, normalright);
                    if(left == -1)
                    {
                        if(dotr <= 0)
                        {
                            SectionIds->SetId(id, right);
                            break;
                        }
                    }
                    else
                    {
                        PlaneOriginals->GetTuple(left, originleft);
                        PlaneNormals->GetTuple(left, normalleft);
                        vtkMath::Subtract(p, originleft, vl);
                        dotl = vtkMath::Dot(vl, normalleft);
                        {
                            if(dotl > 0 && dotr <=0)
                            {
                                SectionIds->SetId(id, right);
                                break;
                            }
                        }
                    }
                } // end of go left
            }
            else if(dotr >0 && dotl > 0) // should go right(the point is on the right side of the current section)
            {

                while(1){
                    left++; right ++;
                    if(right > model->GetNumberOfPoints())
                    {
                        std::cerr<<"failed to find section id (go right to the end)"<<endl;
                        exit(1);
                    }
                    PlaneOriginals->GetTuple(left, originleft);
                    PlaneNormals->GetTuple(left, normalleft);
                    vtkMath::Subtract(p, originleft, vl);
                    dotl = vtkMath::Dot(vl, normalleft);
                    if(right == model->GetNumberOfPoints())
                    {
                        if(dotl > 0)
                        {
                            SectionIds->SetId(id, right);
                            break;
                        }
                    }
                    else
                    {
                        PlaneOriginals->GetTuple(right, originright);
                        PlaneNormals->GetTuple(right, normalright);
                        vtkMath::Subtract(p, originright, vr);
                        dotr = vtkMath::Dot(vr, normalright);
                        {
                            if(dotl > 0 && dotr <=0)
                            {
                                SectionIds->SetId(id, right);
                                break;
                            }
                        }
                    }
                } // end of go right
            }
            else
            {
                std::cerr<<"failed to find section id (unreasonalbe dotl and dotr)"<<endl;
                exit(0);
            }
        }
        // after processing the id-th point, should process its neighbours
        std::cout<<"set "<<id<<" to "<<right<<endl;
        GetSectionIds(t_colon, id, SectionIds, PlaneOriginals, PlaneNormals);
    }
}

void Centerline::GetSectionIds_loop(vtkPolyData *t_colon, vtkIdType seed, vtkIdList *SectionIds,
                               vtkDoubleArray *PlaneOriginals, vtkDoubleArray *PlaneNormals)
{
    int count = 0;
    int difference = 0, tolerance = INFINITY;
    vtkSmartPointer<vtkIdList> currentlevel = vtkSmartPointer<vtkIdList>::New();
    // init
    currentlevel->InsertNextId(seed);

    while(currentlevel->GetNumberOfIds() > 0)
    {
        count += currentlevel->GetNumberOfIds();
        //std::cout<<"processed "<<count<<" points"<<endl;
        vtkSmartPointer<vtkIdList> nextlevel = vtkSmartPointer<vtkIdList>::New();
        for(vtkIdType n = 0; n < currentlevel->GetNumberOfIds(); n++)
        {
            vtkIdType pointid = currentlevel->GetId(n);
            vtkSmartPointer<vtkIdList> connectedVertices = vtkSmartPointer<vtkIdList>::New();
            connectedVertices = GetConnectedVertices(t_colon, pointid);
            vtkIdType currentSectionId = SectionIds->GetId(pointid);

            /*
            std::cout<<"pointid : "<<pointid<<"(current section "<<currentSectionId<<") "<<connectedVertices->GetNumberOfIds()<<": ";
            for(vtkIdType i = 0; i < connectedVertices->GetNumberOfIds(); i++)
            {
                std::cout<<connectedVertices->GetId(i)<<" ";
            }
            std::cout<<endl;
            */

            for(vtkIdType i = 0; i < connectedVertices->GetNumberOfIds(); i++)
            {
                double p[3];
                double originright[3], normalright[3], vr[3];
                double originleft[3], normalleft[3], vl[3];
                double dotl, dotr;
                vtkIdType id = connectedVertices->GetId(i);
                difference = 0;
                t_colon->GetPoint(id, p);
                if(SectionIds->GetId(id) >= 0) // this point has already been assigned
                    continue; // do nothing for this point
                else
                    nextlevel->InsertNextId(id);
                vtkIdType left = currentSectionId-1, right = currentSectionId;
                if(left < 0) // id is at the left end
                {
                    PlaneOriginals->GetTuple(right, originright);
                    PlaneNormals->GetTuple(right, normalright);
                    vtkMath::Subtract(p, originright, vr);
                    if(vtkMath::Dot(vr, normalright) <= 0) // id satisfy the section criteron
                    {
                        SectionIds->SetId(id, right);
                    }
                    // else go right
                    else
                    {
                        while(1){
                            left++; right ++;

                            // check whether the section difference between neighboring points have reached the tolerance
                            if(difference++ >= tolerance)
                            {
                                SectionIds->SetId(id, right);
                                break;
                            }

                            if(right > model->GetNumberOfPoints())
                            {
                                std::cerr<<"failed to find section id (go right to the end)"<<endl;
                                exit(1);
                            }
                            PlaneOriginals->GetTuple(left, originleft);
                            PlaneNormals->GetTuple(left, normalleft);
                            vtkMath::Subtract(p, originleft, vl);
                            dotl = vtkMath::Dot(vl, normalleft);
                            if(right == model->GetNumberOfPoints())
                            {
                                if(dotl > 0)
                                {
                                    SectionIds->SetId(id, right);
                                    break;
                                }
                            }
                            else
                            {
                                PlaneOriginals->GetTuple(right, originright);
                                PlaneNormals->GetTuple(right, normalright);
                                vtkMath::Subtract(p, originright, vr);
                                dotr = vtkMath::Dot(vr, normalright);
                                {
                                    if(dotl > 0 && dotr <=0)
                                    {
                                        SectionIds->SetId(id, right);
                                        break;
                                    }
                                }
                            }
                        } // end of go right
                    }
                }
                else if(right >= model->GetNumberOfPoints()) // id is at the right end
                {
                    PlaneOriginals->GetTuple(left, originleft);
                    PlaneNormals->GetTuple(left, normalleft);
                    vtkMath::Subtract(p, originleft, vl);
                    if(vtkMath::Dot(vl, normalleft) > 0) // id satisfy the section criteron
                    {
                        SectionIds->SetId(id, right);
                    }
                    // else go left
                    else
                    {
                        while(1){
                            left--; right--;

                            // check whether the section difference between neighboring points have reached the tolerance
                            if(difference++ >= tolerance)
                            {
                                SectionIds->SetId(id, right);
                                break;
                            }

                            if(left < -1)
                            {
                                std::cerr<<"failed to find section id (go left to the end)"<<endl;
                                exit(1);
                            }
                            PlaneOriginals->GetTuple(right, originright);
                            PlaneNormals->GetTuple(right, normalright);
                            vtkMath::Subtract(p, originright, vr);
                            dotr = vtkMath::Dot(vr, normalright);
                            if(left == -1)
                            {
                                if(dotr <= 0)
                                {
                                    SectionIds->SetId(id, right);
                                    break;
                                }
                            }
                            else
                            {
                                PlaneOriginals->GetTuple(left, originleft);
                                PlaneNormals->GetTuple(left, normalleft);
                                vtkMath::Subtract(p, originleft, vl);
                                dotl = vtkMath::Dot(vl, normalleft);
                                {
                                    if(dotl > 0 && dotr <=0)
                                    {
                                        SectionIds->SetId(id, right);
                                        break;
                                    }
                                }
                            }
                        } // end of go left
                    }
                }
                else // id is in the middle
                {
                    PlaneOriginals->GetTuple(right, originright);
                    PlaneNormals->GetTuple(right, normalright);
                    vtkMath::Subtract(p, originright, vr);
                    PlaneOriginals->GetTuple(left, originleft);
                    PlaneNormals->GetTuple(left, normalleft);
                    vtkMath::Subtract(p, originleft, vl);
                    dotr = vtkMath::Dot(vr, normalright);
                    dotl = vtkMath::Dot(vl, normalleft);
                    //std::cout<<"dotr: "<<dotr<<" dotl: "<<dotl<<std::endl;
                    if(dotr <= 0 && dotl > 0) // satisfy
                    {
                        SectionIds->SetId(id, right);
                    }
                    else if(dotr <= 0 && dotl <= 0) // should go left(the point is on the left side of the current section)
                    {
                        while(1){
                            left--; right--;

                            // check whether the section difference between neighboring points have reached the tolerance
                            if(difference++ >= tolerance)
                            {
                                SectionIds->SetId(id, right);
                                break;
                            }

                            if(left < -1)
                            {
                                std::cerr<<"failed to find section id (go left to the end)"<<endl;
                                exit(1);
                            }
                            PlaneOriginals->GetTuple(right, originright);
                            PlaneNormals->GetTuple(right, normalright);
                            vtkMath::Subtract(p, originright, vr);
                            dotr = vtkMath::Dot(vr, normalright);
                            if(left == -1)
                            {
                                if(dotr <= 0)
                                {
                                    SectionIds->SetId(id, right);
                                    break;
                                }
                            }
                            else
                            {
                                PlaneOriginals->GetTuple(left, originleft);
                                PlaneNormals->GetTuple(left, normalleft);
                                vtkMath::Subtract(p, originleft, vl);
                                dotl = vtkMath::Dot(vl, normalleft);
                                {
                                    if(dotl > 0 && dotr <=0)
                                    {
                                        SectionIds->SetId(id, right);
                                        break;
                                    }
                                }
                            }
                        } // end of go left
                    }
                    else if(dotr >0 && dotl > 0) // should go right(the point is on the right side of the current section)
                    {

                        while(1){
                            left++; right ++;

                            // check whether the section difference between neighboring points have reached the tolerance
                            if(difference++ >= tolerance)
                            {
                                SectionIds->SetId(id, right);
                                break;
                            }

                            if(right > model->GetNumberOfPoints())
                            {
                                std::cerr<<"failed to find section id (go right to the end)"<<endl;
                                exit(1);
                            }
                            PlaneOriginals->GetTuple(left, originleft);
                            PlaneNormals->GetTuple(left, normalleft);
                            vtkMath::Subtract(p, originleft, vl);
                            dotl = vtkMath::Dot(vl, normalleft);
                            if(right == model->GetNumberOfPoints())
                            {
                                if(dotl > 0)
                                {
                                    SectionIds->SetId(id, right);
                                    break;
                                }
                            }
                            else
                            {
                                PlaneOriginals->GetTuple(right, originright);
                                PlaneNormals->GetTuple(right, normalright);
                                vtkMath::Subtract(p, originright, vr);
                                dotr = vtkMath::Dot(vr, normalright);
                                {
                                    if(dotl > 0 && dotr <=0)
                                    {
                                        SectionIds->SetId(id, right);
                                        break;
                                    }
                                }
                            }
                        } // end of go right
                    }
                    else
                    {
                        std::cerr<<"failed to find section id (unreasonalbe dotl and dotr)"<<endl;
                        exit(0);
                    }
                }
                // after processing the id-th point, should process its neighbours
                //std::cout<<"set "<<id<<" to "<<right<<endl;
            }
        } // for(vtkIdType n = 0; n < nextlevel->GetNumberOfIds(); n++)
        currentlevel->Reset();
        currentlevel->DeepCopy(nextlevel);
    } // while(1)
}
void Centerline::GetSectionIds_loop_v2(vtkPolyData *t_colon, vtkIdType seed, vtkIdList *SectionIds,
                               vtkDoubleArray *PlaneOriginals, vtkDoubleArray *PlaneNormals, vtkPolyDataGroup* CircleGroup)
{
    int count = 0;
    int difference = 0, tolerance = INFINITY;
    double dis2Thres = INFINITY;
    vtkSmartPointer<vtkIdList> currentlevel = vtkSmartPointer<vtkIdList>::New();
    // init
    currentlevel->InsertNextId(seed);

    while(currentlevel->GetNumberOfIds() > 0)
    {
        count += currentlevel->GetNumberOfIds();
        std::cout<<"processed "<<count<<" points"<<endl;
        vtkSmartPointer<vtkIdList> nextlevel = vtkSmartPointer<vtkIdList>::New();
        for(vtkIdType n = 0; n < currentlevel->GetNumberOfIds(); n++)
        {
            vtkIdType pointid = currentlevel->GetId(n);
            vtkSmartPointer<vtkIdList> connectedVertices = vtkSmartPointer<vtkIdList>::New();
            connectedVertices = GetConnectedVertices(t_colon, pointid);
            vtkIdType currentSectionId = SectionIds->GetId(pointid);


            std::cout<<"pointid : "<<pointid<<"(current section "<<currentSectionId<<") "<<connectedVertices->GetNumberOfIds()<<": ";
            for(vtkIdType i = 0; i < connectedVertices->GetNumberOfIds(); i++)
            {
                std::cout<<connectedVertices->GetId(i)<<" ";
            }
            std::cout<<endl;


            for(vtkIdType i = 0; i < connectedVertices->GetNumberOfIds(); i++)
            {
                double p[3];
                double originright[3], normalright[3], vr[3];
                double originleft[3], normalleft[3], vl[3];
                double dotl, dotr;
                vtkIdType id = connectedVertices->GetId(i);
                difference = 0;
                t_colon->GetPoint(id, p);
                if(SectionIds->GetId(id) >= 0) // this point has already been assigned
                    continue; // do nothing for this point
                else
                    nextlevel->InsertNextId(id);
                vtkIdType left = currentSectionId-1, right = currentSectionId;

                // refer to the current circle
                vtkSmartPointer<vtkAppendPolyData> appendFilter = vtkSmartPointer<vtkAppendPolyData>::New();
                vtkSmartPointer<vtkPolyData> currentCircle = vtkSmartPointer<vtkPolyData>::New();
                if(currentSectionId == 0)
                    currentCircle->DeepCopy(CircleGroup->GetMember(0));
                else if(currentSectionId == model->GetNumberOfPoints())
                    currentCircle->DeepCopy(CircleGroup->GetMember(model->GetNumberOfPoints()-1));
                else
                {
                    appendFilter->RemoveAllInputs();
                    appendFilter->AddInputData(CircleGroup->GetMember(left));
                    appendFilter->AddInputData(CircleGroup->GetMember(right));
                    appendFilter->Update();
                    currentCircle->DeepCopy(appendFilter->GetOutput());
                }
                vtkSmartPointer<vtkPointLocator> pointLocator = vtkSmartPointer<vtkPointLocator>::New();
                pointLocator->SetDataSet(currentCircle);
                vtkIdType closestid = pointLocator->FindClosestPoint(p);
                double closestp[3];
                currentCircle->GetPoint(closestid, closestp);
                double dis2 = vtkMath::Distance2BetweenPoints(closestp, p);

                if(left < 0) // id is at the left end
                {
                    PlaneOriginals->GetTuple(right, originright);
                    PlaneNormals->GetTuple(right, normalright);
                    vtkMath::Subtract(p, originright, vr);
                    if(vtkMath::Dot(vr, normalright) <= 0) // id satisfy the section criteron
                    {
                        SectionIds->SetId(id, right);
                    }
                    else if(dis2 > dis2Thres) // try
                    {
                        SectionIds->SetId(id, right);
                    }
                    // else go right
                    else
                    {
                        while(1){
                            left++; right ++;

                            // check whether the section difference between neighboring points have reached the tolerance
                            if(difference++ >= tolerance)
                            {
                                SectionIds->SetId(id, right);
                                break;
                            }

                            if(right > model->GetNumberOfPoints())
                            {
                                std::cerr<<"failed to find section id (go right to the end)"<<endl;
                                exit(1);
                            }
                            PlaneOriginals->GetTuple(left, originleft);
                            PlaneNormals->GetTuple(left, normalleft);
                            vtkMath::Subtract(p, originleft, vl);
                            dotl = vtkMath::Dot(vl, normalleft);
                            if(right == model->GetNumberOfPoints())
                            {
                                if(dotl > 0)
                                {
                                    SectionIds->SetId(id, right);
                                    break;
                                }
                            }
                            else
                            {
                                PlaneOriginals->GetTuple(right, originright);
                                PlaneNormals->GetTuple(right, normalright);
                                vtkMath::Subtract(p, originright, vr);
                                dotr = vtkMath::Dot(vr, normalright);
                                {
                                    if(dotl > 0 && dotr <=0)
                                    {
                                        SectionIds->SetId(id, right);
                                        break;
                                    }
                                }
                            }
                        } // end of go right
                    }
                }
                else if(right >= model->GetNumberOfPoints()) // id is at the right end
                {
                    PlaneOriginals->GetTuple(left, originleft);
                    PlaneNormals->GetTuple(left, normalleft);
                    vtkMath::Subtract(p, originleft, vl);
                    if(vtkMath::Dot(vl, normalleft) > 0) // id satisfy the section criteron
                    {
                        SectionIds->SetId(id, right);
                    }
                    else if(dis2 > dis2Thres) // try
                    {
                        SectionIds->SetId(id, right);
                    }
                    // else go left
                    else
                    {
                        while(1){
                            left--; right--;

                            // check whether the section difference between neighboring points have reached the tolerance
                            if(difference++ >= tolerance)
                            {
                                SectionIds->SetId(id, right);
                                break;
                            }

                            if(left < -1)
                            {
                                std::cerr<<"failed to find section id (go left to the end)"<<endl;
                                exit(1);
                            }
                            PlaneOriginals->GetTuple(right, originright);
                            PlaneNormals->GetTuple(right, normalright);
                            vtkMath::Subtract(p, originright, vr);
                            dotr = vtkMath::Dot(vr, normalright);
                            if(left == -1)
                            {
                                if(dotr <= 0)
                                {
                                    SectionIds->SetId(id, right);
                                    break;
                                }
                            }
                            else
                            {
                                PlaneOriginals->GetTuple(left, originleft);
                                PlaneNormals->GetTuple(left, normalleft);
                                vtkMath::Subtract(p, originleft, vl);
                                dotl = vtkMath::Dot(vl, normalleft);
                                {
                                    if(dotl > 0 && dotr <=0)
                                    {
                                        SectionIds->SetId(id, right);
                                        break;
                                    }
                                }
                            }
                        } // end of go left
                    }
                }
                else // id is in the middle
                {
                    PlaneOriginals->GetTuple(right, originright);
                    PlaneNormals->GetTuple(right, normalright);
                    vtkMath::Subtract(p, originright, vr);
                    PlaneOriginals->GetTuple(left, originleft);
                    PlaneNormals->GetTuple(left, normalleft);
                    vtkMath::Subtract(p, originleft, vl);
                    dotr = vtkMath::Dot(vr, normalright);
                    dotl = vtkMath::Dot(vl, normalleft);
                    //std::cout<<"dotr: "<<dotr<<" dotl: "<<dotl<<std::endl;
                    if(dotr <= 0 && dotl > 0) // satisfy
                    {
                        SectionIds->SetId(id, right);
                    }
                    else if(dis2 > dis2Thres) // try
                    {
                        SectionIds->SetId(id, right);
                    }
                    else if(dotr <= 0 && dotl <= 0) // should go left(the point is on the left side of the current section)
                    {
                        while(1){
                            left--; right--;

                            // check whether the section difference between neighboring points have reached the tolerance
                            if(difference++ >= tolerance)
                            {
                                SectionIds->SetId(id, right);
                                break;
                            }

                            if(left < -1)
                            {
                                std::cerr<<"failed to find section id (go left to the end)"<<endl;
                                exit(1);
                            }
                            PlaneOriginals->GetTuple(right, originright);
                            PlaneNormals->GetTuple(right, normalright);
                            vtkMath::Subtract(p, originright, vr);
                            dotr = vtkMath::Dot(vr, normalright);
                            if(left == -1)
                            {
                                if(dotr <= 0)
                                {
                                    SectionIds->SetId(id, right);
                                    break;
                                }
                            }
                            else
                            {
                                PlaneOriginals->GetTuple(left, originleft);
                                PlaneNormals->GetTuple(left, normalleft);
                                vtkMath::Subtract(p, originleft, vl);
                                dotl = vtkMath::Dot(vl, normalleft);
                                {
                                    if(dotl > 0 && dotr <=0)
                                    {
                                        SectionIds->SetId(id, right);
                                        break;
                                    }
                                }
                            }
                        } // end of go left
                    }
                    else if(dotr >0 && dotl > 0) // should go right(the point is on the right side of the current section)
                    {

                        while(1){
                            left++; right ++;

                            // check whether the section difference between neighboring points have reached the tolerance
                            if(difference++ >= tolerance)
                            {
                                SectionIds->SetId(id, right);
                                break;
                            }

                            if(right > model->GetNumberOfPoints())
                            {
                                std::cerr<<"failed to find section id (go right to the end)"<<endl;
                                exit(1);
                            }
                            PlaneOriginals->GetTuple(left, originleft);
                            PlaneNormals->GetTuple(left, normalleft);
                            vtkMath::Subtract(p, originleft, vl);
                            dotl = vtkMath::Dot(vl, normalleft);
                            if(right == model->GetNumberOfPoints())
                            {
                                if(dotl > 0)
                                {
                                    SectionIds->SetId(id, right);
                                    break;
                                }
                            }
                            else
                            {
                                PlaneOriginals->GetTuple(right, originright);
                                PlaneNormals->GetTuple(right, normalright);
                                vtkMath::Subtract(p, originright, vr);
                                dotr = vtkMath::Dot(vr, normalright);
                                {
                                    if(dotl > 0 && dotr <=0)
                                    {
                                        SectionIds->SetId(id, right);
                                        break;
                                    }
                                }
                            }
                        } // end of go right
                    }
                    else
                    {
                        std::cerr<<"failed to find section id (unreasonalbe dotl and dotr)"<<endl;
                        exit(0);
                    }
                }
                // after processing the id-th point, should process its neighbours
                std::cout<<"set "<<id<<" to "<<right<<endl;
            }
        } // for(vtkIdType n = 0; n < nextlevel->GetNumberOfIds(); n++)
        currentlevel->Reset();
        currentlevel->DeepCopy(nextlevel);
    } // while(1)
}
vtkSmartPointer<vtkIntArray> Centerline::GetSectionHeadAndTail(vtkSmartPointer<vtkDoubleArray> Curvatures, double thres)
{
    vtkSmartPointer<vtkIntArray> SectionHeadAndTail = vtkSmartPointer<vtkIntArray>::New();
    SectionHeadAndTail->SetNumberOfComponents(2);
    int left = -1, right = 0;
    for(int i = 0; i <= model->GetNumberOfPoints(); i++) // at most num+1 section
    {
        if(right < model->GetNumberOfPoints())
        {
            if(Curvatures->GetValue(right) > thres)
            {
                right++;
            }
            else
            {
                int tuple[2];
                tuple[0] = left; tuple[1] = right;
                SectionHeadAndTail->InsertNextTypedTuple(tuple);
                left = right++;
            }
        }
        else
        {
            int tuple[2];
            tuple[0] = left; tuple[1] = right;
            SectionHeadAndTail->InsertNextTypedTuple(tuple);
            break;
        }
    }
    return SectionHeadAndTail;
}
void Centerline::GetSectionIds_loop_combinehighcurvatures(vtkPolyData *t_colon, vtkIdType seed, vtkIdList *SectionIds,
                                                          vtkDoubleArray *PlaneOriginals, vtkDoubleArray *PlaneNormals, vtkDoubleArray *Curvatures)
{
    vtkSmartPointer<vtkIntArray> SectionHeadAndTail = vtkSmartPointer<vtkIntArray>::New();
    SectionHeadAndTail = GetSectionHeadAndTail(Curvatures, 0.01);
    int count = 0;
    int difference = 0, tolerance = INFINITY;
    vtkSmartPointer<vtkIdList> currentlevel = vtkSmartPointer<vtkIdList>::New();
    // init
    currentlevel->InsertNextId(seed);

    while(currentlevel->GetNumberOfIds() > 0)
    {
        count += currentlevel->GetNumberOfIds();
        std::cout<<"processed "<<count<<" points"<<endl;
        vtkSmartPointer<vtkIdList> nextlevel = vtkSmartPointer<vtkIdList>::New();
        for(vtkIdType n = 0; n < currentlevel->GetNumberOfIds(); n++)
        {
            vtkIdType pointid = currentlevel->GetId(n);
            vtkSmartPointer<vtkIdList> connectedVertices = vtkSmartPointer<vtkIdList>::New();
            connectedVertices = GetConnectedVertices(t_colon, pointid);
            vtkIdType currentSectionId = SectionIds->GetId(pointid);

            std::cout<<"pointid : "<<pointid<<"(current section "<<currentSectionId<<") "<<connectedVertices->GetNumberOfIds()<<": ";
            for(vtkIdType i = 0; i < connectedVertices->GetNumberOfIds(); i++)
            {
                std::cout<<connectedVertices->GetId(i)<<" ";
            }
            std::cout<<endl;

            for(vtkIdType i = 0; i < connectedVertices->GetNumberOfIds(); i++)
            {
                double p[3];
                double originright[3], normalright[3], vr[3];
                double originleft[3], normalleft[3], vl[3];
                double dotl, dotr;
                vtkIdType id = connectedVertices->GetId(i);
                difference = 0;
                t_colon->GetPoint(id, p);
                if(SectionIds->GetId(id) >= 0) // this point has already been assigned
                    continue; // do nothing for this point
                else
                    nextlevel->InsertNextId(id);
                //vtkIdType left = currentSectionId-1, right = currentSectionId;
                int tuple[2];
                SectionHeadAndTail->GetTypedTuple(currentSectionId, tuple);
                int left = tuple[0], right = tuple[1];
                if(left < 0) // id is at the left end
                {
                    PlaneOriginals->GetTuple(right, originright);
                    PlaneNormals->GetTuple(right, normalright);
                    vtkMath::Subtract(p, originright, vr);
                    if(vtkMath::Dot(vr, normalright) <= 0) // id satisfy the section criteron
                    {
                        SectionIds->SetId(id, right);
                    }
                    // else go right
                    else
                    {
                        while(1){
                            if(++currentSectionId >= SectionHeadAndTail->GetNumberOfTuples())
                            {
                                std::cerr<<"failed to find section id (go right to the end)"<<endl;
                                exit(1);
                            }
                            SectionHeadAndTail->GetTypedTuple(currentSectionId, tuple);
                            left = tuple[0]; right = tuple[1];
                            // check whether the section difference between neighboring points have reached the tolerance
                            if(difference++ >= tolerance)
                            {
                                SectionIds->SetId(id, right);
                                break;
                            }
                            PlaneOriginals->GetTuple(left, originleft);
                            PlaneNormals->GetTuple(left, normalleft);
                            vtkMath::Subtract(p, originleft, vl);
                            dotl = vtkMath::Dot(vl, normalleft);
                            if(right == model->GetNumberOfPoints())
                            {
                                if(dotl > 0)
                                {
                                    SectionIds->SetId(id, right);
                                    break;
                                }
                            }
                            else
                            {
                                PlaneOriginals->GetTuple(right, originright);
                                PlaneNormals->GetTuple(right, normalright);
                                vtkMath::Subtract(p, originright, vr);
                                dotr = vtkMath::Dot(vr, normalright);
                                {
                                    if(dotl > 0 && dotr <=0)
                                    {
                                        SectionIds->SetId(id, right);
                                        break;
                                    }
                                }
                            }
                        } // end of go right
                    }
                }
                else if(right >= model->GetNumberOfPoints()) // id is at the right end
                {
                    PlaneOriginals->GetTuple(left, originleft);
                    PlaneNormals->GetTuple(left, normalleft);
                    vtkMath::Subtract(p, originleft, vl);
                    if(vtkMath::Dot(vl, normalleft) > 0) // id satisfy the section criteron
                    {
                        SectionIds->SetId(id, right);
                    }
                    // else go left
                    else
                    {
                        while(1){
                            if(--currentSectionId < 0)
                            {
                                std::cerr<<"failed to find section id (go left to the end)"<<endl;
                                exit(1);
                            }
                            SectionHeadAndTail->GetTypedTuple(currentSectionId, tuple);
                            left = tuple[0]; right = tuple[1];
                            // check whether the section difference between neighboring points have reached the tolerance
                            if(difference++ >= tolerance)
                            {
                                SectionIds->SetId(id, right);
                                break;
                            }
                            PlaneOriginals->GetTuple(right, originright);
                            PlaneNormals->GetTuple(right, normalright);
                            vtkMath::Subtract(p, originright, vr);
                            dotr = vtkMath::Dot(vr, normalright);
                            if(left == -1)
                            {
                                if(dotr <= 0)
                                {
                                    SectionIds->SetId(id, right);
                                    break;
                                }
                            }
                            else
                            {
                                PlaneOriginals->GetTuple(left, originleft);
                                PlaneNormals->GetTuple(left, normalleft);
                                vtkMath::Subtract(p, originleft, vl);
                                dotl = vtkMath::Dot(vl, normalleft);
                                {
                                    if(dotl > 0 && dotr <=0)
                                    {
                                        SectionIds->SetId(id, right);
                                        break;
                                    }
                                }
                            }
                        } // end of go left
                    }
                }
                else // id is in the middle
                {
                    PlaneOriginals->GetTuple(right, originright);
                    PlaneNormals->GetTuple(right, normalright);
                    vtkMath::Subtract(p, originright, vr);
                    PlaneOriginals->GetTuple(left, originleft);
                    PlaneNormals->GetTuple(left, normalleft);
                    vtkMath::Subtract(p, originleft, vl);
                    dotr = vtkMath::Dot(vr, normalright);
                    dotl = vtkMath::Dot(vl, normalleft);
                    //std::cout<<"dotr: "<<dotr<<" dotl: "<<dotl<<std::endl;
                    if(dotr <= 0 && dotl > 0) // satisfy
                    {
                        SectionIds->SetId(id, right);
                    }
                    else if(dotr <= 0 && dotl <= 0) // should go left(the point is on the left side of the current section)
                    {
                        while(1){
                            if(--currentSectionId >= SectionHeadAndTail->GetNumberOfTuples())
                            {
                                std::cerr<<"failed to find section id (go left to the end)"<<endl;
                                exit(1);
                            }
                            SectionHeadAndTail->GetTypedTuple(currentSectionId, tuple);
                            left = tuple[0]; right = tuple[1];
                            // check whether the section difference between neighboring points have reached the tolerance
                            if(difference++ >= tolerance)
                            {
                                SectionIds->SetId(id, right);
                                break;
                            }
                            PlaneOriginals->GetTuple(right, originright);
                            PlaneNormals->GetTuple(right, normalright);
                            vtkMath::Subtract(p, originright, vr);
                            dotr = vtkMath::Dot(vr, normalright);
                            if(left == -1)
                            {
                                if(dotr <= 0)
                                {
                                    SectionIds->SetId(id, right);
                                    break;
                                }
                            }
                            else
                            {
                                PlaneOriginals->GetTuple(left, originleft);
                                PlaneNormals->GetTuple(left, normalleft);
                                vtkMath::Subtract(p, originleft, vl);
                                dotl = vtkMath::Dot(vl, normalleft);
                                {
                                    if(dotl > 0 && dotr <=0)
                                    {
                                        SectionIds->SetId(id, right);
                                        break;
                                    }
                                }
                            }
                        } // end of go left
                    }
                    else if(dotr >0 && dotl > 0) // should go right(the point is on the right side of the current section)
                    {

                        while(1){
                            if(++currentSectionId >= SectionHeadAndTail->GetNumberOfTuples())
                            {
                                std::cerr<<"failed to find section id (go right to the end)"<<endl;
                                exit(1);
                            }
                            SectionHeadAndTail->GetTypedTuple(currentSectionId, tuple);
                            left = tuple[0]; right = tuple[1];
                            // check whether the section difference between neighboring points have reached the tolerance
                            if(difference++ >= tolerance)
                            {
                                SectionIds->SetId(id, right);
                                break;
                            }
                            PlaneOriginals->GetTuple(left, originleft);
                            PlaneNormals->GetTuple(left, normalleft);
                            vtkMath::Subtract(p, originleft, vl);
                            dotl = vtkMath::Dot(vl, normalleft);
                            if(right == model->GetNumberOfPoints())
                            {
                                if(dotl > 0)
                                {
                                    SectionIds->SetId(id, right);
                                    break;
                                }
                            }
                            else
                            {
                                PlaneOriginals->GetTuple(right, originright);
                                PlaneNormals->GetTuple(right, normalright);
                                vtkMath::Subtract(p, originright, vr);
                                dotr = vtkMath::Dot(vr, normalright);
                                {
                                    if(dotl > 0 && dotr <=0)
                                    {
                                        SectionIds->SetId(id, right);
                                        break;
                                    }
                                }
                            }
                        } // end of go right
                    }
                    else
                    {
                        std::cerr<<"failed to find section id (unreasonalbe dotl and dotr)"<<endl;
                        exit(0);
                    }
                }
                // after processing the id-th point, should process its neighbours
                std::cout<<"set "<<id<<" to "<<right<<endl;
            }
        } // for(vtkIdType n = 0; n < nextlevel->GetNumberOfIds(); n++)
        currentlevel->Reset();
        currentlevel->DeepCopy(nextlevel);
    } // while(1)
}

double Centerline::LorentzianInterpolationFactor(double l, double angle)
{
    assert(angle >= 0);
    assert(l >= 0 && l <= 1);
    if(angle == 0)
        return l;
    else
        return sin(l * angle) / sin(angle);
}

void Centerline::GetRegionIds(vtkPolyData *t_colon, bool *Is_Fixed, int *RegionIds, int *NumberOfRegions)
{
    NumberOfRegions[0] = 0;
    int region = 0;
    int N = t_colon->GetNumberOfPoints();
    for(int i=0; i<N; i++)
    {
        if(!Is_Fixed[i] && RegionIds[i] == -1)
        {
            vtkSmartPointer<vtkIdList> currentlevel = vtkSmartPointer<vtkIdList>::New();
            currentlevel->InsertNextId(i);
            while(currentlevel->GetNumberOfIds() > 0)
            {
                vtkSmartPointer<vtkIdList> nextlevel = vtkSmartPointer<vtkIdList>::New();
                for(int j=0; j < currentlevel->GetNumberOfIds(); j++)
                {
                    int id = currentlevel->GetId(j);
                    assert(!Is_Fixed[id]);
                    if(RegionIds[id] == -1)
                    {
                        RegionIds[id] = region;
                        vtkSmartPointer<vtkIdList> ptids = vtkSmartPointer<vtkIdList>::New();
                        ptids = GetConnectedVertices(t_colon, id);
                        for(int k = 0; k < ptids->GetNumberOfIds(); k++)
                        {
                            int nbid = ptids->GetId(k);
                            assert(nbid != id);
                            if(!Is_Fixed[nbid] && RegionIds[nbid] == -1)
                            {
                                nextlevel->InsertNextId(nbid);
                            }
                        }
                    }
                    //std::cout<<endl;
                }
                currentlevel->Reset();
                currentlevel->DeepCopy(nextlevel);
            }
            region++;
        }
    }
    NumberOfRegions[0] = region;
}
