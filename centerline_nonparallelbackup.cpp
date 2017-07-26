#include "centerline.h"
#include <iostream>
#include <vtkMath.h>
#include <unistd.h>
Centerline::Centerline()
{
    actor->GetProperty()->SetColor(1,1,0);
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
void Centerline::SmoothCenterline(double std)
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
    std::cout<<NumSmoothedPoints<<endl;
    double gp[3];
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
    for(vtkIdType i = 0; i <NumPoints; i++)
    {
        gp[0] = gx[i + GaussianKernelLen -1];
        gp[1] = gy[i + GaussianKernelLen -1];
        gp[2] = gz[i + GaussianKernelLen -1];
        //std::cout<<i<<" "<<gp[0]<<" "<<gp[1]<<" "<<gp[2]<<" "<<std::endl;
        //std::cout<<i<<" "<<x[i + (GaussianKernelLen-1)/2]<<" "<<y[i + (GaussianKernelLen-1)/2]<<" "<<z[i + (GaussianKernelLen-1)/2]<<" "<<std::endl;
        //std::cout<<i<<" "<<model->GetPoint(i)[0]<<" "<<model->GetPoint(i)[1]<<" "<<model->GetPoint(i)[2]<<" "<<std::endl;
        points->InsertNextPoint(gp);

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
    std::cout<<"centerline smoothed:"<<model->GetNumberOfPoints()<<std::endl;

    free(x);
    free(y);
    free(z);
    free(gx);
    free(gy);
    free(gz);
    free(G);
}
void Centerline::UniformSample(int resolution)
{
    vtkSmartPointer<vtkParametricSpline> spline = vtkSmartPointer<vtkParametricSpline>::New();
    spline->SetPoints(model->GetPoints());
    vtkSmartPointer<vtkParametricFunctionSource> functionSource = vtkSmartPointer<vtkParametricFunctionSource>::New();
    functionSource->SetParametricFunction(spline);
    functionSource->Update();
    functionSource->SetUResolution(resolution);
    functionSource->Update();
    model->DeepCopy(functionSource->GetOutput());
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
                                                     vtkSmartPointer<vtkDoubleArray> Tangents, vtkSmartPointer<vtkDoubleArray> Normals, vtkSmartPointer<vtkDoubleArray> Binormals,
                                                     vtkSmartPointer<vtkPolyData> t_colon)
{
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
        dnormal[0] = tangent[0]; dnormal[1] = tangent[1]; dnormal[2] = tangent[2];
        double RelaxationFactor = 2;
        vtkMath::MultiplyScalar(dnormal, -curvature / RelaxationFactor);
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
    //t_rendermanager->renderModel(newActor);

    // Simple Deformation
    // Get the old coordinates
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
            colors->InsertNextTypedTuple(blue);
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
vtkSmartPointer<vtkPolyData> Centerline::EliminateTorsion(RenderManager* t_rendermanager, vtkSmartPointer<vtkPolyData> t_colon)
{
    bool use_spline = true;  // whether use spline
    double stepSize = 0.0005;
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
    vtkSmartPointer<vtkPoints> ViolationPoints = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkDoubleArray> Directions = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkIntArray> ViolationNums = vtkSmartPointer<vtkIntArray>::New();
    vtkSmartPointer<vtkAppendPolyData> appendFilter = vtkSmartPointer<vtkAppendPolyData>::New();
    vtkSmartPointer<vtkCleanPolyData> cleanFilter = vtkSmartPointer<vtkCleanPolyData>::New();
    vtkSmartPointer<vtkPolyData> IllCutCircles = vtkSmartPointer<vtkPolyData>::New();

    int MaxIter = 1; int modify = 0;
    for(int iter = 0; iter < MaxIter; iter++)
    {
        std::cout<<"Iteration: "<<iter<<endl;
        Tangents->Reset();Tangents->SetNumberOfComponents(3);
        Normals->Reset();Normals->SetNumberOfComponents(3);
        Binormals->Reset();Binormals->SetNumberOfComponents(3);
        Curvatures->Reset();
        Torsions->Reset();
        S->Reset();
        U->Reset();
        CurvaturePoints->Reset();
        ViolationPoints->Reset();
        Directions->Reset();
        ViolationNums->Reset();
        spline->RemoveAllObservers();
        spline->SetPoints(model->GetPoints());

        // Get the length as parameter on each point
        double cumS = 0;
        S->InsertNextValue(cumS);
        for(vtkIdType i = 1; i< model->GetNumberOfPoints(); i++)
        {
            double p[3], lastp[3], distance;
            model->GetPoint(i-1, lastp);
            model->GetPoint(i, p);
            distance = sqrt(vtkMath::Distance2BetweenPoints(p, lastp));
            cumS += distance;
            S->InsertNextValue(cumS);
            //std::cout<<i<<" distance "<<distance<<endl;
        }
        std::cout<<"S max:"<<S->GetValue(model->GetNumberOfPoints()-1)<<endl;
        // Get the normalized parameter, used in spline
        for(vtkIdType i = 0; i< model->GetNumberOfPoints(); i++)
        {
            double u;
            u = S->GetValue(i) / cumS;
            U->InsertNextValue(u);
        }
        std::cout<<"U max:"<<U->GetValue(model->GetNumberOfPoints()-1)<<endl;
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
            Tangents->InsertNextTuple(tangent);
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
            Normals->InsertNextTuple(normal);
            Curvatures->InsertNextValue(curvature);
        }

        //PutNormalsOnSameSide(Normals, Curvatures);

        // Get the binormal on each point
        for(vtkIdType i = 0; i<model->GetNumberOfPoints(); i++)
        {
            double tangent[3], normal[3], binormal[3];
            Tangents->GetTuple(i, tangent);
            Normals->GetTuple(i, normal);
            vtkMath::Cross(tangent, normal, binormal);
            vtkMath::Normalize(binormal);
            Binormals->InsertNextTuple(binormal);
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
            Torsions->InsertNextValue(torsion);
            //std::cout<<"torsion: "<<i<<" "<<torsion<<endl;
        }


        // cut circle
        vtkSmartPointer<vtkCutter> cutter = vtkSmartPointer<vtkCutter>::New();
        cutter->SetInputData(t_colon);
        vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();
        vtkSmartPointer<vtkPolyData> cutline = vtkSmartPointer<vtkPolyData>::New();

        vtkSmartPointer<vtkPolyDataConnectivityFilter> connectivityFilter = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
        connectivityFilter->SetInputConnection(cutter->GetOutputPort());
        connectivityFilter->SetExtractionModeToClosestPointRegion();

        vtkSmartPointer<vtkDoubleArray> Radius = vtkSmartPointer<vtkDoubleArray>::New();

        // calculate the cut circles and curvature points
        double sumr = 0;
        for(vtkIdType i=0; i<model->GetNumberOfPoints(); i += 1)
        {
            double centerPoint[3];
            model->GetPoint(i, centerPoint);
            connectivityFilter->SetClosestPoint(centerPoint);
            plane->SetOrigin(model->GetPoint(i));
            plane->SetNormal(Tangents->GetTuple(i));

            cutter->SetCutFunction(plane);
            cutter->Update();
            connectivityFilter->Update();

            cutline = connectivityFilter->GetOutput();
            //std::cout<<i<<" "<<"cutline points:"<<cutline->GetNumberOfPoints()<<endl;

            double curvaturePoint[3], p[3], v[3], normal[3], angleCos, maxScore = -INFINITY;
            int violationNum = 0, innerNum = 0;
            double direction[3] = {0,0,0};
            double averageCircleR = 0;
            for(vtkIdType j=0; j<cutline->GetNumberOfPoints(); j++)
            {
                cutline->GetPoint(j, p);
                vtkMath::Subtract(p, centerPoint, v);
                averageCircleR += vtkMath::Norm(v);
            }
            averageCircleR /= cutline->GetNumberOfPoints();

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
                    if(r*Krel > 1)
                    {
                        violationNum++;
                        sign = 1;
                        ViolationPoints->InsertNextPoint(p);
                    }
                    double rr = (r < averageCircleR)? r : averageCircleR;
                    double step = 0.05 * pow(rr , 3) * (sign * pow((Krel - 1/r),2) + 0.2 * Krel * Krel);
                    vtkMath::MultiplyScalar(v, step);
                    vtkMath::Add(direction, v, direction);
                }
                if(angleCos > maxScore)
                {
                    maxScore = angleCos;
                    curvaturePoint[0] = p[0];
                    curvaturePoint[1] = p[1];
                    curvaturePoint[2] = p[2];
                }
            }
            if(innerNum > 0)
                vtkMath::MultiplyScalar(direction, 1/innerNum);
            Directions->InsertNextTuple(direction);
            ViolationNums->InsertNextValue(violationNum);

            CurvaturePoints->InsertNextPoint(curvaturePoint);

            double r = (sqrt(vtkMath::Distance2BetweenPoints(curvaturePoint, centerPoint)));
            Radius->InsertNextValue(r);
            if(violationNum > 0)
            {
                appendFilter->RemoveAllInputs();
                appendFilter->AddInputData(cutline);
                appendFilter->AddInputData(IllCutCircles);
                appendFilter->Update();
                cleanFilter->SetInputConnection(appendFilter->GetOutputPort());
                cleanFilter->Update();
                IllCutCircles->DeepCopy(cleanFilter->GetOutput());
                //std::cout<<i<<" "<<IllCutCircles->GetNumberOfPoints()<<"("<<cutline->GetNumberOfPoints()<<")"<<endl;
            }
            sumr+=r;
        }
        if(iter == MaxIter - 1)
        {
            VisualizeTNB(S, Curvatures, Tangents, Normals, Binormals, t_rendermanager);
            VisualizeSpoke(CurvaturePoints, ViolationNums, t_rendermanager);
        }

        // modify the points on curve along the Normal direction
        /*
    if(modify)
    {
        double aver = sumr / model->GetNumberOfPoints();
        vtkSmartPointer<vtkPoints> NewPoints = vtkSmartPointer<vtkPoints>::New();
        for(vtkIdType i = 0; i < model->GetNumberOfPoints(); i += 1)
        {
            double n[3], p[3], np[3], step, k, r, sign = 0, weight = 0.2;
            model->GetPoint(i, p);
            k = Curvatures->GetValue(i);
            r = Radius->GetValue(i);
            if(k > 1/r) sign = 1;
            if(r < aver)
                step =  0.05 * pow(r, (double)3) * (sign * (k - 1/r)*(k - 1/r) + weight * k * k);
            else
                step =  0.05 * pow(aver, (double)3) * (sign * (k - 1/r)*(k - 1/r) + weight * k * k);
            Normals->GetTuple(i, n);
            vtkMath::MultiplyScalar(n, step);
            vtkMath::Add(p, n, np);
            NewPoints->InsertNextPoint(np);
        }
        model->SetPoints(NewPoints);
        SmoothCenterline(3);
    }
    */
        // modify the points on curve along the weighted sum direction
        if(modify)
        {
            vtkSmartPointer<vtkPoints> NewPoints = vtkSmartPointer<vtkPoints>::New();
            for(vtkIdType i = 0; i < model->GetNumberOfPoints(); i += 1)
            {
                double direction[3];
                Directions->GetTuple(i, direction);
                double p[3], newp[3];
                model->GetPoint(i, p);
                vtkMath::Add(p, direction, newp);
                NewPoints->InsertNextPoint(newp);
            }
            model->SetPoints(NewPoints);
            SmoothCenterline(3);
        }
        std::cout<<"Iteration Ends: "<<iter<<" "<<model->GetNumberOfPoints()<<" "<<cumS<<endl;
    }
    // visualize the ill cut circles
    vtkSmartPointer<vtkPolyDataMapper> IllCutCirclesMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    IllCutCirclesMapper->SetInputData(IllCutCircles);
    IllCutCirclesMapper->Update();
    vtkSmartPointer<vtkActor> IllCutCirclesActor = vtkSmartPointer<vtkActor>::New();
    IllCutCirclesActor->SetMapper(IllCutCirclesMapper);
    IllCutCirclesActor->GetProperty()->SetColor(0, 1, 1);
    t_rendermanager->renderModel(IllCutCirclesActor);

    // visualize the violation points
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

    //return Deformation(S, Curvatures, Tangents, Normals, Binormals, t_rendermanager, t_colon);
    return NULL;
}

