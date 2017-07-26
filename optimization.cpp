#include "optimization.h"

vtkSmartPointer<vtkPolyData> Optimize(vtkSmartPointer<vtkPolyData> t_colon, vtkSmartPointer<vtkPolyData> SurfaceLineUp, bool *Is_Fixed, RenderManager *t_rendermanager)
{
    vtkSmartPointer<vtkPolyData> OptimizedSurface = vtkSmartPointer<vtkPolyData>::New();
    OptimizedSurface->DeepCopy(SurfaceLineUp);
    // maintain a list to record the correspondence between original ids and the variable ids
    vtkSmartPointer<vtkIdList> Ids = vtkSmartPointer<vtkIdList>::New();
    int* InvertIds = (int*)malloc(t_colon->GetNumberOfPoints()*sizeof(int));
    for(int i=0; i < t_colon->GetNumberOfPoints(); i++)
    {
        InvertIds[i] = -1;
    }
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for(vtkIdType i=0; i < t_colon->GetNumberOfPoints(); i++)
    {
        if(!Is_Fixed[i])
        {
            double p[3];
            t_colon->GetPoint(i, p);
            Ids->InsertNextId(i);
            points->InsertNextPoint(p);
            InvertIds[i] = Ids->GetNumberOfIds() - 1;
        }
    }

    // check whether there's nan
    /*
    for(int i=0; i<SurfaceLineUp->GetNumberOfPoints(); i++)
    {
        double p[3];
        SurfaceLineUp->GetPoint(i, p);
        std::cout<<i<<" "<<p[0]<<" "<<p[1]<<" "<<p[2]<<endl;
    }
    */
    // visualize the unfixed points
    vtkSmartPointer<vtkPoints> unfixedpoints = vtkSmartPointer<vtkPoints>::New();
    for(int i=0; i < t_colon->GetNumberOfPoints(); i++)
    {
        if(InvertIds[i] > -0.5)
        {
            int id = i;
            double p[3];
            t_colon->GetPoint(id, p);
            unfixedpoints->InsertNextPoint(p);
        }
    }
    std::cout<<"unfixed points: "<<unfixedpoints->GetNumberOfPoints()<<endl;


    int m = Ids->GetNumberOfIds() * 3;
    std::cout<<"total number of variables = "<<m<<endl;
    std::map<vtkIdType, double> coefficientMap;
    coefficientMap.clear();

    double* b_c = (double *)malloc(m*sizeof(double));
    memset(b_c, 0, sizeof(double)*m);

    for(int i=0; i < m; i++)
    {
        if(isnan(b_c[i]))
        std::cout<<Ids->GetId(i/3)<<"()"<<b_c[i]<<endl;
    }
    constructAandb(coefficientMap, b_c, t_colon, SurfaceLineUp, Is_Fixed, Ids, InvertIds, t_rendermanager);

    /*
    for(int i=0; i < m; i++)
    {
        std::cout<<i<<" "<<b_c[i]<<endl;
    }
    */


    // transform coefficientMap into a data structure that eigen recognizes
    std::vector<T> coefficients;
    coefficients.clear();
    for(std::map<vtkIdType, double>::iterator it = coefficientMap.begin(); it!=coefficientMap.end(); it++)
    {
        int idx1 = (vtkIdType)it->first % (vtkIdType)m;
        int idx2 = floor((double(it->first) + 0.01) / m);

        coefficients.push_back(T(idx1, idx2, it->second));
        //if(Ids->GetNumberOfIds() < 300)
          //  std::cout<<it->first<<" "<<idx1<<" "<<idx2<<" "<<it->second<<endl;
    }
    SpMat A(m,m);
    A.setFromTriplets(coefficients.begin(), coefficients.end());

    Eigen::VectorXd b(m);
    b.setZero();
    for(int k = 0; k < m; k++)
    {
        b(k) = b_c[k];
    }

    //Eigen::BiCGSTAB<SpMat, Eigen::IncompleteLUT<double> > solver;
    //solver.preconditioner().setDroptol(0.001);
    //Eigen::BiCGSTAB<SpMat >  solver;
    Eigen::SimplicialCholesky<SpMat> solver;
    solver.compute(A);
    Eigen::VectorXd solved_x = solver.solve(b);

    vtkSmartPointer<vtkPoints> newpoints = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkPoints> checkpoints = vtkSmartPointer<vtkPoints>::New();
    newpoints = SurfaceLineUp->GetPoints();
    for(int vidx = 0; vidx < Ids->GetNumberOfIds(); vidx++)
    {
        double v[3], pold[3], pnew[3];
        int idx = Ids->GetId(vidx);
        t_colon->GetPoint(idx, pold);
        v[0] = solved_x(vidx*3);
        v[1] = solved_x(vidx*3+1);
        v[2] = solved_x(vidx*3+2);
        //std::cout<<vidx<<" "<<pold[0]<<" "<<pold[1]<<" "<<pold[2]<<" - "<<v[0]<<" "<<v[1]<<" "<<v[2]<<endl;
        vtkMath::Add(pold, v, pnew);
        newpoints->SetPoint(idx, pnew);

        if(isnan(v[0]) || isnan(v[1]) || isnan(v[2]))
        {
            checkpoints->InsertNextPoint(pold);
            std::cout<<checkpoints->GetNumberOfPoints()<<" "<<Ids->GetId(vidx)<<endl;
        }
    }
    VisualizePoints(checkpoints, (double)1, (double)0, (double)0, (float)5, t_rendermanager);
    OptimizedSurface->SetPoints(newpoints);

    free(b_c);
    free(InvertIds);
    return OptimizedSurface;
}

void updateA(int m, int i1, int i2, double weight,std::map<vtkIdType, double> &coefficientMap)
{
    m = (vtkIdType)m;
    i1 = (vtkIdType)i1;
    i2 = (vtkIdType)i2;
    //if(i1 == 17080 && i2 == 17080)
    //std::cout<<"updatea "<<i1<<"-"<<i2<<endl;
    std::map<vtkIdType, double>::iterator iter;
    vtkIdType k = (i1*(vtkIdType)3) + (i2*(vtkIdType)3) * m;
    //std::cout<<k<<" ";

    iter = coefficientMap.find(k);
    if (iter == coefficientMap.end()){
        coefficientMap.insert(std::pair<vtkIdType,double>(k,weight));
    }else{
        coefficientMap[k] = iter->second + weight;
    }

    k = (i1*(vtkIdType)3+(vtkIdType)1) + (i2*(vtkIdType)3+(vtkIdType)1) * m;
    //std::cout<<k<<" ";

    iter = coefficientMap.find(k);
    if (iter == coefficientMap.end()){
        coefficientMap.insert(std::pair<vtkIdType,double>(k,weight));
    }else{
        coefficientMap[k] = iter->second + weight;
    }

    k = (i1*(vtkIdType)3+(vtkIdType)2) + (i2*(vtkIdType)3+(vtkIdType)2) * m;
    //if(i1 == 17080 && i2 == 17080)
    //std::cout<<k<<"- "<<i1*(vtkIdType)3+(vtkIdType)2<<" "<<i2*(vtkIdType)+(vtkIdType)2<<endl;

    iter = coefficientMap.find(k);
    if (iter == coefficientMap.end()){
        coefficientMap.insert(std::pair<vtkIdType,double>(k,weight));
    }else{
        coefficientMap[k] = iter->second + weight;
    }
}
vtkSmartPointer<vtkIdList> GetFacetsOfEdge(vtkSmartPointer<vtkPolyData> mesh, int idx1, int idx2)
{
    vtkSmartPointer<vtkIdList> list = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkIdList> list1 = vtkSmartPointer<vtkIdList>::New();
    mesh->GetPointCells(idx1, list1);
    vtkSmartPointer<vtkIdList> list2 = vtkSmartPointer<vtkIdList>::New();
    mesh->GetPointCells(idx2, list2);

    for(int j = 0; j < list1->GetNumberOfIds(); j++)
    {
        for(int k = 0; k < list2->GetNumberOfIds(); k++)
        {
            if(list1->GetId(j) == list2->GetId(k))
            {
                list->InsertNextId(list1->GetId(j));
            }
        }
    }
    return list;
}

double ComputeStretchWeight(vtkSmartPointer<vtkPolyData> t_colon, int idx1, int idx2)
{
    vtkSmartPointer<vtkIdList> list = vtkSmartPointer<vtkIdList>::New();
    list = GetFacetsOfEdge(t_colon, idx1, idx2);

    double area = 0;
    if(list->GetNumberOfIds() == 1 || list->GetNumberOfIds() == 2)
    {
        for(int i=0; i < list->GetNumberOfIds(); i++)
        {
            area += ComputeFacetArea(t_colon, list->GetId(i));
        }
    }
    else
    {
        std::cerr<<"error, list number not correct: "<<list->GetNumberOfIds()<<endl;
        exit(1);
    }
    double p1[3], p2[3];
    t_colon->GetPoint(idx1, p1);
    t_colon->GetPoint(idx2, p2);
    double l = sqrt(vtkMath::Distance2BetweenPoints(p1, p2));
    //std::cout<<list->GetNumberOfIds()<<" - "<<area<<" "<<l<<endl;
    double rv = area / l / l;
    if(isnan(rv))
    {
        std::cout<<idx1<<" "<<p1[0]<<" "<<p1[1]<<" "<<p1[2]<<endl;
        std::cout<<idx2<<" "<<p2[0]<<" "<<p2[1]<<" "<<p2[2]<<endl;
    }
    return area / l / l;
}

void constructAandb(std::map<vtkIdType, double> &coefficientMap, double *b,
                vtkSmartPointer<vtkPolyData> t_colon, vtkSmartPointer<vtkPolyData> SurfaceLineUp,
                bool *Is_Fixed, vtkSmartPointer<vtkIdList> Ids, int* InvertIds, RenderManager *t_rendermanager)
{
    int m = Ids->GetNumberOfIds() * 3;

    // stretch
    // extract the edges
    bool * marked = (bool*)malloc(t_colon->GetNumberOfPoints()*sizeof(bool));
    memset(marked, 0, t_colon->GetNumberOfPoints()*sizeof(bool));
    bool * boundary = (bool*)malloc(t_colon->GetNumberOfPoints()*sizeof(bool));
    memset(boundary, 0, t_colon->GetNumberOfPoints()*sizeof(bool));

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

    /* //check whether the lines are correct
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points = t_colon->GetPoints();
    vtkSmartPointer<vtkCellArray> liness = vtkSmartPointer<vtkCellArray>::New();
    for(int i=0; i < lines->GetNumberOfTuples(); i++)
    {
        int tuple[2];
        lines->GetTypedTuple(i, tuple);
        vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
        line->GetPointIds()->SetId(0, tuple[0]);
        line->GetPointIds()->SetId(1, tuple[1]);
        liness->InsertNextCell(line);
    }
    vtkSmartPointer<vtkPolyData> poly = vtkSmartPointer<vtkPolyData>::New();
    poly->SetPoints(points);
    //poly->SetPolys(liness);
    poly->SetLines(liness);
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(poly);
    mapper->Update();
    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(1,0,0);
    t_rendermanager->renderModel(actor);
    */




    for(int i=0; i < lines->GetNumberOfTuples(); i++)
    {
        int tuple[2];
        lines->GetTypedTuple(i, tuple);
        int idx1 = tuple[0];
        int idx2 = tuple[1];
        //std::cout<<InvertIds[idx1]<<" "<<InvertIds[idx2]<<endl;


        if(GetFacetsOfEdge(t_colon, idx1, idx2)->GetNumberOfIds() <= 1)
        {
            //std::cout<<idx1<<" "<<idx2<<" "<<GetFacetsOfEdge(t_colon, idx1, idx2)->GetNumberOfIds()<<endl;
            boundary[idx1] = true;
            boundary[idx2] = true;
        }

        if(!Is_Fixed[idx1] && !Is_Fixed[idx2]) // if the two points are all flexible, we only need to update A
        {
            int vidx1 = InvertIds[idx1];
            int vidx2 = InvertIds[idx2];
            //std::cout<<m<<" "<<vidx1<<" "<<vidx2<<endl;
            assert(vidx1 >= -0.5 && vidx2 >= -0.5);
            //std::cout<<vidx1<<","<<vidx2<<" - "<<idx1<<","<<idx2<<endl;
            double stretchWeight = ComputeStretchWeight(t_colon, idx1, idx2);
            //std::cout<<"stretch weight = "<<stretchWeight<<endl;
            double weight = REGWEIGHT * stretchWeight * 4;
            if(isnan(weight))
                std::cout<<"nan weight A "<<idx1<<" "<<idx2<<endl;
            updateA(m, vidx1, vidx1, weight, coefficientMap);
            updateA(m, vidx2, vidx2, weight, coefficientMap);
            updateA(m, vidx1, vidx2, -weight, coefficientMap);
            updateA(m, vidx2, vidx1, -weight, coefficientMap);
        }
        else if(!Is_Fixed[idx1] && Is_Fixed[idx2]) // if one point is fixed, we need to update A for the flexible point and add a linear term in b
        {
            int vidx1 = InvertIds[idx1];
            //std::cout<<m<<" "<<vidx1<<endl;
            assert(vidx1 > -0.5);
            double stretchWeight = ComputeStretchWeight(t_colon, idx1, idx2);
            double weight = REGWEIGHT * stretchWeight *4;
            if(isnan(weight))
                std::cout<<"nan weight b1 "<<idx1<<" "<<idx2<<endl;
            updateA(m, vidx1, vidx1, weight, coefficientMap);
            double v[3], pold[3], pnew[3];
            t_colon->GetPoint(idx2, pold);
            SurfaceLineUp->GetPoint(idx2, pnew);
            vtkMath::Subtract(pnew, pold, v);
            //std::cout<<"1- "<<v[0]<<" "<<v[1]<<" "<<v[2]<<endl;
            b[vidx1*3] +=   2 * v[0] * weight/2;
            b[vidx1*3+1] += 2 * v[1] * weight/2;
            b[vidx1*3+2] += 2 * v[2] * weight/2;
            if(isnan(b[vidx1*3]))
                std::cout<<idx1<<" - "<<b[vidx1*3]<<endl;
        }
        else if(Is_Fixed[idx1] && !Is_Fixed[idx2])
        {
            int vidx2 = InvertIds[idx2];
            //std::cout<<m<<" "<<vidx2<<endl;
            assert(vidx2 > -0.5);
            double stretchWeight = ComputeStretchWeight(t_colon, idx1, idx2);
            double weight = REGWEIGHT * stretchWeight * 4;
            if(isnan(weight))
                std::cout<<"nan weight b2 "<<idx1<<" "<<idx2<<endl;
            updateA(m, vidx2, vidx2, weight, coefficientMap);
            double v[3], pold[3], pnew[3];
            t_colon->GetPoint(idx1, pold);
            SurfaceLineUp->GetPoint(idx1, pnew);
            vtkMath::Subtract(pnew, pold, v);
            //std::cout<<"2- "<<v[0]<<" "<<v[1]<<" "<<v[2]<<endl;
            b[vidx2*3] +=   2 * v[0] * weight/2;
            b[vidx2*3+1] += 2 * v[1] * weight/2;
            b[vidx2*3+2] += 2 * v[2] * weight/2;
            if(isnan(b[vidx2*3]))
                std::cout<<idx2<<" - "<<b[vidx2*3]<<endl;
        }
    }

    // bend
    // compute surround areas
    vtkSmartPointer<vtkDoubleArray> SurroundAreas = vtkSmartPointer<vtkDoubleArray>::New();
    for(int i = 0; i < t_colon->GetNumberOfPoints(); i++)
    {
        double surroundarea = 0;
        vtkSmartPointer<vtkIdList> cellids = vtkSmartPointer<vtkIdList>::New();
        t_colon->GetPointCells(i, cellids);
        for(int j = 0; j < cellids->GetNumberOfIds(); j++)
        {
            surroundarea += ComputeFacetArea(t_colon, cellids->GetId(j));
        }
        SurroundAreas->InsertNextValue(surroundarea);
    }

    vtkSmartPointer<vtkPoints> boundarypoints = vtkSmartPointer<vtkPoints>::New();
    for(int i = 0; i < t_colon->GetNumberOfPoints(); i++)
    {
        int idx1 = i;
        if(boundary[idx1]) // don't process boundary edges
        {
            double p[3];
            t_colon->GetPoint(i, p);
            boundarypoints->InsertNextPoint(p);
            continue;
        }
        if(!Is_Fixed[idx1])
        {
            //std::cout<<i<<" - "<<idx1<<endl;
            //double weight = 2* BENDWEIGHT / (SurroundAreas->GetValue(idx1) * SurroundAreas->GetValue(idx1));
            double weight = BENDWEIGHT / (SurroundAreas->GetValue(idx1));
            //std::cout<<"bend weight "<<weight<<endl;
            //weight = 1;
            double center_weight = 0;
            vtkSmartPointer<vtkIdList> neighbours = vtkSmartPointer<vtkIdList>::New();
            neighbours = GetConnectedVertices(t_colon, idx1);
            //if(neighbours->GetNumberOfIds() <= 3)
            //std::cout<<"neighbours:"<<neighbours->GetNumberOfIds()<<endl;
            vtkSmartPointer<vtkDoubleArray> neighbor_weights = vtkSmartPointer<vtkDoubleArray>::New();
            for(int j = 0; j < neighbours->GetNumberOfIds(); j++)
            {
                int idx2 = neighbours->GetId(j);
                double neighbor_weight = 0;
                // compute the weight[j] for the j-th neighbour
                vtkSmartPointer<vtkIdList> cellids = vtkSmartPointer<vtkIdList>::New();
                cellids = GetFacetsOfEdge(t_colon, idx1, idx2);
                assert(cellids->GetNumberOfIds() == 2); // we have already excluded the boundary edges, so this edge must be inside the surface

                for(int k = 0; k < 2; k++)
                {
                    int idx12;
                    vtkSmartPointer<vtkIdList> ptids = vtkSmartPointer<vtkIdList>::New();
                    t_colon->GetCellPoints(cellids->GetId(k), ptids);
                    assert(ptids->GetNumberOfIds() == 3);
                    for(int l = 0; l < ptids->GetNumberOfIds(); l++)
                    {
                        if(ptids->GetId(l)!=idx1 && ptids->GetId(l)!=idx2)
                        {
                            idx12 = ptids->GetId(l);
                            break;
                        }
                    }

                    double p1[3], p2[3], p12[3];
                    t_colon->GetPoint(idx1, p1);
                    t_colon->GetPoint(idx2, p2);
                    t_colon->GetPoint(idx12, p12);
                    double v1[3], v2[3];
                    vtkMath::Subtract(p1, p12, v1);
                    vtkMath::Subtract(p2, p12, v2);
                    double angle = vtkMath::AngleBetweenVectors(v1, v2);

                    neighbor_weight += cos(angle)/sin(angle);
                    //std::cout<<idx1<<" "<<idx2<<" "<<idx12<<" : "<<cos(angle)/sin(angle)<<endl;
                }
                //neighbor_weight = (neighbor_weight > 0.1)? neighbor_weight : 0.1;
                //neighbor_weight = 1;
                //std::cout<<idx2<<" "<<neighbor_weight<<" | ";
                center_weight += neighbor_weight;
                neighbor_weights->InsertNextValue(neighbor_weight);
            }
            //std::cout<<center_weight<<endl;
            int vidx1 = InvertIds[idx1];
            updateA(m, vidx1, vidx1, weight* center_weight * center_weight, coefficientMap);
            //std::cout<<vidx1*3<<" "<<vidx1*3<<" "<<weight* center_weight * center_weight<<endl;

            for(int j = 0; j < neighbours->GetNumberOfIds(); j++)
            {
                int idx2 = neighbours->GetId(j);
                assert(idx2!=idx1);
                int vidx2 = InvertIds[idx2];
                if(!Is_Fixed[idx2])
                {
                    updateA(m, vidx1, vidx2, -weight*center_weight*neighbor_weights->GetValue(j), coefficientMap);
                    //std::cout<<vidx1<<" "<<vidx2<<" "<<-weight*center_weight*neighbor_weights->GetValue(j)<<endl;
                    updateA(m, vidx2, vidx1, -weight*center_weight*neighbor_weights->GetValue(j), coefficientMap);
                    //std::cout<<vidx2<<" "<<vidx1<<" "<<-weight*center_weight*neighbor_weights->GetValue(j)<<endl;
                }
                else
                {
                    double v[3], pold[3], pnew[3];
                    t_colon->GetPoint(idx2, pold);
                    SurfaceLineUp->GetPoint(idx2, pnew);
                    vtkMath::Subtract(pnew, pold, v);

                    // update the b of idx1

                    b[vidx1*3] += 2 * weight/2* center_weight * neighbor_weights->GetValue(j) * v[0];
                    b[vidx1*3+1] += 2 * weight/2 * center_weight * neighbor_weights->GetValue(j) * v[1];
                    b[vidx1*3+2] += 2 * weight/2 * center_weight * neighbor_weights->GetValue(j) * v[2];

                    //testb
                    /*
                    b[vidx1*3] += 2 * weight/2* center_weight * neighbor_weights->GetValue(j) * -50;
                    b[vidx1*3+1] += 2 * weight/2 * center_weight * neighbor_weights->GetValue(j) * 50;
                    b[vidx1*3+2] += 2 * weight/2 * center_weight * neighbor_weights->GetValue(j) * 50;
                    */
                    //std::cout<<vidx1<<" "<<2 * weight/2* center_weight * neighbor_weights->GetValue(j)<<endl;

                    for(int k = 0; k < neighbours->GetNumberOfIds(); k++)
                    {
                        int currentidx = neighbours->GetId(k);
                        if(!Is_Fixed[currentidx])
                        {
                            int currentvidx = InvertIds[currentidx];
                            assert(currentvidx > -0.5);

                            b[currentvidx*3] += -2 * weight/2* neighbor_weights->GetValue(k) * neighbor_weights->GetValue(j) * v[0];
                            b[currentvidx*3+1] += -2 * weight/2* neighbor_weights->GetValue(k) * neighbor_weights->GetValue(j) * v[1];
                            b[currentvidx*3+2] += -2 * weight/2* neighbor_weights->GetValue(k) * neighbor_weights->GetValue(j) * v[2];

                            //testb
                            /*
                            b[currentvidx*3] += -2 * weight/2* neighbor_weights->GetValue(k) * neighbor_weights->GetValue(j) * -50;
                            b[currentvidx*3+1] += -2 * weight/2* neighbor_weights->GetValue(k) * neighbor_weights->GetValue(j) * 50;
                            b[currentvidx*3+2] += -2 * weight/2* neighbor_weights->GetValue(k) * neighbor_weights->GetValue(j) * 50;
                            */
                            //std::cout<<currentvidx<<" "<<2 * weight/2* center_weight * neighbor_weights->GetValue(j)<<endl;
                        }
                    }
                }
            }
            for(int j = 0; j < neighbours->GetNumberOfIds(); j++)
            {
                for(int k = 0; k < neighbours->GetNumberOfIds(); k++)
                {
                    int idx1_ = neighbours->GetId(j);
                    int idx2_ = neighbours->GetId(k);
                    int vidx1_ = InvertIds[idx1_];
                    int vidx2_ = InvertIds[idx2_];
                    if(!Is_Fixed[idx1_] && !Is_Fixed[idx2_])
                    {
                        updateA(m, vidx1_, vidx2_, weight*neighbor_weights->GetValue(j)*neighbor_weights->GetValue(k), coefficientMap);
                        //std::cout<<vidx1_<<" "<<vidx2_<<" "<<weight*neighbor_weights->GetValue(j)*neighbor_weights->GetValue(k)<<endl;
                    }
                }
            }
        }
        else
        {
            vtkSmartPointer<vtkIdList> neighbours = vtkSmartPointer<vtkIdList>::New();
            neighbours = GetConnectedVertices(t_colon, idx1);
            vtkSmartPointer<vtkIdList> variable_neighbours = vtkSmartPointer<vtkIdList>::New();
            double weight = BENDWEIGHT / (SurroundAreas->GetValue(idx1));
            for(int j = 0; j < neighbours->GetNumberOfIds(); j++)
            {
                int currentidx = neighbours->GetId(j);
                if(!Is_Fixed[currentidx])
                    variable_neighbours->InsertNextId(currentidx);
            }
            if(variable_neighbours->GetNumberOfIds() != 0)
            {
                double center_weight = 0;
                vtkSmartPointer<vtkDoubleArray> neighbor_weights = vtkSmartPointer<vtkDoubleArray>::New();
                for(int j = 0; j < neighbours->GetNumberOfIds(); j++)
                {
                    int idx2 = neighbours->GetId(j);
                    double neighbor_weight = 0;
                    // compute the weight[j] for the j-th neighbour
                    vtkSmartPointer<vtkIdList> cellids = vtkSmartPointer<vtkIdList>::New();
                    cellids = GetFacetsOfEdge(t_colon, idx1, idx2);
                    assert(cellids->GetNumberOfIds() == 2); // we have already excluded the boundary edges, so this edge must be inside the surface

                    for(int k = 0; k < 2; k++)
                    {
                        int idx12;
                        vtkSmartPointer<vtkIdList> ptids = vtkSmartPointer<vtkIdList>::New();
                        t_colon->GetCellPoints(cellids->GetId(k), ptids);
                        for(int l = 0; l < ptids->GetNumberOfIds(); l++)
                        {
                            if(ptids->GetId(l)!=idx1 && ptids->GetId(l)!=idx2)
                            {
                                idx12 = ptids->GetId(l);
                                break;
                            }
                        }

                        double p1[3], p2[3], p12[3];
                        t_colon->GetPoint(idx1, p1);
                        t_colon->GetPoint(idx2, p2);
                        t_colon->GetPoint(idx12, p12);
                        double v1[3], v2[3];
                        vtkMath::Subtract(p1, p12, v1);
                        vtkMath::Subtract(p2, p12, v2);
                        double angle = vtkMath::AngleBetweenVectors(v1, v2);

                        neighbor_weight += cos(angle)/sin(angle);
                        //std::cout<<idx1<<" "<<idx2<<" "<<idx12<<" : "<<cos(angle)/sin(angle)<<endl;
                    }
                    //neighbor_weight = (neighbor_weight > 0.1)? neighbor_weight : 0.1;
                    //neighbor_weight = 1;
                    //std::cout<<idx2<<" "<<neighbor_weight<<" | ";
                    center_weight += neighbor_weight;
                    neighbor_weights->InsertNextValue(neighbor_weight);
                }

                for(int j = 0; j < neighbours->GetNumberOfIds(); j++)
                {
                    for(int k = 0; k < neighbours->GetNumberOfIds(); k++)
                    {
                        int idx1_ = neighbours->GetId(j);
                        int idx2_ = neighbours->GetId(k);
                        if(!Is_Fixed[idx1_] && !Is_Fixed[idx2_])
                        {
                            int vidx1_ = InvertIds[idx1_];
                            int vidx2_ = InvertIds[idx2_];
                            updateA(m, vidx1_, vidx2_, weight*neighbor_weights->GetValue(j)*neighbor_weights->GetValue(k), coefficientMap);
                        }
                    }
                }

                for(int j = 0; j < neighbours->GetNumberOfIds(); j++)
                {
                    int idx2 = neighbours->GetId(j);
                    if(!Is_Fixed[idx2])
                    {
                        double pnew[3], pold[3], v[3];
                        t_colon->GetPoint(idx1, pold);
                        SurfaceLineUp->GetPoint(idx1, pnew);
                        vtkMath::Subtract(pnew, pold, v);
                        int vidx2 = InvertIds[idx2];

                        b[vidx2*3] += 2*weight/2 *center_weight*neighbor_weights->GetValue(j) * v[0];
                        b[vidx2*3+1] += 2*weight/2 *center_weight*neighbor_weights->GetValue(j) * v[1];
                        b[vidx2*3+2] += 2*weight/2 *center_weight*neighbor_weights->GetValue(j) * v[2];

                        //testb
                        /*
                        b[vidx2*3] += 2*weight/2 *center_weight*neighbor_weights->GetValue(j) * -50;
                        b[vidx2*3+1] += 2*weight/2 *center_weight*neighbor_weights->GetValue(j) * 50;
                        b[vidx2*3+2] += 2*weight/2 *center_weight*neighbor_weights->GetValue(j) * 50;
                        */
                    }
                    else
                    {
                        for(int k = 0; k < neighbours->GetNumberOfIds(); k++)
                        {
                            int currentidx = neighbours->GetId(k);
                            if(!Is_Fixed[currentidx])
                            {
                                assert(currentidx != idx2);
                                double pnew[3], pold[3], vc[3];
                                t_colon->GetPoint(idx2, pold);
                                SurfaceLineUp->GetPoint(idx2, pnew);
                                vtkMath::Subtract(pnew, pold, vc);
                                int currentvidx = InvertIds[currentidx];

                                b[currentvidx*3] += -2*weight/2 *neighbor_weights->GetValue(j) *neighbor_weights->GetValue(k) * vc[0];
                                b[currentvidx*3+1] += -2*weight/2 *neighbor_weights->GetValue(j) *neighbor_weights->GetValue(k) * vc[1];
                                b[currentvidx*3+2] += -2*weight/2 *neighbor_weights->GetValue(j) *neighbor_weights->GetValue(k) * vc[2];

                                //testb
                                /*
                                b[currentvidx*3] += -2*weight/2 *neighbor_weights->GetValue(j) *neighbor_weights->GetValue(k) * -50;
                                b[currentvidx*3+1] += -2*weight/2 *neighbor_weights->GetValue(j) *neighbor_weights->GetValue(k) * 50;
                                b[currentvidx*3+2] += -2*weight/2 *neighbor_weights->GetValue(j) *neighbor_weights->GetValue(k) * 50;
                                */
                            }
                        }
                    }
                }
            }
        }
    }

    //VisualizePoints(boundarypoints, 1, 0, 0, 5, t_rendermanager);

    free(marked);
    free(boundary);
}

vtkSmartPointer<vtkIdList> GetConnectedVertices(vtkSmartPointer<vtkPolyData> mesh, vtkIdType id)
{
    vtkSmartPointer<vtkIdList> connectedVertices =
            vtkSmartPointer<vtkIdList>::New();

    //get all cells that vertex 'id' is a part of
    vtkSmartPointer<vtkIdList> cellIdList =
            vtkSmartPointer<vtkIdList>::New();
    mesh->GetPointCells(id, cellIdList);

    for(vtkIdType i = 0; i < cellIdList->GetNumberOfIds(); i++)
    {
        vtkSmartPointer<vtkIdList> pointIdList =
                vtkSmartPointer<vtkIdList>::New();
        mesh->GetCellPoints(cellIdList->GetId(i), pointIdList);

        for(vtkIdType j = 0; j < pointIdList->GetNumberOfIds(); j++)
        {
            if(pointIdList->GetId(j) != id)
            {
                bool marked = false;
                for(int k = 0; k < connectedVertices->GetNumberOfIds(); k++)
                {
                    if(connectedVertices->GetId(k) == pointIdList->GetId(j))
                    {
                        marked = true;
                        break;
                    }
                }
                if(!marked)
                    connectedVertices->InsertNextId(pointIdList->GetId(j));
            }
        }
    }
    return connectedVertices;
}
double ComputeFacetArea(vtkSmartPointer<vtkPolyData> t_colon, vtkIdType cellid)
{
    vtkSmartPointer<vtkIdList> ptids = vtkSmartPointer<vtkIdList>::New();
    t_colon->GetCellPoints(cellid, ptids);
    assert(ptids->GetNumberOfIds() == 3);
    double p0[3], p1[3], p2[3];
    t_colon->GetPoint(ptids->GetId(0), p0);
    t_colon->GetPoint(ptids->GetId(1), p1);
    t_colon->GetPoint(ptids->GetId(2), p2);

    double a = sqrt(vtkMath::Distance2BetweenPoints(p0, p1));
    double b = sqrt(vtkMath::Distance2BetweenPoints(p0, p2));
    double c = sqrt(vtkMath::Distance2BetweenPoints(p1, p2));
    double s = (a+b+c)/2;
    return sqrt(s*(s-a)*(s-b)*(s-c));
}
