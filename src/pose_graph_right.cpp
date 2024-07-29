#include <pose_graph_right.h>

int main(int argc, char **argv) {
    if (argc != 2) {
        cout << "Usage: pose_graph_g2o_SE3_lie sphere.g2o" << endl;
        return 1;
    }
    ifstream fin(argv[1]);
    if (!fin) {
        cout << "file " << argv[1] << " does not exist." << endl;
        return 1;
    }

 
    typedef g2o::BlockSolver<g2o::BlockSolverTraits<6, 6>> BlockSolverType;
    typedef g2o::LinearSolverEigen<BlockSolverType::PoseMatrixType> LinearSolverType;
    // auto solver = new g2o::OptimizationAlgorithmLevenberg(g2o::make_unique<BlockSolverType>(g2o::make_unique<LinearSolverType>()));
    auto solver = new g2o::OptimizationAlgorithmGaussNewton(g2o::make_unique<BlockSolverType>(g2o::make_unique<LinearSolverType>()));
    g2o::SparseOptimizer optimizer;   
    optimizer.setAlgorithm(solver);   
    optimizer.setVerbose(true);       

    int vertexCnt = 0, edgeCnt = 0; 

    vector<VertexSE3LieAlgebra *> vectices;
    vector<EdgeSE3LieAlgebra *> edges;
    while (!fin.eof()) {
        string name;
        fin >> name;
        if (name == "VERTEX_SE3:QUAT") {
        
            VertexSE3LieAlgebra *v = new VertexSE3LieAlgebra();
            int index = 0;
            fin >> index;
            v->setId(index);
            v->read(fin);
            optimizer.addVertex(v);
            vertexCnt++;
            vectices.push_back(v);
            if (index == 0)
                v->setFixed(true);
        } else if (name == "EDGE_SE3:QUAT") {
         
            EdgeSE3LieAlgebra *e = new EdgeSE3LieAlgebra();
            int idx1, idx2;    
            fin >> idx1 >> idx2;
            e->setId(edgeCnt++);
            e->setVertex(0, optimizer.vertices()[idx1]);
            e->setVertex(1, optimizer.vertices()[idx2]);
            e->read(fin);
            optimizer.addEdge(e);
            edges.push_back(e);
        }
        if (!fin.good()) break;
    }

    cout << "read total " << vertexCnt << " vertices, " << edgeCnt << " edges." << endl;
    cout << "initial chi2 is " << optimizer.activeChi2() << endl;

    cout << "optimizing ..." << endl;
    optimizer.initializeOptimization();
    optimizer.optimize(30);

    cout << "saving optimization results ..." << endl;


    ofstream fout("evo_right.txt");
    for (VertexSE3LieAlgebra *v:vectices) {
        // fout << "VERTEX_SE3:QUAT ";
        v->write(fout);
    }
    // for (EdgeSE3LieAlgebra *e:edges) {
    //     fout << "EDGE_SE3:QUAT ";
    //     e->write(fout);
    // }
    fout.close();
    return 0;
}
