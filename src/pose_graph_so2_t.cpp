#include <pose_graph_so2_t.h>

int main()
{
    typedef g2o::BlockSolver<g2o::BlockSolverTraits<3, 3>> BlockSolverType;
    typedef g2o::LinearSolverEigen<BlockSolverType::PoseMatrixType> LinearSolverType;
    // typedef g2o::LinearSolverDense<BlockSolverType::PoseMatrixType> LinearSolverType; 
    // typedef g2o::LinearSolverPCG<BlockSolverType::PoseMatrixType> LinearSolverType; 
    auto solver = new g2o::OptimizationAlgorithmLevenberg(g2o::make_unique<BlockSolverType>(g2o::make_unique<LinearSolverType>()));
    // auto solver = new g2o::OptimizationAlgorithmGaussNewton(g2o::make_unique<BlockSolverType>(g2o::make_unique<LinearSolverType>()));
    // auto solver = new g2o::OptimizationAlgorithmDogleg(g2o::make_unique<BlockSolverType>(g2o::make_unique<LinearSolverType>()));
    g2o::SparseOptimizer optimizer; 
    optimizer.setAlgorithm(solver); 
    optimizer.setVerbose(true);     

    MyVertex *v0, *v1, *v2;
    v0 = new MyVertex();
    v1 = new MyVertex();
    v2 = new MyVertex();
    SO2_t v0_initial(0.2, Vector2d(0.5, 0.4));
    SO2_t v1_initial(-0.4, Vector2d(2.3, 0.7));
    SO2_t v2_initial(0.5, Vector2d(4.1, -0.8));

    v0->setId(0);
    v0->setEstimate(v0_initial);
    v1->setId(1);
    v1->setEstimate(v1_initial);
    v2->setId(2);
    v2->setEstimate(v2_initial);
    optimizer.addVertex(v0);
    optimizer.addVertex(v1);
    optimizer.addVertex(v2);

    SO2_t measurement(0.0, Vector2d(2.0, 0.0));

    Matrix3d odom_cov = Matrix3d::Identity();
    odom_cov(0, 0) = 0.2;
    odom_cov(1, 1) = 0.2;
    odom_cov(2, 2) = 0.1;

    Matrix3d prior_cov = Matrix3d::Identity();
    prior_cov(0, 0) = 0.3;
    prior_cov(1, 1) = 0.3;
    prior_cov(2, 2) = 0.1;

    PriorEdge *e0;
    e0 = new PriorEdge();
    e0->setId(0);
    e0->setVertex(0, optimizer.vertices()[0]);
    e0->setMeasurement(SO2_t());
    e0->setInformation(prior_cov.inverse());

    MyEdge *e1, *e2;
    e1 = new MyEdge();
    e2 = new MyEdge();
    e1->setId(1);
    e1->setVertex(0, optimizer.vertices()[0]);
    e1->setVertex(1, optimizer.vertices()[1]);
    e1->setMeasurement(measurement);
    e1->setInformation(odom_cov.inverse());

    e2->setId(2);
    e2->setVertex(0, optimizer.vertices()[1]);
    e2->setVertex(1, optimizer.vertices()[2]);
    e2->setMeasurement(measurement);
    e2->setInformation(odom_cov.inverse());

    optimizer.addEdge(e0);
    optimizer.addEdge(e1);
    optimizer.addEdge(e2);

    cout << "optimizing ..." << endl;
    optimizer.initializeOptimization();
    optimizer.optimize(10);

    SO2_t opt_value_0 = v0->estimate();
    SO2_t opt_value_1 = v1->estimate();
    SO2_t opt_value_2 = v2->estimate();

    std::cout << "vertex 0 = " << opt_value_0.t_.transpose() << " " << opt_value_0.SO2_.log() << std::endl;
    std::cout << "vertex 1 = " << opt_value_1.t_.transpose() << " " << opt_value_1.SO2_.log() << std::endl;
    std::cout << "vertex 2 = " << opt_value_2.t_.transpose() << " " << opt_value_2.SO2_.log() << std::endl;

    return 0;
}