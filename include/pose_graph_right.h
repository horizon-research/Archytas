#include <iostream>
#include <fstream>
#include <string>
#include <Eigen/Core>

#include <g2o/core/base_vertex.h>
#include <g2o/core/base_binary_edge.h>
#include <g2o/core/block_solver.h>
#include <g2o/core/optimization_algorithm_levenberg.h>
#include <g2o/core/optimization_algorithm_gauss_newton.h>
#include <g2o/core/optimization_algorithm_dogleg.h>
#include <g2o/solvers/eigen/linear_solver_eigen.h>
#include <g2o/solvers/dense/linear_solver_dense.h>
// #include <g2o/solvers/cholmod/linear_solver_cholmod.h>
// #include <g2o/solvers/csparse/linear_solver_csparse.h>
#include <g2o/solvers/pcg/linear_solver_pcg.h>

#include <sophus/se3.hpp>

using namespace std;
using namespace Eigen;
using Sophus::SE3d;
using Sophus::SO3d;



typedef Matrix<double, 6, 6> Matrix6d;


Matrix6d JRInv(const SE3d &e)
{
    Matrix6d J;
    J.block(0, 0, 3, 3) = SO3d::hat(e.so3().log());
    J.block(0, 3, 3, 3) = SO3d::hat(e.translation());
    J.block(3, 0, 3, 3) = Matrix3d::Zero(3, 3);
    J.block(3, 3, 3, 3) = SO3d::hat(e.so3().log());
    // J = J * 0.5 + Matrix6d::Identity();
    J = Matrix6d::Identity(); // try Identity if you want
    return J;
}


typedef Matrix<double, 6, 1> Vector6d;

class VertexSE3LieAlgebra : public g2o::BaseVertex<6, SE3d>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    virtual bool read(istream &is) override
    {
        double data[7];
        for (int i = 0; i < 7; i++)
            is >> data[i];
        setEstimate(SE3d(
            Quaterniond(data[6], data[3], data[4], data[5]),
            Vector3d(data[0], data[1], data[2])));
    }

    virtual bool write(ostream &os) const override
    {
        // os << id() << " ";
        // Quaterniond q = _estimate.unit_quaternion();
        // os << _estimate.translation().transpose() << " ";
        // os << q.coeffs()[0] << " " << q.coeffs()[1] << " " << q.coeffs()[2] << " " << q.coeffs()[3] << endl;
        // return true;

        Matrix4d pm = _estimate.matrix();
        os << pm(0, 0) << " " << pm(0, 1) << " " << pm(0, 2) << " " << pm(0, 3) << " " << pm(1, 0) << " " << pm(1, 1) << " " << pm(1, 2) << " " << pm(1, 3) << " " << pm(2, 0) << " " << pm(2, 1) << " " << pm(2, 2) << " " << pm(2, 3) << endl;
        return true;
    }

    virtual void setToOriginImpl() override
    {
        _estimate = SE3d();
    }

 
    virtual void oplusImpl(const double *update) override
    {
        Vector6d upd;
        upd << update[0], update[1], update[2], update[3], update[4], update[5];
        _estimate = _estimate * SE3d::exp(upd);
    }
};


class EdgeSE3LieAlgebra : public g2o::BaseBinaryEdge<6, SE3d, VertexSE3LieAlgebra, VertexSE3LieAlgebra>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    virtual bool read(istream &is) override
    {
        double data[7];
        for (int i = 0; i < 7; i++)
            is >> data[i];
        Quaterniond q(data[6], data[3], data[4], data[5]);
        q.normalize();
        setMeasurement(SE3d(q, Vector3d(data[0], data[1], data[2])));
        for (int i = 0; i < information().rows() && is.good(); i++)
            for (int j = i; j < information().cols() && is.good(); j++)
            {
                is >> information()(i, j);
                if (i != j)
                    information()(j, i) = information()(i, j);
            }
        return true;
    }

    virtual bool write(ostream &os) const override
    {
        VertexSE3LieAlgebra *v1 = static_cast<VertexSE3LieAlgebra *>(_vertices[0]);
        VertexSE3LieAlgebra *v2 = static_cast<VertexSE3LieAlgebra *>(_vertices[1]);
        os << v1->id() << " " << v2->id() << " ";
        SE3d m = _measurement;
        Eigen::Quaterniond q = m.unit_quaternion();
        os << m.translation().transpose() << " ";
        os << q.coeffs()[0] << " " << q.coeffs()[1] << " " << q.coeffs()[2] << " " << q.coeffs()[3] << " ";

        // information matrix
        for (int i = 0; i < information().rows(); i++)
            for (int j = i; j < information().cols(); j++)
            {
                os << information()(i, j) << " ";
            }
        os << endl;
        return true;
    }

 
    virtual void computeError() override
    {
        SE3d v1 = (static_cast<VertexSE3LieAlgebra *>(_vertices[0]))->estimate();
        SE3d v2 = (static_cast<VertexSE3LieAlgebra *>(_vertices[1]))->estimate();
        _error = (_measurement.inverse() * v1.inverse() * v2).log();
    }

 
    virtual void linearizeOplus() override
    {
        SE3d v1 = (static_cast<VertexSE3LieAlgebra *>(_vertices[0]))->estimate();
        SE3d v2 = (static_cast<VertexSE3LieAlgebra *>(_vertices[1]))->estimate();
        Matrix6d J = JRInv(SE3d::exp(_error));

        // _jacobianOplusXi = -J * v2.inverse().Adj() * v1.Adj();
        _jacobianOplusXi = -J * (v2.inverse() * v1).Adj();
        _jacobianOplusXj = J;
    }
};