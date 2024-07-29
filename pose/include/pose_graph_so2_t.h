#include <iostream>
#include <fstream>
#include <string>
#include <Eigen/Core>
#include <g2o/core/base_vertex.h>
#include <g2o/core/base_binary_edge.h>
#include <g2o/core/base_unary_edge.h>
#include <g2o/core/block_solver.h>
#include <g2o/core/optimization_algorithm_levenberg.h>
#include <g2o/core/optimization_algorithm_gauss_newton.h>
#include <g2o/core/optimization_algorithm_dogleg.h>
#include <g2o/solvers/eigen/linear_solver_eigen.h>
#include <g2o/solvers/dense/linear_solver_dense.h>
// #include <g2o/solvers/cholmod/linear_solver_cholmod.h>
// #include <g2o/solvers/csparse/linear_solver_csparse.h>
#include <g2o/solvers/pcg/linear_solver_pcg.h>
#include <sophus/so2.hpp>

using namespace std;
using namespace Eigen;
using Sophus::SO2d;

namespace mygroup
{
    class SO2_t
    {
    public:
        SO2_t()
        {
            SO2_ = SO2d(0.0);
            t_ = Vector2d::Zero();
        }

        SO2_t(Matrix2d r, Vector2d t)
        {
            SO2_ = SO2d(r);
            t_ = t;
        }

        SO2_t(double theta, Vector2d t)
        {
            SO2_ = SO2d(theta);
            t_ = t;
        }

        SO2_t(SO2d g, Vector2d t)
        {
            SO2_ = g;
            t_ = t;
        }

    public:
        Vector2d t_;
        SO2d SO2_;
    };
}

using mygroup::SO2_t;


class MyVertex : public g2o::BaseVertex<3, SO2_t>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    virtual bool read(istream &is) override
    {
    }

    virtual bool write(ostream &os) const override
    {
    }

    virtual void setToOriginImpl() override
    {
        _estimate = SO2_t();
    }

    virtual void oplusImpl(const double *update) override
    {
//        Vector2d t_update(update[0], update[1]);
//        SO2d R_update = SO2d::exp(update[2]);
//        _estimate = SO2_t(_estimate.SO2_ * R_update, _estimate.t_ + _estimate.SO2_.matrix() * t_update);


        Vector2d t_update(update[0], update[1]);
        Vector2d t = _estimate.t_ + t_update;

        SO2d R = _estimate.SO2_ * SO2d::exp(update[2]);

        _estimate = SO2_t(R, t);
    }
};


class MyEdge : public g2o::BaseBinaryEdge<3, SO2_t, MyVertex, MyVertex>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    virtual bool read(istream &is) override
    {
    }

    virtual bool write(ostream &os) const override
    {
    }

    virtual void computeError() override
    {
        SO2_t v1 = (dynamic_cast<MyVertex *>(_vertices[0]))->estimate();
        SO2_t v2 = (dynamic_cast<MyVertex *>(_vertices[1]))->estimate();

        Matrix2d rz = _measurement.SO2_.matrix();
        Vector2d tz = _measurement.t_;

        Matrix2d r1 = v1.SO2_.matrix();
        Vector2d t1 = v1.t_;

        Matrix2d r2 = v2.SO2_.matrix();
        Vector2d t2 = v2.t_;

        Matrix2d re = rz.transpose() * r1.transpose() * r2;
        Vector2d te = -rz.transpose() * tz - rz.transpose() * r1.transpose() * t1 + rz.transpose() * r1.transpose() * t2;

        double phie = SO2d(re).log();
        _error << te, phie;
    }

   virtual void linearizeOplus() override
   {
       SO2_t v1 = (dynamic_cast<MyVertex *>(_vertices[0]))->estimate();
       SO2_t v2 = (dynamic_cast<MyVertex *>(_vertices[1]))->estimate();

       Matrix2d rz = _measurement.SO2_.matrix();
       Vector2d tz = _measurement.t_;

       Matrix2d r1 = v1.SO2_.matrix();
       Vector2d t1 = v1.t_;

       Matrix2d r2 = v2.SO2_.matrix();
       Vector2d t2 = v2.t_;

       // Eigen::Matrix2d dRiT_dtheta;       //  derivative of Ri^T over theta
       // dRiT_dtheta(0, 0) = -1 * r1(1, 0); //  cosX -> -sinX
       // dRiT_dtheta(0, 1) = 1 * r1(0, 0);  //  sinX ->  cosX
       // dRiT_dtheta(1, 0) = -1 * r1(0, 0); // -sinX -> -cosX
       // dRiT_dtheta(1, 1) = -1 * r1(1, 0); //  cosX -> -sinX

       _jacobianOplusXi = Matrix3d::Zero();
       _jacobianOplusXj = Matrix3d::Zero();

       _jacobianOplusXi.block(0, 0, 2, 2) = -rz.transpose() * r1.transpose();
       _jacobianOplusXi.block(0, 2, 2, 1) = rz.transpose() * r1.transpose() * SO2d::hat(1) * (t1 - t2);
       // _jacobianOplusXi.block(0, 2, 2, 1) = rz.transpose() * dRiT_dtheta * (t2 - t1); // 两种方式等价
       _jacobianOplusXi(2, 2) = -1.0;

       _jacobianOplusXj.block(0, 0, 2, 2) = rz.transpose() * r1.transpose();
       _jacobianOplusXj(2, 2) = 1.0;

    //    using std::cout;
    //    using std::endl;

    //    static bool once = true;
    //    static int cnt = 0;
    //    if (once)
    //    {    
    //        cout << ++cnt << endl;
    //        cout << "Analysis _jacobianOplusXi" << endl;
    //        cout << _jacobianOplusXi << endl;
    //        cout << "Analysis _jacobianOplusXj" << endl;
    //        cout << _jacobianOplusXj << endl;

    //        if (cnt > 11)
    //            once = false;
    //    }
   }
};

class PriorEdge : public g2o::BaseUnaryEdge<3, SO2_t, MyVertex>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    virtual bool read(istream &is) override
    {
    }

    virtual bool write(ostream &os) const override
    {
    }

    virtual void computeError() override
    {
        SO2_t v1 = (static_cast<MyVertex *>(_vertices[0]))->estimate();

        Matrix2d rz = _measurement.SO2_.matrix();
        Vector2d tz = _measurement.t_;

        Matrix2d r1 = v1.SO2_.matrix();
        Vector2d t1 = v1.t_;

        Matrix2d re = rz.transpose() * r1;
        Vector2d te = rz.transpose() * (t1 - tz);

        double phie = SO2d(re).log();
        _error << te, phie;
    }

    virtual void linearizeOplus() override
    {
        SO2_t v1 = (static_cast<MyVertex *>(_vertices[0]))->estimate();

        Matrix2d rz = _measurement.SO2_.matrix();
        Vector2d tz = _measurement.t_;

        Matrix2d r1 = v1.SO2_.matrix();
        Vector2d t1 = v1.t_;

        _jacobianOplusXi = Matrix3d::Zero();
        _jacobianOplusXi.block(0, 0, 2, 2) = rz.transpose();
        _jacobianOplusXi(2, 2) = 1.0;
    }
};
