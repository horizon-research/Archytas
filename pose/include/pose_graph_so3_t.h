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
#include <sophus/so3.hpp>
#include <sophus/se3.hpp>


using namespace std;
using namespace Eigen;
using Sophus::SO3d;
typedef Matrix<double, 6, 1> Vector6d;
typedef Matrix<double, 6, 6> Matrix6d;
typedef Matrix<double, 4, 4> Matrix4d;
using Sophus::SE3d;

namespace mygroup
{
    class SO3_t
    {
    public:
        SO3_t()
        {
            SO3_ = SO3d(Matrix3d::Identity());
            t_ = Vector3d::Zero();
        }

        SO3_t(Matrix3d r, Vector3d t)
        {
            SO3_ = SO3d(r);
            t_ = t;
        }

        SO3_t(AngleAxisd phi, Vector3d t)
        {
            SO3_ = SO3d(phi.matrix());
            t_ = t;
        }

        SO3_t(Quaterniond q, Vector3d t)
        {
            SO3_ = SO3d(q);
            t_ = t;
        }

        SO3_t(SO3d g, Vector3d t)
        {
            SO3_ = g;
            t_ = t;
        }

    public:
        Vector3d t_;
        SO3d SO3_;
    };
}

Matrix3d JRInv(const SO3d &e)
{
    Matrix3d J;
    J = SO3d::hat(e.log());
    // J = J * 0.5 + Matrix3d::Identity();
    J = Matrix3d::Identity(); // try Identity if you want
    return J;
}

using mygroup::SO3_t;


class MyVertex : public g2o::BaseVertex<6, SO3_t>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    virtual bool read(istream &is) override
    {
        double data[7];
        for (int i = 0; i < 7; i++)
            is >> data[i];

        setEstimate(SO3_t(Quaterniond(data[6], data[3], data[4], data[5]), Vector3d(data[0], data[1], data[2])));
    }

    virtual bool write(ostream &os) const override
    {
        // os << id() << " ";
        // Quaterniond q(_estimate.SO3_.matrix());
        // os << _estimate.t_.transpose() << " ";
        // os << q.coeffs()[0] << " " << q.coeffs()[1] << " " << q.coeffs()[2] << " " << q.coeffs()[3] << endl;
        // return true;
        SE3d pose(_estimate.SO3_.matrix(), _estimate.t_);
        Matrix4d pm = pose.matrix();
        os << pm(0, 0) << " " << pm(0, 1) << " " << pm(0, 2) << " " << pm(0, 3) << " " << pm(1, 0) << " " << pm(1, 1) << " " << pm(1, 2) << " " << pm(1, 3) << " " << pm(2, 0) << " " << pm(2, 1) << " " << pm(2, 2) << " " << pm(2, 3) << endl;
        return true;

    }

    virtual void setToOriginImpl() override
    {
        _estimate = SO3_t();
    }

    virtual void oplusImpl(const double *update) override
    {
        Vector3d t_update(update[0], update[1], update[2]);
        Vector3d phi_update(update[3], update[4], update[5]);

        // SO3d R_update = SO3d::exp(phi_update);
        // _estimate = SO3_t(_estimate.SO3_ * R_update, _estimate.t_ + _estimate.SO3_.matrix() * t_update);

        Vector3d t = _estimate.t_ + t_update;

        SO3d R = _estimate.SO3_ * SO3d::exp(phi_update);

        _estimate = SO3_t(R, t);
    }
};

class MyEdge : public g2o::BaseBinaryEdge<6, SO3_t, MyVertex, MyVertex>
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
        setMeasurement(SO3_t(q, Vector3d(data[0], data[1], data[2])));
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
        MyVertex *v1 = static_cast<MyVertex *>(_vertices[0]);
        MyVertex *v2 = static_cast<MyVertex *>(_vertices[1]);
        os << v1->id() << " " << v2->id() << " ";
        SO3_t m = _measurement;
        Eigen::Quaterniond q(m.SO3_.matrix());
        os << m.t_.transpose() << " ";
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
        SO3_t v1 = (dynamic_cast<MyVertex *>(_vertices[0]))->estimate();
        SO3_t v2 = (dynamic_cast<MyVertex *>(_vertices[1]))->estimate();

        Matrix3d rz = _measurement.SO3_.matrix();
        Vector3d tz = _measurement.t_;

        Matrix3d r1 = v1.SO3_.matrix();
        Vector3d t1 = v1.t_;

        Matrix3d r2 = v2.SO3_.matrix();
        Vector3d t2 = v2.t_;

        Matrix3d re = rz.transpose() * r1.transpose() * r2;
        // Vector3d te = -rz.transpose() * tz - rz.transpose() * r1.transpose() * t1 + rz.transpose() * r1.transpose() * t2;
        // Vector3d te = rz.transpose() * (r1.transpose() * (t2 - t1) - tz);
        Vector3d te = r1.transpose() * (t2 - t1) - tz;

        Vector3d phie = SO3d(re).log();
        _error << te, phie;
    }

    virtual void linearizeOplus() override
    {
        SO3_t v1 = (dynamic_cast<MyVertex *>(_vertices[0]))->estimate();
        SO3_t v2 = (dynamic_cast<MyVertex *>(_vertices[1]))->estimate();

        Matrix3d rz = _measurement.SO3_.matrix();
        Vector3d tz = _measurement.t_;

        Matrix3d r1 = v1.SO3_.matrix();
        Vector3d t1 = v1.t_;

        Matrix3d r2 = v2.SO3_.matrix();
        Vector3d t2 = v2.t_;

        _jacobianOplusXi = Matrix6d::Zero();
        _jacobianOplusXj = Matrix6d::Zero();

        // Matrix3d J = JRInv(SO3d(rz.transpose() * r1.transpose() * r2));

        // _jacobianOplusXi.block(0, 0, 3, 3) = -rz.transpose() * r1.transpose();
        _jacobianOplusXi.block(0, 0, 3, 3) = -r1.transpose();
        // _jacobianOplusXi.block(0, 3, 3, 3) = rz.transpose() * r1.transpose() * SO3d::hat((t2 - t1)) * r1;
        _jacobianOplusXi.block(0, 3, 3, 3) = r1.transpose() * SO3d::hat((t2 - t1)) * r1;
        _jacobianOplusXi.block(3, 3, 3, 3) = -r2.transpose() * r1;

        // _jacobianOplusXj.block(0, 0, 3, 3) = rz.transpose() * r1.transpose();
        _jacobianOplusXj.block(0, 0, 3, 3) = r1.transpose();
        _jacobianOplusXj.block(3, 3, 3, 3) = Matrix3d::Identity();

// #define DEBUG
#ifdef DEBUG
        using std::cout;
        using std::endl;

        static bool once = true;
        static int cnt = 0;
        if (once)
        {
            cout << ++cnt << endl;
            cout << "Analysis _jacobianOplusXi" << endl;
            cout << _jacobianOplusXi << endl;
            cout << "Analysis _jacobianOplusXj" << endl;
            cout << _jacobianOplusXj << endl;

            if (cnt > 11)
                once = false;
        }
#endif
    }
};
