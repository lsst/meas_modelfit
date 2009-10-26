#include <ndarray/eigen.hpp>

#include <Eigen/Array>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE ndarray-eigen
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(MatrixXd) {
    Eigen::MatrixXd m1 = Eigen::MatrixXd::Random(5,6);
    ndarray::Array<double,2> a1(ndarray::viewMatrixAsArray(m1));
    BOOST_CHECK_EQUAL(m1.data(),a1.getData());
    BOOST_CHECK_EQUAL(m1.rows(),a1.getSize<0>());
    BOOST_CHECK_EQUAL(m1.cols(),a1.getSize<1>());
    ndarray::Array<double,2,2> a2(ndarray::copy(a1));
    Eigen::MatrixXd m2 = ndarray::viewArrayAs<Eigen::MatrixXd>(a2);
    for (int i=0; i<m1.rows(); ++i) {
        for (int j=0; j<m1.cols(); ++j) {
            BOOST_CHECK_EQUAL(m1(i,j),a1[i][j]);
            BOOST_CHECK_EQUAL(m2(i,j),a1[i][j]);
        }
    }
    Eigen::Block<Eigen::MatrixXd> m3 = m1.block(0,1,3,3);
    ndarray::Array<double,2> a3(ndarray::viewMatrixAsArray(m3));
    BOOST_CHECK_EQUAL(m3.data(),a3.getData());
    BOOST_CHECK_EQUAL(m3.rows(),a3.getSize<0>());
    BOOST_CHECK_EQUAL(m3.cols(),a3.getSize<1>());
    ndarray::Array<double,2,2> a4(ndarray::copy(a3));
    Eigen::MatrixXd m4 = Eigen::MatrixXd::Zero(5,6);
    m4.block(0,1,3,3) = ndarray::viewArrayAs< Eigen::Block<Eigen::MatrixXd> >(a4);
    for (int i=0; i<m3.rows(); ++i) {
        for (int j=0; j<m3.cols(); ++j) {
            BOOST_CHECK_EQUAL(m3(i,j),a3[i][j]);
            BOOST_CHECK_EQUAL(m4.block(0,1,3,3)(i,j),a3[i][j]);
        }
    }

}

BOOST_AUTO_TEST_CASE(Matrix3d) {
    Eigen::Matrix3d m1 = Eigen::Matrix3d::Random();
    ndarray::Array<double,2> a1(ndarray::viewMatrixAsArray(m1));
    BOOST_CHECK_EQUAL(m1.data(),a1.getData());
    BOOST_CHECK_EQUAL(m1.rows(),a1.getSize<0>());
    BOOST_CHECK_EQUAL(m1.cols(),a1.getSize<1>());
    ndarray::Array<double,2,2> a2(ndarray::copy(a1));
    Eigen::Matrix3d m2 = ndarray::viewArrayAs<Eigen::Matrix3d>(a2);
    for (int i=0; i<m1.rows(); ++i) {
        for (int j=0; j<m1.cols(); ++j) {
            BOOST_CHECK_EQUAL(m1(i,j),a1[i][j]);
            BOOST_CHECK_EQUAL(m2(i,j),a1[i][j]);
        }
    }
    {
        Eigen::Block<Eigen::Matrix3d,2,2> m3 = m1.block<2,2>(0,0);
        ndarray::Array<double,2> a3(ndarray::viewMatrixAsArray(m3));
        ndarray::Array<double,2,2> a4(ndarray::copy(a3));
        BOOST_CHECK_EQUAL(m3.data(),a3.getData());
        BOOST_CHECK_EQUAL(m3.rows(),a3.getSize<0>());
        BOOST_CHECK_EQUAL(m3.cols(),a3.getSize<1>());
        Eigen::Matrix3d m4 = Eigen::Matrix3d::Zero();
        Eigen::Block<Eigen::Matrix3d,2,2> m5 = m4.block<2,2>(0,0);
        m5 = ndarray::viewArrayAs< Eigen::Block<Eigen::Matrix3d,2,2> >(a4);
        for (int i=0; i<m3.rows(); ++i) {
            for (int j=0; j<m3.cols(); ++j) {
                BOOST_CHECK_EQUAL(m3(i,j),a3[i][j]);
                BOOST_CHECK_EQUAL(m5(i,j),a3[i][j]);
            }
        }
    }
    {
        Eigen::Block<Eigen::Matrix3d,2,1> m3 = m1.block<2,1>(0,2);
        ndarray::Array<double,1> a3(ndarray::viewVectorAsArray(m3));
        ndarray::Array<double,1,1> a4(ndarray::copy(a3));
        BOOST_CHECK_EQUAL(m3.data(),a3.getData());
        BOOST_CHECK_EQUAL(m3.rows(),a3.getSize<0>());
        BOOST_CHECK_EQUAL(m3.cols(),a3.getSize<1>());
        Eigen::Matrix3d m4 = Eigen::Matrix3d::Zero();
        Eigen::Block<Eigen::Matrix3d,2,1> m5 = m4.block<2,1>(0,2);
        m5 = ndarray::viewArrayAs< Eigen::Block<Eigen::Matrix3d,2,1> >(a4);
        for (int i=0; i<m3.size(); ++i) {
                BOOST_CHECK_EQUAL(m3(i),a3[i]);
                BOOST_CHECK_EQUAL(m5(i),a3[i]);
        }
    }
    
}

BOOST_AUTO_TEST_CASE(VectorXd) {
    Eigen::VectorXd m1 = Eigen::VectorXd::Random(5);
    ndarray::Array<double,2> a1(ndarray::viewMatrixAsArray(m1));
    BOOST_CHECK_EQUAL(m1.data(),a1.getData());
    BOOST_CHECK_EQUAL(m1.rows(),a1.getSize<0>());
    BOOST_CHECK_EQUAL(m1.cols(),a1.getSize<1>());
    ndarray::Array<double,2,2> a2(ndarray::copy(a1));
    Eigen::VectorXd m2 = ndarray::viewArrayAs<Eigen::VectorXd>(a2);
    for (int i=0; i<m1.rows(); ++i) {
        for (int j=0; j<m1.cols(); ++j) {
            BOOST_CHECK_EQUAL(m1(i,j),a1[i][j]);
            BOOST_CHECK_EQUAL(m2(i,j),a1[i][j]);
        }
    }
    ndarray::Array<double,1> a3(ndarray::viewVectorAsArray(m1));
    BOOST_CHECK_EQUAL(m1.data(),a1.getData());
    BOOST_CHECK_EQUAL(m1.rows(),a1.getSize<0>());
    for (int i=0; i<m1.rows(); ++i) {
        BOOST_CHECK_EQUAL(m1[i],a3[i]);
    }
}

BOOST_AUTO_TEST_CASE(Vector3d) {
    Eigen::Vector3d m1 = Eigen::Vector3d::Random();
    ndarray::Array<double,2> a1(ndarray::viewMatrixAsArray(m1));
    BOOST_CHECK_EQUAL(m1.data(),a1.getData());
    BOOST_CHECK_EQUAL(m1.rows(),a1.getSize<0>());
    BOOST_CHECK_EQUAL(m1.cols(),a1.getSize<1>());
    ndarray::Array<double,2,2> a2(ndarray::copy(a1));
    Eigen::Vector3d m2 = ndarray::viewArrayAs<Eigen::Vector3d>(a2);
    for (int i=0; i<m1.rows(); ++i) {
        for (int j=0; j<m1.cols(); ++j) {
            BOOST_CHECK_EQUAL(m1(i,j),a1[i][j]);
            BOOST_CHECK_EQUAL(m2(i,j),a1[i][j]);
        }
    }
    ndarray::Array<double,1> a3(ndarray::viewVectorAsArray(m1));
    BOOST_CHECK_EQUAL(m1.data(),a1.getData());
    BOOST_CHECK_EQUAL(m1.rows(),a1.getSize<0>());
    for (int i=0; i<m1.rows(); ++i) {
        BOOST_CHECK_EQUAL(m1[i],a3[i]);
    }
}
