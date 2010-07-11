/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
 * 
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the LSST License Statement and 
 * the GNU General Public License along with this program.  If not, 
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */
 
#include <ndarray.hpp>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE ndarray
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(vectors) {
    ndarray::Vector<int,3> a = ndarray::makeVector(5,6,7);
    BOOST_CHECK_EQUAL(a[0],5);
    BOOST_CHECK_EQUAL(a[1],6);
    BOOST_CHECK_EQUAL(a[2],7);
    ndarray::Vector<int,3> b(a);
    BOOST_CHECK_EQUAL(a,b);
    ndarray::Vector<double,3> c(a);
    BOOST_CHECK_EQUAL(c[0],5.0);
    BOOST_CHECK_EQUAL(c[1],6.0);
    BOOST_CHECK_EQUAL(c[2],7.0);
    ndarray::Vector<double,3> d(5.0);
    BOOST_CHECK_EQUAL(d[0],5.0);
    BOOST_CHECK_EQUAL(d[1],5.0);
    BOOST_CHECK_EQUAL(d[2],5.0);
}

BOOST_AUTO_TEST_CASE(cores) {
    typedef ndarray::detail::Core<double,3> Core;
    ndarray::Vector<int,3> shape = ndarray::makeVector(4,3,2);
    ndarray::Vector<int,3> strides = ndarray::makeVector(6,2,1);
    Core::Ptr core = Core::create(shape,strides);
    BOOST_CHECK_EQUAL(core->getRC(),1);
    Core::Ptr copy = core;
    BOOST_CHECK_EQUAL(core->getRC(),2);
    copy.reset();
    BOOST_CHECK_EQUAL(core->getRC(),1);
}

BOOST_AUTO_TEST_CASE(allocation) {
    ndarray::Vector<int,3> shape = ndarray::makeVector(5,6,7);
    ndarray::Array<float,3,3> a = ndarray::allocate(shape);
    BOOST_CHECK_EQUAL(a.getShape(),shape);
    BOOST_CHECK_EQUAL(a.getData(),a.getOwner().get());
    ndarray::Array<float,3> b = ndarray::allocate(shape);
    BOOST_CHECK_EQUAL(b.getShape(),shape);
    BOOST_CHECK_EQUAL(b.getData(),b.getOwner().get());
    BOOST_CHECK_EQUAL(b.getSize<0>(),shape[0]);
    BOOST_CHECK_EQUAL(b.getSize<1>(),shape[1]);
    BOOST_CHECK_EQUAL(b.getSize<2>(),shape[2]);
    BOOST_CHECK_EQUAL(b.getStride<0>(),6*7);
    BOOST_CHECK_EQUAL(b.getStride<1>(),7);
    BOOST_CHECK_EQUAL(b.getStride<2>(),1);
    BOOST_CHECK_EQUAL(b.getStrides(),ndarray::makeVector(6*7,7,1));
    ndarray::Array<int,1> c = ndarray::allocate(ndarray::makeVector(5));
    BOOST_CHECK_EQUAL(b.size(),5);
    
}

BOOST_AUTO_TEST_CASE(external) {
    double data[3*4*2] = {0};
    ndarray::Vector<int,3> shape = ndarray::makeVector(3,4,2);
    ndarray::Vector<int,3> strides = ndarray::makeVector(8,2,1);
    ndarray::Array<double,3,3> a = ndarray::external(data,shape,strides);
    BOOST_CHECK_EQUAL(a.getData(),data);
    BOOST_CHECK(!a.getOwner());
    BOOST_CHECK_EQUAL(a.getShape(),shape);
    BOOST_CHECK_EQUAL(a.getStrides(),strides);
}

BOOST_AUTO_TEST_CASE(conversion) {
    double data[3*4*2] = {0};
    ndarray::Vector<int,3> shape = ndarray::makeVector(3,4,2);
    ndarray::Vector<int,3> strides = ndarray::makeVector(8,2,1);
    ndarray::Array<double,3,3> a = ndarray::external(data,shape,strides);
    ndarray::Array<double const,3> b = a;
}

BOOST_AUTO_TEST_CASE(shallow) {
    double data[3*4*2] = {0};
    ndarray::Vector<int,3> shape = ndarray::makeVector(3,4,2);
    ndarray::Vector<int,3> strides = ndarray::makeVector(8,2,1);
    ndarray::Array<double,3,3> a = ndarray::external(data,shape,strides);
    ndarray::Array<double,3,1> b = ndarray::external(data,shape,strides);
    BOOST_CHECK(ndarray::shallow(a) == ndarray::shallow(b));
    BOOST_CHECK(ndarray::shallow(a[2]) == ndarray::shallow(b[2]));
    BOOST_CHECK(ndarray::shallow(a[0][1]) == ndarray::shallow(b[0][1]));
    BOOST_CHECK(ndarray::shallow(a[0][1]) != ndarray::shallow(b[1][2]));
    ndarray::Array<double,3,3> c;
    ndarray::shallow(c) = a;
    BOOST_CHECK_EQUAL(a.getData(),c.getData());
    BOOST_CHECK_EQUAL(a.getShape(),c.getShape());
    BOOST_CHECK_EQUAL(a.getStrides(),c.getStrides());
    BOOST_CHECK(ndarray::shallow(a) == ndarray::shallow(c));
    ndarray::Array<double,2> d = c[1];
    BOOST_CHECK_EQUAL(d.getData(),c[1].getData());
    BOOST_CHECK_EQUAL(d.getShape(),c[1].getShape());
    BOOST_CHECK_EQUAL(d.getStrides(),c[1].getStrides());
    BOOST_CHECK(ndarray::shallow(d) == ndarray::shallow(c[1]));
}

BOOST_AUTO_TEST_CASE(casts) {
    double data[3*4*2] = {0};
    ndarray::Vector<int,3> shape = ndarray::makeVector(3,4,2);
    ndarray::Vector<int,3> strides = ndarray::makeVector(8,2,1);
    ndarray::Array<double const,3,1> a = ndarray::external(data,shape,strides);
    ndarray::Array<double const,3,2> b = ndarray::static_dimension_cast<2>(a);
    BOOST_CHECK(ndarray::shallow(a) == ndarray::shallow(b));
    ndarray::Array<double,3,1> c = ndarray::const_array_cast<double>(a);
    BOOST_CHECK(ndarray::shallow(a) == ndarray::shallow(c));
    ndarray::Array<double const,3,3> d = ndarray::dynamic_dimension_cast<3>(a);
    BOOST_CHECK(ndarray::shallow(a) == ndarray::shallow(d));
    ndarray::Array<double const,3,1> e = d[ndarray::view()(0,4,2)()];
    ndarray::Array<double const,3,3> f = ndarray::dynamic_dimension_cast<3>(e);
    BOOST_CHECK(f.empty());
}

BOOST_AUTO_TEST_CASE(indexing) {
    double data[3*4*2] = { 
         0, 1, 2, 3, 4, 5, 6, 7,
         8, 9,10,11,12,13,14,15,
        16,17,18,19,20,21,22,23,
    };
    ndarray::Vector<int,3> a_shape = ndarray::makeVector(4,3,2);
    ndarray::Vector<int,3> a_strides = ndarray::makeVector(6,2,1);
    ndarray::Array<double,3,3> a = ndarray::external(data,a_shape,a_strides);
    BOOST_CHECK(ndarray::shallow(a.front()) == ndarray::shallow(a[0]));
    BOOST_CHECK(ndarray::shallow(a.back()) == ndarray::shallow(a[a_shape[0]-1]));
    int n = 0;
    for (int i=0; i<a_shape[0]; ++i) {
        for (int j=0; j<a_shape[1]; ++j) {
            for (int k=0; k<a_shape[2]; ++k) {
                BOOST_CHECK_EQUAL(a[i][j][k],n);
                ++n;
            }
        }
    }
    ndarray::Vector<int,2> b_shape = ndarray::makeVector(8,3);
    ndarray::Vector<int,2> b_strides = ndarray::makeVector(1,8);
    ndarray::Array<double,2> b = ndarray::external(data,b_shape,b_strides);
    for (int i=0; i<b_shape[0]; ++i) {
        for (int j=0; j<b_shape[1]; ++j) {
            BOOST_CHECK_EQUAL(b[i][j],i+8*j);
        }
    }
    ndarray::Vector<int,2> c_shape = ndarray::makeVector(4,3);
    ndarray::Vector<int,2> c_strides = ndarray::makeVector(1,8);
    ndarray::Array<double,2> c = ndarray::external(data,c_shape,c_strides);
    for (int i=0; i<c_shape[0]; ++i) {
        for (int j=0; j<c_shape[1]; ++j) {
            BOOST_CHECK_EQUAL(c[i][j],i+8*j);
        }
    }
}

BOOST_AUTO_TEST_CASE(iterators) {
    double data[3*4*2] = { 
         0, 1, 2, 3, 4, 5, 6, 7,
         8, 9,10,11,12,13,14,15,
        16,17,18,19,20,21,22,23,
    };
    ndarray::Vector<int,3> a_shape = ndarray::makeVector(4,3,2);
    ndarray::Vector<int,3> a_strides = ndarray::makeVector(6,2,1);
    ndarray::Array<double,3,3> a = ndarray::external(data,a_shape,a_strides);
    ndarray::Array<double,3,3>::Iterator ai_iter = a.begin();
    ndarray::Array<double,3,3>::Iterator const ai_end = a.end();
    for (int i=0; ai_iter != ai_end; ++i, ++ai_iter) {
        ndarray::Array<double,3,3>::Reference::Iterator aj_iter = ai_iter->begin();
        ndarray::Array<double,3,3>::Reference::Iterator const aj_end = ai_iter->end();
        for (int j=0; aj_iter != aj_end; ++j, ++aj_iter) {
            ndarray::Array<double,3,3>::Reference::Reference::Iterator ak_iter = aj_iter->begin();
            ndarray::Array<double,3,3>::Reference::Reference::Iterator const ak_end = aj_iter->end();
            for (int k=0; ak_iter != ak_end; ++k, ++ak_iter) {
                BOOST_CHECK_EQUAL(a[i][j][k],*ak_iter);
            }
        }
    }
    ndarray::Vector<int,2> b_shape = ndarray::makeVector(4,3);
    ndarray::Vector<int,2> b_strides = ndarray::makeVector(1,8);
    ndarray::Array<double,2> b = ndarray::external(data,b_shape,b_strides);
    ndarray::Array<double,2>::Iterator bi_iter = b.begin();
    ndarray::Array<double,2>::Iterator const bi_end = b.end();
    for (int i=0; bi_iter != bi_end; ++i, ++bi_iter) {
        ndarray::Array<double,2>::Reference::Iterator bj_iter = bi_iter->begin();
        ndarray::Array<double,2>::Reference::Iterator const bj_end = bi_iter->end();
        for (int j=0; bj_iter != bj_end; ++j, ++bj_iter) {
            BOOST_CHECK_EQUAL(b[i][j],*bj_iter);
        }
    }

}

BOOST_AUTO_TEST_CASE(views) {
    ndarray::Vector<int,3> shape = ndarray::makeVector(4,3,2);
    ndarray::Array<double,3,3> a = ndarray::allocate(shape);
    BOOST_CHECK(ndarray::shallow(a) == ndarray::shallow(a[ndarray::view()()()]));
    BOOST_CHECK(ndarray::shallow(a) == ndarray::shallow(a[ndarray::view()]));
    BOOST_CHECK(ndarray::shallow(a[1]) == ndarray::shallow(a[ndarray::view(1)]));
    BOOST_CHECK(ndarray::shallow(a[1][2]) == ndarray::shallow(a[ndarray::view(1)(2)]));
    BOOST_CHECK(ndarray::shallow(a) != ndarray::shallow(a[ndarray::view(0,3)]));
    ndarray::Array<double const,2> b = a[ndarray::view()(1,3)(0)];
    BOOST_CHECK(b.getShape() == ndarray::makeVector(4,2));
    BOOST_CHECK(b.getStrides() == ndarray::makeVector(6,2));
    BOOST_CHECK(b.getData() == a.getData() + 2);
    ndarray::Array<double const,2> c = b[ndarray::view(0,4,2)()];
    BOOST_CHECK(c.getShape() == ndarray::makeVector(2,2));
    BOOST_CHECK(c.getStrides() == ndarray::makeVector(12,2));
    BOOST_CHECK(c.getData() == b.getData());
}

BOOST_AUTO_TEST_CASE(predicates) {
    double data1[3*4*2] = { 
         0, 1, 2, 3, 4, 5, 6, 7,
         8, 9,10,11,12,13,14,15,
        16,17,18,19,20,21,22,23,
    };
    double data2[3*4*2] = { 
         0, 1, 2, 3, 4, 5, 6, 7,
         8, 9,10,11,12,13,14,15,
        16,17,18,19,20,21,22,23,
    };
    ndarray::Vector<int,3> shape = ndarray::makeVector(4,3,2);
    ndarray::Vector<int,3> strides = ndarray::makeVector(6,2,1);
    ndarray::Array<double const,3,3> a = ndarray::external(data1,shape,strides);
    ndarray::Array<bool,3,2> b = ndarray::allocate(shape);
    ndarray::Array<bool,3> c = ndarray::allocate(shape);
    ndarray::Array<double,3,1> d = ndarray::external(data2,shape,strides);
    b = (a == 3.0);
    c = !b;
    BOOST_CHECK(ndarray::shallow(a) != ndarray::shallow(d));
    BOOST_CHECK(ndarray::all(a==d));
    BOOST_CHECK(ndarray::any(a==d));
    d[3][1][0] = 5.0;
    BOOST_CHECK(!ndarray::all(a == d));
    BOOST_CHECK(ndarray::any(a == d));
    BOOST_CHECK(ndarray::any(a != d));
    d = -5.0;
    BOOST_CHECK(!ndarray::all(a == d));
    BOOST_CHECK(ndarray::all(a != d));
    BOOST_CHECK(ndarray::all(a > d));
    BOOST_CHECK(!ndarray::any(a == d));
    for (int i=0; i<shape[0]; ++i) {
        for (int j=0; j<shape[1]; ++j) {
            for (int k=0; k<shape[2]; ++k) {
                if (a[i][j][k] == 3) {
                    BOOST_CHECK_EQUAL(b[i][j][k],true);
                    BOOST_CHECK_EQUAL(c[i][j][k],false);
                } else {
                    BOOST_CHECK_EQUAL(b[i][j][k],false);
                    BOOST_CHECK_EQUAL(c[i][j][k],true);
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(binary_ops) {
    float data[3*4*2] = { 
         0, 1, 2, 3, 4, 5, 6, 7,
         8, 9,10,11,12,13,14,15,
        16,17,18,19,20,21,22,23,
    };
    ndarray::Vector<int,3> shape = ndarray::makeVector(4,3,2);
    ndarray::Vector<int,3> strides = ndarray::makeVector(6,2,1);
    ndarray::Array<float,3,2> a = ndarray::external(data,shape,strides);
    ndarray::Array<double,3,1> b = ndarray::allocate(shape);
    ndarray::Array<double,3,3> c = ndarray::allocate(shape);
    c = 0.0;
    double q = 1.2;
    b = a + q;
    c -= (a * b - q);
    for (int i=0; i<shape[0]; ++i) {
        for (int j=0; j<shape[1]; ++j) {
            for (int k=0; k<shape[2]; ++k) {
                BOOST_CHECK_CLOSE(b[i][j][k],a[i][j][k]+q,1E-8);
                BOOST_CHECK_CLOSE(- (b[i][j][k] * a[i][j][k] - q), c[i][j][k],1E-8);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(assignment) {
    double data[3*4*2] = { 
         0, 1, 2, 3, 4, 5, 6, 7,
         8, 9,10,11,12,13,14,15,
        16,17,18,19,20,21,22,23,
    };
    ndarray::Vector<int,3> shape = ndarray::makeVector(3,4,2);
    ndarray::Vector<int,3> strides = ndarray::makeVector(8,2,1);
    ndarray::Array<double,3,3> a = ndarray::external(data,shape,strides);
    ndarray::Array<double,3,3> b = ndarray::allocate(shape);
    b = a;
    ndarray::Array<double,3,2> c = ndarray::allocate(shape);
    ndarray::Array<double const,3,1> d = c;
    c = a;
    BOOST_CHECK(ndarray::shallow(a) != ndarray::shallow(b));
    BOOST_CHECK(ndarray::shallow(a) != ndarray::shallow(c));
    int n = 0;
    for (int i=0; i<shape[0]; ++i) {
        for (int j=0; j<shape[1]; ++j) {
            for (int k=0; k<shape[2]; ++k) {
                BOOST_CHECK_EQUAL(b[i][j][k],n);
                BOOST_CHECK_EQUAL(c[i][j][k],n);
                BOOST_CHECK_EQUAL(d[i][j][k],n);
                ++n;
            }
        }
    }
    double q = 5.3;
    double p = 4.2;
    b = q;
    c = p;
    for (int i=0; i<shape[0]; ++i) {
        for (int j=0; j<shape[1]; ++j) {
            for (int k=0; k<shape[2]; ++k) {
                BOOST_CHECK_EQUAL(b[i][j][k],q);
                BOOST_CHECK_EQUAL(c[i][j][k],p);
            }
        }
    }
    b += a;
    c -= a;
    for (int i=0; i<shape[0]; ++i) {
        for (int j=0; j<shape[1]; ++j) {
            for (int k=0; k<shape[2]; ++k) {
                BOOST_CHECK_CLOSE(b[i][j][k],q,1E-8);
                BOOST_CHECK_CLOSE(c[i][j][k],p,1E-8);
                ++q;
                --p;
            }
        }
    }
}


BOOST_AUTO_TEST_CASE(transpose) {
    double data[3*4*2] = { 
         0, 1, 2, 3, 4, 5, 6, 7,
         8, 9,10,11,12,13,14,15,
        16,17,18,19,20,21,22,23,
    };
    ndarray::Vector<int,3> shape = ndarray::makeVector(3,4,2);
    ndarray::Vector<int,3> strides = ndarray::makeVector(8,2,1);
    ndarray::Array<double,3,3> a = ndarray::external(data,shape,strides);
    ndarray::Array<double const,3> b = a.transpose();
    ndarray::Array<double const,3> c = a.transpose(ndarray::makeVector(1,0,2));
    for (int i=0; i<shape[0]; ++i) {
        for (int j=0; j<shape[1]; ++j) {
            for (int k=0; k<shape[2]; ++k) {
                BOOST_CHECK_EQUAL(a[i][j][k],b[k][j][i]);
                BOOST_CHECK_EQUAL(a[i][j][k],c[j][i][k]);
            }
        }
    }
    BOOST_CHECK(ndarray::shallow(a[ndarray::view()(1)(1)]) == ndarray::shallow(b[ndarray::view(1)(1)()]));
    BOOST_CHECK(ndarray::shallow(a[ndarray::view(0)()(1)]) == ndarray::shallow(b[ndarray::view(1)()(0)]));
    BOOST_CHECK(ndarray::shallow(a[ndarray::view(0)(0)()]) == ndarray::shallow(b[ndarray::view()(0)(0)]));
    BOOST_CHECK(ndarray::shallow(a[ndarray::view()(1)(1)]) == ndarray::shallow(c[ndarray::view(1)()(1)]));
    BOOST_CHECK(ndarray::shallow(a[ndarray::view(0)()(1)]) == ndarray::shallow(c[ndarray::view()(0)(1)]));
    BOOST_CHECK(ndarray::shallow(a[ndarray::view(0)(0)()]) == ndarray::shallow(c[ndarray::view(0)(0)()]));
}
