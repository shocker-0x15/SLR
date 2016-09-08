//
//  Matrix3x3.cpp
//
//  Created by 渡部 心 on 2016/09/08.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#include "Matrix3x3.h"

namespace SLR {
    template struct SLR_API Matrix3x3Template<float>;
    //    template SLR_API const Matrix3x3Template<float> Matrix3x3Template<float>::Identity;
    //    template SLR_API const Matrix3x3Template<float> Matrix3x3Template<float>::Zero;
    //    template SLR_API const Matrix3x3Template<float> Matrix3x3Template<float>::One;
    //    template SLR_API const Matrix3x3Template<float> Matrix3x3Template<float>::Inf;
    //    template SLR_API const Matrix3x3Template<float> Matrix3x3Template<float>::NaN;
    
    template struct SLR_API Matrix3x3Template<double>;
    //    template SLR_API const Matrix3x3Template<double> Matrix3x3Template<double>::Identity;
    //    template SLR_API const Matrix3x3Template<double> Matrix3x3Template<double>::Zero;
    //    template SLR_API const Matrix3x3Template<double> Matrix3x3Template<double>::One;
    //    template SLR_API const Matrix3x3Template<double> Matrix3x3Template<double>::Inf;
    //    template SLR_API const Matrix3x3Template<double> Matrix3x3Template<double>::NaN;
    
    template <typename RealType>
    Matrix3x3Template<RealType> transpose(const Matrix3x3Template<RealType> &m) {
        return Matrix3x3Template<RealType>(Vector3Template<RealType>(m.m00, m.m01, m.m02),
                                           Vector3Template<RealType>(m.m10, m.m11, m.m12),
                                           Vector3Template<RealType>(m.m20, m.m21, m.m22));
    }
    template SLR_API Matrix3x3Template<float> transpose(const Matrix3x3Template<float> &m);
    template SLR_API Matrix3x3Template<double> transpose(const Matrix3x3Template<double> &m);
    
    template <typename RealType>
    Matrix3x3Template<RealType> invert(const Matrix3x3Template<RealType> &m) {
        SLRAssert_NotImplemented();
        Matrix3x3Template<RealType> mat;
        
        return mat;
    }
    template SLR_API Matrix3x3Template<float> invert(const Matrix3x3Template<float> &m);
    template SLR_API Matrix3x3Template<double> invert(const Matrix3x3Template<double> &m);
    
    template <typename RealType>
    Matrix3x3Template<RealType> rotate3x3(RealType angle, const Vector3Template<RealType> &axis) {
        Matrix3x3Template<RealType> matrix;
        Vector3Template<RealType> nAxis = normalize(axis);
        RealType c = std::cos(angle);
        RealType s = std::sin(angle);
        RealType oneMinusC = 1 - c;
        
        matrix.m00 = nAxis.x * nAxis.x * oneMinusC + c;
        matrix.m10 = nAxis.x * nAxis.y * oneMinusC + nAxis.z * s;
        matrix.m20 = nAxis.z * nAxis.x * oneMinusC - nAxis.y * s;
        matrix.m01 = nAxis.x * nAxis.y * oneMinusC - nAxis.z * s;
        matrix.m11 = nAxis.y * nAxis.y * oneMinusC + c;
        matrix.m21 = nAxis.y * nAxis.z * oneMinusC + nAxis.x * s;
        matrix.m02 = nAxis.z * nAxis.x * oneMinusC + nAxis.y * s;
        matrix.m12 = nAxis.y * nAxis.z * oneMinusC - nAxis.x * s;
        matrix.m22 = nAxis.z * nAxis.z * oneMinusC + c;
        
        return matrix;
    }
    template SLR_API Matrix3x3Template<float> rotate3x3(float angle, const Vector3Template<float> &axis);
    template SLR_API Matrix3x3Template<double> rotate3x3(double angle, const Vector3Template<double> &axis);
}
