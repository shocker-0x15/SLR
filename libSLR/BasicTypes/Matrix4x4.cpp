//
//  Matrix4x4.cpp
//
//  Created by 渡部 心 on 2016/01/26.
//  Copyright c 2016年 渡部 心. All rights reserved.
//

#include "Matrix4x4.h"

namespace SLR {
    template struct SLR_API Matrix4x4Template<float>;
//    template SLR_API const Matrix4x4Template<float> Matrix4x4Template<float>::Identity;
//    template SLR_API const Matrix4x4Template<float> Matrix4x4Template<float>::Zero;
//    template SLR_API const Matrix4x4Template<float> Matrix4x4Template<float>::One;
//    template SLR_API const Matrix4x4Template<float> Matrix4x4Template<float>::Inf;
//    template SLR_API const Matrix4x4Template<float> Matrix4x4Template<float>::NaN;

    template struct SLR_API Matrix4x4Template<double>;
//    template SLR_API const Matrix4x4Template<double> Matrix4x4Template<double>::Identity;
//    template SLR_API const Matrix4x4Template<double> Matrix4x4Template<double>::Zero;
//    template SLR_API const Matrix4x4Template<double> Matrix4x4Template<double>::One;
//    template SLR_API const Matrix4x4Template<double> Matrix4x4Template<double>::Inf;
//    template SLR_API const Matrix4x4Template<double> Matrix4x4Template<double>::NaN;
    
    template <typename RealType>
    Matrix4x4Template<RealType> transpose(const Matrix4x4Template<RealType> &m) {
        return Matrix4x4Template<RealType>(Vector4Template<RealType>(m.m00, m.m01, m.m02, m.m03),
                                           Vector4Template<RealType>(m.m10, m.m11, m.m12, m.m13),
                                           Vector4Template<RealType>(m.m20, m.m21, m.m22, m.m23),
                                           Vector4Template<RealType>(m.m30, m.m31, m.m32, m.m33));
    }
    template SLR_API Matrix4x4Template<float> transpose(const Matrix4x4Template<float> &m);
    template SLR_API Matrix4x4Template<double> transpose(const Matrix4x4Template<double> &m);
    
    template <typename RealType>
    Matrix4x4Template<RealType> invert(const Matrix4x4Template<RealType> &m) {
        Matrix4x4Template<RealType> mat = m;
        
        bool colDone[] = {false, false, false, false};
        typedef std::pair<int, int> SwapPair;
        SwapPair swapPairs[] = {SwapPair(0, 0), SwapPair(0, 0), SwapPair(0, 0), SwapPair(0, 0)};
        for (int pass = 0; pass < 4; ++pass) {
            int pvCol = 0;
            int pvRow = 0;
            RealType maxPivot = -1.0f;
            for (int c = 0; c < 4; ++c) {
                if (colDone[c])
                    continue;
                for (int r = 0; r < 4; ++r) {
                    if (colDone[r])
                        continue;
                    
                    RealType absValue = std::fabs(mat[c][r]);
                    if (absValue > maxPivot) {
                        pvCol = c;
                        pvRow = r;
                        maxPivot = absValue;
                    }
                }
            }
            
            mat.swapRows(pvRow, pvCol);
            swapPairs[pass] = SwapPair(pvRow, pvCol);
            
            RealType pivot = mat[pvCol][pvCol];
            if (pivot == 0.0f)
                return Matrix4x4Template<RealType>::NaN;
            
            mat[pvCol][pvCol] = 1.0f;
            mat.scaleRow(pvCol, 1.0f / pivot);
            Vector4Template<RealType> addendRow = mat.row(pvCol);
            for (int r = 0; r < 4; ++r) {
                if (r != pvCol) {
                    RealType s = mat[pvCol][r];
                    mat[pvCol][r] = 0.0f;
                    mat.addRow(r, -s * addendRow);
                }
            }
            
            colDone[pvCol] = true;
        }
        
        for (int pass = 3; pass >= 0; --pass) {
            const SwapPair &pair = swapPairs[pass];
            mat.swapColumns(pair.first, pair.second);
        }
        
        return mat;
    }
    template SLR_API Matrix4x4Template<float> invert(const Matrix4x4Template<float> &m);
    template SLR_API Matrix4x4Template<double> invert(const Matrix4x4Template<double> &m);
    
    template <typename RealType>
    Matrix4x4Template<RealType> lookAt(const Point3Template<RealType> &eye, const Point3Template<RealType> &tgt, const Vector3Template<RealType> &up) {
        Vector3Template<RealType> z = normalize(eye - tgt);
        Vector3Template<RealType> x = normalize(cross(up, z));
        Vector3Template<RealType> y = cross(z, x);
        Vector4Template<RealType> t = Vector4Template<RealType>(-dot(Vector3Template<RealType>(eye), x),
                                                                -dot(Vector3Template<RealType>(eye), y),
                                                                -dot(Vector3Template<RealType>(eye), z), 1.0f);
        
        return Matrix4x4Template<RealType>(Vector4Template<RealType>(x.x, y.x, z.x, 0.0f),
                                           Vector4Template<RealType>(x.y, y.y, z.y, 0.0f),
                                           Vector4Template<RealType>(x.z, y.z, z.z, 0.0f),
                                           t);
    }
    template SLR_API Matrix4x4Template<float> lookAt(const Point3Template<float> &eye, const Point3Template<float> &tgt, const Vector3Template<float> &up);
    template SLR_API Matrix4x4Template<double> lookAt(const Point3Template<double> &eye, const Point3Template<double> &tgt, const Vector3Template<double> &up);
    
    template <typename RealType>
    Matrix4x4Template<RealType> rotate(RealType angle, const Vector3Template<RealType> &axis) {
        Matrix4x4Template<RealType> matrix;
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
        
        matrix.m30 = matrix.m31 = matrix.m32 =
        matrix.m03 = matrix.m13 = matrix.m23 = 0.0f;
        matrix.m33 = 1.0f;
        
        return matrix;
    }
    template SLR_API Matrix4x4Template<float> rotate(float angle, const Vector3Template<float> &axis);
    template SLR_API Matrix4x4Template<double> rotate(double angle, const Vector3Template<double> &axis);
    
    template <typename RealType>
    Matrix4x4Template<RealType> camera(RealType aspect, RealType fovY, RealType near, RealType far) {
        Matrix4x4Template<RealType> matrix;
        RealType f = 1.0f / std::tan(fovY / 2);
        RealType dz = far - near;
        
        matrix.m00 = f / aspect;
        matrix.m11 = f;
        matrix.m22 = -(near + far) / dz;
        matrix.m32 = -1.0f;
        matrix.m23 = -2.0f * far * near / dz;
        matrix.m10 = matrix.m20 = matrix.m30 =
        matrix.m01 = matrix.m21 = matrix.m31 =
        matrix.m02 = matrix.m12 =
        matrix.m03 = matrix.m13 = matrix.m33 = 0.0f;
        
        return matrix;
    }
    template SLR_API Matrix4x4Template<float> camera(float aspect, float fovY, float near, float far);
    template SLR_API Matrix4x4Template<double> camera(double aspect, double fovY, double near, double far);
}
