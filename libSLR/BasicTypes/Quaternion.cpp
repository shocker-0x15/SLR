//
//  Quaternion.cpp
//
//  Created by 渡部 心 on 2016/01/26.
//  Copyright c 2016年 渡部 心. All rights reserved.
//

#include "Quaternion.h"

namespace SLR {
    template struct SLR_API QuaternionTemplate<float>;
    template struct SLR_API QuaternionTemplate<double>;
    
    template <typename RealType>
    void decompose(const Matrix4x4Template<RealType> &mat, Vector3Template<RealType>* T, QuaternionTemplate<RealType>* R, Matrix4x4Template<RealType>* S) {
        T->x = mat[3][0];
        T->y = mat[3][1];
        T->z = mat[3][2];
        
        Matrix4x4Template<RealType> matRS = mat;
        for (int i = 0; i < 3; ++i)
            matRS[3][i] = matRS[i][3] = 0.0f;
        matRS[3][3] = 1.0f;
        
        RealType norm;
        int count = 0;
        Matrix4x4Template<RealType> curR = matRS;
        do {
            Matrix4x4Template<RealType> itR = invert(transpose(curR));
            Matrix4x4Template<RealType> nextR = 0.5f * (curR + itR);
            
            norm = 0;
            for (int i = 0; i < 3; ++i) {
                using std::abs;
                RealType n = abs(curR[0][i] - nextR[0][i]) + abs(curR[1][i] - nextR[1][i]) + abs(curR[2][i] - nextR[2][i]);
                norm = std::max(norm, n);
            }
            curR = nextR;
        } while (++count < 100 && norm > 0.0001);
        *R = QuaternionTemplate<RealType>(curR);
        
        *S = invert(curR) * matRS;
    }
    template SLR_API void decompose(const Matrix4x4Template<float> &mat, Vector3Template<float>* T, QuaternionTemplate<float>* R, Matrix4x4Template<float>* S);
    template SLR_API void decompose(const Matrix4x4Template<double> &mat, Vector3Template<double>* T, QuaternionTemplate<double>* R, Matrix4x4Template<double>* S);
}
