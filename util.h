#ifndef COMMON_UTIL_H
#define COMMON_UTIL_H

#include <Eigen/Dense>
#include "types.h"

namespace common {
static const float eps = 1.1e-6;

int barycentric(const Eigen::Vector3f &a, const Eigen::Vector3f &b,
                const Eigen::Vector3f &c, const Eigen::Vector3f &p,
                Eigen::Vector3f &bary);

bool isBaryValid(Eigen::Vector3d &bary);

inline void putInUnitBox(Mat3Xd &V, Point &minCoeff, double &maxCoeff) {
    if (V.cols() == 0) return;
    minCoeff = V.rowwise().minCoeff();
    V.colwise() -= minCoeff;
    maxCoeff = V.maxCoeff();
    V /= maxCoeff;
}

inline void  restoreFromUnitBox(const Point& minCoeff, double  maxCoeff, Mat3Xd &V) {
    if (V.cols() == 0) return;
    V *= maxCoeff;
    V.colwise() += minCoeff;
}

} // namespace common

#endif //COMMON_UTIL_H
