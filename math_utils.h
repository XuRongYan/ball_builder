//
// Created by 徐溶延 on 2020/11/16.
//

#ifndef MESH_CUTTING_MATH_UTILS_H
#define MESH_CUTTING_MATH_UTILS_H

#define _USE_MATH_DEFINES
#include <cmath>


#include <Eigen/Dense>
#include "types.h"

namespace common::math {
    bool nearNum(double a, double b, double eps = 1e-6);

    double angle(const Point &A, const Point &B);

    double angle(const Point &A, const Point &B, const Point &C);

    double cos(const Point &A, const Point &B);

    bool intersect(const Point &A,
                   const Point &B,
                   const Point &C,
                   const Point &D);

    double distanceToEdge(const Point &p, const Point &A, const Point &B);

    long mod(long i, long size);

    inline size_t prev(size_t i, size_t n) {
        return i == 0 ? n - 1 : i - 1;
    }

    inline size_t next(size_t i, size_t n) {
        return (i + 1) % n;
    }

    inline size_t prev(size_t i, size_t step, size_t n) {
        return step == 1 ? prev(i, n) : prev(prev(i, n), step - 1, n);
    }

    inline size_t next(size_t i, size_t step, size_t n) {
        return step == 1 ? next(i, n) : next(next(i, n), step - 1, n);
    }

    Point bc2Point(const Eigen::Vector3d &A,
                   const Eigen::Vector3d &B,
                   const Eigen::Vector3d &C,
                   const Eigen::VectorXd &bc);

    Point bc2Point(const Eigen::RowVector3d &A,
                   const Eigen::RowVector3d &B,
                   const Eigen::RowVector3d &C,
                   double a, double b);

    Point bc2Point(const Point &A,
                   const Point &B,
                   const Point &C,
                   const Point &D,
                   double a, double b, double c, double d);

    Point bc2Point(const Point &A,
                   const Point &B,
                   const Point &C,
                   const Point &D,
                   const VecXd &bc);

    double faceArea(const Point &A,
                    const Point &B,
                    const Point &C);


    bool bcValid(double a, double b, double c);

    bool bcValid(double a, double b, double c, double d);

    bool bcValid(const Vec3d &bc);

    SphericalCoord cartesian2Spherical(const Point &x);

    Point spherical2cartesian(const SphericalCoord &sc);

    void localBasis3D(const Mat3Xd &V, const Mat3Xi &F, std::vector<Mat3d> &basis);

    std::vector<Mat3d> localBasis3D(const std::vector<Mat3d> &triangles);

    Mat3d localBasis3D(const Mat3d &triangle);

    Mat4d viewTransform(const Point &eye, const Point &ux, const Point &uz);

    Mat4d localTransform(const Point &a, const Point &b, const Point &c);


} // namespace common




#endif //MESH_CUTTING_MATH_UTILS_H
