//
// Created by 徐溶延 on 2020/11/16.
//

#include "math_utils.h"
#include <cmath>
#ifdef USE_OPENMP

#include <omp.h>

#endif

namespace common {
    namespace math {

        double angle(const Point &A, const Point &B) {
            return atan2(A.cross(B).norm(), A.dot(B));
        }

        double angle(const Point &A, const Point &B, const Point &C) {
            return angle(B - A, C - A);
        }

        double cos(const Point &A, const Point &B) {
            return A.dot(B) / (A.norm() * B.norm());
        }

        bool nearNum(double a, double b, double eps) {
            double x = abs(a - b);
            return (x <= eps);
        }

        int compare(const Point &A, const Point &B) {
            if (A.x() != B.x()) {
                return A.x() < B.x();
            }
            return A.y() < B.y();
        }


        bool
        intersect(const Point &A,
                  const Point &B,
                  const Point &C,
                  const Point &D) {
            const auto ab = B - A;
            const auto ad = D - A;
            const auto ac = C - A;
            const auto cd = D - C;
            const auto ca = A - C;
            const auto cb = B - C;
            const auto acab = ac.cross(ab);
            const auto adab = ad.cross(ab);
            const auto cacd = ca.cross(cd);
            const auto cbcd = cb.cross(cd);
            const auto abcd = ab.cross(cd);
            float side_cd = acab.dot(adab);
            float side_ab = cacd.dot(cbcd);

            // normal condition
            if (side_ab < 0 && side_cd < 0) {
                return true;
            }
            // only one point intersected
            if ((side_ab < 0 && side_cd == 0) || (side_cd < 0 && side_ab == 0)) {
                return true;
            }
            // overlap
            if (abcd.norm() == 0) {
                const auto abac = ab.cross(ac);
                if (abac.norm() == 0) {
                    Eigen::Vector3d s1, s2, e1, e2;
                    s1 = A, e1 = B, s2 = C, e2 = D;
                    if (!compare(A, B))
                        std::swap(s1, e1);
                    if (!compare(C, D))
                        std::swap(s2, e2);
                    if (compare(s1, s2) && compare(e2, e1)) {
                        return true;
                    }
                }
            }
            return false;
        }

        double distanceToEdge(const Point &p,
                              const Point &A,
                              const Point &B) {
            Point ap = p - A;
            Point ab = B - A;
            double d = ap.cross(ab).norm() / ab.norm();
            return d;
        }

        long mod(long i, long size) {
            if (i >= 0) return i % size;
            return size + i;
        }

        Point bc2Point(const Eigen::Vector3d &A,
                       const Eigen::Vector3d &B,
                       const Eigen::Vector3d &C,
                       const Eigen::VectorXd &bc) {
            Point pd{0, 0, 0};
            if (bc.size() == 2) {
                pd = (1 - bc.x() - bc.y()) * A + bc.x() * B + bc.y() * C;
            } else if (bc.size() == 3) {
                pd = bc.x() * A + bc.y() * B + bc.z() * C;
            } else {
                std::cerr << "bc2Point: invalid barycentric coordinate" << std::endl;
            }
            return Point(pd.x(), pd.y(), pd.z());
        }


        Point
        bc2Point(const Eigen::RowVector3d &A,
                 const Eigen::RowVector3d &B,
                 const Eigen::RowVector3d &C,
                 double a, double b) {
            const auto &p = a * A + b * B + (1 - a - b) * C;
            return Point(p.x(), p.y(), p.z());
        }

        Point math::bc2Point(const Point &A, const Point &B, const Point &C,
                             const Point &D, double a, double b, double c, double d) {
            assert(bcValid(a, b, c, d));
            return a * A + b * B + c * C + d * D;
        }

        Point bc2Point(const Point &A,
                       const Point &B,
                       const Point &C,
                       const Point &D,
                       const VecXd &bc) {
            assert(bc.size() == 4);
            return bc2Point(A, B, C, D, bc.x(), bc.y(), bc.z(), bc.w());
        }

        bool bcValid(const Vec3d &bc) {
            return bcValid(bc.x(), bc.y(), bc.z());
        }

        bool bcValid(double a, double b, double c) {
            const double eps = 1e-5;
            double sum = a + b + c;
            if (!nearNum(sum, 1, 1e-5)) return false;
            if (a > 1 + eps || a < -eps) return false;
            if (b > 1 + eps || b < -eps) return false;
            if (c > 1 + eps || c < -eps) return false;
            return true;
        }

        bool bcValid(double a, double b, double c, double d) {
            const double eps = 1e-5;
            double sum = a + b + c + d;
            if (!nearNum(sum, 1, 1e-5)) return false;
            if (a > 1 + eps || a < -eps) return false;
            if (b > 1 + eps || b < -eps) return false;
            if (c > 1 + eps || c < -eps) return false;
            if (d > 1 + eps || d < -eps) return false;
            return true;
        }

        double faceArea(const Point &A, const Point &B, const Point &C) {
            const auto &AB = B - A;
            const auto &AC = C - A;
            return 0.5 * fabs((AB.cross(AC)).norm());
        }

        SphericalCoord cartesian2Spherical(const Point &x) {
            double r = std::sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
            double theta = std::acos(x[2] / r);
            double phi = std::atan2(x[1], x[0]);
            return SphericalCoord(r, theta, phi);
        }

        Point spherical2cartesian(const SphericalCoord &sc) {
            double x = sc.r * std::sin(sc.theta) * std::cos(sc.phi);
            double y = sc.r * std::sin(sc.theta) * std::sin(sc.phi);
            double z = sc.r * std::cos(sc.theta);
            return Point(x, y, z);
        }



        void localBasis3D(const Mat3Xd &V, const Mat3Xi &F, std::vector<Mat3d> &basis) {
            basis.resize(F.cols());
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
            for (int i = 0; i < F.cols(); i++) {
                Mat3d triangle;
                for (int j = 0; j < 3; j++) {
                    triangle.col(j) = V.col(F(j, i));
                }
                basis[i] = localBasis3D(triangle);
            }
        }

        std::vector<Mat3d> math::localBasis3D(const std::vector<Mat3d> &triangles) {
            std::vector<Mat3d> results(triangles.size());
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
            for (int i = 0; i < triangles.size(); i++) {
                results[i] = localBasis3D(triangles[i]);
            }
            return results;
        }

        Mat3d localBasis3D(const Mat3d &triangle) {
            Mat3d basis;
            const double l10 = (triangle.col(1) - triangle.col(0)).norm();
            const Point v20 = triangle.col(2) - triangle.col(0);
            const Point v10 = triangle.col(1) - triangle.col(0);
            const double u3 = v20.dot(v10.normalized());
            const double v3 = v20.cross(v10.normalized()).norm();
            basis << 0, l10, u3,
                    0, 0, v3,
                    0, 0, 0;
            return basis;
        }

        Mat4d viewTransform(const Point &eye, const Point &ux, const Point &uz) {
            Mat4d result;
            Point normUx = ux.normalized();
            Point normUz = uz.normalized();
            const Point normUy = normUz.cross(normUx);
            result.row(0) << normUx[0], normUx[1], normUx[2], -eye.dot(normUx);
            result.row(1) << normUy[0], normUy[1], normUy[2], -eye.dot(normUy);
            result.row(2) << normUz[0], normUz[1], normUz[2], -eye.dot(normUz);
            result.row(3) << 0, 0, 0, 1;
            return result;
        }

        Mat4d localTransform(const Point &a, const Point &b, const Point &c) {
            const Point ux = b - a;
            const Point uz = ux.cross(c - a);
            return viewTransform(a, ux, uz);
        }


    } // namespace math
} // namespace common

