#pragma once

#include <vector>
#include <set>
#include <ostream>
#include <Eigen/Core>
#include <array>
#include <map>
#include <set>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <SurfaceMesh/SurfaceMesh.h>



//using namespace Surface_Mesh;
using Surface_Mesh::SurfaceMesh;
using Surface_Mesh::Point;
using SparseMatrixf = Eigen::SparseMatrix<float>;
using SparseMatrixd = Eigen::SparseMatrix<double, Eigen::ColMajor>;
using Tripletf = Eigen::Triplet<float>;
using Tripletd = Eigen::Triplet<double>;
using RowMatf = Eigen::Matrix<float, -1, -1, Eigen::RowMajor>;
using RowMati = Eigen::Matrix<int, -1, -1, Eigen::RowMajor>;
using RowMat3f = Eigen::Matrix<float, 3, 3, Eigen::RowMajor>;
using RowMat2f = Eigen::Matrix<float, 2, 2, Eigen::RowMajor>;
using RowMat32f = Eigen::Matrix<float, 3, 2, Eigen::RowMajor>;
using RowMat34i = Eigen::Matrix<size_t, 3, 4, Eigen::RowMajor>;
using RowMat34f = Eigen::Matrix<float, 3, 4, Eigen::RowMajor>;
using RowMat43f = Eigen::Matrix<float, 4, 3, Eigen::RowMajor>;
using Matd = Eigen::MatrixXd;
using Matf = Eigen::Matrix<float, -1, -1, Eigen::ColMajor>;
using Mati = Eigen::MatrixXi;
using Mat3f = Eigen::Matrix<float, 3, 3, Eigen::ColMajor>;
using Mat3d = Eigen::Matrix<double, 3, 3>;
using Mat4d = Eigen::Matrix<double, 4, 4>;
using Matf = Eigen::Matrix<float, -1, -1>;
using Mat3Xf = Eigen::Matrix3Xf;
using Mat3Xd = Eigen::Matrix3Xd;
using Mat2Xd = Eigen::Matrix2Xd;
using Mat2Xi = Eigen::Matrix2Xi;
using Mat3Xi = Eigen::Matrix3Xi;
using Mat4Xi = Eigen::Matrix4Xi;
using MatX3f = Eigen::MatrixX3f;
using MatX3d = Eigen::MatrixX3d;
using MatX3i = Eigen::MatrixX3i;
using Mat34i = Eigen::Matrix<size_t, 3, 4, Eigen::ColMajor>;
using Mat34f = Eigen::Matrix<float, 3, 4, Eigen::ColMajor>;
using Mat43f = Eigen::Matrix<float, 4, 3, Eigen::ColMajor>;
using RowVec2f = Eigen::RowVector2f;
using RowVec3f = Eigen::RowVector3f;
using RowVec3i = Eigen::RowVector3i;
using RowVec4i = Eigen::Matrix<int, 4, 1, Eigen::RowMajor>;
using RowVecXf = Eigen::RowVectorXf;
using RowVecXi = Eigen::RowVectorXi;
using Vec2f = Eigen::Vector2f;
using Vec3f = Eigen::Vector3f;
using Vec3d = Eigen::Vector3d;
using Vec4d = Eigen::Vector4d;
using Vec2d = Eigen::Vector2d;
using Vec3i = Eigen::Vector3i;
using Vec4i = Eigen::Vector4i;
using VecXf = Eigen::VectorXf;
using VecXd = Eigen::VectorXd;
using VecXi = Eigen::VectorXi;
using RowVec2d = Eigen::RowVector2d;
using RowVec3d = Eigen::RowVector3d;
using RowMatd = Eigen::Matrix<double, -1, -1, Eigen::RowMajor>;
using RowMati = Eigen::Matrix<int, -1, -1, Eigen::RowMajor>;
using RowMat3d = Eigen::Matrix<double, 3, 3, Eigen::RowMajor>;
using RowMat34i = Eigen::Matrix<size_t, 3, 4, Eigen::RowMajor>;
using RowMat34d = Eigen::Matrix<double, 3, 4, Eigen::RowMajor>;
using RowMat43d = Eigen::Matrix<double, 4, 3, Eigen::RowMajor>;
// homogeneous point
using HomoPoint = Eigen::Vector4d;

struct SphericalCoord {
    double r, theta, phi;

    SphericalCoord() {}

    SphericalCoord(double r, double theta, double phi) : r(r), theta(theta), phi(phi) {}
};

namespace common {
    const double EPS3 = 1e-2;

    struct Tet {
        std::vector<Vec3d> V{};
        std::vector<Vec4i> T{};
    };

    enum DIRECTION {
        BILATERAL, // bi-direction
        CW,  // clockwise
        CCW //anti-clockwise
    };


    enum LABEL_TYPE {
        TYPE_NONE = -1,             // inter-point
        TYPE_OUT = 0,               // boundary 1
        TYPE_IN = 1,                    // boundary 0
        TYPE_BOTH = 999,            // boundary belongs to both side of curve
        TYPE_BOUND = 998,           // face or vertex on curve segment
        TYPE_MULTI_LINE = 997,      // two segments in one face
        TYPE_UNCONNECT = 996,       // not belong
    };


    struct argument {
        std::string in_path;     // input mesh's file path
        std::string out_path;    // output polyline's file path
        std::string out_init_path;  // output polyline's initial path
        std::string data_path;    // control point's file path
        std::string topo_path;      // script of controlling topological operations
        bool loop;         // is closed curve
    };

    enum point_type {
        UNDEFINED_POINT = -1,
        VERTEX_POINT,
        EDGE_POINT,
        FACE_POINT
    };

    enum SpacePointType {
        SPACE_POINT,
        SPACE_EDGE,
        SPACE_FACE,
        SPACE_SUSPEND,
        SPACE_INVALID,
    };

    Matd toHomoPointMat(const Matd &mat);

    Matd toHomoVectorMat(const Matd &mat);

    Matd toPointMat(const Matd &mat);

    inline HomoPoint toHomoPoint(const Point &p) {
        return HomoPoint(p.x(), p.y(), p.z(), 1);
    }

    inline HomoPoint toHomoVector(const Point &p) {
        return HomoPoint(p.x(), p.y(), p.z(), 0);
    }

    inline Point toPoint(const HomoPoint &p) {
        return Point(p.x(), p.y(), p.z());
    }
} // namespace common
