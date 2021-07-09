#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include "types.h"
#include <vector>
#include<math.h>
#include <utility>
//#include <spdlog/spdlog.h>
#include "math_utils.h"
#include "util.h"

#define EPS2 1e-6

namespace common {

    Matd toHomoPointMat(const Matd &mat) {
        assert(mat.rows() == 3);
        Matd result(4, mat.cols());
        for (int i = 0; i < mat.cols(); i++) {
            result.col(i) = toHomoPoint(mat.col(i));
        }
        return result;
    }

    Matd toHomoVectorMat(const Matd &mat) {
        assert(mat.rows() == 3);
        Matd result(4, mat.cols());
        for (int i = 0; i < mat.cols(); i++) {
            result.col(i) = toHomoVector(mat.col(i));
        }
        return result;
    }

    Matd toPointMat(const Matd &mat) {
        assert(mat.rows() == 4);
        Matd result(3, mat.cols());
        for (int i = 0; i < mat.cols(); i++) {
            result.col(i) = toPoint(mat.col(i));
        }
        return result;
    }
} // namespace common
