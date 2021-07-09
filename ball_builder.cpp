//
// Created by 徐溶延local on 2021/7/8.
//

#include "ball_builder.h"
#include "math_utils.h"
#define _USE_MATH_DEFINES
#include <math.h>

void buildBall(double r,
               int longitudeNum,
               int altitudeNum,
               MatX3d &V, MatX3i &F) {
    assert(altitudeNum >= 3);
    assert(longitudeNum >= 3);
    double longitudeStep = 2 * M_PI / longitudeNum;
    double altitudeStep = M_PI / (altitudeNum - 1);
    std::vector<SphericalCoord> sphericalCoords;
    sphericalCoords.reserve(2 + (altitudeNum - 2) * longitudeNum);
    std::vector<double> phis(longitudeNum), thetas(altitudeNum);
    double initTheta = -0.5 * M_PI, initPhi = 0;
    for (int i = 0; i < altitudeNum; i++) {
        phis[i] = (double)i * longitudeStep;
    }
    for (int i = 0; i < longitudeNum; i++) {
        thetas[i] = (double)i * altitudeStep;
    }
    sphericalCoords.emplace_back(r, 0, 0);
    for (int i = 1; i < altitudeNum - 1; i++) {
        for (int j = 0; j < longitudeNum; j++) {
            sphericalCoords.emplace_back(r, thetas[i], phis[j]);
        }
    }
    sphericalCoords.emplace_back(r, M_PI, 0);

    V.resize(sphericalCoords.size(), 3);
    for (int i = 0; i < sphericalCoords.size(); i++) {
        V.row(i) = common::math::spherical2cartesian(sphericalCoords[i]);
    }
    F.resize(2 * longitudeNum + (2 * longitudeNum) * (altitudeNum - 3), 3);
    int fid = 0;
    int start = 1;
    for (int i = 0; i < longitudeNum; i++) {
        int a1 = start + i;
        int a2 = start + (i + 1) % longitudeNum;
        F.row(fid++) << a1, a2, 0;
    }
    for (int i = 0; i < altitudeNum - 3; i++) {
        for (int j = 0; j < longitudeNum; j++) {
            int a1 = start + j, a2 = start + (j + 1) % longitudeNum;
            int b1 = start + longitudeNum + j, b2 = start + longitudeNum + (j + 1) % longitudeNum;
            F.row(fid++) << a1, b1, b2;
            F.row(fid++) << b2, a2, a1;
        }
        start += longitudeNum;
    }
    int bottom = V.rows() - 1;
    for (int i = 0; i < longitudeNum; i++) {
        int a1 = start + i;
        int a2 = start + (i + 1) % longitudeNum;
        F.row(fid++) << a2, a1, bottom;
    }
}
