#include <iostream>
#include <string>

#include <igl/writeOBJ.h>

#include "ball_builder.h"

int main(int argc, char** argv) {
    std::string fileName = argv[1];
    double r = std::stod(argv[2]);
    int l = std::stoi(argv[3]);
    int a = std::stoi(argv[4]);

    MatX3d V;
    MatX3i F;
    buildBall(r, l, a, V, F);
    igl::writeOBJ(fileName, V, F);
    return 0;
}
