// Converted with https://www.codeconvert.ai/fortran-to-c++-converter
#include <iostream>
#include <cmath>
#include <functional>
extern "C" void SCALE_AWARE(float dx, float PBL1, float& Psig_bl, float& Psig_shcu);
using std::bind;

void SCALE_AWARE(float dx, float PBL1, float& Psig_bl, float& Psig_shcu) {
    float dxdh;
    Psig_bl = 1.0;
    Psig_shcu = 1.0;
    dxdh = std::max<unsigned long>(2.5 * dx, 10.0) / std::min<unsigned long>(PBL1, 3000.0);
    Psig_bl = ((std::pow(dxdh, 2) + 0.106 * std::pow(dxdh, 0.667)) / (std::pow(dxdh, 2) + 0.066 * std::pow(dxdh, 0.667) + 0.071));
    dxdh = std::max<unsigned long>(2.5 * dx, 10.0) / std::min<unsigned long>(PBL1 + 500.0, 3500.0);
    Psig_shcu = ((std::pow(dxdh, 2) + 0.145 * std::pow(dxdh, 0.667)) / (std::pow(dxdh, 2) + 0.172 * std::pow(dxdh, 0.667) + 0.170));
    if (Psig_bl > 1.0) Psig_bl = 1.0;
    if (Psig_bl < 0.0) Psig_bl = 0.0;
    if (Psig_shcu > 1.0) Psig_shcu = 1.0;
    if (Psig_shcu < 0.0) Psig_shcu = 0.0;
}

int main() {
    float dx, PBL1, Psig_bl, Psig_shcu;
    // Get input values for dx and PBL1
    std::cout << "Enter the value of dx: ";
    std::cin >> dx;
    std::cout << "Enter the value of PBL1: ";
    std::cin >> PBL1;

    SCALE_AWARE(dx, PBL1, Psig_bl, Psig_shcu);

    std::cout << "Psig_bl: " << Psig_bl << std::endl;
    std::cout << "Psig_shcu: " << Psig_shcu << std::endl;

    return 0;
}

