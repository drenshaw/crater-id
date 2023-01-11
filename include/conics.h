#pragma once

#include <tuple>
#include <vector>

#include "vector_math.h"

// template <typename T>
class Conic {
    public:
        Conic(const float=0, const float=0, const float=0, const float=0, const float=0);
        Conic(const std::tuple<float, float, float, float, float>);
        Conic(const std::vector<float>);
        void setGeometricParameters(const std::vector<float>);
        std::vector<float> Geom2Implicit();
        void NormalizeImplicitParameters(std::vector<float>&);
        std::vector<float> getGeom();
        std::vector<float> impl_to_geom(const std::vector<float>);
        void setFromImplicit(const std::vector<float> impl_params);

    private:
        float semimajor_axis_;
        float semiminor_axis_;
        float x_center_;
        float y_center_;
        float angle_;
        // std::tuple<float, float, float, float, float> geom;
};

class ConicImplicit {
    ConicImplicit();
};

class ConicGeometry {
    ConicGeometry();
};

class ConicMatrix {
    ConicMatrix();
};