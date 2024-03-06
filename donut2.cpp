#include<vector>
#include<cmath>
#include<chrono>
#include<iostream>

using namespace std;

// MODIFY
const int FPS = 10;
const double R = 4.0;
const double r = 2.2;
const int WINDOW_WIDTH_CHAR = 61;
const int WINDOW_HEIGHT_CHAR = 31;

// DO NOT MODIFY
const int RAYS_PER_UNIT_X = 4;
const int RAYS_PER_UNIT_Y = 2;
const double WINDOW_HEIGHT = 2.5 * (R + r);
const double WINDOW_WIDTH = 2.5 * (R + r);

struct coord3D {
    double x;
    double y;
    double z;

    bool operator==(const coord3D other) const {
        const double tolerance = 1e-9;
        return abs(x - other.x) < tolerance && abs(y - other.y) < tolerance && abs(z - other.z) < tolerance;
    }
};

vector<coord3D> vector_intersects_donut(coord3D& loc, coord3D& dir) {
    // x = loc.x + dir.x * t
    // y = loc.y + dir.y * t
    // z = loc.z + dir.z * t

    



    return {};
}

double distance(coord3D& a, coord3D& b) {
    return sqrt(a.x * b.x + a.y * b.y + a.z * b.z);
}

// https://web.cs.ucdavis.edu/~amenta/s12/findnorm.pdf
coord3D get_normal(coord3D& point) {

    return {};
}

// cos(theta) = a.dot(b) / (magnitude(a) * magnitude(b))
float angle_between_vectors(coord3D& a, coord3D& b) {

    return 0.0;
}

void render_ascii_donut(coord3D& light, coord3D& view) {

}

// quarternion or rotation matrix
coord3D rotate_point(coord3D& point) {

    return {};
}

int main() {

    // set up window
    // set up light
    // set up view
    // while true
        // render_ascii_donut
        // light = rotate_point(light)
        // view = rorate_point(view)
        // wait till end of frame

    return 0;
}