#include<vector>
#include<cmath>
#include<chrono>
#include<iostream>
#include<limits>

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
        const double tolerance = 1e-2;
        return abs(x - other.x) < tolerance && abs(y - other.y) < tolerance && abs(z - other.z) < tolerance;
    }

    void normalize() {
        double magnitude = x * x + y * y + z * z;
        x = x / magnitude;
        y = y / magnitude;
        z = z / magnitude;
    }
};

// (x^2 + y^2 + z^2 + R^2 + r^2)^2 -4R^2 * (x^2 + y^2) < 0
bool point_inside_donut(coord3D p) {
    static double _4R2 = 4 * R * R;
    static double R2_minus_r2 = R * R - r * r;
    double x2_plus_y2 = p.x * p.x + p.y + p.y;
    double inside_sq = (x2_plus_y2 + p.z * p.z + R2_minus_r2);
    return inside_sq * inside_sq - _4R2 * x2_plus_y2 < 0;
}

coord3D vector_intersects_donut(coord3D& loc, coord3D& dir) {
    // r(t) = loc + dir * t
    // x = loc.x + dir.x * t
    // y = loc.y + dir.y * t
    // z = loc.z + dir.z * t

    static double jump = 0.1;
    static double jump_jump = 0.01;
    
    // the intersection of the vector with the bounds tours
    double t_x1 = -R / (loc.x + dir.x);
    double t_x2 = R / (loc.x + dir.x);
    double t_y1 = -R / (loc.y + dir.y);
    double t_y2 = R / (loc.y + dir.y);
    double t_z1 = -R / (loc.z + dir.z);
    double t_z2 = R / (loc.z + dir.z);
    double t = min(min(min(t_x1, t_x2), min(t_y1, t_y2)), min(t_z1, t_z2));
    double t_max = max(max(max(t_x1, t_x2), min(t_y1, t_y2)), min(t_z1, t_z2));

    while (t < t_max && !point_inside_donut({loc.x + dir.x * t, loc.y + dir.y * t, loc.z + dir.z * t})) {
        t += jump;
    } if (t >= t_max) {
        return {numeric_limits<double>::infinity(), numeric_limits<double>::infinity(), numeric_limits<double>::infinity()};
    } while (point_inside_donut({loc.x + dir.x * t, loc.y + dir.y * t, loc.z + dir.z * t})) {
        t -= jump_jump;
    } 
    
    // estimate theta and phi and estimate x, y, z from those estimations
    
    return {loc.x + dir.x * t, loc.y + dir.y * t, loc.z + dir.z * t};
}

double distance(coord3D& a, coord3D& b) {
    return sqrt(a.x * b.x + a.y * b.y + a.z * b.z);
}

// https://web.cs.ucdavis.edu/~amenta/s12/findnorm.pdf
coord3D get_normal(coord3D& point) {
    // estimate theta and phi

    // get normal from theta and phi
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