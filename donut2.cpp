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
        double magnitude = sqrt(x * x + y * y + z * z);
        x /= magnitude;
        y /= magnitude;
        z /= magnitude;
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

    double t1 = numeric_limits<double>::infinity(), t2 = numeric_limits<double>::infinity();

    static auto min2 = [](double t, double& t1, double& t2) {
        if (t < t1) {
            t2 = t1;
            t1 = t;
        } else if (t < t2) {
            t2 = t;
        }
    };

    if (dir.x != 0) {
        min2((R - loc.x) / dir.x, t1, t2);
        min2((-R - loc.x) / dir.x, t1, t2);
    } if (dir.y != 0) {
        min2((R - loc.y) / dir.y, t1, t2);
        min2((-R - loc.y) / dir.y, t1, t2);
    } if (dir.z != 0) {
        min2((r - loc.z) / dir.z, t1, t2);
        min2((-r - loc.z) / dir.z, t1, t2);
    }

    double t = t1;

    while (t < t2 && !point_inside_donut({loc.x + dir.x * t, loc.y + dir.y * t, loc.z + dir.z * t})) {
        t += jump;
    } if (t >= t2) {
        return {numeric_limits<double>::infinity(), numeric_limits<double>::infinity(), numeric_limits<double>::infinity()};
    } while (point_inside_donut({loc.x + dir.x * t, loc.y + dir.y * t, loc.z + dir.z * t})) {
        t -= jump_jump;
    }
    
    return {loc.x + dir.x * t, loc.y + dir.y * t, loc.z + dir.z * t};
}

double distance(coord3D& a, coord3D& b) {
    return sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2) + pow(a.z - b.z, 2));
}

// https://web.cs.ucdavis.edu/~amenta/s12/findnorm.pdf
coord3D get_normal(coord3D& point) {
    // estimate theta and phi
    static double R2 = R * R;
    static double _2PI = 2 * M_PI;
    static double _3_over_2_PI = 3 * M_PI_2;
    double phi = atan2(point.y, point.x);
    double theta = asin(point.z / r);
    bool upper = point.z > 0;
    bool right = point.x * point.x + point.y * point.y > R2;
    if (point.z == 0.0 && right) {
        theta = 0;
    } else if (point.z == r) {
        theta = M_PI_2;
    } else if (point.z == 0.0 && !right) {
        theta = M_PI;
    } else if (point.z == -r) {
        theta = _3_over_2_PI;
    } else if (upper && right) {
        theta = theta;
    } else if (upper && !right) {
        theta = M_PI - theta;
    } else if (!upper && !right) {
        theta = M_PI - theta;
    } else {
        theta = _2PI + theta;
    }

    // get normal from theta and phi
    double tx = -sin(phi);
    double ty = cos(phi);
    double tz = 0;
    double sx = cos(phi) * -sin(theta);
    double sy = sin(phi) * -sin(theta);
    double sz = cos(theta);
    double nx = ty * sz - tz * sy;
    double ny = tz * sx - tx * sz;
    double nz = tx * sy - ty * sx;

    coord3D normal = {nx, ny, nz};
    normal.normalize();
    return normal;
}

// cos(theta) = a.dot(b) / (magnitude(a) * magnitude(b))
float angle_between_vectors(coord3D& a, coord3D& b) {
    double dot = -a.x * b.x - a.y * b.y - a.z * b.z;
    return acos(dot);
}

void render_ascii_donut(coord3D& light, coord3D& view) {

}

// quarternion or rotation matrix
coord3D rotate_point(coord3D& point) {

    return {};
}

// torus has x max 5 and y max 5
// donut along x and y axis
// theta is around circle -> r
// phi is around donut -> R
int main() {

    return 0;
}