#include<vector>
#include<cmath>
#include<chrono>
#include<iostream>
#include<string>
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

    coord3D(double x, double y, double z) : x(x), y(y), z(z) {}

    bool operator==(const coord3D other) const {
        const double tolerance = 1e-2;
        return abs(x - other.x) < tolerance && abs(y - other.y) < tolerance && abs(z - other.z) < tolerance;
    }

    coord3D operator-(const coord3D other) const {
        return {x - other.x, y - other.y, z - other.z};
    }

    coord3D operator+(const coord3D other) const {
        return {x + other.x, y + other.y, z + other.z};
    }

    double operator*(const coord3D other) const {
        return x * other.x + y * other.y + z * other.z;
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
// requires two normalized vectors
float angle_between_vectors(coord3D& a, coord3D& b) {
    double dot = -a.x * b.x - a.y * b.y - a.z * b.z;
    return acos(dot);
}

char light_to_char(double angle_degrees) {
    static const string light_chars = ".,-~:;=!*#$@";
    static const int num_light_chars = light_chars.size();
    static const double scale = 90.0 / 11.9;
    if (angle_degrees >= 90.0) return '.';
    int index = int((90.0 - angle_degrees) / scale);
    return light_chars[index];
}

void render_ascii_donut(coord3D& light, coord3D& view, double A, double B) {
    static const double x_min = -WINDOW_WIDTH_CHAR / (2 * RAYS_PER_UNIT_X);
    static const double x_max = -x_min;
    static const double y_min = -WINDOW_HEIGHT_CHAR / (2 * RAYS_PER_UNIT_Y);
    static const double y_max = -y_min;
    static const double x_inc = 1.0 / RAYS_PER_UNIT_X;
    static const double y_inc = 1.0 / RAYS_PER_UNIT_Y;
    static const double to_degrees = 180.0 / M_PI;

    const double cosA = cos(A);
    const double sinA = sin(A);
    const double nsinA = -sinA;
    const double cosB = cos(B);
    const double sinB = sin(B);
    const double nsinB = -sinB;

    coord3D light_rotated = {
        light.x * cosB - sinB * (light.y * cosA - light.z * sinA),
        light.x * sinB + cosB * (light.y * cosA - light.z * sinA),
        light.y * sinA + light.z * cosA
    };

    coord3D view_rotated = {
        view.x * cosB - sinB * (view.y * cosA - view.z * sinA),
        view.x * sinB + cosB * (view.y * cosA - view.z * sinA),
        view.y * sinA + view.z * cosA
    };

    coord3D view_dir_rotated = coord3D(0, 0, 0) - view_rotated;
    view_dir_rotated.normalize();

    char buffer[WINDOW_HEIGHT_CHAR * (WINDOW_WIDTH_CHAR + 1) + 1];
    int i = 0;
    for (double y = y_max; y >= y_min; y -= y_inc) {
        for (double x = x_min; x <= x_max; x += x_inc) {
            double inside_paren = y * cosA - view.z * sinA;
            coord3D view_loc_rotated = {x * cosB - sinB * inside_paren, x * sinB - cosB * inside_paren, y * sinA + view.z * cosA};

            coord3D view_intersect_point = vector_intersects_donut(view_loc_rotated, view_dir_rotated);
            bool vector_intersecting_donut = view_intersect_point.z != numeric_limits<double>::infinity();
            if (!vector_intersecting_donut) {
                buffer[i++] = ' ';
                continue;
            }

            coord3D light_dir_to_donut = view_intersect_point - light_rotated;
            light_dir_to_donut.normalize();
            coord3D light_intersect_point = vector_intersects_donut(light_rotated, light_dir_to_donut);
            if (light_intersect_point.x == numeric_limits<double>::infinity() || !(view_intersect_point == light_intersect_point)) {
                buffer[i++] = '.';
                continue;
            }

            coord3D normal = get_normal(view_intersect_point);
            double angle_degrees = angle_between_vectors(normal, light_dir_to_donut) * to_degrees;
            buffer[i++] = light_to_char(angle_degrees);
        } buffer[i++] = '\n';
    } buffer[i] = '\0';

    cout << buffer;
}

// torus has x max 5 and y max 5
// donut along x and y axis
// theta is around circle -> r
// phi is around donut -> R
int main() {

    coord3D light = {0, 5, -15};
    coord3D view = {0, 0, -10};
    double A = 0.0;
    double B = 0.0;

    render_ascii_donut(light, view, A, B);

    return 0;
}