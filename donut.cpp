#include<vector>
#include<cmath>
#include<chrono>
#include<iostream>
#include<string>
#include<limits>
#include<thread>
#include<stdlib.h>

using namespace std;

// Can represent a 3D coordinate or a 3D vector
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

    coord3D rotate(const double (&matrix)[3][3]) const {
        return {x * matrix[0][0] + y * matrix[1][0] + z * matrix[2][0],
            x * matrix[0][1] + y * matrix[1][1] + z * matrix[2][1],
            x * matrix[0][2] + y * matrix[1][2] + z * matrix[2][2]
        };
    }


};

// DONUT DEFINITION ---------------------------------------------
const double R = 4.0;
const double r = 2.2;
// ((x^2 + y^2)^0.5 - R) + z^2 = r^2
// x = (R + r * cos(theta)) * cos(phi)
// y = (R + r * cos(theta)) * sin(phi)
// z = r * sin(theta)
// --------------------------------------------------------------




// MODIFY -------------------------------------------------------
const int FPS = 20;
const int WINDOW_WIDTH_CHAR = 61;
const int WINDOW_HEIGHT_CHAR = 31;
// light and view position to donut facing you at (0, 0, 0)
const coord3D light = {0, 5, -15};
const coord3D view = {0, 0, -10};
// rotataion factors of donut
const double inc_yaw = 6.0 * M_PI / 180.0;
const double inc_pitch = 2.5 * M_PI / 180.0;
const double inc_roll = 0.0 * M_PI / 180.0;
// --------------------------------------------------------------




// DO NOT MODIFY ------------------------------------------------
const int RAYS_PER_UNIT_X = 4;
const int RAYS_PER_UNIT_Y = 2;
const double WINDOW_HEIGHT = 2.5 * (R + r);
const double WINDOW_WIDTH = 2.5 * (R + r);
// --------------------------------------------------------------




// Determines if point is inside donut according to
// (x^2 + y^2 + z^2 + R^2 + r^2)^2 -4R^2 * (x^2 + y^2) < 0
bool point_inside_donut(coord3D p) {
    static double _4R2 = 4 * R * R;
    static double R2_minus_r2 = R * R - r * r;
    double x2_plus_y2 = p.x * p.x + p.y * p.y;
    double inside_sq = (x2_plus_y2 + p.z * p.z + R2_minus_r2);
    return inside_sq * inside_sq - _4R2 * x2_plus_y2 <= 0.0;
}

// Find the point where the vector intersects the donut, if it does
coord3D vector_intersects_donut(coord3D& loc, coord3D& dir) {
    /*
     * This can be modelled by many different equations and
     * systems of equations, but none that I experimented with
     * can be solved closed form. This means that numerical
     * methods and approximation would have to be used. However,
     * I didn't want to use outside libraries in this code.
     * Instead I literally trace the vectors ever so slightly
     * until it hits the donuts to find the intersection point.
     */
    static const double jump = 0.1;
    static const double jump_jump = 0.01;
    static const double R_plus_r = R + r;
    static const double _R_minus_r = -R_plus_r;
    static const double _r = -r;

    double t1 = numeric_limits<double>::infinity(), t2 = numeric_limits<double>::infinity();

    // find the mininum 2 of three doubles and assigns the most min
    // to t1 and second most min to t2
    static auto min2 = [](double t, double& t1, double& t2) {
        if (t < t1) {
            t2 = t1;
            t1 = t;
        } else if (t < t2) {
            t2 = t;
        }
    };

    double t;
    double x, y, z;
    /*
     * Instead of tracing the point all the way from the view location,
     * optimze by only tracing the vector from a rectangular prism that
     * contains the donut. The rectangular prism goes from x: -R-r, R+r,
     * y: -R-r, R+r, z: -r, r. If our vector intersects this prism, it
     * can intersect at most two points. We want to get these two points
     * with the smaller one stored in t1, and the larger stored in t2
     */
    if (dir.x != 0) {
        t = (R_plus_r - loc.x) / dir.x;
        y = loc.y + dir.y * t;
        z = loc.z + dir.z * t;
        bool vector_intersects_positive_x_face = y <= R_plus_r && y >= _R_minus_r && z <= r && z >= _r;
        if (vector_intersects_positive_x_face) min2(t, t1, t2);
        t = (_R_minus_r - loc.x) / dir.x;
        y = loc.y + dir.y * t;
        z = loc.z + dir.z * t;
        bool vector_intersects_negative_x_face = y <= R_plus_r && y >= _R_minus_r && z <= r && z >= _r;
        if (vector_intersects_negative_x_face) min2(t, t1, t2);
    } if (dir.y != 0) {
        t = (R_plus_r - loc.y) / dir.y;
        x = loc.x + dir.x * t;
        z = loc.z + dir.z * t;
        bool vector_intersects_positive_y_face = x <= R_plus_r && x >= _R_minus_r && z <= r && z >= _r;
        if (vector_intersects_positive_y_face) min2(t, t1, t2);
        t = (_R_minus_r - loc.y) / dir.y;
        x = loc.x + dir.x * t;
        z = loc.z + dir.z * t;
        bool vector_intersects_negative_y_face = x <= R_plus_r && x >= _R_minus_r && z <= r && z >= _r;
        if (vector_intersects_negative_y_face) min2(t, t1, t2);
    } if (dir.z != 0) {
        t = (r - loc.z) / dir.z;
        x = loc.x + dir.x * t;
        y = loc.y + dir.y * t;
        bool vector_intersects_positive_z_face = x <= R_plus_r && x >= _R_minus_r && y <= R_plus_r && y >= _R_minus_r;
        if (vector_intersects_positive_z_face) min2(t, t1, t2);
        t = (_r - loc.z) / dir.z;
        x = loc.x + dir.x * t;
        y = loc.y + dir.y * t;
        bool vector_intersects_negative_z_face = x <= R_plus_r && x >= _R_minus_r && y <= R_plus_r && y >= _R_minus_r;
        if (vector_intersects_negative_z_face) min2(t, t1, t2);
    }

    // start the ray tracing at t1 and go until it reaches t2 or the donut
    t = t1;
    while (t < t2 && !point_inside_donut({loc.x + dir.x * t, loc.y + dir.y * t, loc.z + dir.z * t})) {
        t += jump;
    } 
    // if ray doesn';t hit donut, return inf vector
    if (t >= t2) {
        return {
            numeric_limits<double>::infinity(),
            numeric_limits<double>::infinity(),
            numeric_limits<double>::infinity()
        };
    }
    
    // if ray hits donut, backtrack to get more accurate estimate
    while (point_inside_donut({loc.x + dir.x * t, loc.y + dir.y * t, loc.z + dir.z * t})) {
        t -= jump_jump;
    }
    
    return {loc.x + dir.x * t, loc.y + dir.y * t, loc.z + dir.z * t};
}

// Gets the normal of a point on the donut
// https://web.cs.ucdavis.edu/~amenta/s12/findnorm.pdf
coord3D get_normal(coord3D& point) {
    // estimate theta and phi
    static double R2 = R * R;
    static double _2PI = 2 * M_PI;
    static double _3_over_2_PI = 3 * M_PI_2;
    double phi = atan2(point.y, point.x);
    // div can be slightly out of range [-1, 1];
    double div = point.z / r;
    if (div > 1.0) div = 1.0;
    else if (div < -1.0) div = -1.0;
    double theta = asin(div);

    // theta in quad I, II
    bool upper = point.z > 0;
    // theta in quad I, IV
    bool right = point.x * point.x + point.y * point.y > R2;
    if (point.z == 0.0 && right) theta = 0;
    else if (point.z == r) theta = M_PI_2;
    else if (point.z == 0.0 && !right) theta = M_PI;
    else if (point.z == -r) theta = _3_over_2_PI;
    else if (upper && right) theta = theta;
    else if (upper && !right) theta = M_PI - theta;
    else if (!upper && !right) theta = M_PI - theta;
    else theta = _2PI + theta;

    // get normal from theta and phi
    // https://web.cs.ucdavis.edu/~amenta/s12/findnorm.pdf
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

// Calculates the angle between vectors according to
// cos(theta) = a.dot(b) / (magnitude(a) * magnitude(b))
// REQUIRES TWO NORMALIZED VECTORS!!!!!!!!!!
float angle_between_vectors(coord3D& a, coord3D& b) {
    double dot = -a.x * b.x - a.y * b.y - a.z * b.z;
    return acos(dot);
}

// transform an angle to a character
char light_to_char(double angle_degrees) {
    static const string light_chars = ".,-~:;=!*#$@";

    // (num_light_chars - 0.1) to avoid number instability issues which push
    // scale out of [0, num_light_chars - 1] range
    static const double num_light_chars = light_chars.size();
    static const double scale = 90.0 / (num_light_chars - 0.1);

    if (angle_degrees >= 90.0) return '.';
    int index = int((90.0 - angle_degrees) / scale);
    return light_chars[index];
}

// Renders the view of the ascii donut for a view with the donut rotated
// at a given yaw, pitch, and roll
char* render_ascii_donut(double yaw, double pitch, double roll) {
    static const double x_min = -WINDOW_WIDTH_CHAR / (2 * RAYS_PER_UNIT_X);
    static const double x_max = -x_min;
    static const double y_min = -WINDOW_HEIGHT_CHAR / (2 * RAYS_PER_UNIT_Y);
    static const double y_max = -y_min;
    static const double x_inc = 1.0 / RAYS_PER_UNIT_X;
    static const double y_inc = 1.0 / RAYS_PER_UNIT_Y;
    static const double to_degrees = 180.0 / M_PI;
    static const int total_chars = WINDOW_HEIGHT_CHAR * (WINDOW_WIDTH_CHAR + 1) + 1;

    const double cosA = cos(yaw);
    const double sinA = sin(yaw);
    const double cosB = cos(pitch);
    const double sinB = sin(pitch);
    const double cosC = cos(roll);
    const double sinC = sin(roll);

    // Rotation matrix which will rotate the donut with the given yaw, ptich, and roll
    double rotation_matrix[3][3] = {
        {cosB * cosC, sinA * sinB * cosC - cosA * sinC, cosA * sinB * cosC + sinA * sinC},
        {cosB * sinC, sinA * sinB * sinC + cosA * cosC, cosA * sinB * sinC - sinA * cosC},
        {-sinB, sinA * cosB, cosA * cosB}
    };

    /*
     * Note, we do not actually rotate the donut itself, but instead we rotate the light
     * source and position of the viewer with this rotation matrix. This was done to make
     * other calculations easier and cheaper.
     */
    coord3D light_rotated = light.rotate(rotation_matrix);
    coord3D view_rotated = view.rotate(rotation_matrix);
    coord3D view_dir_rotated = coord3D(0, 0, 0) - view_rotated;
    view_dir_rotated.normalize();

    char* buffer = new char[total_chars];
    int i = 0;

    /*
     * We imagine that the view of the donut is actually an array of many parallel
     * beams from the direction of the view point. 
     */
    for (double y = y_max; y >= y_min; y -= y_inc) {
        for (double x = x_min; x <= x_max; x += x_inc) {
            // we take each of these view points and rotate it according to yaw, pitch, and roll
            coord3D view_loc_rotated = coord3D(x, y, view.z).rotate(rotation_matrix);

            // Does this view vector intersect the donut?
            coord3D view_intersect_point = vector_intersects_donut(view_loc_rotated, view_dir_rotated);
            bool vector_intersecting_donut = view_intersect_point.z != numeric_limits<double>::infinity();
            // If not print a space
            if (!vector_intersecting_donut) {
                buffer[i++] = ' ';
                continue;
            }

            // Does light reach that point on the donut, or is it obstructed by another
            // part of the donut?
            coord3D light_dir_to_donut = view_intersect_point - light_rotated;
            light_dir_to_donut.normalize();
            coord3D light_intersect_point = vector_intersects_donut(light_rotated, light_dir_to_donut);
            // If light does not reach this point, print '.'
            if (light_intersect_point.x == numeric_limits<double>::infinity() ||
                !(view_intersect_point == light_intersect_point)) {
                buffer[i++] = '.';
                continue;
            }

            // If the view vector reached the donut and the light hits that point,
            // what is the angle between the surface of the donut and the light vector
            coord3D normal = get_normal(view_intersect_point);
            double angle_degrees = angle_between_vectors(normal, light_dir_to_donut) * to_degrees;
            buffer[i++] = light_to_char(angle_degrees);
        } buffer[i++] = '\n';
    } buffer[i] = '\0'; // add null terminating character for string output

    return buffer;
}

// Get UNIX time
chrono::milliseconds get_time() {
    return chrono::duration_cast<chrono::milliseconds>(
        chrono::system_clock::now().time_since_epoch()
    );
}

// Animates the ascii donut at with modifiable yaw, pitch, and roll increments
int main() {
    static const auto ms_in_frame = chrono::milliseconds(1000 / FPS);
    static const double _2PI = 2.0 * M_PI;

    double yaw = 0.0;
    double pitch = 0.0;
    double roll = 0.0;
    while (true) {
        chrono::milliseconds start_frame_ms = get_time();
        char* donut = render_ascii_donut(yaw, pitch, roll);
        chrono::milliseconds end_render_ms = get_time();

        // sleep from after the frame is rendered until ms_in_frame ms after previous frame
        chrono::milliseconds render_duration_ms = end_render_ms - start_frame_ms;
        chrono::milliseconds sleep_duration_ms = ms_in_frame - render_duration_ms;
        if (sleep_duration_ms.count() > 0) this_thread::sleep_for(sleep_duration_ms);
        
        system("clear");
        cout << donut;

        yaw += inc_yaw;
        pitch += inc_pitch;
        roll += inc_roll;
        if (yaw > _2PI) yaw -= _2PI;
        if (pitch > _2PI) pitch -= _2PI;
        if (roll > _2PI) roll -= _2PI;
    }
    

    return 0;
}