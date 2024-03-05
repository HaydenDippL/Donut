#include<iostream>
#include<vector>
#include<cmath>
#include<unistd.h>

using namespace std;

const double R = 4.0;
const double r = 2.2;
const int WINDOW_WIDTH_CHAR = 61;
const int WINDOW_HEIGHT_CHAR = 31;

const int RAYS_PER_UNIT_X = 4;
const int RAYS_PER_UNIT_Y = 2;

const double WINDOW_HEIGHT = 2.5 * (R + r);
const double WINDOW_WIDTH = 2.5 * (R + r);

const double X_MAX = WINDOW_WIDTH_CHAR / (2 * RAYS_PER_UNIT_X);
const double X_MIN = -1 * X_MAX;
const double Y_MAX = WINDOW_HEIGHT_CHAR / (2 * RAYS_PER_UNIT_Y);
const double Y_MIN = -1 * Y_MAX;

struct coord3D {
    double x;
    double y;
    double z;
};



vector<double> calc_z(double x, double y) {
    double pos_d = sqrt(pow(x, 2) + pow(y, 2));
    double neg_d = -1 * pos_d;

    vector<double> ans;
    double base;
    if ((base = pow(r, 2) - pow(pos_d - R, 2)) >= 0) {
        ans.push_back(sqrt(base));
        ans.push_back(-sqrt(base));
    } if ((base = pow(r, 2) - pow(neg_d - R, 2)) >= 0) {
        ans.push_back(sqrt(base));
        ans.push_back(-sqrt(base));
    } return ans;
}

coord3D get_normal(coord3D coord) {
    double op = sqrt(coord.x * coord.x + coord.z * coord.z);
    double x = 2 * coord.x * (op - R) / op;
    double y = 2 * coord.y;
    double z = 2 * coord.z * (op - R) / op;
    return {x, y, z};
}

int main() {

    coord3D light_source = {8.0, -10.0, 15.0};
    vector<vector<double>> light_matrix(WINDOW_HEIGHT_CHAR, vector<double>(WINDOW_WIDTH_CHAR));

    double y = Y_MIN;
    for (int h = 0; h < WINDOW_HEIGHT_CHAR; ++h) {
        double x = X_MIN;
        for (int w = 0; w < WINDOW_WIDTH_CHAR; ++w) {
            vector<double> zs = calc_z(x, y);
            if (zs.size() == 0) {
                light_matrix[h][w] = 0;
                continue;
            }
            double z = numeric_limits<double>::lowest();
            for (const double potential_z : zs) z = max(z, potential_z);

            light_matrix[h][w] = calc_z(x, y).size() != 0;
            x += 1.0 / RAYS_PER_UNIT_X;
        } y += 1.0 / RAYS_PER_UNIT_Y;
    }

    for (int h = 0; h < WINDOW_HEIGHT_CHAR; ++h) {
        for (int w = 0; w < WINDOW_WIDTH_CHAR; ++w) {
            cout << (light_matrix[h][w] ? '#' : ' ');
        } cout << '\n';
    }

    for (int h = 0; h < WINDOW_HEIGHT_CHAR; ++h) {
        for (int w = 0; w < WINDOW_WIDTH_CHAR; ++w) {
            cout << (light_matrix[h][w]);
        } cout << '\n';
    }

    // render the donut with print statements

    return 0;
}