#include <bits/stdc++.h>
using namespace std;

const double L = 1.0;
const double EPS_COORD = 1e-9;

struct ROW {
    int id;
    vector<tuple<double, double, double, double>> points;
    vector<int> up_neighbors;    
    vector<int> down_neighbors;
    vector<int> left_neighbors;
    vector<int> right_neighbors;
    double area; 
};

static thread_local std::mt19937_64 rng((unsigned)chrono::high_resolution_clock::now().time_since_epoch().count());
static inline double my_random() {
    uniform_real_distribution<double> d(0.0, 1.0);
    return d(rng);
}

int neighbor_check(const vector<ROW>& rows, int p, int l) {
    if (p < 0 || p >= (int)rows.size() || l < 0 || l >= (int)rows.size()) return 0;
    const ROW& row_p = rows[p];
    const ROW& row_l = rows[l];

    double a = get<0>(row_p.points[0]);
    double b = get<1>(row_p.points[0]);
    double c = get<2>(row_p.points[0]);
    double d = get<3>(row_p.points[0]);

    double a1 = get<0>(row_l.points[0]);
    double b1 = get<1>(row_l.points[0]);
    double c1 = get<2>(row_l.points[0]);
    double d1 = get<3>(row_l.points[0]);

    bool x_overlap = (c < d1 - EPS_COORD) && (c1 < d - EPS_COORD);
    bool y_overlap = (a < b1 - EPS_COORD) && (a1 < b - EPS_COORD);

    if ((fabs(b - a1) < EPS_COORD && x_overlap) ||
        (fabs(a - b1) < EPS_COORD && x_overlap) ||
        (fabs(c - d1) < EPS_COORD && y_overlap) ||
        (fabs(d - c1) < EPS_COORD && y_overlap)) {
        return 1;
    }
    return 0;
}

int findroot(int x, vector<int>& parent) {
    if (parent[x] < 0) return x;
    return parent[x] = findroot(parent[x], parent);
}

void unite(int a, int b, vector<int>& parent, vector<int>& boundary_mask) {
    a = findroot(a, parent);
    b = findroot(b, parent);
    if (a == b) return;

    if (parent[a] < parent[b]) {  
        parent[a] += parent[b];
        parent[b] = a;
        boundary_mask[a] |= boundary_mask[b];
    } else {
        parent[b] += parent[a];
        parent[a] = b;
        boundary_mask[b] |= boundary_mask[a];
    }
}

int percolator(int cut, const vector<ROW>& rows) {
    int total_cells = 1 + 3 * cut;

    vector<int> parent(total_cells, 0);
    vector<int> boundary_mask(total_cells, 0);
    vector<bool> occupied(total_cells, false);
    vector<vector<int>> nn(total_cells);

    const int TOUCH_TOP    = 1 << 0;
    const int TOUCH_BOTTOM = 1 << 1;
    const int TOUCH_LEFT   = 1 << 2;
    const int TOUCH_RIGHT  = 1 << 3;

    for (int i = 0; i < total_cells; ++i) {
        const ROW& cell = rows[i];
        nn[i].reserve(cell.up_neighbors.size() + cell.down_neighbors.size() + cell.left_neighbors.size() + cell.right_neighbors.size());
        for (int n : cell.up_neighbors)    nn[i].push_back(n);
        for (int n : cell.down_neighbors)  nn[i].push_back(n);
        for (int n : cell.left_neighbors)  nn[i].push_back(n);
        for (int n : cell.right_neighbors) nn[i].push_back(n);
    }

    vector<int> order(total_cells);
    iota(order.begin(), order.end(), 0);
    shuffle(order.begin(), order.end(), rng);

    int spanning_cluster_size = 0;

    for (int t = 0; t < total_cells; ++t) {
        int s = order[t];
        occupied[s] = true;
        parent[s] = -1;

        auto& p = rows[s].points[0];
        double y_min = get<0>(p);
        double y_max = get<1>(p);
        double x_min = get<2>(p);
        double x_max = get<3>(p);

        if (fabs(y_min - 0.0) < EPS_COORD)  boundary_mask[s] |= TOUCH_BOTTOM;
        if (fabs(y_max - L) < EPS_COORD)    boundary_mask[s] |= TOUCH_TOP;
        if (fabs(x_min - 0.0) < EPS_COORD)  boundary_mask[s] |= TOUCH_LEFT;
        if (fabs(x_max - L) < EPS_COORD)    boundary_mask[s] |= TOUCH_RIGHT;

        for (int nb : nn[s]) {
            if (occupied[nb]) {
                unite(s, nb, parent, boundary_mask);
            }
        }

        int root = findroot(s, parent);
        int mask = boundary_mask[root];

        if (((mask & TOUCH_LEFT) && (mask & TOUCH_RIGHT)) || 
            ((mask & TOUCH_TOP) && (mask & TOUCH_BOTTOM))) {
            spanning_cluster_size = t + 1;
            break;
        }
    }

    return spanning_cluster_size;
}

int cellno(double x, double y, const vector<ROW>& rows) {
    for (int i = 0; i < (int)rows.size(); i++) {
        const auto& bounds = rows[i].points[0];
        if ((get<2>(bounds) <= x) && (x <= get<3>(bounds)) && 
            (get<0>(bounds) <= y) && (y <= get<1>(bounds))) {
            return i;
        }
    }
    return (int)rows.size() - 1; 
}

vector<ROW> lattice(int cut) {
    vector<ROW> rows;
    rows.reserve(3 * cut + 2);
    rows.push_back({0, {make_tuple(0.0, L, 0.0, L)}, {}, {}, {}, {}, L * L});

    for (int e = 0; e < cut; ++e) {
        int selected_index = cellno(my_random(), my_random(), rows);
        ROW saved_cell = rows[selected_index];
        int selected_id = saved_cell.id;

        auto old_point = saved_cell.points[0];
        double y_min = get<0>(old_point);
        double y_max = get<1>(old_point);
        double x_min = get<2>(old_point);
        double x_max = get<3>(old_point);

        double cut_x = x_min + (x_max - x_min) * my_random();
        double cut_y = y_min + (y_max - y_min) * my_random();

        for (int up_id : saved_cell.up_neighbors) {
            auto& dn = rows[up_id].down_neighbors;
            dn.erase(remove(dn.begin(), dn.end(), selected_id), dn.end());
        }
        for (int down_id : saved_cell.down_neighbors) {
            auto& upn = rows[down_id].up_neighbors;
            upn.erase(remove(upn.begin(), upn.end(), selected_id), upn.end());
        }
        for (int left_id : saved_cell.left_neighbors) {
            auto& right_nbrs = rows[left_id].right_neighbors;
            right_nbrs.erase(remove(right_nbrs.begin(), right_nbrs.end(), selected_id), right_nbrs.end());
        }
        for (int right_id : saved_cell.right_neighbors) {
            auto& left_nbrs = rows[right_id].left_neighbors;
            left_nbrs.erase(remove(left_nbrs.begin(), left_nbrs.end(), selected_id), left_nbrs.end());
        }

        double area_bottom_left  = (cut_x - x_min) * (cut_y - y_min);
        double area_bottom_right = (x_max - cut_x) * (cut_y - y_min);
        double area_top_left     = (cut_x - x_min) * (y_max - cut_y);
        double area_top_right    = (x_max - cut_x) * (y_max - cut_y);

        rows[selected_index].points[0] = make_tuple(y_min, cut_y, x_min, cut_x);
        rows[selected_index].area = area_bottom_left;
        rows[selected_index].right_neighbors = {3 * e + 1};
        rows[selected_index].up_neighbors = {3 * e + 3};
        rows[selected_index].left_neighbors.clear();
        rows[selected_index].down_neighbors.clear();

        ROW bottom_right;
        bottom_right.id = 3 * e + 1;
        bottom_right.points = {make_tuple(y_min, cut_y, cut_x, x_max)};
        bottom_right.area = area_bottom_right;
        bottom_right.left_neighbors = {selected_id};
        bottom_right.up_neighbors = {3 * e + 2};

        ROW top_right;
        top_right.id = 3 * e + 2;
        top_right.points = {make_tuple(cut_y, y_max, cut_x, x_max)};
        top_right.area = area_top_right;
        top_right.left_neighbors = {3 * e + 3};
        top_right.down_neighbors = {3 * e + 1};

        ROW top_left;
        top_left.id = 3 * e + 3;
        top_left.points = {make_tuple(cut_y, y_max, x_min, cut_x)};
        top_left.area = area_top_left;
        top_left.right_neighbors = {3 * e + 2};
        top_left.down_neighbors = {selected_id};

        rows.push_back(bottom_right);
        rows.push_back(top_right);
        rows.push_back(top_left);

        vector<int> new_cells = {selected_id, 3 * e + 1, 3 * e + 2, 3 * e + 3};
        
        for (int cell_index : new_cells) {
            for (int left_id : saved_cell.left_neighbors) {
                if (neighbor_check(rows, cell_index, left_id)) {
                    rows[cell_index].left_neighbors.push_back(left_id);
                    rows[left_id].right_neighbors.push_back(cell_index);
                }
            }
            for (int right_id : saved_cell.right_neighbors) {
                if (neighbor_check(rows, cell_index, right_id)) {
                    rows[cell_index].right_neighbors.push_back(right_id);
                    rows[right_id].left_neighbors.push_back(cell_index);
                }
            }
            for (int up_id : saved_cell.up_neighbors) {
                if (neighbor_check(rows, cell_index, up_id)) {
                    rows[cell_index].up_neighbors.push_back(up_id);
                    rows[up_id].down_neighbors.push_back(cell_index);
                }
            }
            for (int down_id : saved_cell.down_neighbors) {
                if (neighbor_check(rows, cell_index, down_id)) {
                    rows[cell_index].down_neighbors.push_back(down_id);
                    rows[down_id].up_neighbors.push_back(cell_index);
                }
            }
        }
    }
    return rows;
}

int main() {
    int cuts = 1000;          
    int total_cells = 1 + 3 * cuts;  
    int iterations = 1000;   

    vector<double> percolation_stats(total_cells + 2, 0.0);

    for (int trial = 1; trial <= iterations; ++trial) {
        vector<ROW> rows = lattice(cuts);  
        int cells_to_span = percolator(cuts, rows);  

        if (cells_to_span > 0 && cells_to_span <= total_cells) {
            for (int g = cells_to_span; g <= total_cells; ++g) {
                percolation_stats[g] += 1.0;
            }
        }

        if (trial % max(1, iterations / 10) == 0) {
            cout << "Completed " << trial << " / " << iterations << endl;
        }
    }

    string filename = "percolation_results.txt";
    ofstream fout(filename);
    fout << "FractionOccupied SpanningProbability" << endl;
    for (int i = 1; i <= total_cells; ++i) {
        fout << (double(i) / total_cells) << " " << (percolation_stats[i] / iterations) << "\n";
    }
    fout.close();

    cout << "\nPercolation Statistics Summary:" << endl;
    cout << "Total trials: " << iterations << endl;
    cout << "Number of cuts: " << cuts << endl;
    cout << "Total cells created: " << total_cells << endl;
    cout << "Results saved to " << filename << endl;

    return 0;
}
