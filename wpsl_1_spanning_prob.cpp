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

int percolator(int cut, const vector<ROW>& rows, bool save_files) {
    int total_cells = 1 + cut;

    ofstream tout;
    ofstream chosenFile;
    if (save_files) {
        tout.open("percolation_n.txt");
        chosenFile.open("chosen_cells.txt");
    }

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

        if (save_files) {
            tout << i << "\t";
            for (size_t j = 0; j < nn[i].size(); ++j) {
                tout << nn[i][j];
                if (j != nn[i].size() - 1) tout << ",";
            }
            tout << "\n";
        }
    }

    vector<int> order(total_cells);
    iota(order.begin(), order.end(), 0);
    shuffle(order.begin(), order.end(), rng);

    int spanning_cluster_size = 0;

    for (int t = 0; t < total_cells; ++t) {
        int s = order[t];
        occupied[s] = true;
        parent[s] = -1;

        if (save_files) {
            chosenFile << rows[s].id << endl;
            for (const auto& point : rows[s].points) {
                double y_min, y_max, x_min, x_max;
                tie(y_min, y_max, x_min, x_max) = point;
                chosenFile << x_min << " " << y_min << "\n" << x_max << " " << y_min << "\n"
                           << x_max << " " << y_max << "\n" << x_min << " " << y_max << "\n\n";
            }
        }

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

    if (save_files) {
        tout.close();
        chosenFile.close();
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
    rows.reserve(cut + 2);
    rows.push_back({0, {make_tuple(0.0, L, 0.0, L)}, {}, {}, {}, {}, L * L});

    for (int step = 0; step < cut; ++step) {
        int selected_index = cellno(my_random(), my_random(), rows);
        ROW saved_cell = rows[selected_index];
        int selected_id = saved_cell.id;

        auto old_point = saved_cell.points[0];
        double y_min = get<0>(old_point);
        double y_max = get<1>(old_point);
        double x_min = get<2>(old_point);
        double x_max = get<3>(old_point);

        for (int up_id : saved_cell.up_neighbors) {
            auto& dn = rows[up_id].down_neighbors;
            dn.erase(remove(dn.begin(), dn.end(), selected_id), dn.end());
        }
        for (int down_id : saved_cell.down_neighbors) {
            auto& upn = rows[down_id].up_neighbors;
            upn.erase(remove(upn.begin(), upn.end(), selected_id), upn.end());
        }
        for (int left_id : saved_cell.left_neighbors) {
            auto& rn = rows[left_id].right_neighbors;
            rn.erase(remove(rn.begin(), rn.end(), selected_id), rn.end());
        }
        for (int right_id : saved_cell.right_neighbors) {
            auto& ln = rows[right_id].left_neighbors;
            ln.erase(remove(ln.begin(), ln.end(), selected_id), ln.end());
        }

        bool vertical = (my_random() < 0.5);

        if (vertical) {
            double cut_x = x_min + (x_max - x_min) * my_random();
            double area_left  = (y_max - y_min) * (cut_x - x_min);
            double area_right = (y_max - y_min) * (x_max - cut_x);
            int new_id = rows.size();

            rows[selected_index].points[0] = make_tuple(y_min, y_max, x_min, cut_x);
            rows[selected_index].area = area_left;
            rows[selected_index].left_neighbors.clear();
            rows[selected_index].right_neighbors = {new_id};
            rows[selected_index].up_neighbors.clear();
            rows[selected_index].down_neighbors.clear();

            ROW right_cell;
            right_cell.id = new_id;
            right_cell.points = {make_tuple(y_min, y_max, cut_x, x_max)};
            right_cell.area = area_right;
            right_cell.left_neighbors = {selected_id};
            rows.push_back(right_cell);

            vector<int> new_cells = {selected_id, new_id};
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
        } else {
            double cut_y = y_min + (y_max - y_min) * my_random();
            double area_bottom = (x_max - x_min) * (cut_y - y_min);
            double area_top    = (x_max - x_min) * (y_max - cut_y);
            int new_id = rows.size();

            rows[selected_index].points[0] = make_tuple(y_min, cut_y, x_min, x_max);
            rows[selected_index].area = area_bottom;
            rows[selected_index].down_neighbors.clear();
            rows[selected_index].up_neighbors = {new_id};
            rows[selected_index].left_neighbors.clear();
            rows[selected_index].right_neighbors.clear();

            ROW top_cell;
            top_cell.id = new_id;
            top_cell.points = {make_tuple(cut_y, y_max, x_min, x_max)};
            top_cell.area = area_top;
            top_cell.down_neighbors = {selected_id};
            rows.push_back(top_cell);

            vector<int> new_cells = {selected_id, new_id};
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
    }
    return rows;
}

int main() {
    int cuts = 5000;           
    int total_cells = 1 + cuts;  
    int iterations = 10000;   

    vector<double> percolation_stats(total_cells + 2, 0.0);

    for (int trial = 1; trial <= iterations; ++trial) {
        vector<ROW> rows = lattice(cuts);  
        
        bool save_files = (trial == iterations); 
        int cells_to_span = percolator(cuts, rows, save_files);  

        if (save_files) {
            ofstream latticeFile("lattice.txt");
            for (const auto& row : rows) {
                latticeFile << row.id << "\n";
                for (const auto& point : row.points) {
                    double y_min, y_max, x_min, x_max;
                    tie(y_min, y_max, x_min, x_max) = point;
                    latticeFile << x_min << " " << y_min << "\n" << x_max << " " << y_min << "\n"
                                << x_max << " " << y_max << "\n" << x_min << " " << y_max << "\n"
                                << x_min << " " << y_min << "\n";
                }
                latticeFile << "\n";
            }
            latticeFile.close();
        }

        if (cells_to_span > 0 && cells_to_span <= total_cells) {
            for (int g = cells_to_span; g <= total_cells; ++g) {
                percolation_stats[g] += 1.0;
            }
        }

        if (trial % max(1, iterations / 10) == 0) {
            cout << "Completed " << trial << " / " << iterations << endl;
        }
    }

    ofstream fout("percolation_results8000.txt");
    fout << "FractionOccupied Frequency" << endl;
    for (int i = 1; i <= total_cells; ++i) {
        fout << (double(i) / total_cells) << " " << (percolation_stats[i] / iterations) << "\n";
    }
    fout.close();

    return 0;
}
