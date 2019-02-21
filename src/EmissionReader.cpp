/*
 * EmissionReader.cpp
 *
 *  Created on: 27.9.2018
 *      Author: kyritsak, Veske
 */

#include "EmissionReader.h"

using namespace std;
namespace femocs {

EmissionReader::EmissionReader(const FieldReader *fr, const HeatReader *hr, Interpolator* i) :
        fields(fr), heat(hr), mesh(NULL), phis_on_line(i)
{
    const int interpolation_rank = 3;
    phis_on_line.set_preferences(false, 3, interpolation_rank);
}

void EmissionReader::initialize(const TetgenMesh* m, bool reinit) {
    mesh = m;

    int n_nodes = fields->size();
    require(n_nodes > 0, "EmissionReader can't use empty fields!");

    // deallocate and allocate currents data
    current_densities.resize(n_nodes);
    nottingham.resize(n_nodes);
    currents.resize(n_nodes);
    thetas_SC.resize(n_nodes);
    markers.resize(n_nodes);

    //deallocate and allocate lines
    rline.resize(n_lines);
    Vline.resize(n_lines);

    //Initialise data
    global_data.Jmax = 0.;
    global_data.Frep = 0.;
    global_data.Jrep = 0.;
    if (reinit) global_data.multiplier = 1.;

    stats.N_calls = 0;
    stats.I_tot.resize(0);
    stats.Jrep.resize(0);
    stats.Jmax.resize(0);
    stats.Frep.resize(0);
    stats.Fmax.resize(0);
}

void EmissionReader::emission_line(const Point3& point, const Vec3& direction, const double rmax) {
    const double rmin = 1.e-5 * rmax;

    // calculate potentials on line starting from face centroid
    // and moving in direction of its norm
    phis_on_line.reserve(n_lines);
    for (int i = 0; i < n_lines; i++){
        rline[i] = rmin + ((rmax - rmin) * i) / (n_lines - 1);
        phis_on_line.append(point + direction * rline[i]);
    }
    phis_on_line.calc_interpolation();

    for (int i = 0; i < n_lines; i++){
        Vline[i] = global_data.multiplier * phis_on_line.get_potential(i);
        rline[i] *= nm_per_angstrom;
    }

    for (int i = 0; i < n_lines; i++){
        Vline[i] -= Vline[0];
        rline[i] -= rline[0];
    }

    // Check data condition (should be monotonous)
    for (int i = 1; i < n_lines; ++i) { // go through points
        if (Vline[i] < Vline[i-1]) { // if decreasing at a point
            double dVdx = 0.0;
            int j;
            for (j = i + 1; j < n_lines; ++j) {
                if (Vline[j] > Vline[i-1]) {
                    dVdx = (Vline[j] - Vline[i-1]) / (rline[j] - rline[i-1]);
                    break;
                }
            }

            if (dVdx == 0.0) {
                if (i > 1)
                    dVdx = (Vline[i-1] - Vline[i-2]) / (rline[i-1] - rline[i-2]);
                else
                    write_verbose_msg("Non-monotonous Vline could not be recovered at i = "
                            + d2s(i));
            }
            for (int k = 0; k <= j; ++k)
                Vline[k] =  Vline[i-1] + (rline[k] - rline[i-1]) * dVdx;
        }
    }
}

void EmissionReader::calculate_globals() {
    global_data.I_tot = 0;
    global_data.I_eff = 0;
    global_data.area = 0;

    if (is_effective.size() != current_densities.size()){
        is_effective.resize(current_densities.size());
        std::fill (is_effective.begin(), is_effective.end(), true);
    }

    double Fsum = 0.;
    for (unsigned int i = 0; i < currents.size(); ++i){ // go through face centroids
        int tri = mesh->quads.to_tri(abs(fields->get_marker(i)));
        // quadrangle area is 1/3 of corresponding triangle area
        double face_area = mesh->tris.get_area(tri) / 3.;
        double Floc = fields->get_elfield_norm(i) * thetas_SC[i] * global_data.multiplier;

        currents[i] = face_area * current_densities[i];
        global_data.I_tot += currents[i];

        if (is_effective[i]){ //if point eligible
            global_data.area += face_area; // increase total area
            global_data.I_eff += currents[i]; // increase total current
            Fsum += Floc * face_area;
        }
    }

    global_data.Jrep = global_data.I_eff / global_data.area;
    global_data.Frep = Fsum / global_data.area;

    stats.N_calls++;
    stats.I_tot.push_back(global_data.I_tot);
    stats.Jrep.push_back(global_data.Jrep);
    stats.Jmax.push_back(global_data.Jmax);
    stats.Fmax.push_back(global_data.Fmax);
    stats.Frep.push_back(global_data.Frep);

}

void EmissionReader::calc_effective_region(double threshold, string mode) {
    is_effective.resize(current_densities.size());

    if (mode == "current"){
        for (int i = 0; i < current_densities.size(); ++i)
            is_effective[i] = current_densities[i] > global_data.Jmax * threshold;
    } else if(mode == "field") {
        for (int i = 0; i < current_densities.size(); ++i)
            is_effective[i] = fields->get_elfield_norm(i) > global_data.Fmax * threshold;
    } else {
        require(false, "calc_effective_region called with wrong mode");
    }
}

int EmissionReader::calc_emission(const Config::Emission &conf, double Veff,
        bool update_eff_region)
{
    constexpr int J_max_error = -10;  // error code of current density is outside the limits
    constexpr int sign_error = -11;   // error code of field in wrong direction
    constexpr int fitting_error = -2; // error code of model fitting failed
    const int n_faces = fields->size();

    struct emission gt;
    gt.W = conf.work_function;    // set workfuntion, must be set in conf. script
    gt.R = 1000.0;       // radius of curvature (overrided by femocs potential distribution)
    gt.gamma = 10;       // enhancement factor (overrided by femocs potential distribution)
    double F, J;         // Local field and current density in femocs units
    bool fitting_failed = false;  // flag of non-fatal error

    global_data.Fmax = 0;
    global_data.Jmax = 0;

    for (int i = 0; i < n_faces; ++i) {
        int quad = abs(fields->get_marker(i));
        int tri = mesh->quads.to_tri(quad);
        Vec3 normal = mesh->tris.get_norm(tri);
        double elfield = fields->get_elfield_norm(i);

        // avoid using dot product as field norm, as field and normal are not precisely perpendicular
        // check, however, that the field and normal are in opposite direction
        if (normal.dotProduct(fields->get_elfield(i)) > 0)
            return sign_error;

        F = global_data.multiplier * elfield;
        gt.mode = 0;
        gt.F = angstrom_per_nm * F;
        gt.Temp = heat->get_temperature(i);
        markers[i] = 0; // marker==0: no full calculation

        // Full calculation with line only for high field points
        if (F > 0.6 * global_data.Fmax && !conf.blunt) {
            // get emission line data
            emission_line(fields->get_point(i), normal, 1.6 * conf.work_function / F);

            gt.Nr = n_lines;
            gt.xr = &rline[0];
            gt.Vr = &Vline[0];
            gt.mode = -21;  // set mode to potential input data
            markers[i] = 1; // marker==1: emission calculated with line
        }

        gt.approx = 0; // simple GTF approximation
        if (Veff <= 0)
            cur_dens_c(&gt); // calculate emission
        else
            cur_dens_SC(&gt, Veff);

        J = gt.Jem * nm2_per_angstrom2; // current density in femocs units
        if (J > conf.J_max)
            return J_max_error;

        // If J is worth it, calculate with full energy integration
        if ((J > 0.1 * global_data.Jmax || gt.ierr) && !conf.cold) {
            gt.approx = 1;
            if (Veff <= 0)
                cur_dens_c(&gt); // calculate emission
            else {
                gt.voltage = Veff;
                cur_dens_SC(&gt, Veff);
            }

            if (gt.ierr != 0) {
                if (gt.ierr == fitting_error)
                    fitting_failed = true;
                else
                    return gt.ierr;
            }

            J = gt.Jem * nm2_per_angstrom2;
            markers[i] = 2;
        }

        current_densities[i] = J;
        nottingham[i] = nm2_per_angstrom2 * gt.heat;
        thetas_SC[i] = gt.theta;

        global_data.Fmax = max(global_data.Fmax, F * gt.theta);
        global_data.Jmax = max(global_data.Jmax, J); // output data
    }

    if (update_eff_region)
        calc_effective_region(0.9, "field");

    calculate_globals();

    if (fitting_failed)
        write_verbose_msg("Model fitting failed. Applied rough FN approximation.");
    return 0;
}

void EmissionReader::write_xyz(ofstream &out) const {
    // start xyz header
    FileWriter::write_xyz(out);

    // write Ovito header
    out << "properties=id:I:1:pos:R:3:marker:I:1:force:R:3:"
         << "ln(rho_norm):R:1:nottingham_heat:R:1\n";

    // write data
    const int n_atoms = fields->size();
    for (int i = 0; i < n_atoms; ++i)
        out << i << ' ' << fields->get_point(i) << ' ' << markers[i] << ' ' << fields->get_elfield(i)
            << ' ' << log(current_densities[i]) << ' ' << nottingham[i] << "\n";
}

void EmissionReader::write_global_data(ofstream &out, const bool first_line) const {
    // specify data header
    if (first_line) out << "time      Itot        Imean        I_fwhm        Area        Jrep"
                        << "         Frep         Jmax        Fmax         multiplier\n";

    else {
        double I_mean = 0.;
        for (auto x : stats.I_tot)
            I_mean += x;

        // specify data
        out << fixed << setprecision(2) << GLOBALS.TIME;
        out << scientific << setprecision(6)
                << " " << global_data.I_tot << " " << I_mean / stats.N_calls
                << " " << global_data.I_eff << " " << global_data.area
                << " " << global_data.Jrep << " " << global_data.Frep
                << " " << global_data.Jmax << " " << global_data.Fmax
                << " " << global_data.multiplier << "\n";
    }
}

void EmissionReader::write_stats(ofstream &out, const bool first_line) const {
    // specify data header
    if (first_line)
        out << "   Itot_mean     Itot_std     Jrep_mean    Jrep_std"
            << "    Jmax_mean    Jmax_std     Frep_mean    Frep_std"
            << "    Fmax_mean    Fmax_std\n";

    else {
        // specify data
        out << " " << stats.Itot_mean << " " << stats.Itot_std
            << " " << stats.Jrep_mean << " " << stats.Jrep_std
            << " " << stats.Jmax_mean << " " << stats.Jmax_std
            << " " << stats.Frep_mean << " " << stats.Frep_std
            << " " << stats.Fmax_mean << " " << stats.Fmax_std
            << "\n";
    }
}

void EmissionReader::write_dat(ofstream &out) const {
    write_global_data(out, first_line(out));
    write_global_data(out, false);
}

void EmissionReader::calc_global_stats(){
    //initialise statistics
    stats.Itot_mean = 0; stats.Itot_std = 0;
    stats.Jmax_mean = 0; stats.Jmax_std = 0;
    stats.Fmax_mean = 0; stats.Fmax_std = 0;
    stats.Jrep_mean = 0; stats.Jrep_std = 0;
    stats.Frep_mean = 0; stats.Frep_std = 0;

    //calculate mean values
    for (int i = 0; i < stats.N_calls; ++i){
        stats.Itot_mean += stats.I_tot[i] / stats.N_calls;
        stats.Jmax_mean += stats.Jmax[i] / stats.N_calls;
        stats.Jrep_mean += stats.Jrep[i] / stats.N_calls;
        stats.Fmax_mean += stats.Fmax[i] / stats.N_calls;
        stats.Frep_mean += stats.Frep[i] / stats.N_calls;
    }

    //calculate standard deviations
    for (int i = 0; i < stats.N_calls; ++i){
        stats.Itot_std += pow(stats.I_tot[i] - stats.Itot_mean, 2);
        stats.Jmax_std += pow(stats.Jmax[i] - stats.Jmax_mean, 2);
        stats.Jrep_std += pow(stats.Jrep[i] - stats.Jrep_mean, 2);
        stats.Fmax_std += pow(stats.Fmax[i] - stats.Fmax_mean, 2);
        stats.Frep_std += pow(stats.Frep[i] - stats.Frep_mean, 2);
    }
    stats.Itot_std = sqrt(stats.Itot_std / stats.N_calls);
    stats.Jmax_std = sqrt(stats.Jmax_std / stats.N_calls);
    stats.Jrep_std = sqrt(stats.Jrep_std / stats.N_calls);
    stats.Fmax_std = sqrt(stats.Fmax_std / stats.N_calls);
    stats.Frep_std = sqrt(stats.Frep_std / stats.N_calls);

    //re-initialise statistics
    stats.I_tot.resize(0);
    stats.Jmax.resize(0);
    stats.Jrep.resize(0);
    stats.Fmax.resize(0);
    stats.Frep.resize(0);
}

string EmissionReader::get_error_codes(vector<int> &errors) const {
    std::sort(errors.begin(), errors.end());

    string retval;
    int error_cntr = 0;
    int prev_error = errors[0];
    for (int e : errors) {
        if (e == prev_error)
            error_cntr++;
        else {
            retval += d2s(error_cntr) + "x of " + d2s(prev_error) + ", ";
            prev_error = e;
            error_cntr = 1;
        }
    }

    retval += d2s(error_cntr) + "x of " + d2s(prev_error);
    return retval;
}

}
