/*
 * FileWriter.cpp
 *
 *  Created on: 26.11.2018
 *      Author: veske
 */

#include "FileWriter.h"
#include "Globals.h"
#include "Macros.h"

#include <fstream>

using namespace std;
namespace femocs {

FileWriter::FileWriter() :
        delta_time(0), last_write_time(0)
{}

FileWriter::FileWriter(const double dt) :
        delta_time(dt), last_write_time(-dt)
{}

void FileWriter::set_write_period(const double dt) {
    require(dt >= 0, "Invalid write period: " + d2s(dt));
    delta_time = dt;
}

bool FileWriter::not_write_time() const {
    return delta_time < 0 || (GLOBALS.TIME - last_write_time) < delta_time;
}

void FileWriter::write(const string &file_name, bool force) {
    if (!force && not_write_time())
        return;

    string ftype = get_file_type(file_name);
    if (!valid_extension(ftype)) {
        write_verbose_msg("Unimplemented file type: " + ftype);
        return;
    }

    ofstream outfile;
    if (ftype == "movie" || ftype == "restart" || ftype == "dat")
        outfile.open(file_name, ios_base::app);
    else
        outfile.open(file_name);

    require(outfile.is_open(), "Can't open a file " + file_name);

    if (ftype == "xyz" || ftype == "movie")
        write_xyz(outfile);
    else if (ftype == "vtk")
        write_vtk(outfile);
    else if (ftype == "msh")
        write_msh(outfile);
    else if (ftype == "dat")
        write_dat(outfile);
    else if (ftype == "bin")
        write_bin(outfile);
    else if (ftype == "restart")
        write_restart(outfile);
    else if (ftype == "ckx")
        write_ckx(outfile);

    outfile.close();
    last_write_time = GLOBALS.TIME;
}

void FileWriter::write_xyz(ofstream &out) const {
    int n_data = size();
    expect(n_data > 0, "Zero data detected!");
    out << n_data << "\n"
            << "Time=" << GLOBALS.TIME << ", Timestep=" << GLOBALS.TIMESTEP << ", " << class_name() << " ";

    out.setf(std::ios::scientific);
    out.precision(12);
}

string FileWriter::get_restart_label() const {
    return "";
}

void FileWriter::write_restart(ofstream &out) const {
    const string label = get_restart_label();
    if (label != "") {
        const int n_data = size();
        expect(n_data > 0, "Zero data detected!");
        out << "$" << label << "\n"
                << n_data << " " << GLOBALS.TIME << " " << GLOBALS.TIMESTEP << "\n";
    }

    write_bin(out);

    if (label != "") out << "\n$End"<< label << "\n";
}

void FileWriter::write_bin(ofstream &out) const {
    require(false, "FileWriter::write_bin not implemented!");
}

void FileWriter::write_ckx(ofstream &out) const {
    require(false, "FileWriter::write_ckx not implemented!");
}

void FileWriter::write_dat(ofstream &out) const {
    // In case of empty file, write data header
    if (out.tellp() == 0)
        out << "Time Timestep";

    // otherwise append data
    else
        out << d2s(GLOBALS.TIME) + " " + d2s(GLOBALS.TIMESTEP);
}

void FileWriter::write_vtk(ofstream &out) const {
    expect(size() > 0, "Zero data detected!");
    out.setf(std::ios::scientific);
    out.precision(12);

    out << "# vtk DataFile Version 3.0\n";
    out << "# " + class_name() + "\n";
    out << "ASCII\nDATASET UNSTRUCTURED_GRID\n";

    write_vtk_points_and_cells(out);
    write_vtk_point_data(out);
    write_vtk_cell_data(out);
}

void FileWriter::write_msh(ofstream &out) const {
    out << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n";
}

} // namespace femocs
