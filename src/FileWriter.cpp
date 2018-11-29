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

FileWriter::FileWriter() : last_write_time(-MODES.WRITE_PERIOD) {}

bool FileWriter::write_time() const {
    constexpr double epsilon = 1e-10;
    return MODES.WRITE_PERIOD >= 0 && (GLOBALS.TIME - last_write_time + epsilon) >= MODES.WRITE_PERIOD;
}

bool FileWriter::first_line(ofstream &out) const {
    return out.tellp() == 0;
}

void FileWriter::write(const string &file_name, unsigned int flags) {
    bool force = flags & FileIO::force;
    if (!force && !write_time())
        return;

    bool update = ! (flags & FileIO::no_update);
    bool append = flags & FileIO::append;

    string ftype = get_file_type(file_name);
    if (!valid_extension(ftype)) {
        write_verbose_msg("Unimplemented file type: " + ftype);
        return;
    }

    ofstream outfile;
    if (ftype == "movie" || ftype == "dat" || append)
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
    if (update) last_write_time = GLOBALS.TIME;
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
    // In case of empty file, write first data header
    if (first_line(out))
        out << "Time Timestep";

    // append data
    out << d2s(GLOBALS.TIME) + ' ' + d2s(GLOBALS.TIMESTEP);
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
