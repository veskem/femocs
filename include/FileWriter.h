/*
 * FileWriter.h
 *
 *  Created on: 26.11.2018
 *      Author: veske
 */

#ifndef FILEWRITER_H_
#define FILEWRITER_H_

#include <string>

using namespace std;
namespace femocs {

/** Some flags to control the file output */
enum FileIO {
    force     = 1 << 0,  ///< omit the control for last write time
    no_update = 1 << 1,  ///< prevent update last_write_time variable after successful write
    append    = 1 << 2   ///< force appending to already excisting file
};

/** General class for holding common routines for writing files */
class FileWriter {
public:
    FileWriter();
    virtual ~FileWriter() {}

    /** By reading the file extension, pick the routine to write the data into file.
     * @param file   path to output file
     * @param force  omit the control for last write time
     * @param update update last write time variable after successful write
     * @param append force appending to already excisting file
     */
    void write(const string &file, unsigned int flags=0);

    /** Size of the data vector */
    virtual int size() const { return 0; }

    /** Determine if enough time has passed since the last file write */
    bool write_time() const;

protected:
    bool first_line(ofstream &out) const;

    /** Check if given file type is implemented for given class */
    virtual bool valid_extension(const string &extension) const = 0;

    // General routines for writing file with given extension
    virtual void write_xyz(ofstream &out) const;
    virtual void write_vtk(ofstream &out) const;
    virtual void write_msh(ofstream &out) const;
    virtual void write_restart(ofstream &out) const;
    virtual void write_bin(ofstream &out) const;
    virtual void write_dat(ofstream &out) const;
    virtual void write_ckx(ofstream &out) const;

    /** Label for restart data */
    virtual string get_restart_label() const;

    /* Routines to write data in VTK format */
    virtual void write_vtk_points_and_cells(ofstream &out) const {}
    virtual void write_vtk_point_data(ofstream &out) const {}
    virtual void write_vtk_cell_data(ofstream &out) const {}

private:
    double last_write_time;  ///< Last time file was written [fs]
};

} // namespace femocs

#endif /* FILEWRITER_H_ */
