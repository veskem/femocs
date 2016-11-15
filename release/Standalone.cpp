/*
 * Standalone.cpp
 *
 *  Created on: 19.04.2016
 *      Author: veske
 */

#include "Femocs.h"
#include "Macros.h"

using namespace std;

//int main(int argc, char* argv[]) {
int main() {
    cout << "Standalone started!" << endl;
    femocs::Femocs femocs("");
    femocs.import_atoms("");
    femocs.run(0.1, "");

    return 0;
}
