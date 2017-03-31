#include "Femocs_wrap.h"
#include <stdio.h>
#include <stdlib.h>

/* Code to simulate Femocs call from Kimocs */
int femocs_interface(int N_atoms, double *x, double *y, double *z, double E_appl,
                                            double *Ex, double *Ey, double *Ez)
{

    int success, tot_success;
    FEMOCS *fem;
    double *Enorm; 
    
    /* Create the femocs object */
    fem = create_femocs("input/md.in");
    tot_success = 0;
    
    /* Import the atoms to femocs */
    femocs_import_file(fem, &success, "");
    tot_success += success;
    

    /* Run Laplace solver */
    femocs_run(fem, &success, E_appl, "");
    tot_success += success;
    
    Enorm = malloc(1000 * sizeof(double));
    
    /* Export electric field on atoms */
    femocs_export_elfield(fem, &success, 1000, Ex, Ey, Ez, Enorm);
    tot_success += success;
    free(Enorm);
    
    return tot_success;  
}



int main() {
    const int n_atoms = 1;
    double *x = NULL, *y = NULL, *z = NULL, Ex[2000], Ey[2000], Ez[2000];
    const double Eappl = 0.2;

    int success = femocs_interface(n_atoms, x, y, z, Eappl, Ex, Ey, Ez);
    
    if (success) 
        printf("Error occured. Success  =  %d\n",success);
    else
        printf("Femocs run successfully\n");
    
    return 0;
}
