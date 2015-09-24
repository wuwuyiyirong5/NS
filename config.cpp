#include "ISOP2P1.h"

void ISOP2P1::config(std::string _config_file)
{
    std::string trash;
    std::ifstream input(_config_file.c_str());
    input >> trash >> trash >> mesh_file;
    input >> trash >> trash >> l_tol;
    input >> trash >> trash >> l_Euler_tol;
    input >> trash >> trash >> n_tol;
    input >> trash >> trash >> viscosity;
    input >> trash >> trash >> t0;
    input >> trash >> trash >> t1;
    input >> trash >> trash >> t;
    input >> trash >> trash >> dt;
    input >> trash >> trash >> CFL;
    input >> trash >> trash >> n_method;
    input >> trash >> trash >> time_step_control;
    input >> trash >> trash >> Stokes_init;
    input >> trash >> trash >> NS_init;
    input >> trash >> trash >> scheme;
    input >> trash >> trash >> body_force;
    input >> trash >> trash >> angle;
    input.close();
    std::cout << "mesh_file = " << mesh_file << std::endl;
    std::cout <<  "l_tol = " << l_tol << std::endl;
    std::cout <<  "l_Euler_tol = " << l_Euler_tol << std::endl;
    std::cout <<  "n_tol = " << n_tol << std::endl;
    std::cout <<  "viscosity = " << viscosity << std::endl;
    std::cout <<  "t0 = " << t0 << std::endl;
    std::cout <<  "t1 = " << t1 << std::endl;
    std::cout <<  "t = " << t << std::endl;
    std::cout <<  "dt = " << dt << std::endl;
    std::cout <<  "CFL = " << CFL << std::endl;
    std::cout <<  "n_method = " << n_method << std::endl;
    std::cout <<  "time_step_control = " << time_step_control << std::endl;
    std::cout <<  "Stokes_init = " << Stokes_init << std::endl;
    std::cout <<  "NS_init = " << NS_init << std::endl;
    std::cout <<  "scheme = " << scheme << std::endl;
    std::cout <<  "body_force = " << body_force << std::endl;
    std::cout <<  "angle = " << angle << std::endl;
    std::cout << "Check the input values." << std::endl;
    getchar();
};

