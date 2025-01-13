#include "LBM.hpp"

int main()
{

    //=======================================================================
    //=======================================================================

    cout << "--- LBM D3Q19 for single component multiphase ---" << endl;
    cout << "--- Computational domain size ---" << endl;
    cout << "--- X = " << lx << " --- " << endl;
    cout << "--- Y = " << ly << " --- " << endl;
    cout << "--- Z = " << lz << " --- " << endl;
    read_parameters(t_max);
    read_obstacles(obst);
    init_density(obst, u_x, u_y, rho, ff);
    cout << "--- Loop Start ---" << endl;

    //=======================================================================
    //=======================================================================

    for (int t_step = 1; t_step <= t_max; ++t_step)
    {
        if (t_step % Nwri == 0 || t_step == 1)
        {
            cout << t_step << endl;
             write_results(obst, rho, upx, upy, upz, p, t_step);
            sum(rho);
        }
        stream(obst, ff);
        getuv(obst, u_x, u_y, u_z, rho, ff);
        calcu_upr(obst, u_x, u_y, u_z, Fx, Fy, Fz, Sx, Sy, Sz, rho, upx, upy, upz);
        calcu_Fxy(obst, rho, Fx, Fy, Fz, Sx, Sy, Sz, p);
        collision(tau, obst, u_x, u_y, u_z, rho, ff, Fx, Fy, Fz, Sx, Sy, Sz);
    }

    cout << "--- Loop End ---" << endl;

    return 0;
}
