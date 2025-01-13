//The task is to translate the code from Chapter 2.12 of the book "Multiphase Lattice Boltzmann Methods: Theory and Application" from Fortran to C.
//D3Q19-LBM-SCMP-SRT

#ifndef LBM_hpp
#define LBM_hpp

#include <iostream>     
#include <cmath>       
#include <string>       
#include <fstream>     
#include <filesystem>   
#include <cstdlib> 

using namespace std;

// This array defines which lattice positions are occupied by fluid nodes (obst=0)
// or solid nodes (obst=1)
int obst[lx][ly][lz];

// Velocity components
double u_x[lx][ly][lz], u_y[lx][ly][lz], u_z[lx][ly][lz];

// Pressure and density
double p[lx][ly][lz], rho[lx][ly][lz];

double psx[lx][ly][lz];

// The real fluid density
// which may differ from the velocity componets in the above; refer to SC model
double upx[lx][ly][lz], upy[lx][ly][lz], upz[lx][ly][lz];

// The force components: Fx, Fy, Fz for the interaction between fluid nodes.
// Sx, Sy, Sz are the interaction (components) between the fluid nodes and solid nodes
// ff is the distribution function
double ff[19][lx][ly][lz], Fx[lx][ly][lz], Fy[lx][ly][lz], Fz[lx][ly][lz];
double Sx[lx][ly][lz], Sy[lx][ly][lz], Sz[lx][ly][lz];

// TT0W is the value of T/T0; RHW and RLW are the coexisting densities
// in the sepcified T/T0.
// For initialization, \rho_l (lower density)
// and \rho_h (higher density) are supposed to be known.
//The below data define the D3Q19 velocity model, xc(ex), yc(ey), zc(ez)
//are the components of e_{ix}, e_{iy}, and e_{iz}, respectively.
double xc[19] = { 0,1,-1,0,0,0,0,1,1,-1,-1,1,-1,1,-1,0,0,0,0 };
double yc[19] = { 0,0,0,1,-1,0,0,1,-1,1,-1,0,0,0,0,1,1,-1,-1 };
double zc[19] = { 0,0,0,0,0,1,-1,0,0,0,0,1,1,-1,-1,1,-1,1,-1 };
int ex[19] = { 0,1,-1,0,0,0,0,1,1,-1,-1,1,-1,1,-1,0,0,0,0 };
int ey[19] = { 0,0,0,1,-1,0,0,1,-1,1,-1,0,0,0,0,1,1,-1,-1 };
int ez[19] = { 0,0,0,0,0,1,-1,0,0,0,0,1,1,-1,-1,1,-1,1,-1 };
// This array gives the opposite direction for e_1, e_2, e_3, .....e_18
// It implements the simple bounce-back rule we use in the collision step
// for solid nodes (obst=1)
int opp[19] = { 0,2,1,4,3,6,5,10,9,8,7,14,13,12,11,18,17,16,15 };

//C-S EOS
//RHW and RLW are the coexisting densities in the corresponding sepcified T/T0.

double TT0W[12] = { 0.975, 0.95, 0.925, 0.9, 0.875, 0.85, 0.825, 0.8, 0.775, 0.75, 0.7, 0.65 };
double RHW[12] = { 0.16, 0.21, 0.23, 0.247, 0.265, 0.279, 0.29, 0.314, 0.30, 0.33, 0.36, 0.38 };
double RLW[12] = { 0.08, 0.067, 0.05, 0.0405, 0.038, 0.032, 0.025, 0.0245, 0.02, 0.015, 0.009, 0.006 };

//Speeds and weighting factors

double cc = 1.0;
double c_squ = cc * cc / 3.0;

//Weighting coefficient in the equilibrium distribution function

double t_k[19] = { 1.0 / 3.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0 };

//Please specify which temperature
//and corresponding \rho_h, \rho_l in above �data� are chosen.
//Initial T/T0, rho_h, and rho_l for the C-S EOS are listed in above �data� section

int k = 5;//important
double TT0 = TT0W[k];
double rho_h = RHW[k];
double rho_l = RLW[k];

//installization

void read_parameters(int& t_max) {
    double visc;
    //Initial radius of the droplet.
    //read(1,*) RR
    RR = 15;
    // \rho_w in calculation of fluid-wall interaction
    //read(1,*) rho_w
    rho_w = 0.12;
    //Relaxation parameter, which is related to viscosity
    //read(1,*) tau
    tau = 1.0;
    //Maximum iteration specified
    //read(1,*) t_max
    t_max = 10000;
    //Output data frequency (can be viewed with TECPLOT)
    //read(1,*) Nwri
    Nwri = 100;
    //close(1)
    visc = c_squ * (tau - 0.5);
    cout << "kinematic viscosity=" << visc << "lu^2/ts" << " tau=" << tau << endl;
}


void read_obstacles(int obst[lx][ly][lz])
{
    for (int z = 0; z < lz; z++)
    {
        for (int y = 0; y < ly; y++)
        {
            for (int x = 0; x < lx; x++)
            {
                obst[x][y][z] = 0;
                if (z == 0)
                {
                    obst[x][y][0] = 1;
                }

            }
        }
    }
}

void init_density(int obst[lx][ly][lz], double u_x[lx][ly][lz], double u_y[lx][ly][lz], double rho[lx][ly][lz], double ff[19][lx][ly][lz])
{
    double u_squ, u_n[19], fequi[19];
    for (int z = 0; z < lz; z++)
    {
        for (int y = 0; y < ly; y++)
        {
            for (int x = 0; x < lx; x++)
            {
                u_x[x][y][z] = 0;
                u_y[x][y][z] = 0;
                u_z[x][y][z] = 0;

                rho[x][y][z] = RLW[k];
                if (pow(x - lx / 2, 2) + pow(y - ly / 2, 2) + pow(z - 25, 2) < pow(RR, 2))
                {
                    rho[x][y][z] = RHW[k];
                }

            }
        }
    }
    for (int z = 0; z < lz; z++)
    {
        for (int y = 0; y < ly; y++)
        {
            for (int x = 0; x < lx; x++)
            {
                u_squ = u_x[x][y][z] * u_x[x][y][z] + u_y[x][y][z] * u_y[x][y][z] + u_z[x][y][z] * u_z[x][y][z];
                for (int k = 0; k < 19; ++k)
                {
                    u_n[k] = xc[k] * u_x[x][y][z] + yc[k] * u_y[x][y][z] + zc[k] * u_z[x][y][z];
                    fequi[k] = t_k[k] * rho[x][y][z] * (cc * u_n[k] / c_squ + (u_n[k] * cc) * (u_n[k] * cc) / (2 * c_squ * c_squ) - u_squ / (2 * c_squ)) + t_k[k] * rho[x][y][z];
                    ff[k][x][y][z] = fequi[k];
                }
            }
        }
    }
}

void calcu_Fxy(int obst[lx][ly][lz], double rho[lx][ly][lz], double Fx[lx][ly][lz], double Fy[lx][ly][lz], double Fz[lx][ly][lz], double Sx[lx][ly][lz], double Sy[lx][ly][lz], double Sz[lx][ly][lz], double p[lx][ly][lz])
{
    const double R = 1.0;
    const double b = 4.0;
    const double a = 1.0;
    const double Tc = 0.3773 * a / (b * R);
    const double TT = TT0W[k] * Tc;
    double G1 = -1.0;
    double G12;
    double lambda = 0.886;
    double Fx1[lx][ly][lz], Fy1[lx][ly][lz], Fz1[lx][ly][lz];


    for (int z = 0; z < lz; z++)
    {
        for (int y = 0; y < ly; y++)
        {
            for (int x = 0; x < lx; x++)
            {
                if (obst[x][y][z] == 0 && rho[x][y][z] != 0.0)
                {

                    p[x][y][z] = rho[x][y][z] * R * TT * (1 + b * rho[x][y][z] / 4.0 + pow(b * rho[x][y][z] / 4.0, 2) - pow(b * rho[x][y][z] / 4.0, 3)) / pow((1 - b * rho[x][y][z] / 4.0), 3) - a * rho[x][y][z] * rho[x][y][z];//CS EOS
                    psx[x][y][z] = sqrt(fabs(2.0 * (p[x][y][z] - rho[x][y][z] * c_squ) / (G1 * c_squ)));

                }
            }
        }
    }

    double psx_w = sqrt(fabs(2.0 * ((rho_w * R * TT * (1 + b * rho_w / 4.0 + pow(b * rho_w / 4.0, 2) - pow(b * rho_w / 4.0, 3)) / pow((1 - b * rho_w / 4.0), 3) - a * rho_w * rho_w) - rho_w * c_squ) / (G1 * c_squ)));

    for (int z = 0; z < lz; z++)
    {
        for (int y = 0; y < ly; y++)
        {
            for (int x = 0; x < lx; x++)
            {
                Fx[x][y][z] = 0.0;
                Fy[x][y][z] = 0.0;
                Fz[x][y][z] = 0.0;

                if (obst[x][y][z] == 0)
                {
                    double sum_x = 0.0, sum_y = 0.0, sum_z = 0.0;
                    for (int k = 0; k < 19; ++k)
                    {
                        int xp = x + ex[k];
                        int yp = y + ey[k];
                        int zp = z + ez[k];
                        if (xp < 0)
                        {
                            xp = lx - 1;
                        }
                        if (xp >= lx)
                        {
                            xp = 0;
                        }
                        if (yp < 0)
                        {
                            yp = ly - 1;
                        }
                        if (yp >= ly)
                        {
                            yp = 0;
                        }
                        if (zp < 0)
                        {
                            zp = lz - 1;
                        }
                        if (zp >= lz)
                        {
                            zp = 0;
                        }
                        if (obst[xp][yp][zp] == 1)
                        {
                            sum_x += t_k[k] * xc[k];
                            sum_y += t_k[k] * yc[k];
                            sum_z += t_k[k] * zc[k];
                        }
                        else
                        {
                            Fx[x][y][z] += t_k[k] * xc[k] * psx[xp][yp][zp];
                            Fy[x][y][z] += t_k[k] * yc[k] * psx[xp][yp][zp];
                            Fz[x][y][z] += t_k[k] * zc[k] * psx[xp][yp][zp];
                        }
                    }
                    Sx[x][y][z] = -G1 * sum_x * psx[x][y][z] * psx_w;
                    Sy[x][y][z] = -G1 * sum_y * psx[x][y][z] * psx_w;
                    Sz[x][y][z] = -G1 * sum_z * psx[x][y][z] * psx_w;

                    Fx[x][y][z] = -G1 * psx[x][y][z] * Fx[x][y][z];
                    Fy[x][y][z] = -G1 * psx[x][y][z] * Fy[x][y][z];
                    Fz[x][y][z] = -G1 * psx[x][y][z] * Fz[x][y][z];
                }
            }
        }
    }
}

void stream(int obst[lx][ly][lz], double ff[19][lx][ly][lz]) {
    int x, y, z, x_e, x_w, y_n, y_s, z_n, z_s;
    double f_hlp[19][lx][ly][lz];

    for (z = 0; z < lz; z++)
    {
        for (y = 0; y < ly; y++)
        {
            for (x = 0; x < lx; x++)
            {
                z_n = ((z + 1) % lz);
                y_n = ((y + 1) % ly);
                x_e = ((x + 1) % lx);
                z_s = (lz - 1) - ((lz - z) % lz);
                y_s = (ly - 1) - ((ly - y) % ly);
                x_w = (lx - 1) - ((lx - x) % lx);

                // stream
                f_hlp[1][x_e][y][z] = ff[1][x][y][z];
                f_hlp[2][x_w][y][z] = ff[2][x][y][z];
                f_hlp[3][x][y_n][z] = ff[3][x][y][z];
                f_hlp[4][x][y_s][z] = ff[4][x][y][z];
                f_hlp[5][x][y][z_n] = ff[5][x][y][z];
                f_hlp[6][x][y][z_s] = ff[6][x][y][z];
                f_hlp[7][x_e][y_n][z] = ff[7][x][y][z];
                f_hlp[8][x_e][y_s][z] = ff[8][x][y][z];
                f_hlp[9][x_w][y_n][z] = ff[9][x][y][z];
                f_hlp[10][x_w][y_s][z] = ff[10][x][y][z];
                f_hlp[11][x_e][y][z_n] = ff[11][x][y][z];
                f_hlp[12][x_w][y][z_n] = ff[12][x][y][z];
                f_hlp[13][x_e][y][z_s] = ff[13][x][y][z];
                f_hlp[14][x_w][y][z_s] = ff[14][x][y][z];
                f_hlp[15][x][y_n][z_n] = ff[15][x][y][z];
                f_hlp[16][x][y_n][z_s] = ff[16][x][y][z];
                f_hlp[17][x][y_s][z_n] = ff[17][x][y][z];
                f_hlp[18][x][y_s][z_s] = ff[18][x][y][z];
            }
        }
    }


    for (int z = 0; z < lz; z++)
    {
        for (int y = 0; y < ly; y++)
        {
            for (int x = 0; x < lx; x++)
            {
                for (int k = 1; k < 19; k++)
                {
                    ff[k][x][y][z] = f_hlp[k][x][y][z];
                }
            }
        }
    }

}


void getuv(int obst[lx][ly][lz], double u_x[lx][ly][lz], double u_y[lx][ly][lz], double u_z[lx][ly][lz], double rho[lx][ly][lz], double ff[19][lx][ly][lz])
{
    int x, y, z;

    for (z = 0; z < lz; z++)
    {
        for (y = 0; y < ly; y++)
        {
            for (x = 0; x < lx; x++)
            {
                rho[x][y][z] = 0.0;

                if (obst[x][y][z] == 0)
                {
                    for (int k = 0; k < 19; k++)
                    {
                        rho[x][y][z] += ff[k][x][y][z];
                    }

                    if (rho[x][y][z] != 0.0)
                    {
                        u_x[x][y][z] = (ff[1][x][y][z] + ff[7][x][y][z] + ff[8][x][y][z] + ff[11][x][y][z] + ff[13][x][y][z] - (ff[2][x][y][z] + ff[9][x][y][z] + ff[10][x][y][z] + ff[12][x][y][z] + ff[14][x][y][z])) / rho[x][y][z];
                        u_y[x][y][z] = (ff[3][x][y][z] + ff[7][x][y][z] + ff[9][x][y][z] + ff[15][x][y][z] + ff[16][x][y][z] - (ff[4][x][y][z] + ff[8][x][y][z] + ff[10][x][y][z] + ff[17][x][y][z] + ff[18][x][y][z])) / rho[x][y][z];
                        u_z[x][y][z] = (ff[5][x][y][z] + ff[11][x][y][z] + ff[12][x][y][z] + ff[15][x][y][z] + ff[17][x][y][z] - (ff[6][x][y][z] + ff[13][x][y][z] + ff[14][x][y][z] + ff[16][x][y][z] + ff[18][x][y][z])) / rho[x][y][z];
                    }
                }
            }
        }
    }
}



void calcu_upr(int obst[lx][ly][lz], double u_x[lx][ly][lz], double u_y[lx][ly][lz], double u_z[lx][ly][lz], double Fx[lx][ly][lz], double Fy[lx][ly][lz], double Fz[lx][ly][lz], double Sx[lx][ly][lz], double Sy[lx][ly][lz], double Sz[lx][ly][lz], double rho[lx][ly][lz], double upx[lx][ly][lz], double upy[lx][ly][lz], double upz[lx][ly][lz])
{
    int x, y, z;

    for (z = 0; z < lz; z++)
    {
        for (y = 0; y < ly; y++)
        {
            for (x = 0; x < lx; x++)
            {
                if (obst[x][y][z] == 0)
                {
                    upx[x][y][z] = u_x[x][y][z] + (Fx[x][y][z] + Sx[x][y][z]) / 2.0 / rho[x][y][z];
                    upy[x][y][z] = u_y[x][y][z] + (Fy[x][y][z] + Sy[x][y][z]) / 2.0 / rho[x][y][z];
                    upz[x][y][z] = u_z[x][y][z] + (Fz[x][y][z] + Sz[x][y][z]) / 2.0 / rho[x][y][z];
                }
                else
                {
                    upx[x][y][z] = u_x[x][y][z];
                    upy[x][y][z] = u_y[x][y][z];
                    upz[x][y][z] = u_z[x][y][z];
                }
            }
        }
    }
}



void collision(double tauc, int obst[lx][ly][lz], double u_x[lx][ly][lz], double u_y[lx][ly][lz], double u_z[lx][ly][lz], double rho[lx][ly][lz], double ff[19][lx][ly][lz], double Fx[lx][ly][lz], double Fy[lx][ly][lz], double Fz[lx][ly][lz], double Sx[lx][ly][lz], double Sy[lx][ly][lz], double Sz[lx][ly][lz])
{
    int x, y, z, k, j;
    double ux, uy, uz, u_squ, u_n[19], u_n1[19], fequ[19], fequ1[19], fequ2[19];
    double temp[19];
    double ux1, uy1, uz1, u_squ1;
    double u_n2[19];

    double meq[19][lx][ly][lz];
    double m[19];
    double m1[19];
    double meq1[19];
    double A[19];
    double S[19];
    int tau = 1;
    int visc = c_squ * (tau - 0.5);
    double miu;

    for (z = 0; z < lz; z++)
    {
        for (y = 0; y < ly; y++)
        {
            for (x = 0; x < lx; x++)
            {
                if (obst[x][y][z] == 1)
                {
                    for (k = 1; k < 19; k++)
                    {
                        temp[k] = ff[k][x][y][z];
                    }
                    for (k = 1; k < 19; k++)
                    {
                        ff[opp[k]][x][y][z] = temp[k];
                    }

                }
                if (obst[x][y][z] == 0)
                {
                    ux = u_x[x][y][z] + tauc * (Fx[x][y][z] + Sx[x][y][z]) / rho[x][y][z];
                    uy = u_y[x][y][z] + tauc * (Fy[x][y][z] + Sy[x][y][z]) / rho[x][y][z];
                    uz = u_z[x][y][z] + tauc * (Fz[x][y][z] + Sz[x][y][z]) / rho[x][y][z];
                    u_squ = ux * ux + uy * uy + uz * uz;
                    for (k = 0; k < 19; k++)
                    {
                        u_n[k] = xc[k] * ux + yc[k] * uy + zc[k] * uz;
                        fequ[k] = t_k[k] * rho[x][y][z] * (cc * u_n[k] / c_squ + (u_n[k] * cc) * (u_n[k] * cc) / (2.0 * c_squ * c_squ) - u_squ / (2.0 * c_squ)) + t_k[k] * rho[x][y][z];
                        fequ2[k] = fequ[k] + (1.0 - 1.0 / tauc) * (ff[k][x][y][z] - fequ[k]);
                        ff[k][x][y][z] = fequ2[k];
                    }
                }
            }
        }
    }
}

void sum(double rho[lx][ly][lz])
{
    double sum = 0;
    for (int x = 0; x < lx; x++)
    {
        for (int y = 0; y < ly; y++)
        {
            for (int z = 0; z < lz; z++)
            {
                sum += rho[x][y][z];

            }
        }
    }
    cout << sum << endl;

}


void write_results(int obst[lx][ly][lz], double rho[lx][ly][lz], double upx[lx][ly][lz], double upy[lx][ly][lz], double upz[lx][ly][lz], double p[lx][ly][lz], int t_step)
{
    std::filesystem::path output_folder_path = "D:\\Test\\";//absolute path
    std::filesystem::path output_file_path = output_folder_path / ("3D" + std::to_string(t_step) + ".plt");

    if (!std::filesystem::exists(output_folder_path))
    {
        std::filesystem::create_directories(output_folder_path);
        std::cout << "Directory created successfully.\n";
    }

    std::ofstream output_file(output_file_path);

    output_file << "variables = x, y, z, rho, upx, upy, upz, p, obst\n";
    output_file << "zone i=" << lx << ", j=" << ly << ", k=" << lz << ", f=point\n";

    for (int z = 0; z < lz; z++)
    {
        for (int y = 0; y < ly; y++)
        {
            for (int x = 0; x < lx; x++)
            {
                output_file << x << ' ' << y << ' ' << z << ' '
                    << rho[x][y][z] << ' ' << upx[x][y][z] << ' ' << upy[x][y][z] << ' ' << upz[x][y][z] << ' ' << p[x][y][z] << ' ' << obst[x][y][z] << '\n';
            }
        }
    }

    output_file.close();
}


#endif /* LBM_hpp */
