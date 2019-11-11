using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using FastFluidSolverMT;

/*
 * Domains.cs
 * Copyright 2017 Christoph Waibel <chwaibel@student.ethz.ch>
 * Modified after WindInflow.cs; Lukas Bystricky <lb13f@my.fsu.edu>
 *
 * This work is licensed under the GNU GPL license version 3 or later.
 */

namespace GHWind
{
    //could have many more domains here 
    //public class Domains  :Domain
    //{


    //}


    /// <summary>
    /// Domain with an exponential wind profile on the inflow (x = 0), 
    /// 0 velocity on the ground (z = 0) and all other boundaries marked as outflow.
    /// </summary>
    public class WindInflow : Domain
    {
        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="Nx">Number of cells (including ghost cells) in x direction</param>
        /// <param name="Ny">Number of cells (including ghost cells) in x direction</param>
        /// <param name="Nz">Number of cells (including ghost cells) in x direction</param>
        /// <param name="length_x">Length of domain in x direction (not including ghost cells)</param>
        /// <param name="length_y">Length of domain in y direction (not including ghost cells)</param>
        /// <param name="length_z">Length of domain in z direction (not including ghost cells)</param>
        /// <param name="Vmet">wind speed measured at the meteorological station</param>
        /// <param name="terrain">Terrain Type [0,3]. 0 = Ocean; 1 = Flat, open country; 2 = Rough, wooded country, urban, industrial, forest; 3 = Towns and cities.</param>
        public WindInflow(int Nx, int Ny, int Nz, double length_x,
            double length_y, double length_z, double Vmet, int terrain)
        {
            this.Nx = Nx;
            this.Ny = Ny;
            this.Nz = Nz;

            this.length_x = length_x;
            this.length_y = length_y;
            this.length_z = length_z;

            hx = length_x / (Nx - 2);
            hy = length_y / (Ny - 2);
            hz = length_z / (Nz - 2);

            boundary_cells = new int[Nx, Ny, Nz];
            obstacle_cells = new int[Nx, Ny, Nz];
            boundary_normal_x = new int[Nx, Ny, Nz];
            boundary_normal_y = new int[Nx, Ny, Nz];
            boundary_normal_z = new int[Nx, Ny, Nz];
            boundary_u = new double[Nx - 1, Ny, Nz];
            boundary_v = new double[Nx, Ny - 1, Nz];
            boundary_w = new double[Nx, Ny, Nz - 1];

            outflow_boundary_x = new int[Nx, Ny, Nz];
            outflow_boundary_y = new int[Nx, Ny, Nz];
            outflow_boundary_z = new int[Nx, Ny, Nz];

            //add_obstacle_cube(15, 20, 15, 20, 0, 5);
            //add_obstacle_cube(15, 18, 18, 20, 5, 12);

            //add_obstacle_cube(18, 25, 25, 30, 0, 8);
            //add_obstacle_cube(30, 40, 15, 25, 0, 10);

            //set outflow boundaries
            //x = 0 will be inflow, z = 0 will be solid ground, all others will be outflow            

            //x outflow
            Parallel.For(1, Ny - 1, j =>
            {
                for (int k = 1; k < Nz - 1; k++)
                {
                    outflow_boundary_x[Nx - 2, j, k] = 1;
                }
            });

            //y outflows
            Parallel.For(1, Nx - 1, i =>
            {
                for (int k = 1; k < Nz - 1; k++)
                {
                    outflow_boundary_y[i, 1, k] = 1;
                    outflow_boundary_y[i, Ny - 2, k] = 1;
                }
            });

            //z outflow
            Parallel.For(1, Nx - 1, i =>
            {
                for (int j = 1; j < Ny - 1; j++)
                {
                    outflow_boundary_z[i, j, Nz - 2] = 1;
                }
            });

            double TerrainDelta = 210;
            double TerrainAlpha = 0.1;
            switch (terrain)
            {
                case 0:
                    TerrainAlpha = 0.1;
                    TerrainDelta = 210;
                    break;
                case 1:
                    TerrainAlpha = 0.14;
                    TerrainDelta = 270;
                    break;
                case 2:
                    TerrainAlpha = 0.22;
                    TerrainDelta = 370;
                    break;
                case 3:
                    TerrainAlpha = 0.33;
                    TerrainDelta = 460;
                    break;
            }
            //x = 0, inflow
            Parallel.For(0, boundary_u.GetLength(1), j =>
            {
                for (int k = 0; k < boundary_u.GetLength(2); k++)
                {
                    double z = (k - 0.5) * hz;
                    boundary_u[0, j, k] = Vmet * Math.Pow(TerrainDelta / 10, TerrainAlpha) * Math.Pow(Math.Max((z / TerrainDelta), 0), TerrainAlpha);
                }
            });

            set_ghost_flags();
            set_boundary_flags();
        }

        /// <summary>
        /// Copy constructor
        /// </summary>
        public WindInflow(WindInflow old)
        {
            Nx = old.Nx;
            Ny = old.Ny;
            Nz = old.Nz;

            hx = old.hx;
            hy = old.hy;
            hz = old.hz;

            length_x = old.length_x;
            length_y = old.length_y;
            length_z = old.length_z;

            boundary_cells = old.boundary_cells;
            obstacle_cells = old.obstacle_cells;
            boundary_normal_x = old.boundary_normal_x;
            boundary_normal_y = old.boundary_normal_y;
            boundary_normal_z = old.boundary_normal_z;
            boundary_u = old.boundary_u;
            boundary_v = old.boundary_v;
            boundary_w = old.boundary_w;
        }
    }


    public class WindInflowOpenFoam : Domain
    {
        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="Nx">Number of cells (including ghost cells) in x direction</param>
        /// <param name="Ny">Number of cells (including ghost cells) in x direction</param>
        /// <param name="Nz">Number of cells (including ghost cells) in x direction</param>
        /// <param name="length_x">Length of domain in x direction (not including ghost cells)</param>
        /// <param name="length_y">Length of domain in y direction (not including ghost cells)</param>
        /// <param name="length_z">Length of domain in z direction (not including ghost cells)</param>
        /// <param name="Uref">Reference velocity at 10m height [m/s]</param>
        /// <param name="z0">Surface roughness height [m]</param>
        public WindInflowOpenFoam(int Nx, int Ny, int Nz, double length_x, double length_y, double length_z,
            double Uref, double z0)
        {
            this.Nx = Nx;
            this.Ny = Ny;
            this.Nz = Nz;

            this.length_x = length_x;
            this.length_y = length_y;
            this.length_z = length_z;

            hx = length_x / (Nx - 2);
            hy = length_y / (Ny - 2);
            hz = length_z / (Nz - 2);

            boundary_cells = new int[Nx, Ny, Nz];
            obstacle_cells = new int[Nx, Ny, Nz];
            boundary_normal_x = new int[Nx, Ny, Nz];
            boundary_normal_y = new int[Nx, Ny, Nz];
            boundary_normal_z = new int[Nx, Ny, Nz];
            boundary_u = new double[Nx - 1, Ny, Nz];
            boundary_v = new double[Nx, Ny - 1, Nz];
            boundary_w = new double[Nx, Ny, Nz - 1];

            outflow_boundary_x = new int[Nx, Ny, Nz];
            outflow_boundary_y = new int[Nx, Ny, Nz];
            outflow_boundary_z = new int[Nx, Ny, Nz];

            //x outflow
            Parallel.For(1, Ny - 1, j =>
            {
                for (int k = 1; k < Nz - 1; k++)
                {
                    outflow_boundary_x[Nx - 2, j, k] = 1;
                }
            });

            //y outflows
            Parallel.For(1, Nx - 1, i =>
            {
                for (int k = 1; k < Nz - 1; k++)
                {
                    outflow_boundary_y[i, 1, k] = 1;
                    outflow_boundary_y[i, Ny - 2, k] = 1;
                }
            });

            //z outflow
            Parallel.For(1, Nx - 1, i =>
            {
                for (int j = 1; j < Ny - 1; j++)
                {
                    outflow_boundary_z[i, j, Nz - 2] = 1;
                }
            });


            // https://www.openfoam.com/documentation/guides/latest/api/classFoam_1_1atmBoundaryLayer.html#details
            const double kappa = 0.41;  // von Karman's constant
            const double Zref = 10.0;         // Reference height [m]
            const double zg = 0.0;            // Minimum z-coordinate [m]
            double Ustar = kappa * (Uref / (Math.Log((Zref + z0) / z0)));   // Friction velocity

            //x = 0, inflow
            Parallel.For(0, boundary_u.GetLength(1), j =>
            {
                for (int k = 0; k < boundary_u.GetLength(2); k++)
                {
                    double z = (k - 0.5) * hz;
                    boundary_u[0, j, k] = (Ustar / kappa) * Math.Max(Math.Log(Math.Max((z - zg + z0) / z0, 0)), 0);
                }
            });

            set_ghost_flags();
            set_boundary_flags();
        }

        /// <summary>
        /// Copy constructor
        /// </summary>
        public WindInflowOpenFoam(WindInflow old)
        {
            Nx = old.Nx;
            Ny = old.Ny;
            Nz = old.Nz;

            hx = old.hx;
            hy = old.hy;
            hz = old.hz;

            length_x = old.length_x;
            length_y = old.length_y;
            length_z = old.length_z;

            boundary_cells = old.boundary_cells;
            obstacle_cells = old.obstacle_cells;
            boundary_normal_x = old.boundary_normal_x;
            boundary_normal_y = old.boundary_normal_y;
            boundary_normal_z = old.boundary_normal_z;
            boundary_u = old.boundary_u;
            boundary_v = old.boundary_v;
            boundary_w = old.boundary_w;
        }
    }



    public class WindInflowAytac : Domain
    {
        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="Nx">Number of cells (including ghost cells) in x direction</param>
        /// <param name="Ny">Number of cells (including ghost cells) in x direction</param>
        /// <param name="Nz">Number of cells (including ghost cells) in x direction</param>
        /// <param name="length_x">Length of domain in x direction (not including ghost cells)</param>
        /// <param name="length_y">Length of domain in y direction (not including ghost cells)</param>
        /// <param name="length_z">Length of domain in z direction (not including ghost cells)</param>
        /// <param name="Vmet">wind speed measured at the meteorological station</param>
        /// <param name="terrain">Terrain Type [0,3]. 0 = Ocean; 1 = Flat, open country; 2 = Rough, wooded country, urban, industrial, forest; 3 = Towns and cities.</param>
        public WindInflowAytac(int Nx, int Ny, int Nz, double length_x,
            double length_y, double length_z)
        {
            this.Nx = Nx;
            this.Ny = Ny;
            this.Nz = Nz;

            this.length_x = length_x;
            this.length_y = length_y;
            this.length_z = length_z;

            hx = length_x / (Nx - 2);
            hy = length_y / (Ny - 2);
            hz = length_z / (Nz - 2);

            boundary_cells = new int[Nx, Ny, Nz];
            obstacle_cells = new int[Nx, Ny, Nz];
            boundary_normal_x = new int[Nx, Ny, Nz];
            boundary_normal_y = new int[Nx, Ny, Nz];
            boundary_normal_z = new int[Nx, Ny, Nz];
            boundary_u = new double[Nx - 1, Ny, Nz];
            boundary_v = new double[Nx, Ny - 1, Nz];
            boundary_w = new double[Nx, Ny, Nz - 1];

            outflow_boundary_x = new int[Nx, Ny, Nz];
            outflow_boundary_y = new int[Nx, Ny, Nz];
            outflow_boundary_z = new int[Nx, Ny, Nz];
            
            //x outflow
            Parallel.For(1, Ny - 1, j =>
            {
                for (int k = 1; k < Nz - 1; k++)
                {
                    outflow_boundary_x[Nx - 2, j, k] = 1;
                }
            });

            //y outflows
            Parallel.For(1, Nx - 1, i =>
            {
                for (int k = 1; k < Nz - 1; k++)
                {
                    outflow_boundary_y[i, 1, k] = 1;
                    outflow_boundary_y[i, Ny - 2, k] = 1;
                }
            });

            //z outflow
            Parallel.For(1, Nx - 1, i =>
            {
                for (int j = 1; j < Ny - 1; j++)
                {
                    outflow_boundary_z[i, j, Nz - 2] = 1;
                }
            });


            //Aytac profile
            //uABL=0.9965; // atmospheric boundary layer friction velocity
            //kappa=0.42; // von Karman constant
            //y0=0.15; // aerodynamic roughness length
            //y =0:0.1:50; // height
            //U=(uABL/kappa)*log((y+y0)/y0); // log-law velocity profile

            double uABL = 0.9965;
            double kappa = 0.42;
            double z0 = 0.15;
            
            //x = 0, inflow
            Parallel.For(0, boundary_u.GetLength(1), j =>
            {
                for (int k = 0; k < boundary_u.GetLength(2); k++)
                {
                    double z = (k - 0.5) * hz;
                    boundary_u[0, j, k] = (uABL / kappa) * Math.Max(Math.Log(Math.Max((z + z0) / z0, 0)), 0);
                }
            });

            set_ghost_flags();
            set_boundary_flags();
        }

        /// <summary>
        /// Copy constructor
        /// </summary>
        public WindInflowAytac(WindInflow old)
        {
            Nx = old.Nx;
            Ny = old.Ny;
            Nz = old.Nz;

            hx = old.hx;
            hy = old.hy;
            hz = old.hz;

            length_x = old.length_x;
            length_y = old.length_y;
            length_z = old.length_z;

            boundary_cells = old.boundary_cells;
            obstacle_cells = old.obstacle_cells;
            boundary_normal_x = old.boundary_normal_x;
            boundary_normal_y = old.boundary_normal_y;
            boundary_normal_z = old.boundary_normal_z;
            boundary_u = old.boundary_u;
            boundary_v = old.boundary_v;
            boundary_w = old.boundary_w;
        }
    }


    public class WindInflowWentao : Domain
    {
        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="Nx">Number of cells (including ghost cells) in x direction</param>
        /// <param name="Ny">Number of cells (including ghost cells) in x direction</param>
        /// <param name="Nz">Number of cells (including ghost cells) in x direction</param>
        /// <param name="length_x">Length of domain in x direction (not including ghost cells)</param>
        /// <param name="length_y">Length of domain in y direction (not including ghost cells)</param>
        /// <param name="length_z">Length of domain in z direction (not including ghost cells)</param>
        /// <param name="Vmet">wind speed measured at the meteorological station</param>
        /// <param name="terrain">Terrain Type [0,3]. 0 = Ocean; 1 = Flat, open country; 2 = Rough, wooded country, urban, industrial, forest; 3 = Towns and cities.</param>
        public WindInflowWentao(int Nx, int Ny, int Nz, double length_x,
            double length_y, double length_z)
        {
            this.Nx = Nx;
            this.Ny = Ny;
            this.Nz = Nz;

            this.length_x = length_x;
            this.length_y = length_y;
            this.length_z = length_z;

            hx = length_x / (Nx - 2);
            hy = length_y / (Ny - 2);
            hz = length_z / (Nz - 2);

            boundary_cells = new int[Nx, Ny, Nz];
            obstacle_cells = new int[Nx, Ny, Nz];
            boundary_normal_x = new int[Nx, Ny, Nz];
            boundary_normal_y = new int[Nx, Ny, Nz];
            boundary_normal_z = new int[Nx, Ny, Nz];
            boundary_u = new double[Nx - 1, Ny, Nz];
            boundary_v = new double[Nx, Ny - 1, Nz];
            boundary_w = new double[Nx, Ny, Nz - 1];

            outflow_boundary_x = new int[Nx, Ny, Nz];
            outflow_boundary_y = new int[Nx, Ny, Nz];
            outflow_boundary_z = new int[Nx, Ny, Nz];


            //set outflow boundaries
            //x = 0 will be inflow, z = 0 will be solid ground, all others will be outflow            

            //x outflow
            Parallel.For(1, Ny - 1, j =>
            {
                for (int k = 1; k < Nz - 1; k++)
                {
                    outflow_boundary_x[Nx - 2, j, k] = 1;
                }
            });

            //y outflows
            Parallel.For(1, Nx - 1, i =>
            {
                for (int k = 1; k < Nz - 1; k++)
                {
                    outflow_boundary_y[i, 1, k] = 1;
                    outflow_boundary_y[i, Ny - 2, k] = 1;
                }
            });

            //z outflow
            Parallel.For(1, Nx - 1, i =>
            {
                for (int j = 1; j < Ny - 1; j++)
                {
                    outflow_boundary_z[i, j, Nz - 2] = 1;
                }
            });


            //x = 0, inflow
            double p1, p2, p3, p4, p5, p6, p7, p8;
            p1 = 0.0001469;
            p2 = -0.004079;
            p3 = 0.04558;
            p4 = -0.2643;
            p5 = 0.8549;
            p6 = -1.56;
            p7 = 2.113;
            p8 = 2.671;

            Parallel.For(0, boundary_u.GetLength(1), j =>
            {
                for (int k = 0; k < boundary_u.GetLength(2); k++)
                {
                    double z = (k - 0.5) * hz;
                    //this is from curve fitting tool
                    boundary_u[0, j, k] = p1 * Math.Pow(z, 7) + p2 * Math.Pow(z, 6) + p3 * Math.Pow(z, 5) + p4 * Math.Pow(z, 4)
                        + p5 * Math.Pow(z, 3) + p6 * Math.Pow(z, 2) + p7 * z + p8;
                }
            });

            set_ghost_flags();
            set_boundary_flags();
        }

        /// <summary>
        /// Copy constructor
        /// </summary>
        public WindInflowWentao(WindInflow old)
        {
            Nx = old.Nx;
            Ny = old.Ny;
            Nz = old.Nz;

            hx = old.hx;
            hy = old.hy;
            hz = old.hz;

            length_x = old.length_x;
            length_y = old.length_y;
            length_z = old.length_z;

            boundary_cells = old.boundary_cells;
            obstacle_cells = old.obstacle_cells;
            boundary_normal_x = old.boundary_normal_x;
            boundary_normal_y = old.boundary_normal_y;
            boundary_normal_z = old.boundary_normal_z;
            boundary_u = old.boundary_u;
            boundary_v = old.boundary_v;
            boundary_w = old.boundary_w;
        }
    }


    /// <summary>
    /// Domain for the lid driven cavity
    /// </summary>
    public class CavityDomain : Domain
    {
        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="Nx">Number of cells in x direction (including ghost cells)</param>
        /// <param name="Ny">Number of cells in y direction (including ghost cells)</param>
        /// <param name="Nz">Number of cells in z direction (including ghost cells)</param>
        public CavityDomain(int Nx, int Ny, int Nz)
        {
            this.Nx = Nx;
            this.Ny = Ny;
            this.Nz = Nz;

            length_x = 1;
            length_y = 1;
            length_z = 1;

            hx = length_x / (Nx - 2);
            hy = length_y / (Ny - 2);
            hz = length_z / (Nz - 2);

            boundary_cells = new int[Nx, Ny, Nz];
            obstacle_cells = new int[Nx, Ny, Nz];
            boundary_normal_x = new int[Nx, Ny, Nz];
            boundary_normal_y = new int[Nx, Ny, Nz];
            boundary_normal_z = new int[Nx, Ny, Nz];
            boundary_u = new double[Nx - 1, Ny, Nz];
            boundary_v = new double[Nx, Ny - 1, Nz];
            boundary_w = new double[Nx, Ny, Nz - 1];
            outflow_boundary_x = new int[Nx, Ny, Nz];
            outflow_boundary_y = new int[Nx, Ny, Nz];
            outflow_boundary_z = new int[Nx, Ny, Nz];

            //C# default values for int or double arrays are 0, so we only need to set nonzero fields

            set_ghost_flags();
            set_boundary_flags();

            /**************************************************************************************
             * u boundary conditions
             *************************************************************************************/
            Parallel.For(0, boundary_u.GetLength(0), i =>
            {
                for (int j = 0; j < boundary_u.GetLength(1); j++)
                {
                    boundary_u[i, j, Nz - 1] = 1;
                }
            });
        }

        /// <summary>
        /// Copy constructor
        /// </summary>
        /// <param name="old"></param>
        public CavityDomain(CavityDomain old)
        {
            Nx = old.Nx;
            Ny = old.Ny;
            Nz = old.Nz;

            hx = old.hx;
            hy = old.hy;
            hz = old.hz;

            length_x = old.length_x;
            length_y = old.length_y;
            length_z = old.length_z;

            boundary_cells = old.boundary_cells;
            obstacle_cells = old.obstacle_cells;
            boundary_normal_x = old.boundary_normal_x;
            boundary_normal_y = old.boundary_normal_y;
            boundary_normal_z = old.boundary_normal_z;
            boundary_u = old.boundary_u;
            boundary_v = old.boundary_v;
            boundary_w = old.boundary_w;

            outflow_boundary_x = old.outflow_boundary_x;
            outflow_boundary_y = old.outflow_boundary_y;
            outflow_boundary_z = old.outflow_boundary_z;
        }
    }   


    public class DomainOpen : Domain
    {
        public double windspeed = 0;

        private const double a = 1.25;
        private const double d = 2.25;
        private const double nu = 1;

        /***************************************************************************
         * Constructor
         **************************************************************************/
        public DomainOpen(int Nx, int Ny, int Nz, int length_x, int length_y, int length_z,
            int[, ,] obstacles, double windspeed) 
        {
            this.windspeed = windspeed;

            this.Nx = Nx;
            this.Ny = Ny;
            this.Nz = Nz;

            this.length_x = length_x;
            this.length_y = length_y;
            this.length_z = length_z;

            hx = length_x / (Nx - 2);
            hy = length_y / (Ny - 2);
            hz = length_z / (Nz - 2);

            boundary_cells = new int[Nx, Ny, Nz];
            obstacle_cells = new int[Nx, Ny, Nz];
            boundary_normal_x = new int[Nx, Ny, Nz];
            boundary_normal_y = new int[Nx, Ny, Nz];
            boundary_normal_z = new int[Nx, Ny, Nz];
            boundary_u = new double[Nx - 1, Ny, Nz];
            boundary_v = new double[Nx, Ny - 1, Nz];
            boundary_w = new double[Nx, Ny, Nz - 1];

            outflow_boundary_x = new int[Nx, Ny, Nz];
            outflow_boundary_y = new int[Nx, Ny, Nz];
            outflow_boundary_z = new int[Nx, Ny, Nz];


            //set up obstacle
            for (int i = 0; i < Nx; i++)
            {
                for (int j = 0; j < Ny; j++)
                {
                    for (int k = 0; k < Nz; k++)
                    {
                        obstacle_cells[i, j, k] = obstacles[i, j, k];
                    }
                }
            }


            //set outflow boundaries
            //x = 0 will be inflow, z = 0 will be solid ground, all others will be outflow
            
            //x outflow
            for (int j = 1; j < Ny - 1; j++)
            {
                for (int k = 1; k < Nz - 1; k++)
                {
                    outflow_boundary_x[Nx - 2, j, k] = 1;
                }
            }

            //y outflows
            for (int i = 1; i < Nx - 1; i++)
            {
                for (int k = 1; k < Nz - 1; k++)
                {
                    outflow_boundary_y[i, 1, k] = 1;
                    outflow_boundary_y[i, Ny - 2, k] = 1;
                }
            }

            //z outflow
            for (int i = 1; i < Nx - 1; i++)
            {
                for (int j = 1; j < Ny - 1; j++)
                {
                    outflow_boundary_z[i, j, Nz - 2] = 1;
                }
            }

            set_ghost_flags();
            set_boundary_flags();
            update_boundary_conditions(0);
        }

        public void update_boundary_conditions(double t)
        {
            /**************************************************************************************
            * inflow boundary conditions
            *************************************************************************************/

            for (int j = 0; j < boundary_u.GetLength(1); j++)
            {
                for (int k = 0; k < boundary_u.GetLength(2); k++)
                {
                    boundary_u[0, j, k] = windspeed;
                }
            }

        }
        /***********************************************************************
         * Copy constructor
         ***********************************************************************/
        public DomainOpen(DomainOpen old)
        {
            Nx = old.Nx;
            Ny = old.Ny;
            Nz = old.Nz;

            hx = old.hx;
            hy = old.hy;
            hz = old.hz;

            length_x = old.length_x;
            length_y = old.length_y;
            length_z = old.length_z;

            boundary_cells = old.boundary_cells;
            obstacle_cells = old.obstacle_cells;
            boundary_normal_x = old.boundary_normal_x;
            boundary_normal_y = old.boundary_normal_y;
            boundary_normal_z = old.boundary_normal_z;
            boundary_u = old.boundary_u;
            boundary_v = old.boundary_v;
            boundary_w = old.boundary_w;
        }

    }
}
