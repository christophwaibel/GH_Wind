using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;
using FastFluidSolverrr;
using System.Threading.Tasks;

/*
 * GHFFDliddrivencavity.cs
 * Copyright 2017 Christoph Waibel <chwaibel@student.ethz.ch>
 * 
 * This work is licensed under the GNU GPL license version 3.
*/

namespace GHWind
{
    public class GHFFDliddrivencavity : GH_Component
    {


        Domain omega;
        FluidSolver ffd;
        DataExtractor de;

        double t;
        bool resetFFD;



        /// <summary>
        /// Initializes a new instance of the GHFFDliddrivencavity class.
        /// </summary>
        public GHFFDliddrivencavity()
            : base("FFD lid driven cavity", "FFDldc",
            "Fast Fluid Dynamics Solver",
            "EnergyHubs", "Wind Simulation")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            //#0, #1
            pManager.AddNumberParameter("Domain size", "Domain size", "Domain x,y,z size in [m].", GH_ParamAccess.list);
            pManager.AddIntegerParameter("Domain discretization", "Nx,Ny,Nz", "Domain discretization Nx, Ny, Nz, i.e. how many fluid cells in each direction.", GH_ParamAccess.list);
            //later make another input for defining more precisely the domain. like, internal flow, external flow, inflows, outflows...)

            //#2
            pManager.AddNumberParameter("Time Step", "dt", "Calculation time step dt.", GH_ParamAccess.item);

            //#3
            pManager.AddBooleanParameter("Run?", "Run?", "Run the solver. (Loop via Grasshopper timer component)", GH_ParamAccess.item);

            //#4
            pManager.AddBooleanParameter("Write Results?", "Results?", "Write Results for visualization in GH.", GH_ParamAccess.item);
            pManager[4].Optional = true;

            //#5
            pManager.AddBooleanParameter("Export VTK", "ExpVTK", "Export Results to VTK", GH_ParamAccess.item);
            pManager[5].Optional = true;

            //#6
            pManager.AddBooleanParameter("Reset", "Reset", "Reset domain", GH_ParamAccess.item);
            pManager[6].Optional = true;

            //#7
            pManager.AddNumberParameter("Re", "Re", "Reynolds Number", GH_ParamAccess.item);


        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            //# 0-1
            pManager.AddGenericParameter("v centred", "v centred", "velocities, cell centred", GH_ParamAccess.list);
            pManager.AddGenericParameter("p centred", "p centred", "pressure, cell centred", GH_ParamAccess.item);

            //# 2-3
            pManager.AddGenericParameter("v staggered", "v staggered", "velocities, on staggered grid", GH_ParamAccess.list);
            pManager.AddGenericParameter("p staggered", "p staggered", "pressure, on staggered grid", GH_ParamAccess.item);

            //#4
            pManager.AddNumberParameter("u", "u", "u in [m/s] at vertical center line over z", GH_ParamAccess.list);
            //#5
            pManager.AddNumberParameter("v", "v", "v in [m/s] at vertical center line over x", GH_ParamAccess.list);

        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // *********************************************************************************
            // get info from grasshopper
            // *********************************************************************************
            List<double> xyzsize = new List<double>();
            if (!DA.GetDataList(0, xyzsize)) { return; };

            List<int> Nxyz = new List<int>();
            if (!DA.GetDataList(1, Nxyz)) { return; };
            int Nx = Nxyz[0];
            int Ny = Nxyz[1];
            int Nz = Nxyz[2];


            // time step
            double dt = 0.1;
            if (!DA.GetData(2, ref dt)) { return; }


            bool run = false;
            if (!DA.GetData(3, ref run)) { return; }



            DA.GetData(6, ref resetFFD);

            double re = 400;
            if (!DA.GetData(7, ref re)) { return; }





            // *********************************************************************************
            //from Lukas 
            // *********************************************************************************
            double nu = 1 / re;       // increase viscosity to impose turbulence. the higher velocity, the higher visc., 1e-3

            // Set initial velocity conditions
            double[, ,] u0 = new double[Nx + 1, Ny + 2, Nz + 2];
            double[, ,] v0 = new double[Nx + 2, Ny + 1, Nz + 2];
            double[, ,] w0 = new double[Nx + 2, Ny + 2, Nz + 1];

            // Create empty arrays for body forces
            double[, ,] f_x = new double[Nx + 1, Ny + 2, Nz + 2];
            double[, ,] f_y = new double[Nx + 2, Ny + 1, Nz + 2];
            double[, ,] f_z = new double[Nx + 2, Ny + 2, Nz + 1];

            // Create structure for solver parameters
            FluidSolver.solver_struct solver_prams = new FluidSolver.solver_struct();
            solver_prams.tol = 1e-4;
            solver_prams.min_iter = 1;
            solver_prams.max_iter = 30;
            solver_prams.verbose = false;
            solver_prams.backtrace_order = 2;
            solver_prams.mass_correction = false;
            solver_prams.mass_corr_alpha = 0.7;


            // Create FFD solver and domain
            if (ffd == null)
            {
                // Create domain
                // Pass in the number of cells in x, y and z direction (now including ghost cells)
                omega = new CavityDomain(Nx + 2, Ny + 2, Nz + 2);
                ffd = new FluidSolver(omega, dt, nu, u0, v0, w0, solver_prams);
                de = new DataExtractor(omega, ffd);
                t = 0;
            }

            //reset FFD solver and domain
            if (resetFFD)
            {
                omega = new CavityDomain(Nx + 2, Ny + 2, Nz + 2);
                ffd = new FluidSolver(omega, dt, nu, u0, v0, w0, solver_prams);
                de = new DataExtractor(omega, ffd);
                t = 0;
                resetFFD = false;
            }





            //run solver. the solving-loop (new timestep) is executed in Grasshopper with a timer-component.
            if (run) ffd.time_step(f_x, f_y, f_z);





            // *******************************************************************************************
            // *********************************     Output Results       ********************************
            double[, ,] p = new double[Nx, Ny, Nz];
            double[, ,] vu = new double[Nx, Ny, Nz];
            double[, ,] vv = new double[Nx, Ny, Nz];
            double[, ,] vw = new double[Nx, Ny, Nz];

            double[, ,] pstag = new double[Nx + 1, Ny + 1, Nz + 1];
            double[, ,] vustag = new double[Nx + 1, Ny + 1, Nz + 1];
            double[, ,] vvstag = new double[Nx + 1, Ny + 1, Nz + 1];
            double[, ,] vwstag = new double[Nx + 1, Ny + 1, Nz + 1];

            for (int i = 0; i < Nx; i++)
            {
                for (int j = 0; j < Ny; j++)
                {
                    for (int k = 0; k < Nz; k++)
                    {
                        if (omega.obstacle_cells[i + 1, j + 1, k + 1] != 1)
                        {
                            p[i, j, k] = de.get_pressure(i * omega.hx + 0.5 * omega.hx, j * omega.hy + 0.5 * omega.hy, k * omega.hz + 0.5 * omega.hz);
                            double[] vel = de.get_velocity(i * omega.hx + 0.5 * omega.hx, j * omega.hy + 0.5 * omega.hy, k * omega.hz + 0.5 * omega.hz);
                            vu[i, j, k] = vel[0];
                            vv[i, j, k] = vel[1];
                            vw[i, j, k] = vel[2];

                            pstag[i, j, k] = de.get_pressure(i * omega.hx, j * omega.hy, k * omega.hz);
                            double[] velcen = de.get_velocity(i * omega.hx, j * omega.hy, k * omega.hz);
                            vustag[i, j, k] = velcen[0];
                            vvstag[i, j, k] = velcen[1];
                            vwstag[i, k, k] = velcen[2];

                        }
                        else
                        {
                            p[i, j, k] = 0;
                            vu[i, j, k] = 0;
                            vv[i, j, k] = 0;
                            vw[i, j, k] = 0;

                            //pstag[i, j, k] = 0;
                            //vustag[i, j, k] = 0;
                            //vvstag[i, j, k] = 0;
                            //vwstag[i, k, k] = 0;

                            pstag[i, j, k] = de.get_pressure(i * omega.hx, j * omega.hy, k * omega.hz);
                            double[] velcen = de.get_velocity(i * omega.hx, j * omega.hy, k * omega.hz);
                            vustag[i, j, k] = velcen[0];
                            vvstag[i, j, k] = velcen[1];
                            vwstag[i, k, k] = velcen[2];
                        }
                    }
                }
            }

            //last x slice
            for (int j = 0; j < Ny + 1; j++)
            {
                for (int k = 0; k < Nz + 1; k++)
                {
                    pstag[Nx, j, k] = de.get_pressure((Nx) * omega.hx, j * omega.hy, k * omega.hz);
                    double[] vcen = de.get_velocity((Nx) * omega.hx, j * omega.hy, k * omega.hz);
                    vustag[Nx, j, k] = vcen[0];
                    vvstag[Nx, j, k] = vcen[1];
                    vwstag[Nx, j, k] = vcen[2];
                }
            }

            //last y slice
            for (int i = 0; i < Nx + 1; i++)
            {
                for (int k = 0; k < Nz + 1; k++)
                {
                    pstag[i, Ny, k] = de.get_pressure(i * omega.hx, (Ny) * omega.hy, k * omega.hz);
                    double[] vcen = de.get_velocity(i * omega.hx, (Ny) * omega.hy, k * omega.hz);
                    vustag[i, Ny, k] = vcen[0];
                    vvstag[i, Ny, k] = vcen[1];
                    vwstag[i, Ny, k] = vcen[2];
                }
            }

            //last z slice
            for (int i = 0; i < Nx + 1; i++)
            {
                for (int j = 0; j < Ny + 1; j++)
                {
                    pstag[i, j, Nz] = de.get_pressure(i * omega.hx, j * omega.hy, (Nz) * omega.hz);
                    double[] vcen = de.get_velocity(i * omega.hx, j * omega.hy, (Nz) * omega.hz);
                    vustag[i, j, Nz] = vcen[0];
                    vvstag[i, j, Nz] = vcen[1];
                    vwstag[i, j, Nz] = vcen[2];
                }
            }

            List<double[, ,]> veloutCen = new List<double[, ,]> { };
            veloutCen.Add(vu);
            veloutCen.Add(vv);
            veloutCen.Add(vw);

            List<double[, ,]> veloutStag = new List<double[, ,]> { };
            veloutStag.Add(vustag);
            veloutStag.Add(vvstag);
            veloutStag.Add(vwstag);


            DA.SetDataList(0, veloutCen);
            DA.SetData(1, p);
            DA.SetDataList(2, veloutStag);
            DA.SetData(3, pstag);





            List<double> ulist = new List<double>();
            List<double> vlist = new List<double>();

            //u velocities. (y called) z-section, in the vertical center... x=0.5Nx, y=0.5Ny
            //(Ny+1) / 2
            for (int k = 0; k < Nz + 1; k++)
            {
                double[] vcen = de.get_velocity((Nx / 2) * omega.hx + 0.5 * omega.hx,
                                                (Ny / 2) * omega.hy + 0.5 * omega.hy, 
                                                                        k * omega.hz);
                ulist.Add(vcen[0]);
            }

            for (int i = 0; i < Nx + 1; i++)
            {
                double[] vcen = de.get_velocity(i * omega.hx, 
                                                (Ny / 2) * omega.hy + 0.5 * omega.hy, 
                                                (Nz / 2) * omega.hz + 0.5 * omega.hz);

                vlist.Add(vcen[2]);
            }
            

            DA.SetDataList(4, ulist);
            DA.SetDataList(5, vlist);










        }

        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                // return Resources.IconForThisComponent;
                return GHWind.Properties.Resources.liddrivencavity;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("{28750158-9877-4fb6-9968-c826e5c13855}"); }
        }
    }
}