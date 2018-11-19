using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;

using FastFluidSolverMT;
using System.Threading.Tasks;
using System.IO;

/*
 * GHFFDSolver.cs
 * Copyright 2018 Christoph Waibel <chwaibel@student.ethz.ch>
 * 
 * This work is licensed under the GNU GPL license version 3.
*/


namespace GHWind
{
    public class GHFFDSolverV2 : GH_Component
    {
        Domain omega;
        FluidSolver ffd;
        DataExtractor de;

        double t;
        bool resetFFD;

        /// <summary>
        /// Initializes a new instance of the GHFFDSolverV2 class.
        /// </summary>
        public GHFFDSolverV2()
            : base("FFD Solver V2", "FFD_V2",
            "Fast Fluid Dynamics Solver V2",
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
            pManager.AddGenericParameter("Geometry", "Geometry", "Geometry as list of doubles [6] {xmin, xmax, ymin, ymax, zmin, zmax}, representing the obstacle cubes.", GH_ParamAccess.list);

            //#3
            pManager.AddNumberParameter("Time Step", "dt", "Calculation time step dt.", GH_ParamAccess.item);

            //#4
            pManager.AddNumberParameter("Horizon", "Horizon", "Calculation time horizon, i.e. sum of dt. until termination. Should be sufficient for domain to converge, but shouldn't be too much, otherwise wasted time. Numeric convergence indicators would help here.", GH_ParamAccess.item);

            //#5
            pManager.AddNumberParameter("Wind Speed", "Vmet", "Wind Speed [m/s] at meteorological station at 10 m height above ground.", GH_ParamAccess.item);

            //#6
            pManager.AddIntegerParameter("Terrain", "Terrain", "Terrain coefficients for wind speed. 0 = Ocean; 1 = Flat, open country; 2 = Rough, wooded country, urban, industrial, forest; 3 = Towns and Cities.", GH_ParamAccess.item);

            //#7
            pManager.AddBooleanParameter("Run?", "Run?", "Run the solver. (Loop via Grasshopper timer component)", GH_ParamAccess.item);

            //#8
            pManager.AddBooleanParameter("Results?", "Results?", "Output Data Extractor class? E.g. for Cp calculation, or flow visualization.", GH_ParamAccess.item);
            pManager[8].Optional = true;

            //#9
            pManager.AddBooleanParameter("Export VTK", "ExpVTK", "Export Results to VTK", GH_ParamAccess.item);
            pManager[9].Optional = true;

            //#10
            pManager.AddBooleanParameter("Reset", "Reset", "Reset domain", GH_ParamAccess.item);
            pManager[10].Optional = true;

            //#11
            pManager.AddTextParameter("Solver Parameters", "params", "FFD solver parameters. Provide a semicolon-separated string, e.g. '1.511e-5;1e-4;1;30;2;false;0.7;false'. Items: 'kinematic viscosity (double); tolerance (double); min_iter (int); max_iter (int); backtrace_order (int, 1 or 2); mass_correction (true or false); mass_corr_alpha (double), verbose (true or false).", GH_ParamAccess.item);
            pManager[11].Optional = true;

            //#12
            pManager.AddBooleanParameter("Residuals?", "Residuals?", "Calculate residuals for convergence analysis? Writes text file to 'C:\residuals.txt'.", GH_ParamAccess.item);
            pManager[12].Optional = true;

            //#13
            pManager.AddIntegerParameter("mean_dt", "mean_dt", "m*dt for outputting mean flow field (instead of snapshot). m should be identified by observing the residuals. Default is m=10.", GH_ParamAccess.item);
            pManager[13].Optional = true;
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            //#0,1
            pManager.AddGenericParameter("v centred", "v centred", "velocities, cell centred", GH_ParamAccess.list);
            pManager.AddGenericParameter("p centred", "p centred", "pressure, cell centred", GH_ParamAccess.item);

            //#2,3
            pManager.AddGenericParameter("v staggered", "v staggered", "velocities, on staggered grid", GH_ParamAccess.list);
            pManager.AddGenericParameter("p staggered", "p staggered", "pressure, on staggered grid", GH_ParamAccess.item);

            //#4
            pManager.AddGenericParameter("DE", "DE", "Data Extractor, containing omega and FFD classes", GH_ParamAccess.item);

           // pManager.AddTextParameter("VTK path", "VTK path", "Output path of VTK results file", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // *********************************************************************************
            // Inputs
            // *********************************************************************************
            List<double> xyzsize = new List<double>();
            if (!DA.GetDataList(0, xyzsize)) { return; };

            List<int> Nxyz = new List<int>();
            if (!DA.GetDataList(1, Nxyz)) { return; };
            int Nx = Nxyz[0];
            int Ny = Nxyz[1];
            int Nz = Nxyz[2];


            List<double[]> geom = new List<double[]>();
            if (!DA.GetDataList(2, geom)) { return; };


            // time step
            double dt = 0.1;
            if (!DA.GetData(3, ref dt)) { return; }

            // horizon
            double t_end = 1;
            if(!DA.GetData(4, ref t_end)) { return; }

            // wind speed
            double Vmet = 10;
            if (!DA.GetData(5, ref Vmet)) { return; }

            //terrain type
            int terrain = 0;
            if (!DA.GetData(6, ref terrain)) { return; }


            bool run = false;
            if (!DA.GetData(7, ref run)) { return; }



            //List<Mesh> mshCp = new List<Mesh>();
            //DA.GetDataList(10, mshCp);
            bool writeresults = false;
            DA.GetData(8, ref writeresults);


            DA.GetData(10, ref resetFFD);

            bool calcres = false;
            DA.GetData(12, ref calcres);

            int m = 10;
            DA.GetData(13, ref m);

            string strparam = null;
            DA.GetData(11, ref strparam);

            string[] str_params =null;
            if (strparam != null) str_params = strparam.Split(';');
          
            double nu = 1.511e-5;       // increase viscosity to impose turbulence. the higher velocity, the higher visc., 1e-3
            FluidSolver.solver_struct solver_params = new FluidSolver.solver_struct();
            if (str_params != null)
            {
                nu = Convert.ToDouble(str_params[0]);
                solver_params.tol = Convert.ToDouble(str_params[1]);
                solver_params.min_iter = Convert.ToInt16(str_params[2]);
                solver_params.max_iter = Convert.ToInt16(str_params[3]);
                solver_params.backtrace_order = Convert.ToInt16(str_params[4]);
                solver_params.mass_correction = str_params[5].Equals("false") ? false : true;
                solver_params.mass_corr_alpha = Convert.ToDouble(str_params[6]);
                solver_params.verbose = str_params[7].Equals("false") ? false : true;
            }
            else
            {
                solver_params.tol = 1e-4;
                solver_params.min_iter = 1;
                solver_params.max_iter = 30;
                solver_params.backtrace_order = 2;
                solver_params.mass_correction = false;
                solver_params.mass_corr_alpha = 0.7;
                solver_params.verbose = false;
            }





            // *********************************************************************************
            // Set-up FFD Solver
            // *********************************************************************************
            // Set initial velocity conditions
            double[,,] u0 = new double[Nx + 1, Ny + 2, Nz + 2];
            double[,,] v0 = new double[Nx + 2, Ny + 1, Nz + 2];
            double[,,] w0 = new double[Nx + 2, Ny + 2, Nz + 1];

            // Create empty arrays for body forces
            double[,,] f_x = new double[Nx + 1, Ny + 2, Nz + 2];
            double[,,] f_y = new double[Nx + 2, Ny + 1, Nz + 2];
            double[,,] f_z = new double[Nx + 2, Ny + 2, Nz + 1];




            // Create FFD solver and domain
            if (ffd == null || resetFFD)
            {
                omega = new WindInflow(Nx + 2, Ny + 2, Nz + 2, xyzsize[0], xyzsize[1], xyzsize[2], Vmet, terrain);
                foreach (double[] geo in geom)
                {
                    omega.add_obstacle(geo[0], geo[1], geo[2], geo[3], geo[4], geo[5]);
                }

                ffd = new FluidSolver(omega, dt, nu, u0, v0, w0, solver_params);
                de = new DataExtractor(omega, ffd);
                t = 0;

                if (resetFFD) resetFFD = false;            //reset FFD solver and domain

                Rhino.RhinoApp.WriteLine("GRASSHOPPER FFD Air Flow Simulation.");
                Rhino.RhinoApp.WriteLine("GH Plug-in: https://github.com/christophwaibel/GH_Wind");
                Rhino.RhinoApp.WriteLine("FFD Solver: https://github.com/lukasbystricky/GSoC_FFD");
                Rhino.RhinoApp.WriteLine("________________________________________________________");
                Rhino.RhinoApp.WriteLine("...Domain initialized");
                Rhino.RhinoApp.WriteLine("________________________________________________________");
            }














            // *******************************************************************************************
            // Run each time GH is updated
            // *******************************************************************************************
            //run solver. the solving-loop (new timestep) is executed in Grasshopper with a timer-component.
            //!!!!!!!!!!!!!!! CHANGE
            if (run)
            {
                int counter = 0;
                FluidSolver[] ffd_old = new FluidSolver[m];

                File.AppendAllText(@"C:\\residuals.txt", "pmin; pmax; pavg; umin; umax; uavg; vmin; vmax; vavg; wmin; wmax; wavg;\n");
                while (t < t_end)
                {
                    Rhino.RhinoApp.WriteLine(Convert.ToString(t) + " of " + Convert.ToString(t_end));

                    double[,,] p_t2 = new double[ffd.p.GetLength(0), ffd.p.GetLength(1), ffd.p.GetLength(2)];
                    Array.Copy(ffd.p, 0, p_t2, 0, ffd.p.Length);
                    double[,,] u_t2 = new double[ffd.u.GetLength(0), ffd.u.GetLength(1), ffd.u.GetLength(2)];
                    Array.Copy(ffd.u, 0, u_t2, 0, ffd.u.Length);
                    double[,,] v_t2 = new double[ffd.v.GetLength(0), ffd.v.GetLength(1), ffd.v.GetLength(2)];
                    Array.Copy(ffd.v, 0, v_t2, 0, ffd.v.Length);
                    double[,,] w_t2 = new double[ffd.w.GetLength(0), ffd.w.GetLength(1), ffd.w.GetLength(2)];
                    Array.Copy(ffd.w, 0, w_t2, 0, ffd.w.Length);

                    ffd.time_step(f_x, f_y, f_z);
                    if (t > dt && calcres)
                    {
                        double[] p_residuals;
                        double[,,] p_t1 = ffd.p;
                        FastFluidSolverMT.Utilities.calculate_residuals(p_t1, p_t2, out p_residuals);
                        Rhino.RhinoApp.WriteLine("p residuals: {0};{1};{2}", p_residuals[0], p_residuals[1], p_residuals[2]);
                        double[] u_residuals;
                        double[,,] u_t1 = ffd.u;
                        FastFluidSolverMT.Utilities.calculate_residuals(u_t1, u_t2, out u_residuals);
                        Rhino.RhinoApp.WriteLine("u residuals: {0};{1};{2}", u_residuals[0], u_residuals[1], u_residuals[2]);
                        double[] v_residuals;
                        double[,,] v_t1 = ffd.v;
                        FastFluidSolverMT.Utilities.calculate_residuals(v_t1, v_t2, out v_residuals);
                        Rhino.RhinoApp.WriteLine("v residuals: {0};{1};{2}", v_residuals[0], v_residuals[1], v_residuals[2]);
                        double[] w_residuals;
                        double[,,] w_t1 = ffd.w;
                        FastFluidSolverMT.Utilities.calculate_residuals(w_t1, w_t2, out w_residuals);
                        Rhino.RhinoApp.WriteLine("w residuals: {0};{1};{2}", w_residuals[0], w_residuals[1], w_residuals[2]);

                        File.AppendAllText(@"C:\\residuals.txt", Convert.ToString(p_residuals[0]) + ";" + Convert.ToString(p_residuals[1]) + ";" + Convert.ToString(p_residuals[2]) + ";" +
                            Convert.ToString(u_residuals[0])+";"+Convert.ToString(u_residuals[1])+";"+Convert.ToString(u_residuals[2])+";"+
                            Convert.ToString(v_residuals[0])+";"+Convert.ToString(v_residuals[1])+";"+Convert.ToString(v_residuals[2])+";"+
                            Convert.ToString(w_residuals[0]) + ";" + Convert.ToString(w_residuals[1]) + ";" + Convert.ToString(w_residuals[2])+"\n");
                    }
 
                    if (t >= t_end - m * dt)
                    {
                        ffd_old[counter] = new FluidSolver(ffd);
                        counter++;
                    }

                    t += dt;
                }

                //averaging results 
                FluidSolver ffd_mean = new FluidSolver(ffd);
                ffd_mean.p = new double[ffd.p.GetLength(0), ffd.p.GetLength(1), ffd.p.GetLength(2)];
                ffd_mean.u = new double[ffd.u.GetLength(0), ffd.u.GetLength(1), ffd.u.GetLength(2)];
                ffd_mean.v = new double[ffd.v.GetLength(0), ffd.v.GetLength(1), ffd.v.GetLength(2)];
                ffd_mean.w = new double[ffd.w.GetLength(0), ffd.w.GetLength(1), ffd.w.GetLength(2)];
                for (int i = 0; i < ffd_mean.p.GetLength(0); i++)
                {
                    for(int j=0; j < ffd_mean.p.GetLength(1); j++)
                    {
                        for(int k=0; k < ffd_mean.p.GetLength(2); k++)
                        {
                            for (int u=0; u<counter; u++)
                            {
                                ffd_mean.p[i, j, k] += ffd_old[u].p[i, j, k]; 
                            }
                            ffd_mean.p[i, j, k] /= counter;
                        }
                    }
                }

                for (int i = 0; i < ffd_mean.u.GetLength(0); i++)
                {
                    for (int j = 0; j < ffd_mean.u.GetLength(1); j++)
                    {
                        for (int k = 0; k < ffd_mean.u.GetLength(2); k++)
                        {
                            for (int u = 0; u < counter; u++)
                            {
                                ffd_mean.u[i, j, k] += ffd_old[u].u[i, j, k];
                            }
                            ffd_mean.u[i, j, k] /= counter;
                        }
                    }
                }

                for (int i = 0; i < ffd_mean.v.GetLength(0); i++)
                {
                    for (int j = 0; j < ffd_mean.v.GetLength(1); j++)
                    {
                        for (int k = 0; k < ffd_mean.v.GetLength(2); k++)
                        {
                            for (int u = 0; u < counter; u++)
                            {
                                ffd_mean.v[i, j, k] += ffd_old[u].v[i, j, k];
                            }
                            ffd_mean.v[i, j, k] /= counter;
                        }
                    }
                }

                for (int i = 0; i < ffd_mean.w.GetLength(0); i++)
                {
                    for (int j = 0; j < ffd_mean.w.GetLength(1); j++)
                    {
                        for (int k = 0; k < ffd_mean.w.GetLength(2); k++)
                        {
                            for (int u = 0; u < counter; u++)
                            {
                                ffd_mean.w[i, j, k] += ffd_old[u].w[i, j, k];
                            }
                            ffd_mean.w[i, j, k] /= counter;
                        }
                    }
                }

                de = new DataExtractor(omega, ffd_mean);

            }


  



            // *******************************************************************************************
            // *******************************************************************************************
            // TO DO:   vtk export
            //   pp.export_data_vtk(String.Concat("lid_driven_cavity_", tstep, ".vtk"), Nx, Ny, Nz, tstep );
            //bool run2 = (bool)Component.Params.Input[5].Sources[0].VolatileData;
            //while (true)












            // *******************************************************************************************
            // Redraw on or off
            // *******************************************************************************************
            //return mean over m*dt, instead of only one snapshot
            if (writeresults)
            {
                //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                // I could move all this away, an only output de data extractor
                double[,,] p = new double[Nx, Ny, Nz];
                double[,,] vu = new double[Nx, Ny, Nz];
                double[,,] vv = new double[Nx, Ny, Nz];
                double[,,] vw = new double[Nx, Ny, Nz];

                double[,,] pstag = new double[Nx + 1, Ny + 1, Nz + 1];
                double[,,] vustag = new double[Nx + 1, Ny + 1, Nz + 1];
                double[,,] vvstag = new double[Nx + 1, Ny + 1, Nz + 1];
                double[,,] vwstag = new double[Nx + 1, Ny + 1, Nz + 1];

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
                            }
                            else
                            {
                                p[i, j, k] = 0;
                                vu[i, j, k] = 0;
                                vv[i, j, k] = 0;
                                vw[i, j, k] = 0;

                            }
                            pstag[i, j, k] = de.get_pressure(i * omega.hx, j * omega.hy, k * omega.hz);
                            double[] velcen = de.get_velocity(i * omega.hx, j * omega.hy, k * omega.hz);
                            vustag[i, j, k] = velcen[0];
                            vvstag[i, j, k] = velcen[1];
                            vwstag[i, j, k] = velcen[2];
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

                List<double[,,]> veloutCen = new List<double[,,]> { };
                veloutCen.Add(vu);
                veloutCen.Add(vv);
                veloutCen.Add(vw);

                List<double[,,]> veloutStag = new List<double[,,]> { };
                veloutStag.Add(vustag);
                veloutStag.Add(vvstag);
                veloutStag.Add(vwstag);



                DA.SetDataList(0, veloutCen);
                DA.SetData(1, p);
                DA.SetDataList(2, veloutStag);
                DA.SetData(3, pstag);
                DA.SetData(4, de);
            }
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
                return Properties.Resources.solver;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("b9713e08-6665-4c3b-b1c1-18ac1d970874"); }
        }
    }
}