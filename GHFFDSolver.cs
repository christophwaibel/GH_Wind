using System;
using System.Collections.Generic;
using System.Drawing;

using Grasshopper.Kernel;
using Rhino.Geometry;
using FastFluidSolver;

namespace GHWind
{
    public class GHFFDSolver : GH_Component
    {
        GH_Document doc;
        IGH_Component Component;
        Domain omega;
        FluidSolver ffd;
        DataExtractor de;
        Form1 f;
        int tstep;
        double t;

        Visualizer visualize;

        /// <summary>
        /// Each implementation of GH_Component must provide a public 
        /// constructor without any arguments.
        /// Category represents the Tab in which the component will appear, 
        /// Subcategory the panel. If you use non-existing tab or panel names, 
        /// new tabs/panels will automatically be created.
        /// </summary>
        public GHFFDSolver()
            : base("FFD", "FFD",
            "Fast Fluid Dynamics Solver",
            "EnergyHubs", "Wind Simulation")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddBoxParameter("Domain", "Domain", "Fluid Domain, as 3D Box", GH_ParamAccess.item);
            pManager[0].Optional = true;

            pManager.AddIntegerParameter("Domain discretization", "Nx,Ny,Nz", "Domain discretization Nx, Ny, Nz.", GH_ParamAccess.list);
            //later make another input for defining more precisely the domain. like, internal flow, external flow, inflows, outflows...)



            //write another seperate component, which does the discretization
            pManager.AddIntegerParameter("obstacles", "obstacles", "Add obstacles from Obstacle-Component", GH_ParamAccess.list);



            pManager.AddNumberParameter("Time Step", "dt", "Calculation time step", GH_ParamAccess.item);

            pManager.AddNumberParameter("Wind Speed", "vel", "Wind Speed in m/s", GH_ParamAccess.item);

            pManager.AddNumberParameter("Wind Direction", "WindDir", "Wind Direction, 0-360°", GH_ParamAccess.item);
            pManager[5].Optional = true;

            pManager.AddBooleanParameter("Run?", "Run?", "Run the solver and stop it", GH_ParamAccess.item);

            pManager.AddBooleanParameter("Write GH Results", "Results", "Write Results for visualization in GH", GH_ParamAccess.item);
            pManager[7].Optional = true;

            pManager.AddIntegerParameter("xyz", "xyz", "Draw cells on x(0), y(1) or z(2) section", GH_ParamAccess.item);
            pManager[8].Optional = true;

            pManager.AddIntegerParameter("Quantity", "Quantity", "Velocity (0), Pressure (1), or Species Concentration (2)", GH_ParamAccess.item);
            pManager[8].Optional = true;

            pManager.AddBooleanParameter("Export VTK", "ExpVTK", "Export Results to VTK", GH_ParamAccess.item);
            pManager[9].Optional = true;

            pManager.AddBooleanParameter("Reset", "Reset", "Reset domain", GH_ParamAccess.item);
            pManager[10].Optional = true;

        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddLineParameter("velocities", "velocities", "add velocities", GH_ParamAccess.item);
            pManager.HideParameter(0);


           // pManager.AddTextParameter("VTK path", "VTK path", "Output path of VTK results file", GH_ParamAccess.item);
        }


        private readonly List<Line> _lines = new List<Line>();
        private readonly List<int> _widths = new List<int>();
        private readonly List<Color> _colors = new List<Color>();
        protected override void BeforeSolveInstance()
        {

            _lines.Clear();
            _widths.Clear();
            _colors.Clear();

        }


        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object can be used to retrieve data from input parameters and 
        /// to store data in output parameters.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // *********************************************************************************
            //grasshopper current document
            // *********************************************************************************
            Component = this;
            doc = Component.OnPingDocument();



            DA.SetData(0, new Line(new Point3d(0, 0, 0), new Vector3d(1, 1, 1)));


            // *********************************************************************************
            // get info from grasshopper
            // *********************************************************************************
            List<int> Nxyz = new List<int>();
            if (!DA.GetDataList(1, Nxyz)) { return; };
            int Nx = Nxyz[0];
            int Ny = Nxyz[1];
            int Nz = Nxyz[2];
            //Nx = 50;
            //Ny =20;
            //Nz = 10;


            //input: 1d list of 0,1 for obstacles. should be of size Nx*Ny*Nz
            List<int> obst = new List<int>();
            if (!DA.GetDataList(2, obst)) { return; };
            int[, ,] obstacles = new int[Nx + 2, Ny + 2, Nz + 2];
            for (int x = 0; x < Nx; x++)
            {
                for (int y = 0; y < Ny; y++)
                {
                    for (int z = 0; z < Nz; z++)
                    {
                        //obstacles[x, y, z] = obst[x + y * Nx + z * Ny * Nx];        //from 1d-index to 3d-matrix
                        obstacles[x, y, z] = 0;
                    }
                }
            }





            // time step
            double dt = 0.01;
            if (!DA.GetData(3, ref dt)) { return; }


            // wind speed
            double windvel = 10;
            if (!DA.GetData(4, ref windvel)) { return; }

            // wind direction in degree. 0/360 is wind from north to south. 90 is from south to north
            //double winddir = 0;
            //if (!DA.GetData(3, ref winddir)) { return; }

            //double tf = 1;      //final timestep


            //double[, ,] u0 = new double[Nx + 2, Ny + 2, Nz + 2];
            //double[, ,] v0 = new double[Nx + 2, Ny + 2, Nz + 2];
            //double[, ,] w0 = new double[Nx + 2, Ny + 2, Nz + 2];

            bool resetFFD = false;
            if (DA.GetData(10, ref resetFFD)) { return; };









            // *********************************************************************************
            //from Lukas 
            // *********************************************************************************
            double nu = 0.001;

            double[, ,] u0 = new double[Nx + 1, Ny + 2, Nz + 2];    // initial x velocity
            double[, ,] v0 = new double[Nx + 2, Ny + 1, Nz + 2];    // initial y velocity
            double[, ,] w0 = new double[Nx + 2, Ny + 2, Nz + 1];    // initial z velocity

            double[, ,] f_x = new double[Nx + 1, Ny + 2, Nz + 2];
            double[, ,] f_y = new double[Nx + 2, Ny + 1, Nz + 2];
            double[, ,] f_z = new double[Nx + 2, Ny + 2, Nz + 1];


            FluidSolver.solver_struct solver_prams = new FluidSolver.solver_struct();
            solver_prams.tol = 1e-5;
            solver_prams.min_iter = 1;
            solver_prams.max_iter = 20;
            solver_prams.verbose = false;
            solver_prams.backtrace_order = 2;




            if (omega == null) omega = new DomainOpen(Nx + 2, Ny + 2, Nz + 2, 200, 100, 50, obstacles, windvel);

            if (ffd == null)
            {
                ffd = new FluidSolver(omega, dt, nu, u0, v0, w0, solver_prams);
                de = new DataExtractor(omega, ffd);
                tstep = 0;
                t = 0;

                //for (int i = 0; i < Nx; i++)
                //{
                //    for (int u = 0; u < Ny; u++)
                //    {
                //        //double[] vel = de.get_velocity((i - 0.5) * omega.hx, (u - 0.5) * omega.hy, (5 - 0.5) * omega.hz);
                //        //Line arrowlines = new Rhino.Geometry.Line(new Rhino.Geometry.Point3d(i * omega.hx, u * omega.hy, 5 * omega.hz),
                //        //     new Vector3d(vel[0], vel[1], vel[2]));
                //        Line arrowlines = new Line(new Point3d(i * omega.hx, u * omega.hy, 5 * omega.hz), new Vector3d(10, 10, 10));
                //        _lines.Add(arrowlines);
                //        _widths.Add(Math.Min(10, (int)(arrowlines.Length * 0.1) + 1));
                //        _colors.Add(Color.CornflowerBlue);
                //    }
                //}

                //DA.SetDataList(0, _lines);

            }
            else if (ffd != null && resetFFD == true)
            {
                ffd = new FluidSolver(omega, dt, nu, u0, v0, w0, solver_prams);
                de = new DataExtractor(omega, ffd);
                tstep = 0;
                t = 0;
                resetFFD = false;
            }












            bool run = false;
            if (!DA.GetData(6, ref run)) { return; }


            //if (omega == null)
            //{
            //    omega = new DomainOpen(Nx + 2, Ny + 2, Nz + 2);
            //}

            //if (ffd == null)
            //{
            //    ffd = new FluidSolver(omega, dt, nu, u0, v0, w0, false);
            //    tstep = 0; 
            //    t = 0;    
            //}
            // PostProcessor pp = new PostProcessor(ffd, omega);

            if (f == null)
            {
                f = new Form1();
                f.Show(Grasshopper.Instances.DocumentEditor);
                Grasshopper.Instances.DocumentEditor.FormShepard.RegisterForm(f);
                f.checkBox1.Text = "run the solver";




            }

            if (visualize == null && ffd != null)
            {
                visualize = new Visualizer(ffd, omega, Rhino.RhinoDoc.ActiveDoc);
                visualize.AddObstacleGeometry();
                visualize.AddDomainBoundingBoxGeometry();

                //visualize.AddArrows();
            }




            //_lines.Clear();
            //_widths.Clear();
            double top, low, third;
            top = 15;
            low = 0;
            third = (top - low) / 5;
            double RR, GG, BB;
            Color c = new Color();
            // *******************************************************************************************
            // *******************************************************************************************
            // TO DO:   fix this loop with
            //   pp.export_data_vtk(String.Concat("lid_driven_cavity_", tstep, ".vtk"), Nx, Ny, Nz, tstep );
            //bool run2 = (bool)Component.Params.Input[5].Sources[0].VolatileData;
            //while (true)
            while (tstep < 20)
            {
                //if (f.checkBox1.Checked == true)
                //{

                    //doc.NewSolution(false);
                    //Rhino.RhinoDoc.ActiveDoc.Views.Redraw();
                    //Rhino.RhinoApp.Wait();

                    t += dt;
                    tstep++;


                    Rhino.RhinoApp.WriteLine("Time t = {0}", t);
                    f.setlabel(t.ToString());
                    f.Update();

                    ffd.time_step(f_x, f_y, f_z);
                    //pp.export_data_vtk (String.Concat("lid_driven_cavity_", tstep, ".vtk"), Nx, Ny, Nz, tstep );  

                    if (tstep == 49)
                    {
                        //visualize.AddArrows(new int[] { 1, Nx }, new int[] { 1, Ny }, new int[] { 1, Nz });
                        //adding arrows for velocities
                        for (int i = 0; i < Nx; i++)
                        {
                            for (int u = 0; u < Ny; u++)
                            {
                                double[] vel = de.get_velocity((i - 0.5) * omega.hx, (u - 0.5) * omega.hy, (5 - 0.5) * omega.hz);
                                Line arrowlines = new Line(new Point3d(i * omega.hx, u * omega.hy, 5 * omega.hz),
                                     new Vector3d(vel[0], vel[1], vel[2]));
                                _lines.Add(arrowlines);
                                _widths.Add(Math.Min(10, (int)(arrowlines.Length * 0.1) + 1));
                                if (arrowlines.Length > third && arrowlines.Length <= 2 * third)
                                {
                                    RR = (arrowlines.Length - third) * (255 / third);
                                    GG = 0;
                                    BB = 255 - ((arrowlines.Length - third) * (255 / third));
                                }
                                else if (arrowlines.Length > 2 * third)
                                {
                                    RR = 255;
                                    GG = (arrowlines.Length - 2 * third) * (255 / third);
                                    BB = 0;
                                }
                                else
                                {
                                    RR = 0;
                                    GG = 0;
                                    BB = 255;
                                }
                                if (RR > 255) RR = 255;
                                if (GG > 255) GG = 255;
                                if (BB > 255) BB = 255;
                                c = Color.FromArgb((int)RR, (int)GG, (int)BB);
                                _colors.Add(c);

                            }
                        }
                    } 



                }
            //    else
            //    {
            //        Rhino.RhinoDoc.ActiveDoc.Views.Redraw();
            //        Rhino.RhinoApp.Wait();
            //    }



            //}


            

            //Rhino.RhinoDoc.ActiveDoc.Views.Redraw();
            //Rhino.RhinoApp.Wait();






            // ____________________________________________________________________
            // // THIS SHOWS HOW TO DYNAMICALLY USE GRASSHOPPER SLIDER INPUTS
            //List<Grasshopper.Kernel.Special.GH_NumberSlider> sliders = new List<Grasshopper.Kernel.Special.GH_NumberSlider>();
            //foreach (IGH_Param param in Component.Params.Input)
            //{
            //    Grasshopper.Kernel.Special.GH_NumberSlider slider = param.Sources[0] as Grasshopper.Kernel.Special.GH_NumberSlider;
            //    if (slider != null)
            //    {
            //        sliders.Add(slider);
            //        x = (int)slider.CurrentValue;
            //    }
            //}

            //while (terminate != true)
            //{
            //    doc.NewSolution(false);

            //    fx = x * 10;
            //    maximize(x, fx);
            //    sliders[0].TickValue = xNew;
            //    x = xNew;
            //}


        }



        public override void DrawViewportWires(IGH_PreviewArgs args)
        {
            base.DrawViewportWires(args);

            for (int i = 0; i < _lines.Count; i++)
            {
                args.Display.DrawLine(_lines[i], _colors[i], _widths[i]);
                args.Display.DrawArrowHead(_lines[i].To, _lines[i].Direction, _colors[i], 40, 0);
            }
            Rhino.RhinoDoc.ActiveDoc.Views.Redraw();
            Rhino.RhinoApp.Wait();
        }








        /// <summary>
        /// Provides an Icon for every component that will be visible in the User Interface.
        /// Icons need to be 24x24 pixels.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                // You can add image files to your project resources and access them like this:
                //return Resources.IconForThisComponent;
                return null;
            }
        }

        /// <summary>
        /// Each component must have a unique Guid to identify it. 
        /// It is vital this Guid doesn't change otherwise old ghx files 
        /// that use the old ID will partially fail during loading.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("{2f9d884f-31ed-4afd-921a-f407565448fa}"); }
        }
    }
}
