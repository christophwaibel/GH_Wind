using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;
using System.Drawing;

namespace GHWind
{
    public class GHVisualizerVectors : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the GHTest class.
        /// </summary>
        public GHVisualizerVectors()
            : base("Vector Visualizer", "Vector Visualizer",
                "Dynamic vector visualizer for the FFD solver. Draws and updates velocity vectors (and pressure values) at every timestep.",
                "EnergyHubs", "Wind Simulation")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddPointParameter("origin", "origin", "origin", GH_ParamAccess.item);
            pManager.AddIntegerParameter("xyz section", "xyz", "which section? x=0, y=1, z=2", GH_ParamAccess.item);
            pManager.AddIntegerParameter("section position", "section", "where exactly to draw section", GH_ParamAccess.item); 
            pManager.AddNumberParameter("min", "min", "minimum value. needed for colour gradient", GH_ParamAccess.item);
            pManager.AddNumberParameter("max", "max", "maximum value. needed for colour gradient", GH_ParamAccess.item);

            pManager.AddGenericParameter("velocity", "velocity", "velocity", GH_ParamAccess.list);
            pManager.AddGenericParameter("pressure", "pressure", "pressure", GH_ParamAccess.item);

            pManager.AddNumberParameter("hx", "hx", "hx", GH_ParamAccess.item);
            pManager.AddNumberParameter("hy", "hy", "hy", GH_ParamAccess.item);
            pManager.AddNumberParameter("hz", "hz", "hz", GH_ParamAccess.item);

            pManager.AddBooleanParameter("show pressure", "show pressure", "colour vectors with pressure values? Otherwise its velocity. (or black).", GH_ParamAccess.item);
            pManager[10].Optional = true;

            pManager.AddIntegerParameter("colour sheme", "colours", "Colour sheme. 0: Blue (min) - Red - Yellow (max); 1: Blue (min) - Green - Red (max); 2: Black only.", GH_ParamAccess.item);
            pManager[11].Optional = true;

            pManager.AddNumberParameter("scale length", "scale length", "scale length of vectors", GH_ParamAccess.item);
            pManager[12].Optional = true;
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddLineParameter("draw lines", "draw lines", "draw coloured lines", GH_ParamAccess.item);
            pManager.HideParameter(0);
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
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            //List<Line> lines = new List<Line>();
            //DA.GetDataList(0, lines);
            Point3d origin = Point3d.Unset;
            if (!DA.GetData(0, ref origin)) { origin = new Point3d(0,0,0); }

            int xyz = 1;
            if (!DA.GetData(1, ref xyz)) { xyz = 1; }

            int sectionheight = 1;
            if (!DA.GetData(2, ref sectionheight)) { sectionheight = 1; }

            //double scale = double.NaN;
            //if (!DA.GetData(3, ref scale)) { scale = 1; }

            double low = double.NaN;
            if (!DA.GetData(3, ref low)) { low = 0; }

            double top = double.NaN;
            if (!DA.GetData(4, ref top)) { top = 15; }

            double[, ,] vu, vv, vw;
            List<double[, ,]> vel = new List<double[, ,]> { };
            DA.GetDataList(5, vel);
            vu = vel[0];
            vv = vel[1];
            vw = vel[2];

            double[, ,] p = new double[,,] { };
            DA.GetData(6, ref p);

            double hx = double.NaN;
            double hy = double.NaN;
            double hz =double.NaN;
            if (!DA.GetData(7, ref hx)){return;}
            if (!DA.GetData(8,ref hy)){return;}
            if (!DA.GetData(9, ref hz)){return;}

            bool dispP = false;
            DA.GetData(10, ref dispP);

            int colourSheme = 0;
            DA.GetData(11, ref colourSheme);

            double scale = 1;
            DA.GetData(12, ref scale);



            Color c = new Color();

            double minp = double.MaxValue;
            double maxp = double.MinValue;
            for (int i = 0; i < p.GetLength(0); i++)
            {
                for (int j = 0; j < p.GetLength(1); j++)
                {
                    for (int k = 0; k < p.GetLength(2); k++)
                    {
                        if (minp > p[i, j, k]) minp = p[i, j, k];
                        if (maxp < p[i, j, k]) maxp = p[i, j, k];
                    }
                }
            }


            switch (xyz)
            {
                case 0: //x
                    for (int j = 0; j < vv.GetLength(1); j++)
                    {
                        for (int k = 0; k < vw.GetLength(2); k++)
                        {
                            Line arrowlines = new Line(new Point3d(sectionheight * hx + 0.5 * hx + origin[0], j * hy + 0.5 * hy + origin[1], k * hz + 0.5 * hz + origin[2]),
                                     new Vector3d(vu[sectionheight, j, k] * scale, vv[sectionheight, j, k] * scale, vw[sectionheight, j, k] * scale));
                            _lines.Add(arrowlines);
                            _widths.Add(1);

                            double quantity = arrowlines.Length / scale;
                            if (dispP)
                            {
                                quantity = p[sectionheight, j, k];
                            }
                            c = Utilities.GetRGB(colourSheme, quantity, top, low);
                            _colors.Add(c);
                        }
                    }
                    break;

                case 1: //y
                    for (int i = 0; i < vu.GetLength(0); i++)
                    {
                        for (int k = 0; k < vw.GetLength(2); k++)
                        {
                            Line arrowlines = new Line(new Point3d(i * hx + 0.5 * hx + origin[0], sectionheight * hy + 0.5 * hy + origin[1], k * hz + 0.5 * hz + origin[2]),
                                     new Vector3d(vu[i, sectionheight, k] * scale, vv[i, sectionheight, k] * scale, vw[i, sectionheight, k] * scale));
                            _lines.Add(arrowlines);
                            _widths.Add(1);

                            double quantity = arrowlines.Length/scale;
                            if (dispP)
                            {
                                quantity = p[i, sectionheight, k];
                            }
                            c = Utilities.GetRGB(colourSheme, quantity, top, low);
                            _colors.Add(c);
                        }
                    }
                    break;

                case 2: //z
                    for (int i = 0; i < vu.GetLength(0); i++)
                    {
                        for (int j = 0; j < vv.GetLength(1); j++)
                        {
                            Line arrowlines = new Line(new Point3d(i * hx + 0.5 * hx + origin[0], j * hy + 0.5 * hy + origin[1], sectionheight * hz + 0.5 * hz + origin[2]),
                                     new Vector3d(vu[i, j, sectionheight] * scale, vv[i, j, sectionheight] * scale, vw[i, j, sectionheight] * scale));
                            _lines.Add(arrowlines);
                            _widths.Add(1);

                            double quantity = arrowlines.Length/ scale;
                            if (dispP)
                            {
                                quantity = p[i, j, sectionheight];
                                //if (minp < 0) quantity += Math.Abs(minp);
                                //low = 0;
                                //top = maxp + Math.Abs(minp);
                                //third = (top - low) / 5;
                            }
                            c = Utilities.GetRGB(colourSheme, quantity, top, low);
                            _colors.Add(c);
                        }
                    }
                    break;


            }


            
          
        }

        public override void DrawViewportWires(IGH_PreviewArgs args)
        {
            base.DrawViewportWires(args);

            for (int i = 0; i < _lines.Count; i++)
            {
                args.Display.DrawLine(_lines[i], _colors[i],_widths[i]);
                args.Display.DrawArrowHead(_lines[i].To, _lines[i].Direction, _colors[i], _lines[i].Length * 2.0, 0);
            }
            Rhino.RhinoDoc.ActiveDoc.Views.Redraw();
            Rhino.RhinoApp.Wait();
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
                //return null;
                return GHWind.Properties.Resources.visu_vec;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("{ac7fa7e2-9f81-4f5c-8858-b9d76b08acc1}"); }
        }
    }
}