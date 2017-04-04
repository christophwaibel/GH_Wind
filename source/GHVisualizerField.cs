using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;
using System.Drawing;

/*
 * GHVisualizerField.cs
 * Copyright 2017 Christoph Waibel <chwaibel@student.ethz.ch>
 * 
 * This work is licensed under the GNU GPL license version 3.
*/

namespace GHWind
{
    public class GHVisualizerField : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the GHVisualizerField class.
        /// </summary>
        public GHVisualizerField()
            : base("Field Visualizer", "Field Visualizer",
                "Dynamic field visualizer for the FFD solver. Draws and updates velocity and pressure values on a field at every timestep.",
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
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddMeshParameter("value field", "value field", "section showing pressure or velocity values", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Point3d origin = Point3d.Unset;
            if (!DA.GetData(0, ref origin)) { origin = new Point3d(0, 0, 0); }

            int xyz = 1;
            if (!DA.GetData(1, ref xyz)) { xyz = 1; }

            int sectionheight = 1;
            if (!DA.GetData(2, ref sectionheight)) { sectionheight = 1; }

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
            double hz = double.NaN;
            if (!DA.GetData(7, ref hx)) { return; }
            if (!DA.GetData(8, ref hy)) { return; }
            if (!DA.GetData(9, ref hz)) { return; }

            bool dispP = false;
            DA.GetData(10, ref dispP);

            int colourSheme = 0;
            DA.GetData(11, ref colourSheme);



            //min max pressure values
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



            Point3f[][] MVert = new Point3f[vu.GetLength(0)][];
            Color[][] Cols = new Color[vu.GetLength(0)][];
            int[][] index = new int[vu.GetLength(0)][];
            int counter = 0;

            Mesh MshColSection = new Mesh();
            switch (xyz)
            {
                case 0: //x
                    //8---9---10--11
                    //| x | x | x |
                    //4---5---6---7
                    //| x | x | x |
                    //0---1---2---3
                    //
                    //- loop through and get colours for all vertices
                    //- make vertices as point3d
                    MVert = new Point3f[vu.GetLength(1)][];
                    Cols = new Color[vu.GetLength(1)][];
                    index = new int[vu.GetLength(1)][];
                    counter = 0;

                    for (int j = 0; j < vu.GetLength(1); j++)
                    {
                        MVert[j] = new Point3f[vu.GetLength(2)];
                        Cols[j] = new Color[vu.GetLength(2)];
                        index[j] = new int[vu.GetLength(2)];

                        for (int k = 0; k < vu.GetLength(2); k++)
                        {
                            MVert[j][k] = new Point3f((float)(sectionheight * hx + origin[0]), (float)(j * hy + origin[1]), (float)(k * hz + origin[2]));
                            MshColSection.Vertices.Add(MVert[j][k]);

                            Cols[j][k] = new Color();
                            index[j][k] = counter;
                            counter++;

                            double quantity = 0;
                            if (dispP)
                            {
                                quantity = p[sectionheight, j, k];
                                //if (minp <= 0) quantity += Math.Abs(minp);
                                //low = 0;
                                //top = maxp + Math.Abs(minp);
                                //third = (top - low) / 5;
                            }
                            else
                            {
                                Line arrowlines = new Line(new Point3d(sectionheight * hx + origin[0], j * hy + origin[1], k * hz + origin[2]),
                                        new Vector3d(vu[sectionheight, j, k], vv[sectionheight, j, k], vw[sectionheight, j, k]));
                                quantity = arrowlines.Length;
                            }


                            Cols[j][k] = Utilities.GetRGB(colourSheme, quantity, top, low);
                            MshColSection.VertexColors.SetColor(index[j][k], Cols[j][k]);
                        }
                    }



                    for (int j = 0; j < vu.GetLength(1) - 1; j++)
                    {
                        for (int k = 0; k < vu.GetLength(2) - 1; k++)
                        {
                            MshColSection.Faces.AddFace(index[j][k], index[j + 1][k], index[j + 1][k + 1], index[j][k + 1]);
                        }
                    }
                    break;

                case 1: //y
                    //8---9---10--11
                    //| x | x | x |
                    //4---5---6---7
                    //| x | x | x |
                    //0---1---2---3
                    //
                    //- loop through and get colours for all vertices
                    //- make vertices as point3d
                    MVert = new Point3f[vu.GetLength(0)][];
                    Cols = new Color[vu.GetLength(0)][];
                    index = new int[vu.GetLength(0)][];
                    counter = 0;

                    for (int i = 0; i < vu.GetLength(0); i++)
                    {
                        MVert[i] = new Point3f[vu.GetLength(2)];
                        Cols[i] = new Color[vu.GetLength(2)];
                        index[i] = new int[vu.GetLength(2)];

                        for (int k = 0; k < vu.GetLength(2); k++)
                        {
                            MVert[i][k] = new Point3f((float)(i * hx + origin[0]), (float)(sectionheight * hy + origin[1]), (float)(k * hz + origin[2]));
                            MshColSection.Vertices.Add(MVert[i][k]);

                            Cols[i][k] = new Color();
                            index[i][k] = counter;
                            counter++;

                            double quantity = 0;
                            if (dispP)
                            {
                                quantity = p[i, sectionheight, k];
                                //if (minp <= 0) quantity += Math.Abs(minp);
                                //low = 0;
                                //top = maxp + Math.Abs(minp);
                                //third = (top - low) / 5;
                            }
                            else
                            {
                                Line arrowlines = new Line(new Point3d(i * hx + origin[0], sectionheight * hy + origin[1], k * hz + origin[2]),
                                        new Vector3d(vu[i, sectionheight, k], vv[i, sectionheight, k], vw[i, sectionheight, k]));
                                quantity = arrowlines.Length;
                            }


                            Cols[i][k] = Utilities.GetRGB(colourSheme, quantity, top, low);
                            MshColSection.VertexColors.SetColor(index[i][k], Cols[i][k]);
                        }
                    }



                    for (int i = 0; i < vu.GetLength(0) - 1; i++)
                    {
                        for (int k = 0; k < vu.GetLength(2) - 1; k++)
                        {
                            MshColSection.Faces.AddFace(index[i][k], index[i + 1][k], index[i + 1][k + 1], index[i][k + 1]);
                        }
                    }
                    break;

                case 2: //z
                    //8---9---10--11
                    //| x | x | x |
                    //4---5---6---7
                    //| x | x | x |
                    //0---1---2---3
                    //
                    //- loop through and get colours for all vertices
                    //- make vertices as point3d
                     MVert = new Point3f[vu.GetLength(0)][];
                    Cols = new Color[vu.GetLength(0)][];
                    index = new int[vu.GetLength(0)][];
                    counter = 0;


                    for (int i = 0; i < vu.GetLength(0); i++)
                    {
                        MVert[i] = new Point3f[vv.GetLength(1)];
                        Cols[i] = new Color[vv.GetLength(1)];
                        index[i] = new int[vv.GetLength(1)];

                        for (int j = 0; j < vv.GetLength(1); j++)
                        {
                            MVert[i][j] = new Point3f((float)(i * hx + origin[0]), (float)(j * hy + origin[1]), (float)(sectionheight * hz + origin[2]));
                            MshColSection.Vertices.Add(MVert[i][j]);

                            Cols[i][j] = new Color();
                            index[i][j] = counter;
                            counter++;

                            double quantity = 0;
                            if (dispP)
                            {
                                quantity = p[i, j, sectionheight];
                                //if (minp <= 0) quantity += Math.Abs(minp);
                                //low = 0;
                                //top = maxp + Math.Abs(minp);
                                //third = (top - low) / 5;
                            }
                            else
                            {
                                Line arrowlines = new Line(new Point3d(i * hx + origin[0], j * hy + origin[1], sectionheight * hz + origin[2]),
                                        new Vector3d(vu[i, j, sectionheight], vv[i, j, sectionheight], vw[i, j, sectionheight]));
                                quantity = arrowlines.Length;
                            }


                            Cols[i][j] = Utilities.GetRGB(colourSheme, quantity, top, low);
                            MshColSection.VertexColors.SetColor(index[i][j], Cols[i][j]);
                        }
                    }



                    for (int i = 0; i < vu.GetLength(0) - 1; i++)
                    {
                        for (int j = 0; j < vv.GetLength(1) - 1; j++)
                        {
                            MshColSection.Faces.AddFace(index[i][j], index[i + 1][j], index[i + 1][j + 1], index[i][j + 1]);
                        }
                    }
                    break;
            }

            

            DA.SetData(0, MshColSection);

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
                return GHWind.Properties.Resources.visu_field;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("{4525ceb8-47cc-4f9e-a0f1-e4dfedb84f53}"); }
        }
    }
}