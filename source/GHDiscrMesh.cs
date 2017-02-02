using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;

namespace GHWind
{
    public class GHDiscrMesh : GH_Component
    {

        public GHDiscrMesh()
            : base("Discretize Mesh", "Discretize Meshes",
                "Discretize Meshes into cubes. Geometry for FFD solver.",
                "EnergyHubs", "Wind Simulation")
        {
        }


        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddNumberParameter("xyzsize", "xyzsize", "xyzsize", GH_ParamAccess.list);
            pManager.AddIntegerParameter("xyzcells", "xyzcells", "xyzcells", GH_ParamAccess.list);
            pManager.AddMeshParameter("solid mesh", "solid mesh", "solid mesh, meaning it needs to be closed", GH_ParamAccess.list);
            pManager.AddPointParameter("origin", "origin", "origin", GH_ParamAccess.item);
            //pManager.AddBooleanParameter("export to a text file?", "export csv", "export the geometry into a csv file?", GH_ParamAccess.item);
            //pManager[4].Optional = true;
        }


        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddBoxParameter("boxes", "boxes", "boxes", GH_ParamAccess.list);
            pManager.AddGenericParameter("cubes as doubles", "cubes", "cubes", GH_ParamAccess.list);
            //pManager.AddTextParameter("export path", "export path", "path of the csv file", GH_ParamAccess.item);
        }


        protected override void SolveInstance(IGH_DataAccess DA)
        {
            //Script to discretize a topography-mesh into cubes for the FFD solver
            //  2016-09-07  Christoph Waibel
            //



            //input:
            //  x,y,z grid size
            //  x,y,z cell numbers
            //  origin point (all calculations use this as reference)
            //  mesh geometry
            List<double> xyzsize = new List<double>();
            if (!DA.GetDataList(0, xyzsize)) { return; };
            List<int> xyzcells = new List<int>();
            if (!DA.GetDataList(1, xyzcells)) { return; };
            List<Mesh> meshs = new List<Mesh>();
            if (!DA.GetDataList(2, meshs)) { return; }
            Point3d origin = Point3d.Unset;
            if (!DA.GetData(3, ref origin)) { return; }

            double tolerance = 0.01;

            //output:
            //  list of cubes (as double [6]{xmin, xmax, ymin, ymax, zmin, zmax})

            //1. create list of points on grid cells 0-point
            List<Point3d> pts = new List<Point3d>();
            double xcelldist = xyzsize[0] / xyzcells[0];
            double ycelldist = xyzsize[1] / xyzcells[1];
            double zcelldist = xyzsize[2] / xyzcells[2];
            for (int x = 0; x < xyzcells[0]; x++)
            {
                for (int y = 0; y < xyzcells[1]; y++)
                {
                    pts.Add(new Point3d((x + 0.5) * xcelldist + origin[0], (y + 0.5) * ycelldist + origin[1], origin[2]));
                }
            }

            Vector3d vec = new Vector3d(0, 0, 1);
            //foreach grid cell
            List<Box> box = new List<Box>();
            List<double[]> cubes = new List<double[]>();
            foreach (Point3d pt in pts)
            {
                foreach (Mesh mesh in meshs)
                {
                    //2. shoot rays straight up (from z-origin to zmax)
                    //Ray3d ray = new Ray3d(new Point3d(pt[0] + 0.5 * (double)xcelldist, pt[1] + 0.5 * (double)ycelldist, pt[2]), vec);
                    Ray3d ray = new Ray3d(new Point3d(pt[0], pt[1], pt[2]), vec);
                    //3. check for intersections with mesh
                    Point3d ptX1 = new Point3d(pt[0], pt[1], pt[2]);
                    bool blninters = false;
                    double inters1 = Rhino.Geometry.Intersect.Intersection.MeshRay(mesh, ray);


                    if (inters1 > 0)
                    {
                        Point3d ptat1 = ray.PointAt(inters1);
                        ray = new Ray3d(new Point3d(ptat1[0], ptat1[1], ptat1[2] + tolerance), vec);

                        //flatten point to the next/previous grid step
                        ptat1 = new Point3d(ptat1.X, ptat1.Y, Math.Round((ptat1[2] - origin[2]) / zcelldist, 0) * zcelldist + origin[2]);

                        double inters2 = Rhino.Geometry.Intersect.Intersection.MeshRay(mesh, ray);
                        if (inters2 > 0)
                        {
                            Point3d ptat2 = ray.PointAt(inters2);
                            blninters = true;
                            //repeat until there are no other intersections anymore
                            while (blninters == true)
                            {


                                //5. create cube with xmin xmax ymin ymax zmin (all already known), and zmax (calculated with mesh intersection)
                                //      rounding the intersection-z-value to match the grid-z-values
                                Point3d ptmax = new Point3d(ptat1[0] + 0.5*xcelldist, ptat1[1] +0.5* ycelldist, Math.Round((ptat2[2] - origin[2]) / zcelldist, 0) * zcelldist + origin[2]);
                                box.Add(new Box(new BoundingBox(new Point3d (ptat1[0]-0.5*xcelldist , ptat1 [1]-0.5*ycelldist , ptat1[2]), ptmax)));
                                cubes.Add(new double[] { ptat1[0] - origin[0] - 0.5 * xcelldist, ptmax[0] - origin[0], ptat1[1] - origin[1]-0.5*ycelldist , ptmax[1] - origin[1], ptat1[2] - origin[2], ptmax[2] - origin[2] });

                                //check another intersection couple  */------/*
                                ptat1 = new Point3d(ptat2[0], ptat2[1], ptat2[2] + tolerance);
                                ray = new Ray3d(ptat1, vec);
                                inters1 = Rhino.Geometry.Intersect.Intersection.MeshRay(mesh, ray);
                                if (inters1 > 0)
                                {
                                    ptat1 = ray.PointAt(inters1);
                                    ray = new Ray3d(new Point3d(ptat1[0], ptat1[1], ptat1[2] + tolerance), vec);
                                    inters2 = Rhino.Geometry.Intersect.Intersection.MeshRay(mesh, ray);
                                    if (inters2 > 0)
                                    {
                                        ptat2 = ray.PointAt(inters2);
                                    }
                                    else
                                    {
                                        blninters = false;
                                    }
                                }
                                else
                                {
                                    blninters = false;
                                }
                            }

                        }

                    }


                }


            }










            DA.SetDataList(0, box);
            DA.SetDataList(1, cubes);
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
                return GHWind.Properties.Resources.discr;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("{01f3853b-a2fa-419f-91e9-7258da235449}"); }
        }
    }
}