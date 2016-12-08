using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;
using System.IO;
namespace GHWind
{
    public class GHExportGeometry : GH_Component
    {

        bool export;

        public GHExportGeometry()
            : base("Export Geometry",  "Export geometry",
                "Export geometry into a *.csv file",
                "EnergyHubs", "Wind Simulation")
        {
        }


        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("cubes as doubles", "cubes", "cubes as list of double[xmin, xmax, ymin, ymax, zmin, zmax]", GH_ParamAccess.list);
            pManager.AddTextParameter("path", "path", "path to export geometry data to", GH_ParamAccess.item);
            pManager.AddBooleanParameter("export?", "export?", "export data? use a button", GH_ParamAccess.item);
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
        }


        protected override void SolveInstance(IGH_DataAccess DA)
        {

            List<double[]> geom = new List<double[]>();
            if (!DA.GetDataList(0, geom)) { return; };


            string path = null;
            if (!DA.GetData(1, ref path)) { return; }

            DA.GetData(2, ref export);

            
            
            
            //EXPORT GEOMETRY
            if (export)
            {

                string[] lines;
                var list = new List<string>();
                foreach (double[] geo in geom)
                {
                    string line = geo[0].ToString() + "," + geo[1].ToString() + "," + geo[2].ToString() + "," + geo[3].ToString() + "," + geo[4].ToString() + "," + geo[5].ToString();
                    list.Add(line);
                }
                lines = list.ToArray();
                File.WriteAllLines(path, lines);
                export = false;
            }

        }




        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                // return Resources.IconForThisComponent;
                return GHWind.Properties.Resources.discr_export;
            }
        }



        public override Guid ComponentGuid
        {
            get { return new Guid("{9d1d9ace-f25a-4bf8-8282-660335fd2bd5}"); }
        }
    }
}