using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;

namespace GHWind
{
    public class GHPrintCellValues : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the GHPrintCellValues class.
        /// </summary>
        public GHPrintCellValues()
          : base("Print Cells Values", "Print Cell Values",
                            "Print all cell values of the entire domain as 1D-array in this order: foreach in Nx, foreach in Ny, foreach in Nz, print cell[x,y,z].",
                            "EnergyHubs", "Wind Simulation")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("obstacle array 3D", "obstacle array 3D", "obstacle array 3D", GH_ParamAccess.item);
            pManager.AddGenericParameter("u, v, w -vel cell center 3D", "u, v, w -vel cell center 3D", "u, v, w -vel cell center 3D", GH_ParamAccess.list);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("bool obstacle list 1D", "bool obstacle list 1D", "bool obstacle list 1D", GH_ParamAccess.list);
            pManager.AddGenericParameter("u-vel list 1D", "u-vel list 1D", "u-vel list 1D", GH_ParamAccess.list);
            pManager.AddGenericParameter("v-vel list 1D", "v-vel list 1D", "v-vel list 1D", GH_ParamAccess.list);
            pManager.AddGenericParameter("w-vel list 1D", "w-vel list 1D", "w-vel list 1D", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            int[,,] obst = new int[,,] { };
            DA.GetData(0, ref obst);

            
            List<double[,,]> vel = new List<double[,,]> { };
            DA.GetDataList(1, vel);
            //vu = vel[0];
            //vv = vel[1];
            //vw = vel[2];


            List<int> obstList = new List<int> { };
            List<double> uList = new List<double> { };
            List<double> vList = new List<double> { };
            List<double> wList = new List<double> { };
            for (int x=0; x < obst.GetLength(0)-2; x++)
            {
                for (int y=0; y<obst.GetLength(1)-2; y++)
                {
                    for(int z=0; z < obst.GetLength(2)-2; z++)
                    {
                        obstList.Add(obst[x+1, y+1, z+1]);
                        uList.Add(vel[0][x, y, z]);
                        vList.Add(vel[1][x, y, z]);
                        wList.Add(vel[2][x, y, z]);
                    }
                }
            }
            DA.SetDataList(0, obstList);
            DA.SetDataList(1, uList);
            DA.SetDataList(2, vList);
            DA.SetDataList(3, wList);
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
                return null;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("e02685d2-29e5-4548-b1fb-9416ec9af3fd"); }
        }
    }
}