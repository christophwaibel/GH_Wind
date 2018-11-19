using System;
using System.Drawing;
using Grasshopper.Kernel;

namespace GHWind
{
    public class GHWindInfo : GH_AssemblyInfo
    {
        public override string Name
        {
            get
            {
                return "GHWind";
            }
        }
        public override Bitmap Icon
        {
            get
            {
                //Return a 24x24 pixel bitmap to represent this GHA library.
                return null;
            }
        }
        public override string Description
        {
            get
            {
                //Return a short string describing the purpose of this GHA library.
                return "";
            }
        }
        public override Guid Id
        {
            get
            {
                return new Guid("ebdd9491-bc5f-4223-b0fc-baedcfecad54");
            }
        }

        public override string AuthorName
        {
            get
            {
                //Return a string identifying you or your company.
                return "Empa";
            }
        }
        public override string AuthorContact
        {
            get
            {
                //Return a string representing your preferred contact details.
                return "";
            }
        }
    }
}
