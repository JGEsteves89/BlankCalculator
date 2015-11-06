using MathNet.Spatial.Euclidean;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BlankCalculator {
    public class Mesh {
        public Mesh() {
            //Plane and FixedPoints initializor
            FixedPoints = new List<double[]>();
            oVecPlane = new double[] { 0, 0, 0, 0 };
            Vertices = new List<double[]>();
            FacetNormals = new List<double[]>();
            TrianglesVertices = new List<int[]>();
            Edges = new List<int[]>();
            TrianglesEdges = new List<int[]>();
            IndiceOfFixedPoints = new List<int>();
        }

        internal List<double[]> FixedPoints;
        internal double[] oVecPlane;
        internal List<double[]> Vertices;
        internal List<int[]> TrianglesVertices;
        internal List<int[]> TrianglesEdges;
        internal List<int[]> Edges;
        internal List<double[]> FacetNormals;
        internal Point3D oRoot;
        internal UnitVector3D vDir1;
        internal UnitVector3D vDir2;
        internal List<int> IndiceOfFixedPoints;
        internal Point3D OnePointFix;

        public int OneIndFix { get; internal set; }

        internal void ComputeFixedPoints() {
            double min = 10000000;
            double cur = -10000000;
            int iMin = -1;
            foreach (double[] FixPoint in FixedPoints) {
                min = 10000000;
                cur = -10000000;
                Point3D RefPt = new Point3D(FixPoint);
                for (int i = 0; i < Vertices.Count; i++) {
                    Point3D CurPt = new Point3D(Vertices[i]);
                    cur = RefPt.DistanceTo(CurPt);
                    if (cur < min) {
                        min = cur;
                        iMin = i;
                    }
                }
                IndiceOfFixedPoints.Add(iMin);
            }
            IndiceOfFixedPoints = IndiceOfFixedPoints.Distinct().ToList();
            OneIndFix = IndiceOfFixedPoints[0];
            OnePointFix = new Point3D(new double[] { Vertices[OneIndFix][0], Vertices[OneIndFix][1], Vertices[OneIndFix][2] });
        }

        internal void ComputePlaneofFlattening() {
            oRoot = new Point3D(new double[] { oVecPlane[0], oVecPlane[1], oVecPlane[2] });
            vDir1 = new UnitVector3D(new double[] { oVecPlane[3], oVecPlane[4], oVecPlane[5] });
            vDir2 = new UnitVector3D(new double[] { oVecPlane[6], oVecPlane[7], oVecPlane[8] });
        }
    }
}
