using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Spatial.Euclidean;

namespace BlankCalculator {
    public partial class Form1 : Form {
        public Form1() {
            InitializeComponent();
        }
        public CATInterop CAT = new CATInterop();
        private void Form1_Load(object sender, EventArgs e) {

            //Plane and FixedPoints initializor
            List<double[]> FixedPoints = new List<double[]>();
            double[] oVecPlane = new double[] { 0, 0, 0, 0 };

            //Interop with catia to retrive exported stlPath
            string StlPath = CAT.SaveSurfaceAndFixedPoints(ref FixedPoints, ref oVecPlane);

            List<double[]> Vertices = new List<double[]>();
            List<int[]> TrianglesVertices = new List<int[]>();
            List<int[]> Edges = new List<int[]>();
            List<int[]> TrianglesEdges = new List<int[]>();

            //Parse stlFile to Vertice and triangle matrix
            STLReader.STLRead(StlPath, ref Vertices, ref TrianglesVertices,ref TrianglesEdges, ref Edges);
            if (Vertices.Count == 0 || TrianglesVertices.Count == 0) {
                MessageBox.Show("Não foi encontrada nenhuma informação na superficie exportada. Tente outra vez.");
                Environment.Exit(0);
            }

            //Detection of the most apropriated verticePoints
            double min = 10000000;
            double cur = -10000000;
            int iMin = -1;
            List<int> IndiceOfFixedPoints = new List<int>();
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
            int OneIndFix = IndiceOfFixedPoints[0];
            Point3D OnePointFix = new Point3D(new double[] { Vertices[OneIndFix][0], Vertices[OneIndFix][1], Vertices[OneIndFix][2] });
          
      //Calculation of Plane of flattening
             Point3D oRoot = new Point3D(new double[] { oVecPlane[0], oVecPlane[1], oVecPlane[2] });
            UnitVector3D vDir1 = new UnitVector3D(new double[] { oVecPlane[3], oVecPlane[4], oVecPlane[5] });
            UnitVector3D vDir2 = new UnitVector3D(new double[] { oVecPlane[6], oVecPlane[7], oVecPlane[8] });

            //Solve using the vertice to vertice interpretation of the article
            //VerticeToVerticeSolver vvSolver = new VerticeToVerticeSolver();
            //Vector<double> X1 = vvSolver.Solve(Vertices, TrianglesVertices, IndiceOfFixedPoints, oRoot, vDir1, vDir2);
            ////Print Result of triangles in CATIA;
            //CAT.PrintTriangles(X1, TrianglesVertices, oRoot, vDir1, vDir2);

            //Solve using the edge to edge interpretation of the article
            EdgeToEdgeSolver eeSolver = new EdgeToEdgeSolver();
            Vector<double> X2 = eeSolver.Solve(Vertices, TrianglesEdges,Edges, IndiceOfFixedPoints, oRoot, vDir1, vDir2);
            //Print Result of triangles in CATIA;
            
            CAT.PrintTriangles(X2, TrianglesVertices, oRoot, vDir1, vDir2, OneIndFix, OnePointFix);

            Application.Exit();
        }
    }
}
