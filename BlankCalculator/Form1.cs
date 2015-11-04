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

            List<double[]> FixedPoints = new List<double[]>();
            double[] oVecPlane = new double[] { 0, 0, 0, 0 };
            string StlPath = CAT.SaveSurfaceAndFixedPoints(ref FixedPoints, ref oVecPlane);

            List<double[]> Vertices = new List<double[]>();
            List<int[]> Triangles = new List<int[]>();

            STLReader.STLRead(StlPath, ref Vertices, ref Triangles);
            if (Vertices.Count == 0 || Triangles.Count == 0) {
                MessageBox.Show("Não foi encontrada nenhuma informação na superficie exportada. Tente outra vez.");
                Environment.Exit(0);
            }

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

            double[,] MatrixA = new double[Vertices.Count * 2, Vertices.Count * 2];

            for (int i = 0; i < Vertices.Count; i++) {
                int[] CurTri = new int[] { 0, 0, 0 };
                int indexDown0 = -1;
                foreach (int[] t in Triangles) {
                    for (int j = 0; j < 3; j++) {
                        if (t[j] == i) { indexDown0 = j; break; }
                    }
                    if (indexDown0 != -1) { CurTri = t; break; }
                }
                int indexDown1 = indexDown0 - 1;
                if (indexDown1 < 0) indexDown1 = 2;
                int indexDown2 = indexDown1 - 1;
                if (indexDown2 < 0) indexDown2 = 2;

                Point3D vd0, vd1, vd2;
                vd0 = new Point3D(Vertices[CurTri[indexDown0]]);
                vd1 = new Point3D(Vertices[CurTri[indexDown1]]);
                vd2 = new Point3D(Vertices[CurTri[indexDown2]]);

                Line3D lU, lD;
                lU = new Line3D(vd1, vd0);
                lD = new Line3D(vd1, vd2);

                double Ang = lU.Direction.AngleTo(lD.Direction).Radians;
                double Len = lU.Length / lD.Length;


                MatrixA[CurTri[indexDown0] * 2, CurTri[indexDown0] * 2] = 1;
                MatrixA[CurTri[indexDown0] * 2, CurTri[indexDown0] * 2 + 1] = 0;

                MatrixA[CurTri[indexDown0] * 2, CurTri[indexDown1] * 2] = Len * Math.Cos(Ang) - 1;
                MatrixA[CurTri[indexDown0] * 2, CurTri[indexDown1] * 2 + 1] = -Len * Math.Sin(Ang);

                MatrixA[CurTri[indexDown0] * 2, CurTri[indexDown2] * 2] = -Len * Math.Cos(Ang);
                MatrixA[CurTri[indexDown0] * 2, CurTri[indexDown2] * 2 + 1] = Len * Math.Sin(Ang);


                MatrixA[CurTri[indexDown0] * 2 + 1, CurTri[indexDown0] * 2] = 0;
                MatrixA[CurTri[indexDown0] * 2 + 1, CurTri[indexDown0] * 2 + 1] = 1;

                MatrixA[CurTri[indexDown0] * 2 + 1, CurTri[indexDown1] * 2] = Len * Math.Sin(Ang);
                MatrixA[CurTri[indexDown0] * 2 + 1, CurTri[indexDown1] * 2 + 1] = Len * Math.Cos(Ang) - 1;

                MatrixA[CurTri[indexDown0] * 2 + 1, CurTri[indexDown2] * 2] = -Len * Math.Sin(Ang);
                MatrixA[CurTri[indexDown0] * 2 + 1, CurTri[indexDown2] * 2 + 1] = -Len * Math.Cos(Ang);
            }

            double[,] MatrixCa = new double[IndiceOfFixedPoints.Count * 2, Vertices.Count * 2];
            double[] VectorR = new double[IndiceOfFixedPoints.Count * 2];

            Point3D oRoot = new Point3D(new double[] { oVecPlane[0], oVecPlane[1], oVecPlane[2] });
            UnitVector3D vDir1 = new UnitVector3D(new double[] { oVecPlane[3], oVecPlane[4], oVecPlane[5] });
            UnitVector3D vDir2 = new UnitVector3D(new double[] { oVecPlane[6], oVecPlane[7], oVecPlane[8] });
            Plane oPlane = new Plane(vDir1.CrossProduct(vDir2),  oRoot);

            for (int i = 0; i < IndiceOfFixedPoints.Count; i++) {
                MatrixCa[i * 2, IndiceOfFixedPoints[i]*2] = 1;
                VectorR[i * 2] = new Point3D(Vertices[IndiceOfFixedPoints[i]]).ProjectOn(oPlane).X;
                MatrixCa[i * 2 + 1, IndiceOfFixedPoints[i]*2+1] = 1;
                VectorR[i * 2 + 1] = new Point3D(Vertices[IndiceOfFixedPoints[i]]).ProjectOn(oPlane).Y;
            }

            Matrix<double> Ca = Matrix<double>.Build.DenseOfArray(MatrixCa);
            Console.WriteLine(Ca);
            Vector<double> R = Vector<double>.Build.DenseOfArray(VectorR);
            Console.WriteLine(R);
            Matrix<double> A = Matrix<double>.Build.DenseOfArray(MatrixA);
            Console.WriteLine(A);
            double Penalty = 10;

            Matrix<double> Ak = A.Transpose() * A + Penalty * Ca.Transpose() * Ca;
            Console.WriteLine(Ak);
            Vector<double> X = Ak.Solve(Penalty * Ca.Transpose() * R);
            Console.WriteLine(X);

            CAT.PrintTriangles(X,Triangles, oRoot, vDir1, vDir2);


        }
    }
}
