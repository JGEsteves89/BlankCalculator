﻿using System;
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
            List<double[]> FacetNormals = new List<double[]>();
            List<int[]> TrianglesVertices = new List<int[]>();
            List<int[]> Edges = new List<int[]>();
            List<int[]> TrianglesEdges = new List<int[]>();

            //Parse stlFile to Vertice and triangle matrix
            STLReader.STLRead(StlPath, ref Vertices, ref TrianglesVertices, ref TrianglesEdges, ref Edges, ref FacetNormals);
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
            IndiceOfFixedPoints = IndiceOfFixedPoints.Distinct().ToList();
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
            EdgeToEdgeSolver2 eeSolver = new EdgeToEdgeSolver2();
            Vector<double> X2 = eeSolver.Solve(Vertices, TrianglesEdges, Edges, IndiceOfFixedPoints, oRoot, vDir1, vDir2);
            //Print Result of triangles in CATIA;

            CAT.PrintTriangles(X2, TrianglesVertices, oRoot, vDir1, vDir2, OneIndFix, OnePointFix);





            UnitVector3D dx, dy, dz;
            dx = new UnitVector3D(new double[] { 1, 0, 0 });
            dy = new UnitVector3D(new double[] { 0, 1, 0 });
            dz = new UnitVector3D(new double[] { 0, 0, 1 });


            UnitVector3D PlaneNormal = vDir1.CrossProduct(vDir2);

            double E = 200 * 10 ^ 6;
            double v = 0.3;
            double h = 1;

            double[,] MatrixK = new double[Vertices.Count * 3, Vertices.Count * 3];
            double[] Vectorq = new double[Vertices.Count * 3];

            double[,] MatrixD = new double[3, 3];
            MatrixD[0, 0] = (E / ((1 + v) * (1 - 2 * v))) * (1 - v);
            MatrixD[0, 1] = (E / ((1 + v) * (1 - 2 * v))) * (v);
            MatrixD[1, 1] = (E / ((1 + v) * (1 - 2 * v))) * (1 - v);
            MatrixD[1, 0] = (E / ((1 + v) * (1 - 2 * v))) * (v);
            MatrixD[2, 2] = (E / ((1 + v) * (1 - 2 * v))) * (1 - v / 2);
            Matrix<double> D = Matrix<double>.Build.DenseOfArray(MatrixD);

            for (int j = 0; j < TrianglesVertices.Count; j++) {
                int[] tri = TrianglesVertices[j];

                UnitVector3D dx0, dy0, dz0;
                UnitVector3D dxn, dyn, dzn;

                Point3D p1, p2, p3;
                Point3D p10, p20, p30;
                Point3D p1n, p2n, p3n;
                Point3D p11, p21, p31;

                p1 = new Point3D(Vertices[tri[0]]);
                p2 = new Point3D(Vertices[tri[1]]);
                p3 = new Point3D(Vertices[tri[2]]);

                p1n = new Point3D(new double[] { X2[tri[0] * 2], X2[tri[0] * 2 + 1], 0 });
                p2n = new Point3D(new double[] { X2[tri[1] * 2], X2[tri[1] * 2 + 1], 0 });
                p3n = new Point3D(new double[] { X2[tri[2] * 2], X2[tri[2] * 2 + 1], 0 });

                dx0 = new UnitVector3D(p2.X - p1.X, p2.Y - p1.Y, p2.Z - p1.Z);
                dz0 = new UnitVector3D(FacetNormals[j]);
                dy0 = dz0.CrossProduct(dx0);

                dxn = new UnitVector3D(p2n.X - p1n.X, p2n.Y - p1n.Y, p2n.Z - p1n.Z);
                dzn = PlaneNormal;
                dyn = dxn.CrossProduct(dzn);

                p10 = new Point3D(new double[] { 0, 0, 0 });
                p20 = new Point3D(new double[] { p2.ToVector().DotProduct(dx0.ToVector()) - p1.ToVector().DotProduct(dx0.ToVector()), 0, 0 });
                p30 = new Point3D(new double[] { p3.ToVector().DotProduct(dx0.ToVector()) - p1.ToVector().DotProduct(dx0.ToVector()), p3.ToVector().DotProduct(dy0.ToVector()) - p1.ToVector().DotProduct(dy0.ToVector()), 0 });

                p11 = new Point3D(new double[] { 0, 0, 0 });
                p21 = new Point3D(new double[] { p2n.ToVector().DotProduct(dxn.ToVector()) - p1n.ToVector().DotProduct(dxn.ToVector()), 0, 0 });
                p31 = new Point3D(new double[] { p3n.ToVector().DotProduct(dxn.ToVector()) - p1n.ToVector().DotProduct(dxn.ToVector()), p3n.ToVector().DotProduct(dyn.ToVector()) - p1n.ToVector().DotProduct(dyn.ToVector()), 0 });

                double[,] Transform = new double[3 * 3, 3 * 3];
                for (int i = 0; i < 3; i++) {
                    Transform[i * 3 + 0, i * 3 + 0] = Math.Cos(dx0.AngleTo(dx).Radians);
                    Transform[i * 3 + 0, i * 3 + 1] = Math.Cos(dx0.AngleTo(dy).Radians);
                    Transform[i * 3 + 0, i * 3 + 2] = Math.Cos(dx0.AngleTo(dz).Radians);

                    Transform[i * 3 + 1, i * 3 + 0] = Math.Cos(dy0.AngleTo(dx).Radians);
                    Transform[i * 3 + 1, i * 3 + 1] = Math.Cos(dy0.AngleTo(dy).Radians);
                    Transform[i * 3 + 1, i * 3 + 2] = Math.Cos(dy0.AngleTo(dz).Radians);

                    Transform[i * 3 + 2, i * 3 + 0] = Math.Cos(dyn.AngleTo(dx).Radians);
                    Transform[i * 3 + 2, i * 3 + 1] = Math.Cos(dyn.AngleTo(dy).Radians);
                    Transform[i * 3 + 2, i * 3 + 2] = Math.Cos(dyn.AngleTo(dz).Radians);
                }

                double A = (p10.X * (p20.Y - p30.Y) + p20.X * (p30.Y - p10.Y) + p30.X * (p10.Y - p20.Y)) / 2;

                double[,] MatrixB = new double[3, 3 * 3];
                MatrixB[0, 0] = (1 / (2 * A)) * (p20.Y - p30.Y);
                MatrixB[1, 1] = (1 / (2 * A)) * (p30.X - p20.X);
                MatrixB[2, 0] = (1 / (2 * A)) * (p30.X - p20.X);
                MatrixB[2, 1] = (1 / (2 * A)) * (p20.Y - p30.Y);

                MatrixB[0, 3] = (1 / (2 * A)) * (p30.Y - p10.Y);
                MatrixB[1, 4] = (1 / (2 * A)) * (p10.X - p30.X);
                MatrixB[2, 3] = (1 / (2 * A)) * (p10.X - p30.X);
                MatrixB[2, 4] = (1 / (2 * A)) * (p30.Y - p10.Y);

                MatrixB[0, 6] = (1 / (2 * A)) * (p10.Y - p20.Y);
                MatrixB[1, 7] = (1 / (2 * A)) * (p20.X - p10.X);
                MatrixB[2, 6] = (1 / (2 * A)) * (p20.X - p10.X);
                MatrixB[2, 7] = (1 / (2 * A)) * (p10.Y - p20.Y);

                Matrix<double> T = Matrix<double>.Build.DenseOfArray(Transform);
                Console.WriteLine(T);
                Matrix<double> B = Matrix<double>.Build.DenseOfArray(MatrixB);
                Console.WriteLine(B);
                Matrix<double> Ke = B.Transpose() * D * B * h * A;
                Console.WriteLine(Ke);
                Ke = T.Transpose() * Ke * T;
                Console.WriteLine(Ke);

                for (int i1 = 0; i1 < 3; i1++) {
                    for (int i2 = 0; i2 < 3; i2++) {
                        for (int l = 0; l < 3; l++) {
                            for (int c = 0; c < 3; c++) {
                                MatrixK[tri[i1] * 3+l, tri[i2] * 3+c] += Ke[i1*3+l, i2*3+c];
                            }
                        }
                    }

                }
                Matrix<double> Temp = Matrix<double>.Build.DenseOfArray(MatrixK);
                Console.WriteLine(Temp);
                double[] vectorqt = new double[3 * 3];
                vectorqt[0 * 3] += p11.X - p10.X;
                vectorqt[0 * 3 + 1] += p11.Y - p10.Y;

                vectorqt[1 * 3] += p21.X - p20.X;
                vectorqt[1 * 3 + 1] += p21.Y - p20.Y;

                vectorqt[2 * 3] += p31.X - p30.X;
                vectorqt[2 * 3 + 1] += p31.Y - p30.Y;

                Vector<double> qt = Vector<double>.Build.DenseOfArray(vectorqt);
                Console.WriteLine(qt);
                qt = T.Transpose() * qt * T;
                Console.WriteLine(qt);
                for (int i1 = 0; i1 < 3; i1++) {
                    for (int l = 0; l < 3; l++) {
                        Vectorq[tri[i1] * 3 + l] += qt[i1 * 3 + l];
                    }
                }
                Vector<double> VTemp = Vector<double>.Build.DenseOfArray(Vectorq);
                Console.WriteLine(VTemp);
            }
            Matrix<double> K = Matrix<double>.Build.DenseOfArray(MatrixK);
            Console.WriteLine(K);
            Vector<double> q = Vector<double>.Build.DenseOfArray(Vectorq);
            Console.WriteLine(q);
            Vector<double> Fi = K * q;
            Console.WriteLine(Fi);

            double w = 0.5;
            Vector<double> qi = Vector<double>.Build.DenseOfArray(q.ToArray());
            Vector<double> dqi = Vector<double>.Build.Dense(Vertices.Count);

            for (int i = 0; i < 100; i++) {
               // Console.WriteLine(dqi);
                dqi = K.Solve(Fi);
                Console.WriteLine(dqi);
               // Console.WriteLine(qi);
                qi = qi + w * dqi;
               // Console.WriteLine(qi);
              //  Console.WriteLine(Fi);
                Fi = K * qi;
               // Console.WriteLine(Fi);
            }

            Application.Exit();
        }
    }
}
