using MathNet.Numerics.LinearAlgebra;
using MathNet.Spatial.Euclidean;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BlankCalculator {
    public class EdgeToEdgeSolver2 {

        internal Vector<double> Solve(List<double[]> Vertices, List<int[]>TrianglesEdges, List<int[]> Edges, List<int> IndiceOfFixedPoints, Point3D oRoot, UnitVector3D vDir1, UnitVector3D vDir2) {
            double[,] MatrixA = new double[Edges.Count * 2, Vertices.Count * 2];

            for (int i = 0; i < Edges.Count; i++) {
                int[] CurTri = new int[] { 0, 0, 0 };
                int indexDown0 = -1;
                foreach (int[] t in TrianglesEdges) {
                    for (int j = 0; j < 3; j++) {
                        if (t[j] == i) { indexDown0 = j; break; }
                    }
                    if (indexDown0 != -1) { CurTri = t; break; }
                }

                int v0 = 0, v1 = 0, v2 = 0;
                int indexDown1 = indexDown0 - 1;
                if (indexDown1 < 0) indexDown1 = 2;
                if (Edges[CurTri[indexDown0]][0] == Edges[CurTri[indexDown1]][0]) {
                    v0 = Edges[CurTri[indexDown0]][1];
                    v1 = Edges[CurTri[indexDown0]][0];
                    v2 = Edges[CurTri[indexDown1]][1];
                }else if(Edges[CurTri[indexDown0]][0] == Edges[CurTri[indexDown1]][1]) {
                    v0 = Edges[CurTri[indexDown0]][1];
                    v1 = Edges[CurTri[indexDown0]][0];
                    v2 = Edges[CurTri[indexDown1]][0];
                } else if (Edges[CurTri[indexDown0]][1] == Edges[CurTri[indexDown1]][0]) {
                    v0 = Edges[CurTri[indexDown0]][0];
                    v1 = Edges[CurTri[indexDown0]][1];
                    v2 = Edges[CurTri[indexDown1]][1];
                } else if (Edges[CurTri[indexDown0]][1] == Edges[CurTri[indexDown1]][1]) {
                    v0 = Edges[CurTri[indexDown0]][0];
                    v1 = Edges[CurTri[indexDown0]][1];
                    v2 = Edges[CurTri[indexDown1]][0];
                }

                Point3D vd0, vd1, vd2;
                vd0 = new Point3D(Vertices[v0]);
                vd1 = new Point3D(Vertices[v1]);
                vd2 = new Point3D(Vertices[v2]);

                Line3D lU, lD;
                lU = new Line3D(vd1, vd0);
                lD = new Line3D(vd1, vd2);

                double Ang = lU.Direction.AngleTo(lD.Direction).Radians;
                double Len =  lU.Length/ lD.Length;


                MatrixA[CurTri[indexDown0] * 2, v0 * 2] = 1;
                MatrixA[CurTri[indexDown0] * 2, v0 * 2 + 1] = 0;

                MatrixA[CurTri[indexDown0] * 2, v2 * 2] = -Len * Math.Cos(Ang);
                MatrixA[CurTri[indexDown0] * 2, v2 * 2 + 1] = Len * Math.Sin(Ang);

                MatrixA[CurTri[indexDown0] * 2, v1 * 2] = Len * Math.Cos(Ang) - 1;
                MatrixA[CurTri[indexDown0] * 2, v1 * 2 + 1] = -Len * Math.Sin(Ang);


                MatrixA[CurTri[indexDown0] * 2 + 1, v0 * 2] = 0;
                MatrixA[CurTri[indexDown0] * 2 + 1, v0 * 2 + 1] = 1;

                MatrixA[CurTri[indexDown0] * 2 + 1, v2 * 2] = -Len * Math.Sin(Ang);
                MatrixA[CurTri[indexDown0] * 2 + 1, v2 * 2 + 1] = -Len * Math.Cos(Ang);

                MatrixA[CurTri[indexDown0] * 2 + 1, v1 * 2] = Len * Math.Sin(Ang);
                MatrixA[CurTri[indexDown0] * 2 + 1, v1 * 2 + 1] = Len * Math.Cos(Ang) - 1;
            }

            double[,] MatrixCa = new double[IndiceOfFixedPoints.Count * 2, Vertices.Count * 2];
            double[] VectorR = new double[IndiceOfFixedPoints.Count * 2];

            Plane oPlane = new Plane(vDir1.CrossProduct(vDir2), oRoot);

            for (int i = 0; i < IndiceOfFixedPoints.Count; i++) {
                MatrixCa[i * 2, IndiceOfFixedPoints[i] * 2] = 1;
                VectorR[i * 2] = new Point3D(Vertices[IndiceOfFixedPoints[i]]).ProjectOn(oPlane).X;
                MatrixCa[i * 2 + 1, IndiceOfFixedPoints[i] * 2 + 1] = 1;
                VectorR[i * 2 + 1] = new Point3D(Vertices[IndiceOfFixedPoints[i]]).ProjectOn(oPlane).Y;
            }

            Matrix<double> Ca = Matrix<double>.Build.DenseOfArray(MatrixCa);
            Console.WriteLine(Ca);
            Vector<double> R = Vector<double>.Build.DenseOfArray(VectorR);
            Console.WriteLine(R);
            Matrix<double> A = Matrix<double>.Build.DenseOfArray(MatrixA);
            Console.WriteLine(A);
            double Penalty = 1000;

            Matrix<double> Ak = A.Transpose() * A + Penalty * Ca.Transpose() * Ca;
            Console.WriteLine(Ak);
            Vector<double> X = Ak.Solve(Penalty * Ca.Transpose() * R);
            Console.WriteLine(X);
            return X;
        }
    }
}
