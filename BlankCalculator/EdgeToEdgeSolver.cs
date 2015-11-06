using MathNet.Numerics.LinearAlgebra;
using MathNet.Spatial.Euclidean;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BlankCalculator {
    public static class EdgeToEdgeSolver {
        internal static Vector<double> Solve(Mesh M) {
            return Solve(M.Vertices, M.TrianglesEdges, M.Edges, M.IndiceOfFixedPoints, M.oRoot, M.vDir1, M.vDir2);
        }
        public static Vector<double> Solve(List<double[]> Vertices, List<int[]> TrianglesEdges, List<int[]> Edges, List<int> IndiceOfFixedPoints, Point3D oRoot, UnitVector3D vDir1, UnitVector3D vDir2) {

            double[,] MatrixA = new double[Edges.Count*2, Vertices.Count * 2];

            for (int i = 0; i < Edges.Count; i++) {
                int[] CurTri = new int[] { 0, 0, 0 };
                int indexDown0 = -1;
                foreach (int[] t in TrianglesEdges) {
                    for (int j = 0; j < 3; j++) {
                        if (t[j] == i) { indexDown0 = j; break; }
                    }
                    if (indexDown0 != -1) { CurTri = t; break; }
                }
                int indexDown1 = indexDown0 - 1;
                if (indexDown1 < 0) indexDown1 = 2;

                int vi0 =0, vi1 = 0, vi2 = 0;
                if (Edges[CurTri[indexDown0]][0]== Edges[CurTri[indexDown1]][0]) {
                    vi0 = Edges[CurTri[indexDown0]][1];
                    vi1= Edges[CurTri[indexDown0]][0];
                    vi2 = Edges[CurTri[indexDown1]][1];
                } else if (Edges[CurTri[indexDown0]][0] == Edges[CurTri[indexDown1]][1]) {
                    vi0 = Edges[CurTri[indexDown0]][1];
                    vi1 = Edges[CurTri[indexDown0]][0];
                    vi2 = Edges[CurTri[indexDown1]][0];
                } else if (Edges[CurTri[indexDown0]][1] == Edges[CurTri[indexDown1]][0]) {
                    vi0 = Edges[CurTri[indexDown0]][0];
                    vi2 = Edges[CurTri[indexDown0]][1];
                    vi1 = Edges[CurTri[indexDown1]][1];
                } else if (Edges[CurTri[indexDown0]][1] == Edges[CurTri[indexDown1]][0]) {
                    vi0 = Edges[CurTri[indexDown0]][0];
                    vi1 = Edges[CurTri[indexDown0]][1];
                    vi2 = Edges[CurTri[indexDown1]][0];
                }

                Point3D vd0, vd1, vd2;
                vd0 = new Point3D(Vertices[vi0]);
                vd1 = new Point3D(Vertices[vi1]);
                vd2 = new Point3D(Vertices[vi2]);

                Line3D lU, lD;
                lU = new Line3D(vd1, vd0);
                lD = new Line3D(vd1, vd2);

                double Ang = lU.Direction.AngleTo(lD.Direction).Radians;
                double Len = lU.Length/ lD.Length ;


                MatrixA[i * 2, vi0 * 2] = 1;
                MatrixA[i * 2, vi0 * 2 + 1] = 0;

                MatrixA[i * 2, vi2 * 2] = -Len * Math.Cos(Ang);
                MatrixA[i * 2, vi2 * 2 + 1] = Len * Math.Sin(Ang);

                MatrixA[i * 2, vi1 * 2] = Len * Math.Cos(Ang) - 1;
                MatrixA[i * 2, vi1 * 2 + 1] = -Len * Math.Sin(Ang);


                MatrixA[i * 2 + 1, vi0 * 2] = 0;
                MatrixA[i * 2 + 1, vi0 * 2 + 1] = 1;

                MatrixA[i * 2 + 1, vi2 * 2] = -Len * Math.Sin(Ang);
                MatrixA[i * 2 + 1, vi2 * 2 + 1] = -Len * Math.Cos(Ang);

                MatrixA[i * 2 + 1, vi1 * 2] = Len * Math.Sin(Ang);
                MatrixA[i* 2 + 1, vi1 * 2 + 1] = Len * Math.Cos(Ang) - 1;
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
            double Penalty = 100;

            Matrix<double> Ak = A.Transpose() * A + Penalty * Ca.Transpose() * Ca;
            Console.WriteLine(Ak);
            Vector<double> X = Ak.Solve(Penalty * Ca.Transpose() * R);

            X = Ak.Solve(Penalty * Ca.Transpose() * R);
            Console.WriteLine(X);
            bool valid = true;
            if (!valid) {
                var solver = new MathNet.Numerics.LinearAlgebra.Double.Solvers.BiCgStab();
                var preconditioner = new MathNet.Numerics.LinearAlgebra.Double.Solvers.MILU0Preconditioner();
                var iterator = new MathNet.Numerics.LinearAlgebra.Solvers.Iterator<double>(
                    new MathNet.Numerics.LinearAlgebra.Solvers.ResidualStopCriterion<double>(1.0e-8),
                    new MathNet.Numerics.LinearAlgebra.Solvers.IterationCountStopCriterion<double>(1000));

                var x = Ak.SolveIterative(Penalty * Ca.Transpose() * R, solver, iterator);
                var status = iterator.Status;
                Console.WriteLine(X);
            }
            return X;
        }


    }
}
