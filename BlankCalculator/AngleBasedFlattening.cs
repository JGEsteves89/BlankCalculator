using MathNet.Numerics.LinearAlgebra;
using MathNet.Spatial.Euclidean;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BlankCalculator {
    public static class AngleBasedFlattening {
        internal static Vector<double> Solve(Mesh M) {
            return Solve(M.Vertices, M.TrianglesVertices, M.Edges, M.IndiceOfFixedPoints, M.oRoot, M.vDir1, M.vDir2);
        }
        public static Vector<double> Solve(List<double[]> Vertices, List<int[]> Triangles, List<int[]> Edges, List<int> IndiceOfFixedPoints, Point3D oRoot, UnitVector3D vDir1, UnitVector3D vDir2) {

            Vector<double> Solution = Vector<double>.Build.Dense(Triangles.Count * 3);
            

            List<List<int[]>> Wheels = new List<List<int[]>>();
            for (int i = 0; i < Vertices.Count; i++) {
                List<int> TrianglesWithThisVertice = new List<int>();
                for (int t = 0; t < Triangles.Count; t++) {
                    if (Triangles[t][0] == i || Triangles[t][1] == i || Triangles[t][2] == i) TrianglesWithThisVertice.Add(t);
                }
                int[] SheredEdge = new int[] { 0, 0 };
                int Last = 0;
                int Firt = 0;
                if (Triangles[TrianglesWithThisVertice[0]][0] == i) { SheredEdge = new int[] {i, Triangles[TrianglesWithThisVertice[0]][1] }; }
                else if( Triangles[TrianglesWithThisVertice[0]][1] == i) { SheredEdge = new int[] {i, Triangles[TrianglesWithThisVertice[0]][2] }; } 
                else if (Triangles[TrianglesWithThisVertice[0]][2] == i) { SheredEdge = new int[] { i, Triangles[TrianglesWithThisVertice[0]][0] }; }
                List<int[]> WheelEdges = new List<int[]>();
                WheelEdges.Add(SheredEdge);
                do {
                    int found=-1 ;
                    for (int t = 0; t < TrianglesWithThisVertice.Count; t++) {
                        if (t != Last) {
                            if (SheredEdge.Contains(Triangles[TrianglesWithThisVertice[t]][0]) && SheredEdge.Contains(Triangles[TrianglesWithThisVertice[t]][1])) {
                                SheredEdge = new int[] { i, Triangles[TrianglesWithThisVertice[t]][2] }; found = t;
                            } else if (SheredEdge.Contains(Triangles[TrianglesWithThisVertice[t]][1]) && SheredEdge.Contains(Triangles[TrianglesWithThisVertice[t]][2])) {
                                SheredEdge = new int[] { i, Triangles[TrianglesWithThisVertice[t]][0] }; found = t;
                            } else if (SheredEdge.Contains(Triangles[TrianglesWithThisVertice[t]][0]) && SheredEdge.Contains(Triangles[TrianglesWithThisVertice[t]][2])) {
                                SheredEdge = new int[] { i, Triangles[TrianglesWithThisVertice[t]][1] }; found = t;
                            }
                        }
                    }
                    if (found != -1) {
                        Last = found;
                        if (Last == Firt) goto ClosedWeel;
                        WheelEdges.Add(SheredEdge);
                    } else goto Continue;
                } while (true);
                ClosedWeel:
                Wheels.Add(WheelEdges);
                Continue:
                continue;
            }

            Matrix <double> StiffMatrix = Matrix<double>.Build.Dense(Vertices.Count+ Triangles.Count+ Wheels.Count, Triangles.Count * 3);

            //Vertex consistency ∑ei = 2π −∑αi
            Vector<double> VertexConsistency = Vector<double>.Build.Dense(Vertices.Count + Triangles.Count + Wheels.Count);
            for (int i = 0; i < Vertices.Count; i++) {
                Point3D CurrentVertex = new Point3D(Vertices[i]);
                List<double> Angles = new List<double>();
                foreach (int[] t in Triangles) {
                    int i1 = -1, i2 = 0, i0 = 0;
                    if (t[0] == i) { i1 = 0; } else if (t[1] == i) { i1 = 1; } else if (t[2] == i) { i1 = 2; }
                    if (i1 != -1) {
                        i2 = i1 + 1;
                        if (i2 == 3) i2 = 0;
                        i0 = i1 - 1;
                        if (i0 == -1) i0 = 2;
                        Point3D vd0, vd1, vd2;
                        vd0 = new Point3D(Vertices[t[i0]]);
                        vd1 = new Point3D(Vertices[t[i1]]);
                        vd2 = new Point3D(Vertices[t[i2]]);

                        Line3D lU, lD;
                        lU = new Line3D(vd1, vd0);
                        lD = new Line3D(vd1, vd2);
                        Angles.Add(lD.Direction.AngleTo(lU.Direction).Radians);
                    }
                }
                VertexConsistency[i] = 2 * Math.PI - Angles.Sum();
            }

            //Triangle Consistency eα +eβ +eγ = π − (α + β + γ) 
            Vector<double> TriangleConsistency = Vector<double>.Build.Dense(Vertices.Count + Triangles.Count + Wheels.Count);
            for (int t = 0; t < Triangles.Count; t++) {
                double[] Angles = new double[] { 0, 0, 0 };
                for (int i = 0; i < 3; i++) {
                    int i1 = i, i2 = 0, i0 = 0;
                    i2 = i1 + 1;
                    if (i2 == 3) i2 = 0;
                    i0 = i1 - 1;
                    if (i0 == -1) i0 = 2;
                    Point3D vd0, vd1, vd2;
                    vd0 = new Point3D(Vertices[Triangles[t][i0]]);
                    vd1 = new Point3D(Vertices[Triangles[t][i1]]);
                    vd2 = new Point3D(Vertices[Triangles[t][i2]]);

                    Line3D lU, lD;
                    lU = new Line3D(vd1, vd0);
                    lD = new Line3D(vd1, vd2);
                    Angles[i] = lD.Direction.AngleTo(lU.Direction).Radians;
                    StiffMatrix[Vertices.Count + t, t * 3 + i] = 1;
                }
                TriangleConsistency[Vertices.Count + t] = Math.PI - Angles.Sum();
            }

            //Wheel Consistency ∑cot(βi) eβi − cot(γi) eγi =∑log(sinγi) − log(sinβi)
            Vector<double> WheelConsistency = Vector<double>.Build.Dense(Vertices.Count + Triangles.Count + Wheels.Count);
            for (int t = 0; t < Wheels.Count; t++) {
                double cotβi = 0;
                double cotγi = 0;
                double logSinβi = 0;
                double logSinγi = 0;
                double[] Angles = new double[] { 0, 0, 0 };
                for (int i = 0; i < Wheels[t].Count; i++) {
                    Point3D vd0, vd1, vd2;
                    vd0 = new Point3D(Vertices[Wheels[t][i][0]]);
                    vd1 = new Point3D(Vertices[Wheels[t][i][1]]);
                    int other = i + 1;
                    if (other == Wheels[t].Count) other = 0;
                    vd2 = new Point3D(Vertices[Wheels[t][other][1]]);

                    Line3D lV1, lV2, lOp;
                    lV1 = new Line3D(vd0, vd1);
                    lV2 = new Line3D(vd0, vd2);
                    lOp = new Line3D(vd1, vd2);
                    cotβi += Acot(Math.Abs(lV2.Direction.AngleTo(lOp.Direction).Radians));
                    cotγi += Acot(Math.Abs(lV1.Direction.AngleTo(lOp.Direction).Radians));

                    logSinβi += Math.Log(Math.Sin(Math.Abs(lV2.Direction.AngleTo(lOp.Direction).Radians)));
                    logSinγi += Math.Log(Math.Sin(Math.Abs(lV1.Direction.AngleTo(lOp.Direction).Radians)));

                    StiffMatrix[Vertices.Count + t, t * 3 + i] = 1;
                }
                TriangleConsistency[Vertices.Count + t] = Math.PI - Angles.Sum();
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

                double Ang = -lU.Direction.AngleTo(lD.Direction).Radians;
                double Len = lD.Length / lU.Length;


                MatrixA[CurTri[indexDown0] * 2, CurTri[indexDown0] * 2] = 1;
                MatrixA[CurTri[indexDown0] * 2, CurTri[indexDown0] * 2 + 1] = 0;

                MatrixA[CurTri[indexDown0] * 2, CurTri[indexDown1] * 2] = -Len * Math.Cos(Ang);
                MatrixA[CurTri[indexDown0] * 2, CurTri[indexDown1] * 2 + 1] = Len * Math.Sin(Ang);

                MatrixA[CurTri[indexDown0] * 2, CurTri[indexDown2] * 2] = Len * Math.Cos(Ang) - 1;
                MatrixA[CurTri[indexDown0] * 2, CurTri[indexDown2] * 2 + 1] = -Len * Math.Sin(Ang);


                MatrixA[CurTri[indexDown0] * 2 + 1, CurTri[indexDown0] * 2] = 0;
                MatrixA[CurTri[indexDown0] * 2 + 1, CurTri[indexDown0] * 2 + 1] = 1;

                MatrixA[CurTri[indexDown0] * 2 + 1, CurTri[indexDown1] * 2] = -Len * Math.Sin(Ang);
                MatrixA[CurTri[indexDown0] * 2 + 1, CurTri[indexDown1] * 2 + 1] = -Len * Math.Cos(Ang);

                MatrixA[CurTri[indexDown0] * 2 + 1, CurTri[indexDown2] * 2] = Len * Math.Sin(Ang);
                MatrixA[CurTri[indexDown0] * 2 + 1, CurTri[indexDown2] * 2 + 1] = Len * Math.Cos(Ang) - 1;
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
        public static double Acot(double x) {
            return x == 0 ? 0 : Math.Atan(1 / x);
        }
    }

}
