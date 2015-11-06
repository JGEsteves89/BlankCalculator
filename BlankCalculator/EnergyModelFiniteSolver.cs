using MathNet.Numerics.LinearAlgebra;
using MathNet.Spatial.Euclidean;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BlankCalculator {
    public static class EnergyModelFiniteSolver {
        internal static Vector<double> solve(Vector<double> X, Mesh M) {
            return solve(M.Vertices, M.TrianglesVertices, M.FacetNormals, X, M.oRoot, M.vDir1, M.vDir2,M);
        }
        public static Vector<double> solve(
            List<double[]> Vertices,
            List<int[]> TrianglesVertices,
            List<double[]> FacetNormals,
            Vector<double> X, Point3D oRoot,
            UnitVector3D vDir1,
            UnitVector3D vDir2, Mesh M) {


            //Method from article -Wang,Smith. Surface flattening based on energy model. Computer-Aided Design (2002) 823-833	

            //PseudoCode
            //Input: P - Set of nodes, in the initial position and N is the number of nodes
            //Output: the final positions of P with E(o) minimized

            //  FOR i = 1 TO n
            //        mi = p / 3 Sum(Ak), where Ak is the area

            //  WHILE (Relative Area difference Es > Permissible accuracy or Relative edge length difference Ec > Permissible accuracy)
            //  AND Variation of E(o) > Permissible percentage €
            //  AND the number of iterations < Permissible number N
            //        FOR i = 1 TO n
            //                Compute Tensile force of Node Pi: Fi = Sum(C * (Dist(PiPj) - Dist(QiQj)))nij where(P - 2D - Q - 3D nij - Vector Pi to Pj)
            //                Compute new position of Pi  qi = qi + dtqi. + dt ^ 2 / 2 qi..where qi.= qi.+ dtqi..and qi..= Fi / mi
            //                Compute Penalty force and aplly to Fpi
            //                Compute new position of Pi  qi = qi + dtqi. + dt ^ 2 / 2 qi..where qi.= qi.+ dtqi..and qi..= Fpi / mi
            //        Compute new Es= Sum(TotalAreaNow - TotalAreaBefore) / TotalAreaNow
            //        Compute new Ec= Sum(TotalLenghtNow - TotalLenghtBefore) / TotalLenghtNow
            //        Compute new E(o)Sum(E(pi)) where E(pi) = 0.5 * Sum(C * (Dist(PiPj) - Dist(QiQj))) ^ 2

            double[] Mass = new double [Vertices.Count];
            int[] MassCounter = new int[Vertices.Count];
            double Permissible = 0.00000001;
            Vector<double> Fi = Vector<double>.Build.Dense(2);
            List<Vector<double>> dqi = new List<Vector<double>>();
            List<Vector<double>> qi = new List<Vector<double>>();
            double LastEc = 0, Ec = 0;
            double LastEs = 0, Es = 0;
            double LastEo = 0, Eo = 0;
            double C = 0.5;
            double ro = 1;
            double t = 0.01;
            int N = 100;
            int Iteration = 0;
            double Total = 0;
            double LastTotal = 0;
            foreach (int[] tri in TrianglesVertices) {
                double A = getArea(Vertices[tri[0]], Vertices[tri[1]], Vertices[tri[2]]);
                double P = getPerimeter(Vertices[tri[0]], Vertices[tri[1]], Vertices[tri[2]]);
                Mass[tri[0]] += ((double)1 / (double)3) * A * ro;
                Mass[tri[1]] += ((double)1 / (double)3) * A * ro;
                Mass[tri[2]] += ((double)1 / (double)3) * A * ro;
                MassCounter[tri[0]] += 1;
                MassCounter[tri[1]] += 1;
                MassCounter[tri[2]] += 1;
                LastEc += A;
                LastEs += P;
            }
            for (int i = 0; i < Vertices.Count; i++) {
                Mass[i] = Mass[i] / MassCounter[i];
                qi.Add(Vector<double>.Build.DenseOfArray(new double[] { X[i * 2], X[i * 2 + 1], 0 }));
                dqi.Add(Vector<double>.Build.DenseOfArray(new double[] { 0, 0, 0 }));
            }
            do {
                Iteration += 1;
                LastEo = Eo;
                Eo = 0;
                Total = 0;
                for (int i = 0; i < Vertices.Count; i++) {
                    Point3D Pi = new Point3D(qi[i][0], qi[i][1], qi[i][2]);
                    Point3D Qi = new Point3D(Vertices[i][0], Vertices[i][1], Vertices[i][2]);
                    Fi = Vector<double>.Build.Dense(3);
                    for (int j = i+1; j < Vertices.Count; j++) {
                        if (i == j) continue;
                        Point3D Pj = new Point3D(qi[j][0], qi[j][1], qi[j][2]);
                        Point3D Qj = new Point3D(Vertices[j][0], Vertices[j][1], Vertices[j][2]);
                        UnitVector3D nij = new UnitVector3D((Pi.ToVector()- Pj.ToVector()).ToArray());
                        Fi += C * ((Pi.DistanceTo(Pj) - Qi.DistanceTo(Qj)) * nij).ToVector();
                        Eo += Math.Pow(C * ((Pi.DistanceTo(Pj) - Qi.DistanceTo(Qj))), 2);
                    }
                    Vector<double> ddqi = Fi / Mass[i];
                    Console.WriteLine(Fi);
                    dqi[i] += t * ddqi;
                    qi[i] += dqi[i] * t + 0.5 * Math.Pow(t, 2) * ddqi;

                    Total += Math.Abs(qi[i][0] - X[i * 2]);
                    Total += Math.Abs(qi[i][1] - X[i * 2] + 1);
                }
                Console.WriteLine(Total- LastTotal);
                LastTotal = Total;
                LastEc = Ec;
                Ec = 0;
                LastEs = Es;
                Es = 0;
                foreach (int[] tri in TrianglesVertices) {
                    double A = getArea(Vertices[tri[0]], Vertices[tri[1]], Vertices[tri[2]]);
                    double P = getPerimeter(Vertices[tri[0]], Vertices[tri[1]], Vertices[tri[2]]);
                    Ec += A;
                    Es += P;
                }
                if (((Ec - LastEc) / Ec < Permissible || (Es - LastEs) / Es < Permissible) 
                    && (Eo - LastEo) / Eo < Permissible
                    && Iteration > N) break;
            } while (true);
            Vector<double> XResult = Vector<double>.Build.DenseOfArray(X.ToArray());

            for (int i = 0; i < Vertices.Count; i++) {
                //Console.WriteLine(qi[i][0] - X[i * 2]);
                Total += Math.Abs(qi[i][0] - X[i * 2]);
                XResult[i * 2] = qi[i][0];
                //Console.WriteLine(qi[i][1] - X[i * 2+1]);
                Total += Math.Abs(qi[i][1] - X[i * 2]+1);
                XResult[i * 2+1] = qi[i][1];
            }
            Console.WriteLine(Total);
            return XResult;
        }
        public static double getArea(double[]p1, double[] p2, double[] p3) {
            double[] x = new double[] { p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2] };
            double[] y = new double[] { p3[0] - p1[0], p3[1] - p1[1], p3[2] - p1[2] };
             return  0.5 * Math.Sqrt(Math.Pow(x[1] * y[2] - x[2] * y[1], 2) + Math.Pow(x[2] * y[0] - x[0] * y[2], 2) + Math.Pow(x[0] * y[1] - x[1] * y[0], 2));
        }
        public static double getPerimeter(double[] p1, double[] p2, double[] p3) {
            double[] x = new double[] { p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2] };
            double[] y = new double[] { p3[0] - p1[0], p3[1] - p1[1], p3[2] - p1[2] };
            double[] z = new double[] { p3[0] - p2[0], p3[1] - p1[1], p2[2] - p1[2] };
            double total = 0;
            for (int i = 0; i < 3; i++) {
                total += Math.Pow(x[i], 2);
                total += Math.Pow(y[i], 2);
                total += Math.Pow(z[i], 2);
            }
            return total;
        }
    }
}
