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
        public enum FiniteSolverMethod {
            NewtonRaphson,
            EnergyModel
        }
        public enum PlanarGraphSolver {
            EdgeRatioBased,
            AngleBased
        }
        private void Form1_Load(object sender, EventArgs e) {
            Mesh M = new Mesh();
            FiniteSolverMethod Method = FiniteSolverMethod.EnergyModel;
            PlanarGraphSolver PlanarSolver = PlanarGraphSolver.EdgeRatioBased;
            //Interop with catia to retrive exported stlPath
            string StlPath = CAT.SaveSurfaceAndFixedPoints(ref M.FixedPoints, ref M.oVecPlane);

            //Parse stlFile to Vertice and triangle matrix
            STLReader.STLRead(StlPath, ref M);
            if (M.Vertices.Count == 0 || M.TrianglesVertices.Count == 0) {
                MessageBox.Show("Não foi encontrada nenhuma informação na superficie exportada. Tente outra vez.");
                Environment.Exit(0);
            }

            //Detection of the most apropriated verticePoints
            M.ComputeFixedPoints();

            //Calculation of Plane of flattening
            M.ComputePlaneofFlattening();

            Vector<double> X= Vector<double>.Build.Dense(1);
            if (PlanarSolver == PlanarGraphSolver.EdgeRatioBased) {
                //Solve using the edge to edge interpretation of the article
                X = EdgeToEdgeSolver.Solve(M);
            } else if (PlanarSolver == PlanarGraphSolver.AngleBased) {
                //Solve using Angle Based Flattening Method
                X = AngleBasedFlattening.Solve(M);
            }


            //Print Result of triangles in CATIA;
            CAT.PrintTriangles(X, M);
            Vector<double> Y = Vector<double>.Build.Dense(X.Count);
            if (Method== FiniteSolverMethod.NewtonRaphson) {
                //Try to calculate de displacement vector based on finite-element analysis using NewtonRapson Method
                //ATENTION - IS NOT COMPLETE/WITH ERRORS!
                 Y =NewtonRaphsonFiniteSolver.solve(X, M);
            } else if (Method == FiniteSolverMethod.EnergyModel) {
                //Try to calculate de displacement vector based on energy analysis
                //ATENTION - IS NOT COMPLETE/WITH ERRORS!
                Y = EnergyModelFiniteSolver.solve(X, M);
            }
            CAT.PrintTriangles(Y, M);
            Application.Exit();
        }
    }
}
