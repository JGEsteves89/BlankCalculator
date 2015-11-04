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
            string StlPath = CAT.SaveSurfaceAsStl();

            List<double[]> Vertices= new List<double[]>();
            List<int[]> Triangles= new List<int[]> ();

            STLReader.STLRead(StlPath,ref Vertices,ref Triangles);
            if(Vertices.Count==0 || Triangles.Count == 0) {
                MessageBox.Show("Não foi encontrada nenhuma informação na superficie exportada. Tente outra vez.");
                Environment.Exit(0);
            }

            double[,] Equations = new double[Vertices.Count * 2, Vertices.Count * 2];

            for (int i = 0; i < Vertices.Count; i++) {
                int[] CurTri = new int[] { 0, 0, 0 };
                int indexDown0 = -1;
                foreach (int[] t in Triangles) {
                    for (int j = 0; j < 3; j++) {
                        if (t[j] == i) { indexDown0 = j; break; }
                    }
                    if (indexDown0 != -1) {CurTri = t;break;}
                }
                int indexDown1 = indexDown0 - 1;
                if (indexDown1 < 0) indexDown1 = 2;
                int indexDown2 = indexDown1 - 1;
                if (indexDown2 < 0) indexDown2 = 2;
            
                Point3D vd0, vd1, vd2;
                vd0 = new Point3D(Vertices[CurTri[indexDown0]]);
                vd1 = new Point3D(Vertices[CurTri[indexDown1]]);
                vd2 = new Point3D(Vertices[CurTri[indexDown2]]);

                Line3D  lU, lD;
                lU = new Line3D(vd1, vd0);
                lD = new Line3D(vd1, vd2);

                double Ang = lU.Direction.AngleTo(lD.Direction).Radians;
                double Len = lU.Length / lD.Length;


                Equations[CurTri[indexDown0] * 2, CurTri[indexDown0] * 2] = 1;
                Equations[CurTri[indexDown0] * 2, CurTri[indexDown0] * 2+1] = 0;

                Equations[CurTri[indexDown0] * 2, CurTri[indexDown1] * 2] = Len*Math.Cos(Ang)-1;
                Equations[CurTri[indexDown0] * 2, CurTri[indexDown1] * 2 + 1] = -Len * Math.Sin(Ang);

                Equations[CurTri[indexDown0] * 2, CurTri[indexDown2] * 2] = -Len * Math.Cos(Ang);
                Equations[CurTri[indexDown0] * 2, CurTri[indexDown2] * 2 + 1] = Len * Math.Sin(Ang);


                Equations[CurTri[indexDown0] * 2+1, CurTri[indexDown0] * 2] = 0;
                Equations[CurTri[indexDown0] * 2 + 1, CurTri[indexDown0] * 2 + 1] = 1;

                Equations[CurTri[indexDown0] * 2 + 1, CurTri[indexDown1] * 2] = Len * Math.Sin(Ang);
                Equations[CurTri[indexDown0] * 2 + 1, CurTri[indexDown1] * 2 + 1] = Len * Math.Cos(Ang) - 1;

                Equations[CurTri[indexDown0] * 2 + 1, CurTri[indexDown2] * 2] = -Len * Math.Sin(Ang);
                Equations[CurTri[indexDown0] * 2 + 1, CurTri[indexDown2] * 2 + 1] = -Len * Math.Cos(Ang);
            }
            var A = Matrix<double>.Build.DenseOfArray(Equations);
            Console.WriteLine(A);
            var b = Vector<double>.Build.Dense(Vertices.Count * 2);
            Console.WriteLine(b);
            var x = A.Solve(b);
            Console.WriteLine(x);
        }
    }
}
