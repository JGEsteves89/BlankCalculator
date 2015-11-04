using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;


namespace BlankCalculator {
    public static class STLReader {
    
        public static void STLRead(string CaminhoStl, ref List<double[]> Vertices, ref List<int[]> Triagles) {
            Vertices = new List<double[]>();
            Triagles = new List<int[]>();
            int[] LastTriangle = new int[] { 0, 0, 0 };
            int i = 0;
            List<string> Lines = System.IO.File.ReadAllLines(CaminhoStl).ToList();
            foreach (string line in Lines) {
                if (line.Trim().StartsWith("vertex")) {
                    string newline = line.Replace(".", ",");
                    double[] aux = new double[] { double.Parse(newline.Split(new string[] { " " }, StringSplitOptions.RemoveEmptyEntries)[1]),
                                                double.Parse(newline.Split(new string[]{" "}, StringSplitOptions.RemoveEmptyEntries)[2]),
                                                double.Parse(newline.Split(new string[]{" "}, StringSplitOptions.RemoveEmptyEntries)[3]) };
                    LastTriangle[i] = FindIndex(Vertices, aux);
                    if (LastTriangle[i] == -1) {
                        Vertices.Add(aux);
                        LastTriangle[i] = Vertices.Count - 1;
                    }
                    i += 1;
                } else if (line.Trim().StartsWith("endloop")) {
                    i = 0;
                    Triagles.Add(new int[] { LastTriangle[0], LastTriangle[1], LastTriangle[2] });
                }
            }
        }
        public static int FindIndex( List<double[]> Vertices, double[] aux) {
            for (int i = 0; i < Vertices.Count; i++) {
                if (Vertices[i][0].Equals(aux[0])) {
                    if (Vertices[i][1].Equals(aux[1])) {
                        if (Vertices[i][2].Equals(aux[2])) {
                            return i;
                        }
                    }
                }
            }
            return -1;
        }
    }

}
