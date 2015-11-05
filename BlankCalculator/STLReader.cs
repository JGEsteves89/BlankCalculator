using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading.Tasks;


namespace BlankCalculator {
    public static class STLReader {
    
        public static void STLRead(string CaminhoStl, ref List<double[]> Vertices, ref List<int[]> TrianglesVertices, ref List<int[]> TrianglesEdges, ref List<int[]> Edges, ref List<double[]> FacetsNormal) {
            Vertices = new List<double[]>();
            FacetsNormal = new List<double[]>();
            TrianglesVertices = new List<int[]>();
            Edges = new List<int[]>();
            int[] LastTriangle = new int[] { 0, 0, 0 };
            int[] LastEdges= new int[] { 0, 0, 0 };
            int i = 0;
            List<string> Lines = System.IO.File.ReadAllLines(CaminhoStl).ToList();
            foreach (string line in Lines) {
                if (line.Trim().StartsWith("vertex")) {
                    double[] aux = new double[] { double.Parse(line.Split(new string[] { " " }, StringSplitOptions.RemoveEmptyEntries)[1], CultureInfo.InvariantCulture),
                                                double.Parse(line.Split(new string[]{" "}, StringSplitOptions.RemoveEmptyEntries)[2], CultureInfo.InvariantCulture),
                                                double.Parse(line.Split(new string[]{" "}, StringSplitOptions.RemoveEmptyEntries)[3], CultureInfo.InvariantCulture) };
                    LastTriangle[i] = FindVertice(Vertices, aux);
                    if (LastTriangle[i] == -1) {
                        Vertices.Add(aux);
                        LastTriangle[i] = Vertices.Count - 1;
                    }
                    i += 1;
                } else if (line.Trim().StartsWith("endloop")) {
                    i = 0;
                    TrianglesVertices.Add(new int[] { LastTriangle[0], LastTriangle[1], LastTriangle[2] });

                    LastEdges[0] = FindEdge(Edges, new int[] { LastTriangle[0], LastTriangle[1] });
                    if (LastEdges[0] == -1) {
                        Edges.Add(new int[] { LastTriangle[0], LastTriangle[1] });
                        LastEdges[0] = Edges.Count - 1;
                    }

                    LastEdges[1] = FindEdge(Edges, new int[] { LastTriangle[1], LastTriangle[2] });
                    if (LastEdges[1] == -1) {
                        Edges.Add(new int[] { LastTriangle[1], LastTriangle[2] });
                        LastEdges[1] = Edges.Count - 1;
                    }

                    LastEdges[2] = FindEdge(Edges, new int[] { LastTriangle[2], LastTriangle[0] });
                    if (LastEdges[2] == -1) {
                        Edges.Add(new int[] { LastTriangle[2], LastTriangle[0] });
                        LastEdges[2] = Edges.Count - 1;
                    }
                    TrianglesEdges.Add(new int[] { LastEdges[0], LastEdges[1], LastEdges[2] });
                } else if (line.Trim().StartsWith("facet normal")) {
                    double[] aux = new double[] { double.Parse(line.Split(new string[] { " " }, StringSplitOptions.RemoveEmptyEntries)[2], CultureInfo.InvariantCulture),
                                                double.Parse(line.Split(new string[]{" "}, StringSplitOptions.RemoveEmptyEntries)[3], CultureInfo.InvariantCulture),
                                                double.Parse(line.Split(new string[]{" "}, StringSplitOptions.RemoveEmptyEntries)[4], CultureInfo.InvariantCulture) };
                    FacetsNormal.Add(aux);
                }
            }
        }
        public static int FindVertice( List<double[]> Vertices, double[] aux) {
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
        public static int FindEdge(List<int[]> Edges, int[] aux) {
            for (int i = 0; i < Edges.Count; i++) {
                if (Edges[i][0].Equals(aux[0])) {
                    if (Edges[i][1].Equals(aux[1])) {
                              return i;
                    }
                } else if (Edges[i][0].Equals(aux[1])) {
                    if (Edges[i][1].Equals(aux[0])) {
                        return i;
                    }
                }
            }
            return -1;
        }
    }

}
