using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MECMOD;
using INFITF;
using HybridShapeTypeLib;
using PARTITF;
using SPATypeLib;
using KnowledgewareTypeLib;
using System.Windows.Forms;
using MathNet.Numerics.LinearAlgebra;

namespace BlankCalculator {
    public class CATInterop {
        public INFITF.Application CATIA;
        public PartDocument oPartDoc;
        public Selection oSel;
        public SPAWorkbench oSpa;
        public CATInterop() {
            try {
                CATIA = (INFITF.Application)Microsoft.VisualBasic.Interaction.GetObject(null, "CATIA.Application");
            } catch (Exception) {
                if (System.Windows.Forms.MessageBox.Show("Não foi detectado nenhum CATIA aberto. Pretende Abrir?", "", System.Windows.Forms.MessageBoxButtons.YesNo) == System.Windows.Forms.DialogResult.Yes) {
                    CATIA = (INFITF.Application)Microsoft.VisualBasic.Interaction.CreateObject("CATIA.Application");
                } else {
                    Environment.Exit(0);
                }
            }
        }

        internal string SaveSurfaceAndFixedPoints(ref List<double[]> FixedPoints, ref double[] oPlane) {
            string StlPath= System.IO.Path.GetTempPath() + "STLFileBlankCalculator"; ;
            oPartDoc = (PartDocument)CATIA.ActiveDocument;
            oSel = oPartDoc.Selection;
            oSpa = (SPAWorkbench)oPartDoc.GetWorkbench("SPAWorkbench");

            FixedPoints = new List<double[]>();

            oSel.Clear();
            oSel.Add(SelectSurface("Selecione qual a superficie que pretende planificar. Esc para sair."));
            oSel.Copy();
            oSel.Clear();
            PartDocument NewPart = (PartDocument)CATIA.Documents.Add("Part");
            NewPart.Selection.Clear();
            NewPart.Selection.Add(NewPart.Part);
            NewPart.Selection.PasteSpecial("CATPrtResultWithOutLink");
            if (System.IO.File.Exists(StlPath + ".stl")) {
                System.IO.File.Delete(StlPath + ".stl");
            }
            NewPart.Selection.Clear();
            NewPart.Part.Update();
            CATIA.DisplayFileAlerts = false;
            NewPart.ExportData(StlPath, "stl");
            NewPart.Close();
            CATIA.DisplayFileAlerts = true;

            object[] Vec = new object[3];
            oSel.Clear();
            Reference Ref1 = SelectPoint("Selecione o conjunto de pontos fixos. Esc para sair.");
            oSel.Clear();
            if (Ref1 == null) Environment.Exit(0);
            oSpa.GetMeasurable(Ref1).GetPoint(Vec);
            FixedPoints.Add(new double[] { (double)Vec[0], (double)Vec[1], (double)Vec[2] });
            do {
                Ref1 = SelectPoint("Selecione o conjunto de pontos fixos (" + FixedPoints.Count + " selecionados). Esc para terminar.");
                oSel.Clear();
                if (Ref1 == null) break;
                oSpa.GetMeasurable(Ref1).GetPoint(Vec);
                FixedPoints.Add(new double[] { (double)Vec[0], (double)Vec[1], (double)Vec[2] });
            } while (true);
            if(FixedPoints.Count==0) Environment.Exit(0);
            oSel.Clear();
            oPartDoc.Part.Update();
            System.Windows.Forms.Application.DoEvents();
            System.Threading.Thread.Sleep(500);
            Vec = new object[9];
            Reference Ref2 = SelectPlane("Selecione qual o plano do planificado. Esc para terminar.");
            oSel.Clear();
            oSpa.GetMeasurable(Ref2).GetPlane(Vec);
            oPlane = new double[] { (double)Vec[0], (double)Vec[1], (double)Vec[2], 
                                    (double)Vec[3], (double)Vec[4], (double)Vec[5], 
                                    (double)Vec[6], (double)Vec[7], (double)Vec[8] };
            return StlPath + ".stl";
        }
        private Reference SelectSurface(string msg) {
            object[] InputobjectType = new object[] { "Face" };
            Reference Ref= null;
            oSel.Clear();
            AGAIN:
            object status = oSel.SelectElement2(InputobjectType, msg, true);
            if (status.ToString() == "Cancel" || status.ToString() == "Undo") {
                Environment.Exit(0);
            }

            if (oSel.Item2(1).Reference.DisplayName.ToLower().IndexOf("Selection".ToLower()) !=-1) {
                try {
                    Ref = oSel.Item2(1).Reference;
                } catch (Exception) {
                    try {
                        Ref = (Reference)((AnyObject)oSel.Item2(1).Value).Parent;
                    } catch (Exception) {
                        MessageBox.Show("Certifique-se que a selecção está correcta.");
                        goto AGAIN;
                    }
                }
            }
            return Ref;
        }
        private Reference SelectPlane(string msg) {
            object[] InputobjectType = new object[] { "Plane" };
            Reference Ref = null;
            oSel.Clear();
            AGAIN:
            object status = oSel.SelectElement2(InputobjectType, msg, true);
            if (status.ToString() == "Cancel" || status.ToString() == "Undo") {
                Environment.Exit(0);
            }
            if (oSel.Item2(1).Reference.DisplayName.ToLower().IndexOf("Plane".ToLower()) != -1) {
                try {
                    Ref = oSel.Item2(1).Reference;
                } catch (Exception) {
                    try {
                        Ref = (Reference)((AnyObject)oSel.Item2(1).Value).Parent;
                    } catch (Exception) {
                        MessageBox.Show("Certifique-se que a selecção está correcta.");
                        goto AGAIN;
                    }
                }
            }
            return Ref;
        }
        private Reference SelectPoint(string msg) {
            object[] InputobjectType = new object[] { "Point" };
            Reference Ref = null;
            oSel.Clear();
            AGAIN:
            object status = oSel.SelectElement2(InputobjectType, msg, true);
            if (status.ToString() == "Cancel" || status.ToString() == "Undo") {
               return null;
            }
            if (oSel.Item2(1).Reference.DisplayName.ToLower().IndexOf("Point".ToLower()) != -1) {
                try {
                    Ref = oSel.Item2(1).Reference;
                    oSel.VisProperties.SetRealWidth(4, 1);
                    oSel.VisProperties.SetRealColor(System.Drawing.Color.Blue.R, System.Drawing.Color.Blue.G, System.Drawing.Color.Blue.B,1);
                } catch (Exception) {
                    MessageBox.Show("Certifique-se que a selecção está correcta.");
                    goto AGAIN;
                }
            }
            return Ref;
        }

        internal void PrintTriangles(Vector<double> x,
            List<int[]> triangles,
            MathNet.Spatial.Euclidean.Point3D oRoot,
            MathNet.Spatial.Euclidean.UnitVector3D vDir1,
            MathNet.Spatial.Euclidean.UnitVector3D vDir2,
            int iTrans,
            MathNet.Spatial.Euclidean.Point3D oTrans,
            bool just2D= false) {

            Part oPart = oPartDoc.Part;
            HybridShapeFactory hsf = (HybridShapeFactory)oPart.HybridShapeFactory;
            HybridBody hb1;
            try {
                hb1=(HybridBody)oPart.HybridBodies.GetItem("Blank Calculator Result");
                hb1.get_Name();
                oSel.Clear();
                oSel.Add(hb1);
                oSel.Delete();
            } catch (Exception) {}

            hb1 = (HybridBody)oPart.HybridBodies.Add();
            hb1.set_Name("Blank Calculator Result");

            //CATIA.RefreshDisplay = false;
            //CATIA.Interactive = false;
            List<Reference> RsltPoints = new List<Reference>();
            MathNet.Spatial.Euclidean.CoordinateSystem Axis = new MathNet.Spatial.Euclidean.CoordinateSystem(oRoot, vDir1, vDir2, vDir1.CrossProduct(vDir2));
            if (just2D) {
                for (int i = 0; i < x.Count / 2; i++) {
                    MathNet.Spatial.Euclidean.Point3D  PtMath = new MathNet.Spatial.Euclidean.Point3D(new double[] { x[i * 2], x[i * 2 + 1], 0 });
                    Point PTCat = hsf.AddNewPointCoord(PtMath.X, PtMath.Y, PtMath.Z);
                    PTCat.Compute();
                    RsltPoints.Add(oPart.CreateReferenceFromObject(PTCat));
                }
            } else {
                MathNet.Spatial.Euclidean.Point3D PtMath = new MathNet.Spatial.Euclidean.Point3D(new double[] { x[iTrans * 2], x[iTrans * 2 + 1], 0 });
                PtMath = Axis.TransformToCoordSys(PtMath);
                double[] vecTranslation = new double[] { oTrans.X - PtMath.X, oTrans.Y - PtMath.Y, oTrans.Z - PtMath.Z };
                for (int i = 0; i < x.Count / 2; i++) {
                    PtMath = new MathNet.Spatial.Euclidean.Point3D(new double[] { x[i * 2], x[i * 2 + 1], 0 });
                    PtMath = Axis.TransformToCoordSys(PtMath);
                    Point PTCat = hsf.AddNewPointCoord(PtMath.X + vecTranslation[0], PtMath.Y + vecTranslation[1], PtMath.Z + vecTranslation[2]);
                    PTCat.Compute();
                    RsltPoints.Add(oPart.CreateReferenceFromObject(PTCat));
                }
            }
              

            foreach (int[] item in triangles) {
                Line lUp = hsf.AddNewLinePtPt(RsltPoints[item[0]], RsltPoints[item[1]]);
                lUp.Compute();
                hb1.AppendHybridShape(lUp);
                Line l0 = hsf.AddNewLinePtPt(RsltPoints[item[1]], RsltPoints[item[2]]);
                l0.Compute();
                hb1.AppendHybridShape(l0);
                Line lDown = hsf.AddNewLinePtPt(RsltPoints[item[2]], RsltPoints[item[0]]);
                lDown.Compute();
                hb1.AppendHybridShape(lDown);
            }
            
            //CATIA.RefreshDisplay = true;
            //CATIA.Interactive = true;
        }
    }
}
