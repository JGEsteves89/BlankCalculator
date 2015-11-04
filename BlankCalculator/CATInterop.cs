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

namespace BlankCalculator {
    public class CATInterop {
        public INFITF.Application CATIA;
        public PartDocument oPartDoc;
        public Selection oSel;

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

        internal string SaveSurfaceAsStl() {
            string StlPath= System.IO.Path.GetTempPath() + "STLFileBlankCalculator"; ;
            oPartDoc = (PartDocument)CATIA.ActiveDocument;
            oSel = oPartDoc.Selection;
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

            if (oSel.Item2(1).Reference.DisplayName.IndexOf("Selection")!=-1) {
                try {
                    Ref = oSel.Item2(1).Reference;
                } catch (Exception) {
                    try {
                        Ref = (Reference)((AnyObject)oSel.Item2(1).Value).Parent;
                    } catch (Exception) {
                        MessageBox.Show("Certifique.se que o product se encontra activo.");
                        goto AGAIN;
                    }
                }
            }


            return Ref;

        }

    }
}
