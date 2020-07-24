void generate_ntuple(){
  TNtuple* ntuple = new TNtuple("ntuple", "data from bob azmoun", "position:intensity");
  ntuple->ReadFile("Values_Cleaved_Sanded.csv");//Values_Cleaved_Sanded,Values_Cleaved_Not_Sanded
  int length = 785;// 539 for sculpted tip. 785 for cleaved and sanded tip, 790 for cleaved not sanded
  TH1F* hIntensity = new TH1F("hIntensity", "Intensity vs Position;Position;Intensity", length +1, 0, length);
  ntuple->Draw("position>>hIntensity", "intensity","hist");
  ntuple->SaveAs("Values_Cleaved_Sanded.root");
  return;
}
