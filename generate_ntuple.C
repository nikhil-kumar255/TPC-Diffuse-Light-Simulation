void generate_ntuple(){
  TNtuple* ntuple = new TNtuple("ntuple", "data from bob azmoun", "position:intensity");
  ntuple->ReadFile("test.csv");
  TH1F* hIntensity = new TH1F("hIntensity", "Intensity vs position, I think?;intensity;position", 540, 0, 539);
  ntuple->Draw("position>>hIntensity", "intensity");
  ntuple->SaveAs("ntuple.root");
  return;
}
