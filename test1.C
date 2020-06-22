void test1() {
	/*TNtuple* ntuple = new TNtuple("ntuple", "data from bob azmoun", "position:intensity");
	ntuple->ReadFile("test.csv");*/
   float NBins = 540; 
   float Full_Length = 14;//cm
   float Fiber_Scale = Full_Length/NBins;
	TFile* f = TFile::Open("ntuple.root");
	TNtuple* ntuple = (TNtuple*)(f->Get("ntuple"));
	TH1F* hIntensity = new TH1F("hIntensity", "Intensity vs Position;position;intensity", NBins, 0, NBins-1);
	ntuple->Draw("position>>hIntensity", "intensity");
	TH2F* h2= new TH2F("h2", "Filled Test Plane", NBins, -Full_Length/2, Full_Length / 2, NBins, -Full_Length / 2, Full_Length / 2);
	/*e ^ (-((x - x_0) ^ 2 + (y - y_0 ^ 2) / sigma)
	F(x = x_i, y = y_i) = F(x_i) * F(y_i)*/
TH2F* h4 = new TH2F("h4","Ideal Gaussian",NBins,-Full_Length/2, Full_Length / 2, NBins, -Full_Length / 2, Full_Length / 2);
gRandom = new TRandom3();
double sigma_ideal = 7;
double sigma_y = sigma_ideal;
double sigma_x = sigma_ideal;
	for (int i = 0; i < NBins; i++) {
         
		float x_i = (-Full_Length/2)+i*Fiber_Scale;
		float x_0_i = i;
			for (int j = 0; j < NBins; j++) {
				float y_0_i = j;
				float y_i = (-Full_Length / 2) + j * Fiber_Scale;
				h2->Fill(x_i,y_i,hIntensity->GetBinContent(hIntensity->FindBin(x_0_i))* hIntensity->GetBinContent(hIntensity->FindBin(y_0_i)));
                h4->Fill(x_i,y_i,exp(-2*((x_i*x_i)/(sigma_x * sigma_x)+ (y_i * y_i) / (sigma_y * sigma_y))));

		}
	
	}
	/*
	x_i is a coordinate in h2(runs from -Full_Length/2 to Full_Length/2 in Fiber_Scale
	
	
	*/
	double x,y,z;
	double x_mean=h2->GetMean(1);
	double y_mean = h2->GetMean(2);
	z = 10;// z=10 cm . from Bob's setup
	/*h4->Draw("colz");*/ // draw the ideal angular distribution
	TH1D* px = h4->ProjectionX("px", -7, 7); // where firstYbin = -7 and lastYbin = 7
	px->Draw();

	TH2F* h3 = new TH2F("h3","theta vs. phi;Theta;Phi",500, -1, 1, 500, -2*TMath::Pi(), 2*TMath::Pi());
	for (int i = 0; i < 5000000; i++) {
       h2->GetRandom2(x, y);
	    x = x - x_mean;
	    y = y - y_mean;
	    double phi =atan2(y,x);
	    double theta=atan2(sqrt(x*x +y*y),z);
		h3->Fill(theta,phi);
	
	}
	//h3->Draw("colz");
	// to handle histogram centering issues. Get mean of x and y from TH2 functions. and modify angle formulas to reflect shift of distribution center.  

	
	
	//for (int i = 1; i <= NBins - 270; i++)
	//	hIntensity->SetBinContent(i, hIntensity->GetBinContent(i + 270));
	//for (int i = NBins - 270; i <= NBins; i++) hIntensity->SetBinContent(i, 0);
	/*hIntensity->SetBins(540,-7,7);*/
    
	//h->GetBinContent(h->FindBin(x)); // Find bin coresponding to coordinate and return value of that bin
	//h2->Fill(hIntensity, hIntensity);
	//ntuple->Fit("gaus","area");
	//position=hIntensity->GetRandom()
	//pos_cm=position*scaling_factor
	// get some angle=pos_cm/distance_from_laser_in_cm different eqn obv
	return;
}