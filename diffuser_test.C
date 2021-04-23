const TVector3 DrawRandomRayFrom(TH2F* hXY, float d); //draw a ray (angles right, length arbitrary) from a distribution that makes an intensity map hXY a distance d from a source.
const TVector3 DiffusePhoton(const TVector3 photon_direction, TF1* angleDist); //modify a vector by two orthogonal random draws from the angle (deg) distribution.
TH2F* LoadLaserProfileFromFile(const char* profname, const char* tuplefile, const char* tuplename, float length);//build a 2D intensity profile from a file containing a 1D ntuple, and name the resulting profile 'profname'.  The total number of entries span a length of 'length' in cm.
TH2F* LoadLaserProfileFromFile2(const char* profname, const char* tuplefile, const char* tuplename, float length, float center);//version 2 of the above. monte carlo method
TH2F* LoadLaserProfileFromGaussian(const char* profname, float sigma, float length); //build a 2D intensity profile from a 1D gaussian with sigma 'sigma'cm spanning a length 'length' in cm.
const TVector3 PrismFacet(const TVector3 photon_direction, float IndexBefore, float IndexAfter, int numberprisms, const TVector3 Prism_NVector);

void diffuser_test() {
	float prismIndex = 1.500029;
	float nPrisms = 1;
	int DiffuseBounds = 80;// change to 55 for ed50, 80 for ed40
	float dist = 11.0;// distance to profile card
	float nsteps = 9;
	float prismAngle_theta = 10;
	float prismAngle_phi = 0;
	float laser_tilt = 0;
	float EdTransPercentile = 90;
	float nphoton_satistics = 1000;
	//TH2F* hIntensity;
	TVector3 prism_normal;
	TVector3 laser_nominal;	
	std::vector<float> Angle_in;
	std::vector<float> Angle_out;
	
	TF1* fEdAngle = new TF1("fEdAngle", "([0]+[1]*exp(-(x-[2])**2/(2*[3]**2)))", -DiffuseBounds, DiffuseBounds);
	fEdAngle->SetParameters(0, 0.9483733002234809, 0, 16.1331011690544950);//test of gaussian fit to image of 40 deg//***** 40 degree Edmond Optics Diffuser; Bound -80 to 80
	//hIntensity = LoadLaserProfileFromFile2("CleavedSandedTip", "Values_Cleaved_Sanded.root", "ntuple", 14.0, 400);
    laser_nominal.SetXYZ(0, sin(laser_tilt), cos(laser_tilt)); 
	for (int i = 0; i < nsteps; i++) {
		TVector3 photon_direction = laser_nominal;
		//TVector3 xz_dir_pre = photon_direction;
		//xz_dir_pre.SetY(0);
		float diffuser_tilt_angle = i * prismAngle_theta - 40;
		Angle_in.push_back(diffuser_tilt_angle);printf("**Angle In=%f\n", Angle_in[i]);
		prism_normal.SetXYZ(0, 0, 1);
		prism_normal.RotateX((diffuser_tilt_angle) * TMath::Pi() / 180);
		prism_normal.RotateZ(prismAngle_phi * TMath::Pi() / 180);
		photon_direction = PrismFacet(photon_direction, 1, prismIndex, 1, prism_normal);
		TH2F* hPhotonAtSurface = new TH2F("hPhotonAtSurface", "Photon Position At wall;x(cm);y(cm)", 80, -100, 100, 80, -100, 100);
		for(int j=0; j<nphoton_satistics;j++){
		TVector3 photon_direction_temp = DiffusePhoton(photon_direction, fEdAngle);
		printf("**direction=%f*,%f*,%f\n", photon_direction_temp(0), photon_direction_temp(1), photon_direction_temp(2));
		float distance_scale = 10 / photon_direction_temp(2); //what is photon position at some random wall. say 10 units away
			photon_direction_temp *= distance_scale;
		printf("**direction=%f*,%f*,%f\n", photon_direction_temp(0), photon_direction_temp(1), photon_direction_temp(2));
		TVector3 xzt_dir_post = photon_direction_temp;
		hPhotonAtSurface->Fill(photon_direction_temp(0), photon_direction_temp(1));
		TH2F* hPhotonAtSurface = new TH2F("hPhotonAtSurface", "Photon Position At wall;x(cm);y(cm)", 80, -100, 100, 80, -100, 100);
		}		
		float wall_coord_x = hPhotonAtSurface->GetMean(1);
		float wall_coord_y = hPhotonAtSurface->GetMean(2);
		TVector3 v3(wall_coord_x, wall_coord_y, 10);
		photon_direction = v3;

		TVector3 xz_dir_post = photon_direction;
		float Post_angle_rad = xz_dir_post.Theta();
		//xz_dir_post.SetY(0);
		if(TMath::Pi()/2 -.1<xz_dir_post.Phi()&&xz_dir_post.Phi()< TMath::Pi()/2+.1)
		{
			Post_angle_rad *= -1;// shorthand for multiply by neg 1 and replace
		}
		Angle_out.push_back(Post_angle_rad * 180 / TMath::Pi());
		printf("Angle Out=%f\n", Angle_out[i]);
	}
	TCanvas* c1 = new TCanvas("c1", "A Simple Graph Example", 900, 900);
	c1->cd(1);
	TGraph* gone = new TGraph(Angle_out.size(), &(Angle_in[0]), &(Angle_out[0]));
	gone->SetTitle("Angle in vs angle out;Angle in;Angle out");
	gone->Draw("AC*");
	return;

};




const TVector3 DrawRandomRayFrom(TH2F* hXY, float d) {
	//draw a ray (angles right, length 1) from a distribution that makes an intensity map hXY a distance d from a source.
	//we center it automatically in this code, so it does not need to be centered in advance.
	double x, y;
	hXY->GetRandom2(x, y);
	float x_mean = hXY->GetMean(1);
	float y_mean = hXY->GetMean(2);
	TVector3 ray(x - x_mean, y - y_mean, d);
	return ray.Unit(); //return a unit vector parallel to ray.
}

const TVector3 DiffusePhoton(const TVector3 photon_direction, TF1* angleDist) {
	//assume distribution is 1D, uncorrelated, so that we can apply two independent orthogonal rotations.
	TVector3 output_direction = photon_direction;//diffuse wrt photon direction

	TVector3 localA = output_direction.Orthogonal();//we have no intrinsic guarantee /which/ orthogonal direction this is
	TVector3 localB = output_direction.Cross(localA);//but we know that p x A will be orthogonal to both p and A, so these suffice.
	output_direction.Rotate(TMath::Pi() / 180 * angleDist->GetRandom(), localA);//convert to radians and rotate around A
	output_direction.Rotate(TMath::Pi() / 180 * angleDist->GetRandom(), localB);//then convert to radians and rotate  around B.
	return output_direction;
}

TH2F* LoadLaserProfileFromFile(const char* profname, const char* tuplefile, const char* tuplename, float length) {
	TFile* f = TFile::Open(tuplefile);
	TNtuple* ntuple = (TNtuple*)(f->Get(tuplename));
	float NBins = ntuple->GetEntries();//540;// Number of Bins in Bob's Histogram
	TH1F* hIntensity = new TH1F("hIntensity", "Intensity vs Position;position;intensity", NBins, 0, NBins);
	ntuple->Draw("position>>hIntensity", "intensity", "goff");

	float Full_Length = 14;//cm; length of Bob's histogram
	float Fiber_Scale = Full_Length / NBins; // Length of one bin in cm

	TH2F* hist = new TH2F(Form("h%s", profname), Form("2D Output from %s;x position(cm);y position(cm)", profname), NBins, -Full_Length / 2, Full_Length / 2, NBins, -Full_Length / 2, Full_Length / 2);
	for (int i = 0; i < NBins; i++) {
		float x_i = (-Full_Length / 2) + (i + 0.5) * Fiber_Scale;
		for (int j = 0; j < NBins; j++) {
			float y_i = (-Full_Length / 2) + (j + 0.5) * Fiber_Scale;
			hist->Fill(x_i, y_i, hIntensity->GetBinContent(hIntensity->FindBin(i)) * hIntensity->GetBinContent(hIntensity->FindBin(j)));
		}
	}
	return hist;
}
//*/
//*
TH2F* LoadLaserProfileFromFile2(const char* profname, const char* tuplefile, const char* tuplename, float length, float center) {
	TFile* f = TFile::Open(tuplefile);
	TNtuple* ntuple = (TNtuple*)(f->Get(tuplename));
	float NBins = ntuple->GetEntries();//540;// Number of Bins in Bob's Histogram
	TH1F* hIntensity = new TH1F("hIntensity", "Intensity vs Position(Pixel);position(Pixel);intensity", NBins, 0, NBins);
	ntuple->Draw("position>>hIntensity", "intensity", "goff");
	float Full_Length = 14;//cm; length of Bob's histogram
	float Fiber_Scale = Full_Length / NBins; // Length of one bin in cm
	TH1F* hIntensityPos = new TH1F("hIntensityPos", "Intensity vs Position(cm);position(cm);intensity", NBins, -Full_Length / 2, Full_Length / 2);
	TH2F* hist = new TH2F(Form("h%s", profname), Form("2D Output from %s;x position(cm);y position(cm)", profname), NBins, -Full_Length / 2, Full_Length / 2, NBins, -Full_Length / 2, Full_Length / 2);
	for (int i = 0; i < NBins; i++) {
		float x_i = (-Full_Length / 2) + (i + 0.5) * Fiber_Scale;
		//hIntensityPos->Fill(x_i, hIntensity->GetBinContent(hIntensity->FindBin(i)));
		for (int j = 0; j < NBins; j++) {
			float y_i = (-Full_Length / 2) + (j + 0.5) * Fiber_Scale;
			float r_i = sqrt((x_i * x_i + y_i * y_i));

			hist->Fill(x_i, y_i, (hIntensity->GetBinContent(hIntensity->FindBin(r_i / Fiber_Scale + center))));
		}
	}
	/*
	 TCanvas* c = new TCanvas("TestBin", "TestBin", 2400, 1440);
	 c->Divide(3, 1);
	 c->cd(1);
	 TLine line;
	 hIntensity->Draw("hist");
	 line.DrawLine(400, 0, 400, 260);
	c->cd(2);
	hIntensityPos->Draw("hist");
	line.DrawLine(0, 0, 0, 260);
	c->cd(3);
	hist->Draw("colz");
	c->SaveAs(Form("Distribution_Check_%s.pdf", profname));
	*/
	return hist;
}
//*/

TH2F* LoadLaserProfileFromGaussian(const char* profname, float sigma, float length) {
	double sigma_ideal = sigma;
	double sigma_y = sigma_ideal;
	double sigma_x = sigma_ideal;

	int NBins = 100;
	float Full_Length = length;//cm; length of histogram
	float Fiber_Scale = Full_Length / NBins; // Length of one bin in cm


	TH2F* hist = new TH2F(Form("h%s", profname), Form("2D Output from %s;x position(cm);y position(cm)", profname), NBins, -Full_Length / 2, Full_Length / 2, NBins, -Full_Length / 2, Full_Length / 2);


	for (int i = 0; i < NBins; i++) {
		float x_i = (-Full_Length / 2) + (i + 0.5) * Fiber_Scale;
		for (int j = 0; j < NBins; j++) {
			float y_i = (-Full_Length / 2) + (j + 0.5) * Fiber_Scale;
			hist->Fill(x_i, y_i, exp(-2 * ((x_i * x_i) / (sigma_x * sigma_x) + (y_i * y_i) / (sigma_y * sigma_y))));
		}
	}
	return hist;
}
//*
const TVector3 PrismFacet(const TVector3 photon_direction, float IndexBefore, float IndexAfter, int numberprisms, const TVector3 Prism_NVector) {
	TVector3 output_direction = photon_direction;// Tilt w.r.t photon directiona dn prism facet direction
	float theta_before = output_direction.Angle(Prism_NVector);//angle with respect to the (backward-facing) normal of the kth facet of prism L.
	//******* the negative here(attached to prism normal) is causing errors. I expect it is because it is taking the difference of the actual directions, not the angle used in snell's law(complementary to this).
	float arg_asin = IndexBefore / IndexAfter * sin(theta_before);//test adding neg to theta before
	/*
	if (arg_asin >= 1) { //out of probability range, photon is absorbed}
		PrismLoss++;
		islost = true;
		break; //skip out of this loop
	}
	*/
	float theta_after = asin(arg_asin);//angle we want for the outgoing photon wrt the (forward-facing)normal.

	TVector3 ortho = output_direction.Cross(Prism_NVector);//the rotation axis perpendicular to the prism and photon normals.
	output_direction = Prism_NVector; //to get the new direction, start with the normal vector
	output_direction.Rotate(theta_after, ortho);//then rotate about the orthogonal direction the specified amount.
						//2nd apex angle version for comparison

	return output_direction;

}
