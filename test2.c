const TVector3 DrawRandomRayFrom(TH2F* hXY, float d); //draw a ray (angles right, length arbitrary) from a distribution that makes an intensity map hXY a distance d from a source.
const TVector3 DiffusePhoton(const TVector3 photon_direction, TF1* angleDist); //modify a vector by two orthogonal random draws from the angle (deg) distribution.
TH2F* LoadLaserProfileFromFile(const char* profname, const char* tuplefile, const char* tuplename, float length);//build a 2D intensity profile from a file containing a 1D ntuple, and name the resulting profile 'profname'.  The total number of entries span a length of 'length' in cm.
TH2F* LoadLaserProfileFromFile2(const char* profname, const char* tuplefile, const char* tuplename, float length, float center);//version 2 of the above. monte carlo method
TH2F* LoadLaserProfileFromGaussian(const char* profname, float sigma, float length); //build a 2D intensity profile from a 1D gaussian with sigma 'sigma'cm spanning a length 'length' in cm.
const TVector3 PrismFacet(const TVector3 photon_direction, float IndexBefore, float IndexAfter, int numberprisms, const TVector3 Prism_NVector);

std::vector<float> SimulateLasers(int nLasers, float laser_tilt_angle, TH2F* hIntensity, float dist,
	int nDiffusers, float diffTransmission, TF1* diffAngleProfile,
	int nPrisms, float prismTransmission, float prismIndex, float prismAngle,
	const char* canvasname);

void test2() {

	const int nLasers = 1;
	int nDiffusers = 0;
	float laserthetaparameter = 0.1;// only meaningful when not using source file
	float laser_tilt_angle = 0;
	float prismAngle = -5;
	//float prismIndex = 1.0;
	float prismIndex = 1.500029;
	float nPrisms = 1;
	int DiffuseBounds = 55;// change to 55 for ed50, 80 for ed40
	float dist = 11.0;// distance to profile card
	//set up the diffuser parameters:
	//null angle diffuser returns almost exactly what we take in:
	TF1* fNullAngle = new TF1("fThorAngle", "[0]*exp(-0.5*(x/[1])**2)", -30, 30);
	fNullAngle->SetParameters(1, 0.001); //test with very narrow gaussian to make sure we get back what we put in.

	//thorlabs diffuser with a by-hand fit to the model on their page
	TF1* fThorAngle = new TF1("fThorAngle", "([0]**2/((x-[1])**2+[0]**2))", -30, 30);
	fThorAngle->SetParameters(7.5, 0);//thorlabs model is 7.5, itterating for edmond 50 deg->12 seems fair
	fThorAngle->SetTitle("Thorlabs diffuser angular distribution;angle (deg);arb. prob.");
	float ThorTransPercentile = 68.35;// note 68.35 for 10-220, 57.6 for 10-120, 90 for holographic diffuser from edmond
	  //90;//percent.
	 //Edmond optics diffuser with a by-hand fit to the model on their page
	///*
	TF1* fEdAngle = new TF1("fEdAngle", "([0]+[1]*exp(-(x-[2])**2/(2*[3]**2)))", -DiffuseBounds, DiffuseBounds);
	//***** 40 degree Edmond Optics Diffuser; Bound -80 to 80
	//fEdAngle->SetParameters(0, 0.9483733002234809, 0, 16.1331011690544950);//test of gaussian fit to image of 40 deg
	// note first parameter A=0.06310490730073468
	// note second parameter H=0.9483733002234809// with A=x_0=0->0.9991798327879257
	// note third parameter x_0=-0.40679912703211746
	// note fourth parameter W=16.1331011690544950// with A=x_0=0->17.41119227449015
	//*****  50 Degree Edmond Optics Diffuser; Bound -55 to 55
	fEdAngle->SetParameters(0, 0.58197, 0, 29.0055);//test of gaussian fit to image of 50 deg
  // note first parameter A=0(by force)
  // note second parameter H=0.58197
  // note third parameter x_0=3.81880
  // note fourth parameter W=29.0055
	fEdAngle->SetTitle("Edmond optics diffuser angular distribution;angle (deg);arb. prob.");
	float EdTransPercentile = 90;// note 68.35 for 10-220, 57.6 for 10-120, 90 for holographic diffuser from edmond
	//*/
	TH2F* hIntensity;
	//load Bob's sculpted-tip laser profile from file(stored in root/Bin):
	//hIntensity=LoadLaserProfileFromFile("SculptedTip","ntuple.root","ntuple",14.0);
	hIntensity = LoadLaserProfileFromFile2("CleavedSandedTip", "Values_Cleaved_Sanded.root", "ntuple", 14.0, 400);
	//hIntensity=LoadLaserProfileFromFile2("CleavedNotSandedTip","Values_Cleaved_Not_Sanded.root","ntuple",14.0,390);

	//load a gaussian laser profile:
	//hIntensity=LoadLaserProfileFromGaussian("PureGaussianSig1",1,14.0);// sigma 7 is  spread
	//load an exponetial
	//hIntensity = LoadLaserProfileFromExponential("Exponential0.1", 0.1, 14.0);// sigma 7

	//now you can make loops that run this multiple times, if you like.  
	float laser_tilt_angles[] = { 8,11,14,18 };// for 4 and 5:,19,22     add on s[nDiffusers]
	// tilt angle repository: { 8,12,14,18 } for ed40 both CleavedTips,
	// ed50:for CleavedNotSandedTip, { 9,13,17,18 } for CleavedSandedTip
	float prismAngles[] = { -3.5,-6,-7,-9 };//{ 0,-4.5,-6,-8}for ed40 gauss1;
	// prism angle repository: {-3,-6,-7,-8} for ed40 CleavedNotSandedTip, {-3.5,-6,-7,-9} for ed40 CleavedSandedTip
	// ed50:for CleavedNotSandedTip, for CleavedSandedTip
	////////////////////////////////////code variation to run
	/////////////////1 max:Intensity & Max uniformity
	/*
	int angle_step = 1;
	for (int angle_itt = 0; angle_steps < 11; angle_steps++)
{
			prismAngle = prismAngle - angle_step;
			SimulateLasers(nLasers, laser_tilt_angle, hIntensity, dist,
			nDiffusers, EdTransPercentile, fEdAngle,
			nPrisms, 100, prismIndex,& prismAngle,
			Form("TEST_Diffuser_snsu3_CleavedSanded11_Ed40_%d", nDiffusers));

	}*/
	////////////////end 1
	////////////////2
	/*
	for (int nDiffusers = 0; nDiffusers < 4; nDiffusers++) {//*/
	//*
	std::vector<float> Light_Frac;
	std::vector<float> vecI;
	std::vector<float> Ave_Ratio;
	std::vector<float> AngLe;
	std::vector<float> Peak_x;
	std::vector<float> Peak_y;
	for (int prismAngle_itt = 0; prismAngle_itt < 10; prismAngle_itt++) {//*/
		vecI = SimulateLasers(nLasers, laser_tilt_angle, hIntensity, dist,
		nDiffusers, EdTransPercentile, fEdAngle,
		nPrisms, 100, prismIndex, -prismAngle_itt,
		Form("TEST2_TiltedDiffuser_snsu3_CleavedSanded_Ed50_Angle%d", -prismAngle_itt));
	//*/
	// comment out the brace below if not running loop oversomething
		//*
	AngLe.push_back(prismAngle_itt);
	//printf("ratio, yield=%f,%f\n", vecI[1], vecI[0]);
	Light_Frac.push_back(vecI[1] * 100);
	Ave_Ratio.push_back(vecI[0]);
	Peak_x.push_back(vecI[2]);
	Peak_y.push_back(vecI[3]);
	}//*
	//TCanvas* c1 = new TCanvas("c1", "A Simple Graph Example", 200, 10, 500, 300);
	TCanvas* c1 = new TCanvas("c1", "A Simple Graph Example", 700, 900);
	c1->Divide(2, 3);
	c1->cd(1);
	TGraph* gone = new TGraph(Light_Frac.size(), &(Ave_Ratio[0]), &(Light_Frac[0]));
	gone->SetTitle("Yield vs uniformity;Center to Inner Ratio;Yield");
	gone->Draw("AC*");
	c1->cd(2);
	TGraph* gtwo = new TGraph(Light_Frac.size(), &(AngLe[0]), &(Light_Frac[0]));
	gtwo->SetTitle("Yield vs Angle;Yield;Degrees");
	gtwo->Draw("AC*");
	c1->cd(3);
	TGraph* gthree = new TGraph(Light_Frac.size(), &(AngLe[0]), &(Ave_Ratio[0]));
	gthree->SetTitle("Ratio vs Angle;Angle(deg);Center to Inner Ratio");
	gthree->Draw("AC*");
	c1->cd(4);
	TGraph* gfour = new TGraph(Light_Frac.size(), &(AngLe[0]),&(Peak_y[0]) );
	gfour->SetTitle("Peak y position vs angle;Angle(deg);Peak y position(cm)");
	gfour->Draw("AC*");
	c1->cd(5);
	TGraph* gfive = new TGraph(Light_Frac.size(),&(AngLe[0]) ,&(Peak_x[0]) );
	gfive->SetTitle("Peak x position vs angle;Angle(deg);Peak x position(cm)");
	gfive->Draw("AC*");
  //*/
  /////////////////end 2
	return;
};




std::vector<float> SimulateLasers(int nLasers, float laser_tilt_angle, TH2F* hIntensity, float dist,
	int nDiffusers, float diffTransmission, TF1* diffAngleProfile,
	int nPrisms, float prismTransmission, float prismIndex, float prismAngle,
	const char* canvasname)
{
	//define the basis vectors local to the laser:
	TVector3 laser_nominal[nLasers];
	TVector3 laser_transverse[nLasers];

	//define the position of the lasers:
	TVector3 laser_position[nLasers];
	float laser_position_angle0 = 0;

	float angle_increment = 2 * TMath::Pi() / nLasers; // assume they're equally spaced

	//define the parameters of the laser stack:
	float laser_tilt = laser_tilt_angle * TMath::Pi() / 180; //and have a common tilt, leaning outward in the rz plane.
	//nPrisms=1; //hard-wired for now.
	float prism_normal_theta_deg =  prismAngle ; //one for each facet positive=tilted out
	float prism_normal_phi_deg=  0 ; //one for each facet positive=image rotated CW on the CM
	float prism_index[] = { 1.0,1.0,prismIndex };//one for each region separated by a facet


	TVector3 prism_normal[nLasers];
	//define 
	//
	TH1F* PhotonSurvival = new TH1F("PhotonSurvival", "Photon Survival at various stages", 16, -0.5, 15.5);

	//rotate each laser to the proper position and update the beam nominal and transverse direction.
	for (int i = 0; i < nLasers; i++)
	{
		laser_nominal[i].SetXYZ(0, sin(laser_tilt), cos(laser_tilt));// to alternate tilt incorporate (-1)^i:pow(-1,i)*
		//printf("%f\n",laser_nominal[i].X());
		//laser_nominal[i].SetXYZ(0, 0, 1);
		laser_transverse[i].SetXYZ(1, 0, 0);
		laser_position[i].SetXYZ(0.0, 40.0, 0);
		//rotate the basis vectors of the laser:
		//laser_nominal[i].RotateX(laser_tilt);
		//laser_transverse[i].RotateX(laser_tilt);
		//printf("%f\n", laser_nominal[i].X());
		laser_nominal[i].RotateZ(laser_position_angle0 + angle_increment * i);
		laser_transverse[i].RotateZ(laser_position_angle0 + angle_increment * i);
		//rotate the position of the laser:
		laser_position[i].RotateZ(laser_position_angle0 + angle_increment * i);

		
			//set up the nominal orientation if this were at (x=0,y>0):
			prism_normal[i].SetXYZ(0, 0, 1);
			prism_normal[i].RotateX(prism_normal_theta_deg * TMath::Pi() / 180);
			prism_normal[i].RotateZ(prism_normal_phi_deg * TMath::Pi() / 180);
			//rotate this for the actual position of this laser stack:
			prism_normal[i].RotateZ(laser_position_angle0 + angle_increment * i);
		
	}




	//define the surface we are illuminating:
	TVector3 surface_point(0, 0, 105.2);// surface is actually 105.2
	TVector3 surface_normal(0, 0, 1);


	//laser intensity as a function of distance from the beam center:
	float w0 = 0.1;
	float zsample = 1.0;//all in cm.
	float lambda = 266 * 1e-7;//10^-7 converts from nm in cm.
	float denom = w0 * w0 + lambda * lambda * zsample * zsample / (TMath::Pi() * TMath::Pi() * w0 * w0);
	float thetamax = TMath::Pi() / 4;
	float samplemax = 100; //
	float samplemax2 = TMath::Tan(thetamax) * zsample; // max radial distance = max angle at zsample distance.
	printf("samplemax2=%f\n", lambda);
	if (samplemax2 > samplemax)
		samplemax = samplemax2;
	TF1* laser_r = new TF1("laser_r", "x*exp(-2*x^2/[0])", 0, samplemax); // intensity as function of r (note we have performed the integral over phi!)
	//TF1 *laser_r=new TF1("laser_r","x/[0]",0,TMath::Tan(thetamax)*zsample); // intensity as function of r
	laser_r->SetParameter(0, denom);




	//plane in vector notation:
	//(v-p0).n=0 -- the position of point v relative to the defined point has no component out of the plane 
	//line in vector notation:
	//v=t*L+L0 -- the position of point v is the starting position plus some multiple of its direction vector.
	//sub in:
	//(t*L+L0-p0).n=0
	//t*(L.n)+L0.n-p0.n=0
	//t=(p0-L0).n/(L.n) -- as long as L.n!=0;
	//the intersection point is:
	//v={(p0-L0).n/(L.n)}*L+L0
	//(and so the /starting/ point as a function of intersection point is:
	//L0=v-{(p0-L0).n/(L.n)}*L
	//L0=(v-p0.n/(L.n)*L)+L0.n/(L.n)*L
	//since in our case the laser position has z=0,  L0.n=Lz*zhat=0
	//L0=v-p0.n/(L.n)*L
	// )
	//conceptually, point backward from the intersection point in the negative laser direction until we hit the plane defined by the laser direction and the laser position. This tells us the 

	//when in doubt, google "root [class]" like "root TH2F"
	// TH2F(name, title and axes, number of bins, minimum coordinate, maximum, number of y bins, minimum, maximum)
	TH2F* hPhotonAtSurface = new TH2F("hPhotonAtSurface", "Photon Position At CM;x(cm);y(cm)", 80, -100, 100, 80, -100, 100);
	TH2F* hPhotonAngle = new TH2F("hPhotonAngle", "Photon Angle;#theta (x);#phi (y)", 50, 0, 2, 50, 0, 6.5);
	TH1F* hPreDiffusionAngle = new TH1F("hPreDiffusionAngle", "Photon Y Angle Before Diffusion;#theta (deg)", 50, -45, 45);
	TH1F* hPostDiffusionAngle = new TH1F("hPostDiffusionAngle", "Photon Y Angle After Diffusion;#theta (deg)", 50, -45, 45);
	TH2F* hDiffuserDeflection = new TH2F("hDiffuserDeflection", "Photon Angle before and after diffuser;#theta (pre);#theta (post)", 50, -10, 180, 50, -10, 180);
	TH2F* hPhotonDirection = new TH2F("hPhotonDirection", "hPhotonDirection; (x); (y)", 50, -2, 2, 50, -2, 2);
	TH1F* hPhoton = new TH1F("hPhoton", "hPhoton", 100, 0, TMath::Pi());
	//*//Initialize File
	FILE* pFile;
	int n;
	char name[100];
	pFile = fopen("test_trace.csv", "w");
	//*//
	int nPhotons = 5000000;
	//****Begin Loss Parameters
	//int PreFieldCageSurvival = 1;
	int IFCLoss = 0;
	int OFCLoss = 0;
	//int PrePrismSurvival = 1;
	int PrismLoss = 0;
	//int PreDiffuserSurvival[10] = { };
	int DiffuserAbsorbtion[10] = { };
	int DiffuserGeometric[10] = { };
	//****End Loss Parameters

	for (int i = 0; i < nPhotons; i++) {
		int L = i % nLasers;
		int photonstep = 0;
		//PhotonSurvival->Fill("GenP",1); photonstep++;
		PhotonSurvival->Fill(photonstep); photonstep++;
		//draw a photon direction from our source distribution:
		TVector3 local_direction = DrawRandomRayFrom(hIntensity, dist); //d=10.0cm 
		float photon_theta = local_direction.Theta();
		float photon_phi = local_direction.Phi();

		//draw a photon positional offset from our laser width:
		float photon_r = laser_r->GetRandom();
		//hPhotonAngle->Fill(photon_theta, photon_phi);//fill x and y.  annoying that it's backward from the arguments earlier, but what're you gonna do? ;)
		//*/
		//****** end standard gaussian laser
		//****************Custom Laser fiber section
		/*
		  z = 10;// z=10 cm . from Bob's setup
		  h2->GetRandom2(x, y);
		  x = x - x_mean;
		  y = y - y_mean;
		  double photon_phi= atan2(y,x);
		  float photon_r = 0; // placeholder for bob
		  double photon_theta = atan2(sqrt(x * x + y * y), z);
		//*/
		//***************End Custom laser Fiber section
		//***************Ideal Laser Fiber section
		/*
		  z = 10;
		  h4->GetRandom2(x, y);
		  double photon_phi = atan2(y, x);
		  float photon_r = 0; // placeholder for bob
		  double photon_theta = atan2(sqrt(x * x + y * y), z);
		//*/
		//***************End Ideal laser Fiber section



		//************Lens Section
		float focal_length = 1.0;//focal length in centimeters //float photon_theta=something to do with photon_r
		//photon_theta = photon_r / focal_length + photon_theta; //thin lens equation angle change. f is focal length for DIVERGING lens
		//************End Lens Section

		//photon starts pointing in the laser_nominal direction;
		TVector3 photon_direction = laser_nominal[L];
		//Rotate photon by theta and phi
		photon_direction.Rotate(photon_theta, laser_transverse[L]);
		TVector3 photon_position = laser_position[L];
		//I don't think we handled this right.  I also don't think we need it:
		TVector3 photon_offset = laser_transverse[L] * photon_r;
		photon_offset.Rotate(photon_phi, laser_nominal[L]);
		photon_direction.Rotate(photon_phi, laser_nominal[L]);
		photon_position = photon_position + photon_offset;
		TVector3 xz_dir_pre = photon_direction;
		// initialize apex angle variables 
		float apex_angle;//something seems very wrong with this. 1+0'
		float apex_angle_2;//seems more reasonable. 1'-2
		float theta_zero;
		float theta_zero_prime;
		float theta_one;
		float theta_one_prime;
		float theta_two;
		//****Prism Section
		//PrePrismSurvival++;
		xz_dir_pre = photon_direction;
		int prismfacet = 0;
		bool Geolost = false;
		bool islost = false;
		//const TVector3 PrismFacet(const TVector3 photon_direction, float IndexBefore, float IndexAfter, int numberprisms, const TVector3 Prism_NVector);
		photon_direction = PrismFacet(photon_direction, 1, prismIndex, nPrisms, prism_normal[L]);
//should consider critical angle of total internal reflection, too, just to be thorough...
			if (photon_direction.Theta() * 180 / TMath::Pi() > 90) {
				islost = true;
				PrismLoss++;
				break;
			}
			prismfacet++;
		
		if (islost) continue;
		PhotonSurvival->Fill(photonstep); photonstep++;
		//PhotonSurvival->Fill("PrG", 1); photonstep++;
		//need to consider reflections and other light losses at some point.  

		//end of prism

		//****Diffuser Section
		//the diffuser adds an additional random set of rotations w.r.t. the photon direction
		//rather than rotate theta and phi around the photon axis, we get two
		// orthogonal vectors relative to the photon direction, and do a 'rotate
		// about x axis' followed by a 'rotate about y axis' in that frame:

		//angle before the diffuser, for our records.  by fiat, let's always measure the 1D angular spread as the spread in the XZ plane:  y=0;

		xz_dir_pre.SetY(0);
		hPreDiffusionAngle->Fill(xz_dir_pre.Theta() * 180 / TMath::Pi());


		for (int j = 0; j < nDiffusers; j++) {
			photon_direction = DiffusePhoton(photon_direction, diffAngleProfile);
			float rand1 = rand() % 100; //random
			//PreDiffuserSurvival[j++];
			if (rand1 > diffTransmission) { //out of probability range, photon is absorbed}
				islost = true;
				DiffuserAbsorbtion[j]++;
				break; //skip out of this loop
			}
			//PhotonSurvival->Fill(Form("DA%d",j),1); photonstep++;
			PhotonSurvival->Fill(photonstep); photonstep++;
			if (photon_direction.Theta() * 180 / TMath::Pi() > 90) {
				DiffuserGeometric[j]++;
				Geolost = true;

				break;
			}
			//PhotonSurvival->Fill(Form("DG%d", j),1); photonstep++;
			PhotonSurvival->Fill(photonstep); photonstep++;
		}
		if (Geolost) continue;

		//for now, let's always measure the 1D angular spread as the spread in the XZ plane:  y=0;
		TVector3 xz_dir_post = photon_direction;
		xz_dir_post.SetY(0);
		hPostDiffusionAngle->Fill(xz_dir_post.Theta() * 180 / TMath::Pi());
		if (islost) continue;  //skip to the next photon


	   //end of diffuser
		hDiffuserDeflection->Fill(xz_dir_pre.Theta() * 180 / 3.14, xz_dir_post.Theta() * 180 / 3.14);
		//*****Collisions Section
		 //find the position where this photon intersects the CM:
		float denominator = photon_direction.Dot(surface_normal);
		if (denominator < 0.001) { //parallel to surface.  no solution}
			continue; //skip to the next iteration through the loop
		}

		//PhotonSurvival->Fill("FChk",1); photonstep++;
		PhotonSurvival->Fill(photonstep); photonstep++;

		TVector3 intersect = (surface_normal.Dot(surface_point - photon_position) / denominator) * photon_direction + photon_position;


		hPhotonAtSurface->Fill(intersect.X(), intersect.Y());









	}
		//*//              Ideal gaussian Slice, angular distribution
//TH1D* px = hPhotonAtSurface->ProjectionX("px", -100, 100); // where firstYbin = -100 and lastYbin = 100 //original attempt
		TH1D* px = new TH1D("px", "Position vs Intensity;Position;Intensity", 201, -100, 100);
		TH1D* pax = new TH1D("pax", "Angle vs Intensity;Angle;Intensity", 20, -30, 30);
		double xI, yI;
		double zI = 105.2;//distance of fiber to CM
		double yPORT = 40; //cm; trial method to create slice histogram 
		double R_out = 80;

		for (int j = 0; j < 2 * R_out + 1; j++) {

			float xI_0 = j - R_out;
			//hPhotonAtSurface->GetRandom2(xI, yI);
			double ThetaI = (180 / TMath::Pi()) * atan2(xI_0, zI);// in degrees when including 180/TMath::Pi()
			px->Fill(xI_0, hPhotonAtSurface->GetBinContent(hPhotonAtSurface->FindBin(xI_0, yPORT)));
			//remember, find bin gives bin number based on coordinates, and  get bin content gives the value stored in that bin	
			pax->Fill(ThetaI, px->GetBinContent(px->FindBin(xI_0)));
		}//Draw a lot of useful histograms:
		int nCells = hPhotonAtSurface->GetNcells();//total number of bins in histogram
		TH1F* hIntensityProfile = new TH1F("hIntensityProfile", "Histogram of Intensity per bin;intensity;nbins", 100, 1, 5 * nPhotons / nCells);
		TH1F* hIntensityRadius = new TH1F("hIntensityRadius", "Mean Intensity vs Radius;radius;mean intensity", 100, 0, 100);
		TH1F* hIntensityTheta = new TH1F("hIntensityTheta", "Mean Intensity vs Angle;Theta;mean intensity", 50, -(1 / 2) * TMath::Pi(), (1 / 2) * TMath::Pi());
		TH2F* hRadiusProfile = new TH2F("hRadiusProfile", "Histogram of Intensity vs Radius for all bins, radius,intensity", 100, 0, 100, 100, 1, 5 * nPhotons / nCells);
		TH1F* hRadius = new TH1F("hRadius", "nbins vs radius, for normalization;radius;nbins", 100, 0, 100);
		TH1F* hTheta = new TH1F("hTheta", "nbins vs Angle(Theta), for normalization;Theta;nbins", 50, 0, (1 / 2) * TMath::Pi());

		//read the bins from photonAtSurface and histogram their contents to get the per-radius average intensities
		for (int j = 0; j < nCells; j++) {
			float content = hPhotonAtSurface->GetBinContent(j);
			int binx, biny, binz;
			hPhotonAtSurface->GetBinXYZ(j, binx, biny, binz);//get the per-axis bins for this global bin
			float xpos = hPhotonAtSurface->GetXaxis()->GetBinCenter(binx);
			float ypos = hPhotonAtSurface->GetYaxis()->GetBinCenter(biny);
			float zpos = 105.2;
			float radius = sqrt(xpos * xpos + ypos * ypos);
			float theta = atan2(radius, zpos);
			hIntensityProfile->Fill(content);
			hIntensityRadius->Fill(radius, content);
			hIntensityTheta->Fill(theta, content);
			hRadiusProfile->Fill(radius, content);
			hRadius->Fill(radius);
			hTheta->Fill(theta);
		}
		//average by dividing by the number of bins that apply to each radius:
		hIntensityRadius->Divide(hRadius);
		hIntensityTheta->Divide(hTheta);
		

		float light_fraction = hPhotonAtSurface->GetEntries() / nPhotons;
		int bounds[3][2] = { {25,27},{45,55},{71,75} };
		float diffs[3];
		for (int j = 0; j < 3; j++) diffs[j] = bounds[j][1] - bounds[j][0] + 1;
		//calculate some means from the intensity vs radius (averaging over several bins):
		float central_ave = hIntensityRadius->Integral(bounds[1][0], bounds[1][1]) / diffs[1];
		float low_ave = hIntensityRadius->Integral(bounds[0][0], bounds[0][1]) / diffs[0];
		float high_ave = hIntensityRadius->Integral(bounds[2][0], bounds[2][1]) / diffs[2];
		float ratio_ave = central_ave / low_ave; // how much brighter is central than low?
		float asymmetry = (low_ave - high_ave) / (low_ave + high_ave); //how much brighter is low than high?
		//hIntensity->GetXaxis()->GetBinLowEdge([0][0]);
		//hIntensity->GetXaxis()->GetBinLowEdge([0][0]);
		TCanvas* c = new TCanvas(canvasname, canvasname, 900, 900);
		c->Divide(3, 3);
		c->cd(1);
		hIntensity->Draw("colz");
		c->cd(2);
		diffAngleProfile->Draw();
		c->cd(3);
		hDiffuserDeflection->Draw("colz");
		c->cd(4);
		hPhotonAtSurface->Draw("colz");
		c->cd(5);
		hIntensityRadius->Draw("hist");
		hIntensityRadius->SaveAs("hIntensityRadius.hist.root");
		float x0, y0, x1, y1;
		TLine* line = new TLine(x0, y0, x1, y1);
		line->SetLineColor(kRed);
		//low_ave lines
		line->DrawLine(hIntensityRadius->GetXaxis()->GetBinLowEdge(bounds[0][0]), 0, hIntensityRadius->GetXaxis()->GetBinLowEdge(bounds[0][0]), hIntensityRadius->GetBinContent(hIntensityRadius->GetMaximumBin()));
		line->DrawLine(hIntensityRadius->GetXaxis()->GetBinUpEdge(bounds[0][1]), 0, hIntensityRadius->GetXaxis()->GetBinUpEdge(bounds[0][1]), hIntensityRadius->GetBinContent(hIntensityRadius->GetMaximumBin()));
		//central_ave lines
		line->SetLineColor(kGreen);
		line->DrawLine(hIntensityRadius->GetXaxis()->GetBinLowEdge(bounds[1][0]), 0, hIntensityRadius->GetXaxis()->GetBinLowEdge(bounds[1][0]), hIntensityRadius->GetBinContent(hIntensityRadius->GetMaximumBin()));
		line->DrawLine(hIntensityRadius->GetXaxis()->GetBinUpEdge(bounds[1][1]), 0, hIntensityRadius->GetXaxis()->GetBinUpEdge(bounds[1][1]), hIntensityRadius->GetBinContent(hIntensityRadius->GetMaximumBin()));
		//High_Ave lines
		line->SetLineColor(kViolet);
		line->DrawLine(hIntensityRadius->GetXaxis()->GetBinLowEdge(bounds[2][0]), 0, hIntensityRadius->GetXaxis()->GetBinLowEdge(bounds[2][0]), hIntensityRadius->GetBinContent(hIntensityRadius->GetMaximumBin()));
		line->DrawLine(hIntensityRadius->GetXaxis()->GetBinUpEdge(bounds[2][1]), 0, hIntensityRadius->GetXaxis()->GetBinUpEdge(bounds[2][1]), hIntensityRadius->GetBinContent(hIntensityRadius->GetMaximumBin()));
		c->cd(6);
		float texpos = 0.9; float texshift = 0.08;
		TLatex* tex = new TLatex(0.0, texpos, Form("PhotonSource=%s", hIntensity->GetName()));
		//
		tex->Draw(); texpos -= texshift; //draw this line and shift our position down to be ready for the next line
		tex = new TLatex(0.0, texpos, Form("nDiffs=%d, Trans=%1.1f%%", nDiffusers, diffTransmission));
		tex->Draw(); texpos -= texshift; //draw this line and shift our position down to be ready for the next line
		tex = new TLatex(0.0, texpos, Form("nPrisms=%d, Index=%1.2f, #Theta=%2.1f, Trans=%1.1f%%", nPrisms, prismIndex, prismAngle, prismTransmission));
		tex->Draw(); texpos -= texshift; //draw this line and shift our position down to be ready for the next line
		tex = new TLatex(0.0, texpos, Form("%% Photons on CM: %1.1f", light_fraction * 100));
		tex->Draw(); texpos -= texshift;
		tex = new TLatex(0.0, texpos, Form("Center(%d-%d) to Inner (%d-%d) ratio: %1.3f",
			bounds[1][0], bounds[1][1], bounds[0][0], bounds[0][1], ratio_ave));
		tex->Draw(); texpos -= texshift;
		tex = new TLatex(0.0, texpos, Form("Inner (%d-%d) to Outer (%d-%d) asym: %1.3f",
			bounds[0][0], bounds[0][1], bounds[2][0], bounds[2][1], asymmetry));
		tex->Draw(); texpos -= texshift;
		if (low_ave / high_ave > 1.1) {
			tex = new TLatex(0.0, texpos, Form("Inner %1.3f%% brighter than Outer!  Angle out more", (low_ave / high_ave - 1) * 100));
		}
		else   if (high_ave / low_ave > 1.1) {
			tex = new TLatex(0.0, texpos, Form("Outer %1.3f%% brighter than Inner!  Angle in more", (high_ave / low_ave - 1) * 100));
		}
		else {
			tex = new TLatex(0.0, texpos, Form("Outer and Inner approximately balanced."));
		}
		tex->Draw(); texpos -= texshift; //draw this line and shift our position down to be ready for the next line
		c->cd(7);
		float texpos2 = 0.9; float texshift2 = 0.08;
		tex = new TLatex(0.0, texpos2, Form("Geometric Losses"));
		tex->Draw(); texpos2 -= texshift2; //draw this line and shift our position down to be ready for the next line
		//tex = new TLatex(0.0, texpos2, Form("Note the second number is:"));
		//tex->Draw(); texpos2 -= texshift2; //draw this line and shift our position down to be ready for the next line
		//tex = new TLatex(0.0, texpos2, Form("Lost Photons/Surviving Photons(pre-culling)"));
		//tex->Draw(); texpos2 -= texshift2; //draw this line and shift our position down to be ready for the next line
		//float RatIFC = IFCLoss / PreFieldCageSurvival;
		//float RatOFC = OFCLoss / PreFieldCageSurvival;
		//float RatPrism = PrismLoss / PrePrismSurvival;
		//float RatGeoD[10] = {};
	   // float RatAbsD[10] = {};
		tex = new TLatex(0.0, texpos2, Form("Inner Field Cage Loss=%d", IFCLoss));
		tex->Draw(); texpos2 -= texshift2; //draw this line and shift our position down to be ready for the next line
		tex = new TLatex(0.0, texpos2, Form("Outer Field Cage Loss=%d", OFCLoss));
		tex->Draw(); texpos2 -= texshift2; //draw this line and shift our position down to be ready for the next line
		int TotalGeoLoss = IFCLoss + OFCLoss;
		tex = new TLatex(0.0, texpos2, Form("Total Geometric Loss=%d", TotalGeoLoss));
		tex->Draw(); //texpos -= texshift; //draw this line and shift our position down to be ready for the next line
		c->cd(8);
		float texpos3 = 0.9; float texshift3 = 0.08;
		for (int j = 0; j < nDiffusers; j++) {
			//RatGeoD[j] = DiffuserGeometric[j] / (PreDiffuserSurvival[j]+1);
			//RatAbsD[j] = DiffuserAbsorbtion[j] / (PreDiffuserSurvival[j] + 1);
			tex = new TLatex(0.0, texpos3, Form("Absorption Loss from Diffuser %d =%d", j + 1, DiffuserAbsorbtion[j]));
			tex->Draw(); texpos3 -= texshift3; //draw this line and shift our position down to be ready for the next line
			tex = new TLatex(0.0, texpos3, Form("Geometric Loss from Diffuser %d =%d", j + 1, DiffuserGeometric[j]));
			tex->Draw(); texpos3 -= texshift3; //draw this line and shift our position down to be ready for the next line
		}
		tex = new TLatex(0.0, texpos3, Form("Geometric Loss from Prism =%d", PrismLoss));
		tex->Draw(); texpos3 -= texshift3; //draw this line and shift our position down to be ready for the next line
		//*
		c->cd(9)->SetLogy();
		PhotonSurvival->SetFillColor(38);
		PhotonSurvival->SetBarWidth(0.5);
		//PhotonSurvival->SetBarOffset(0.1);
		PhotonSurvival->Draw("bar");
		//*/
		//
		c->SaveAs(Form("%s.pdf", canvasname));
		fclose(pFile);
		
	    int xp, yp, zp;
		hPhotonAtSurface->GetBinXYZ(hPhotonAtSurface->GetMaximumBin(), xp, yp, zp);
	    std::vector<float> vecD;
		vecD.push_back(ratio_ave);
		vecD.push_back(light_fraction);
		vecD.push_back(xp);
		vecD.push_back(yp);
		return vecD;;
	
}

//****Support Functions:
//and the supporting functions:

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
const TVector3 PrismFacet(const TVector3 photon_direction, float IndexBefore, float IndexAfter, int numberprisms, const TVector3 Prism_NVector){
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

