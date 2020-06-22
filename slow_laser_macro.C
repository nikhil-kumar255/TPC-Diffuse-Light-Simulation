void slow_laser_macro() {

	const int nLasers = 12;
	//incorporate Bob's laser profile
	float NBins = 540;// Number of Bins in Bob's Histogram
	float Full_Length = 14;//cm; length of Bob's histogram
	float Fiber_Scale = Full_Length / NBins; // Length of one bin in cm
	TFile* f = TFile::Open("ntuple.root");
	TNtuple* ntuple = (TNtuple*)(f->Get("ntuple"));
	TH1F* hIntensity = new TH1F("hIntensity", "Intensity vs Position;position;intensity", NBins, 0, NBins - 1);
	ntuple->Draw("position>>hIntensity", "intensity");
	TH2F* h2 = new TH2F("h2", "Filled Test Plane;x position(cm);y position(cm)", NBins, -Full_Length / 2, Full_Length / 2, NBins, -Full_Length / 2, Full_Length / 2);
	// section to setup "ideal" Gaussian Distribution
	TH2F* h4 = new TH2F("h4", "Ideal Gaussian", NBins, -Full_Length / 2, Full_Length / 2, NBins, -Full_Length / 2, Full_Length / 2);
	gRandom = new TRandom3();
	double sigma_ideal = 7 ;
	double sigma_y = sigma_ideal;
	double sigma_x = sigma_ideal;
	//
	for (int i = 0; i < NBins; i++) {
		float x_i = (-Full_Length / 2) + i * Fiber_Scale;
		float x_0_i = i;
		for (int j = 0; j < NBins; j++) {
			float y_0_i = j;
			float y_i = (-Full_Length / 2) + j * Fiber_Scale;
			h2->Fill(x_i, y_i, hIntensity->GetBinContent(hIntensity->FindBin(x_0_i)) * hIntensity->GetBinContent(hIntensity->FindBin(y_0_i)));// Fill Custom Distribution
			h4->Fill(x_i, y_i, exp(-2 * ((x_i * x_i) / (sigma_x * sigma_x) + (y_i * y_i) / (sigma_y * sigma_y))));//Fill Ideal Distribution
		}
	}
		double x,y,z;
		double x_mean = h2->GetMean(1);//Get x mean of Bob's Profile
		double y_mean = h2->GetMean(2);//Get y mean of Bob's Profile
	//define the basis vectors local to the laser:
	TVector3 laser_nominal[nLasers];
	TVector3 laser_transverse[nLasers];

	//define the position of the lasers:
	TVector3 laser_position[nLasers];//(0,40.0,0); //position laser emitter on the y axis.
	float laser_position_angle0 = 0;

	float angle_increment = 2 * TMath::Pi() / nLasers; // assume they're equally spaced
	float laser_tilt = 0 * TMath::Pi() / 180; //and have a common tilt, leaning outward in the rz plane,

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
	}




	//define the surface we are illuminating:
	TVector3 surface_point(0, 0, 100.0);
	TVector3 surface_normal(0, 0, 1);



	//set up the distributions for theta and phi angles relative to laser_nominal:
	//TF1 *laser_theta=new TF1("laser_theta","exp([0]*x)",0,TMath::Pi()/2.0);
	//TF1 *laser_theta=new TF1("laser_test","exp(-[0]*(x-[1])^2)",0,TMath::Pi()/2.0);
	TF1* laser_phi = new TF1("laser_phi", "1", 0, 2 * TMath::Pi());
	//laser_theta->SetParameters(100,0.5);
	//laser_theta->Draw();

	//laser intensity as a function of distance from the beam center:
	float w0 = 0.1;
	float zsample = 1.0;//all in cm.
	float lambda = 600 * 1e-7;//10^-7 converts from nm in cm.
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
	// denom is placeolder for this case
	TF1* laser_theta = new TF1("laser_theta", "exp(-x/[0])", 0, TMath::Pi() / 2);
		laser_theta->SetParameter(0, 1);
	//return;

	  //TODO: find a realistic waist
	  //TODO: define a lens?
	  //TODO: add a legend to the plots so they show the variables used!
	  //could weight by intensity if we wanted.
	  //could invert the calculation to get the beam position as a function of point of intersection


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
	TH2F* hPhotonAtSurface = new TH2F("hPhotonAtSurface", "hPhotonAtSurface",80, -100, 100, 80, -100, 100);
	TH2F* hPhotonAngle = new TH2F("hPhotonAngle", "hPhotonAngle;#theta (x);#phi (y)", 50, 0, 2, 50, 0, 6.5);
	TH2F* hPhotonDirection = new TH2F("hPhotonDirection", "hPhotonDirection; (x); (y)", 50, -2, 2, 50, -2, 2);
	TH1F* hPhoton = new TH1F("hPhoton", "hPhoton", 100, 0, TMath::Pi());

	int nPhotons = 4000000;
	for (int i = 0; i < nPhotons; i++) {





		int L = i % nLasers;


		//hPhoton->Fill(laser_shape->GetRandom(),1.0/nPhotons);
		//double photon_theta = laser_theta->GetRandom();
		
		//float photon_phi = laser_phi->GetRandom();

		//****************Custom Laser fiber section
		///*
		z = 10;// z=10 cm . from Bob's setup
		h2->GetRandom2(x, y);
		x = x - x_mean;
		y = y - y_mean;
		double photon_phi= atan2(y,x);
		float photon_r = 0; // placeholder for bob
        double photon_theta = atan2(sqrt(x * x + y * y), z); 
		//*/
	    //***************End Custom laser Fiber section
		//Ideal Laser Fiber section
		/*
		z = 10;
		h4->GetRandom2(x, y);
		double photon_phi = atan2(y, x);
		float photon_r = 0; // placeholder for bob
		double photon_theta = atan2(sqrt(x * x + y * y), z);
		//*/
		//

		//float photon_r = laser_r->GetRandom();
		//hPhotonAngle->Fill(photon_theta, photon_phi);//fill x and y.  annoying that it's backward from the arguments earlier, but what're you gonna do? ;)

		//************Lens Section
		float focal_length = 1.0;//focal length in centimeters //float photon_theta=something to do with photon_r
		//photon_theta = photon_r / focal_length + photon_theta; //thin lens equation angle change. f is focal length for DIVERGING lens
		//************End Lens Section
		TVector3 photon_direction = laser_nominal[L];
		//for Nikhil:  photon_direction.Rotate(something to do with photon_theta)
		photon_direction.Rotate(photon_theta, laser_transverse[L]);
		TVector3 photon_position = laser_position[L];
		TVector3 photon_offset = laser_transverse[L] * photon_r;
		photon_offset.Rotate(photon_phi, laser_nominal[L]);
		photon_direction.Rotate(photon_phi, laser_nominal[L]);
		photon_position = photon_position + photon_offset;
		//hPhotonAtSurface->Fill(photon_direction.X(),photon_direction.Y());
		//hPhotonAtSurface->Fill(photon_position.X(),photon_position.Y());
		//continue;
		//double z=photon_position.z();
		//double photon_theta = atan2(sqrt(x*x+y*y),z);
		float denominator = photon_direction.Dot(surface_normal);
		if (denominator < 0.001) { //parallel to surface.  no solution}
			continue; //skip to the next iteration through the loop
		}
		TVector3 intersect = (surface_normal.Dot(surface_point - photon_position) / denominator) * photon_direction + photon_position;
		//if (intersect.X() * intersect.X() + intersect.Y() * intersect.Y() < 297.5625) { continue; }  // exclude solutions that end in center\
		//if (laser_position[L].X() < 0){if laser_position[L].Y() < 0} { continue; }  // exclude solutions that end in center
		//TVector3 photon_trajectory = photon_direction * 120;//create a photon track
		////TF2 * Field_Cage_Limits = new TF2("Field_Cage_Limits", "x^2+y^2", 0, 17.25, 0, 17.25);
		//TF3 * Field_Cage_Limits = new TF3("Field_Cage_Limits", "x^2+y^2", 0, 17.25, 0, 17.25, 0, 100);
		//if (photon_trajectory == Field_Cage_Limits) { continue; }
		double R = 17.25;// TPC Inner Field Cage Radius
		double R_out = 80;// TPC outer Field Cage, Inner surface Radius
		double v_x = photon_direction.x();
		double v_y = photon_direction.y();
		double x_0 = photon_position.x();
		double y_0 = photon_position.y();
		double v_z = photon_direction.z();
		double z_0 = photon_position.z();
		/*double b = -16 * v_x * v_x * y_0 * y_0 + 4761 * v_x * v_x + 32 * v_x * v_y * x_0 * y_0 - 16 * v_y * v_y * x_0 * x_0 + 4761*v_y * v_y;*/
        /*double c = v_x * v_x + v_y * v_y;*/
		//Generalized Parameters b & c for TPC of variable inner field cage radius 
		double b = (-2 * v_x * x_0 - 2 * v_y * y_0) * (-2 * v_x * x_0 - 2 * v_y * y_0) - 4 * (-v_x * v_x - v_y * v_y) *(R * R - x_0 * x_0 - y_0 * y_0); 
		double b_out = (-2 * v_x * x_0 - 2 * v_y * y_0) * (-2 * v_x * x_0 - 2 * v_y * y_0) - 4 * (-v_x * v_x - v_y * v_y) * (R_out * R_out - x_0 * x_0 - y_0 * y_0);
		double c = 2*(v_x*v_x + v_y*v_y);
		if (b >= 0 && c > 0) {
			/*double t_1 = (-sqrt(b) - 4*v_x * x_0 - 4*v_y * y_0) / (4*c);*/
			/*double t_2 = (sqrt(b) - 4*v_x * x_0 - 4*v_y * y_0) / (4 * c);*/
			//Generalized Parameter t for TPC of variable inner field cage radius 
			double t_1 = (-sqrt(b) - 2 * v_x * x_0 - 2 * v_y * y_0) / (c);
			double t_2 = (sqrt(b) - 2 * v_x * x_0 - 2 * v_y * y_0) / (c);
			double z_1 = z_0 + v_z * t_1;
			double z_2 = z_0 + v_z * t_2;
			if ((z_1 > z_0&& z_1 < intersect.z()) || (z_2 > z_0&& z_2 < intersect.z()) ) continue;
		
		}
		if (b_out >= 0 && c > 0) {
			double t_3 = (-sqrt(b_out) - 2 * v_x * x_0 - 2 * v_y * y_0) / (c);
			double t_4 = (sqrt(b_out) - 2 * v_x * x_0 - 2 * v_y * y_0) / (c);
			double z_3 = z_0 + v_z * t_3;
			double z_4 = z_0 + v_z * t_4;
	if ((z_3 > z_0&& z_3 < intersect.z()) || (z_4 > z_0&& z_4 < intersect.z()))continue;
		}

		//t = (-sqrt((-2 v_x x_0 - 2 v_y y_0)^2 - 4 (-v_x^2 - v_y^2) (R^2 - x_0^2 - y_0^2)) - 2 v_x x_0 - 2 v_y y_0)/(2 (v_x^2 + v_y^2))
	
	    //Double_t m = photon_direction.Mag();

		hPhotonAtSurface->Fill(intersect.X(), intersect.Y());
		
		/*



	  //hPhoton->Fill(laser_shape->GetRandom(),1.0/nPhotons);
	  //float photon_theta=laser_theta->GetRandom();
	  float photon_phi=laser_phi->GetRandom();
		float photon_r=laser_r->GetRandom();
		float  photon_theta=TMath::ATan2(photon_r,zsample);

	  TVector3 photon_direction=laser_nominal[L];
	  
	  photon_direction.Rotate(photon_theta,laser_transverse[L]);
	  photon_direction.Rotate(photon_phi,laser_nominal[L]);
	  hPhotonDirection->Fill(photon_direction.X(),photon_direction.Y());
	  float denominator=photon_direction.Dot(surface_normal);
	  if (abs(denominator)<0.001) { //parallel to surface.  no solution}
		continue; //skip to the next iteration through the loop
	  }
		TVector3 intersect=(surface_normal.Dot(surface_point-laser_position[L])/denominator)*photon_direction+laser_position[L];
	  hPhotonAtSurface->Fill(intersect.X(),intersect.Y());
		 */
	}
	//hPhotonDirection->Draw("colz");
	//hPhotonAngle->Draw("colz");
	//hPhotonAtSurface->Draw("surf1");
	// Reference h2->Fill(x_i, y_i, hIntensity->GetBinContent(hIntensity->FindBin(x_0_i)) * hIntensity->GetBinContent(hIntensity->FindBin(y_0_i)));
	// rEFERENCE TH1F("hIntensity", "Intensity vs Position;position;intensity", NBins, 0, NBins - 1);
	//********** 
	//              Ideal gaussian Slice, angular distribution
	//TH1D* px = hPhotonAtSurface->ProjectionX("px", -100, 100); // where firstYbin = -100 and lastYbin = 100 //******original attempt
	TH1D* px = new TH1D("px", "Position vs Intensity;Position;Intensity", 201, -100, 100);
	TH1D* pax = new TH1D("pax", "Angle vs Intensity;Angle;Intensity",30, - 30, 30);
	double xI, yI;
	double zI = 100;//distance of fiber to CM
	double yPORT = 60; //cm; trial method to create slice histogram 
	double R_out = 80;
	
	for (int i = 0; i < 2*R_out+1; i++) {
		
        float xI_0 = i- R_out;
		//hPhotonAtSurface->GetRandom2(xI, yI);
		double ThetaI = (180 / TMath::Pi())*atan2(xI_0 , zI);// in degrees when including 180/TMath::Pi()
		px->Fill(xI_0,hPhotonAtSurface->GetBinContent(hPhotonAtSurface->FindBin(xI_0,yPORT)));
		//remember, find bin gives bin number based on coordinates, and  get bin content gives the value stored in that bin	
		pax->Fill(ThetaI, px->GetBinContent(px->FindBin(xI_0)));
	}
	//double PhiI = atan2(yI, xI);
	pax->Draw("hist");
	//**********


	//*****************************
	//Draw photons striking cm with all given previous conditions
	/*
	hPhotonAtSurface->Draw("colz");
	//*/
	//*****************************

	//hPhoton->Draw("same");
	return;
}

/* things you can look up/do:
rebuild in matlab
figure out rotations needed.
get some testable cases to confirm
pick an interesting option or two:
gaussian from last time at (0,0,1) check it looks the same
gaussian from last time at 30 degrees radially outward.
figure out if there is a root package for lookign at where lines cross planes -- if not, write the equations for it.
*/