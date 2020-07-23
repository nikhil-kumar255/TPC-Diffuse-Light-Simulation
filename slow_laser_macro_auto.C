const TVector3 DrawRandomRayFrom(TH2F* hXY, float d); //draw a ray (angles right, length arbitrary) from a distribution that makes an intensity map hXY a distance d from a source.
const TVector3 DiffusePhoton(const TVector3 photon_direction, TF1* angleDist); //modify a vector by two orthogonal random draws from the angle (deg) distribution.
TH2F* LoadLaserProfileFromFile(const char *profname,const char*tuplefile,const char *tuplename, float length);//build a 2D intensity profile from a file containing a 1D ntuple, and name the resulting profile 'profname'.  The total number of entries span a length of 'length' in cm.

																			   
																			   /*

bool HitsCylinder(const TVector3 photon_direction, const TVector3 photon_origin, float radius, float zf ); //return true if a photon traverses a cylinder  of radius r between its starting point and zf
bool HitsIFC(const TVector3 photon_direction,const TVector3 photon_origin){ return HitsCylinder(photon_direction,photon_origin, 17.25,100.0);};
bool HitsOFC(const TVector3 photon_direction,const TVector3 photon_origin){ return HitsCylinder(photon_direction,photon_origin, 80.0,100.0);};
TH2F* LoadLaserProfileFromGaussian(const char *profname, float sigma, float length); //build a 2D intensity profile from a 1D gaussian with sigma 'sigma'cm spanning a length 'length' in cm.
//*/

//and the main feature:  set up n lasers at a definted tilt angle and a specified itnensity on a card at a distace dist away.
//                       run that laser through nDiffusers with a transmission of diffTrans each, and the given profile
//                       then run it through n prisms, with a defined transmission, index of refraction, and wedge angle.
//this will produce a new canvas called 'canvasname' (and we can rig it to save that canvas, too, by adding a line at the end).
//todo:  could also have it return an array of useful parameters if you want.
void SimulateLasers(int nLasers, float laser_tilt_angle,TH2F *hIntensity, float dist,
		      int nDiffusers,float diffTransmission, TF1 *diffAngleProfile,
		    int nPrisms, float prismTransmission, float prismIndex, float prismAngle,
		    const char* canvasname);


void slow_laser_macro_auto() {

  const int nLasers = 1;
  int nDiffuser = 0;
  float laserthetaparameter=0.01;
  float laser_tilt_angle =  10;
  //set up the diffuser parameters:
  //null angle diffuser returns almost exactly what we take in:
  TF1 *fNullAngle=new TF1("fThorAngle","[0]*exp(-0.5*(x/[1])**2)",-30,30);
  fNullAngle->SetParameters(1,0.001); //test with very narrow gaussian to make sure we get back what we put in.

  //thorlabs diffuser with a by-hand fit to the model on their page
  TF1* fThorAngle = new TF1("fThorAngle", "([0]**2/((x-[1])**2+[0]**2))", -30, 30);
  fThorAngle->SetParameters(7.5, 0);//thorlabs model is 7.5,0
  fThorAngle->SetTitle("diffuser angular distribution;angle (deg);arb. prob.");
  float thorTransPercentile= 68;// note 68.35 for 10-220, 57.6 for 10-120
	68;//percent.

  TH2F* hIntensity;
  //load Bob's sculpted-tip laser profile from file:
  hIntensity=LoadLaserProfileFromFile("SculptedTip","ntuple.root","ntuple",14.0);
  float dist=10.0;//
  //load a gaussian laser profile:
  //hIntensity=LoadLaserProfileFromGaussian("PureGaussianSig7",7.0,14.0);

  //now you can make loops that run this multiple times, if you like.
  SimulateLasers(12,10,hIntensity,10.0,
		 4,thorTransPercentile,fThorAngle,
		 0,100,0,0,
		 "canvas");
   return;
};

void SimulateLasers(int nLasers, float laser_tilt_angle,TH2F *hIntensity, float dist,
		      int nDiffusers,float diffTransmission, TF1 *diffAngleProfile,
		    int nPrisms, float prismTransmission, float prismIndex, float prismAngle,
		    const char* canvasname){   
  //define the basis vectors local to the laser:
  TVector3 laser_nominal[nLasers];
  TVector3 laser_transverse[nLasers];

  //define the position of the lasers:
  TVector3 laser_position[nLasers];
  float laser_position_angle0 = 0;

  float angle_increment = 2 * TMath::Pi() / nLasers; // assume they're equally spaced
  float laser_tilt = laser_tilt_angle * TMath::Pi() / 180; //and have a common tilt, leaning outward in the rz plane,

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
  TH2F* hPhotonAtSurface = new TH2F("hPhotonAtSurface", "Photon Position AtSurface;x(cm);y(cm)",80, -100, 100, 80, -100, 100);
  TH2F* hPhotonAngle = new TH2F("hPhotonAngle", "Photon Angle;#theta (x);#phi (y)", 50, 0, 2, 50, 0, 6.5);
  TH1F* hPreDiffusionAngle = new TH1F("hPreDiffusionAngle", "Photon Y Angle Before Diffusion;#theta (deg)", 50, -45, 45);
  TH1F* hPostDiffusionAngle = new TH1F("hPostDiffusionAngle", "Photon Y Angle After Diffusion;#theta (deg)", 50, -45, 45);
  TH2F* hPhotonDirection = new TH2F("hPhotonDirection", "hPhotonDirection; (x); (y)", 50, -2, 2, 50, -2, 2);
  TH1F* hPhoton = new TH1F("hPhoton", "hPhoton", 100, 0, TMath::Pi());
  //*//Initialize File
  FILE* pFile;
  int n;
  char name[100];
  pFile = fopen("test_trace.csv", "w");
  //*//
  int nPhotons = 1000000;

  
  for (int i = 0; i < nPhotons; i++) {
    int L = i % nLasers;

    //draw a photon direction from our source distribution:
    TVector3 local_direction=DrawRandomRayFrom(hIntensity,dist); //d=10.0cm 
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

    //****Diffuser Section
    //the diffuser adds an additional random set of rotations w.r.t. the photon direction
    //rather than rotate theta and phi around the photon axis, we get two
    // orthogonal vectors relative to the photon direction, and do a 'rotate
    // about x axis' followed by a 'rotate about y axis' in that frame:

    //angle before the diffuser, for our records.  by fiat, let's always measure the 1D angular spread as the spread in the XZ plane:  y=0;
    TVector3 xz_dir_pre_diffusion=photon_direction;
    xz_dir_pre_diffusion.SetY(0);
    hPreDiffusionAngle->Fill(xz_dir_pre_diffusion.Theta()*180/TMath::Pi());
		
		
    bool islost = false;
    for (int i = 0; i < nDiffusers; i++) {
      photon_direction = DiffusePhoton(photon_direction, diffAngleProfile);
      float rand1 = rand() % 100; //random
      if (rand1 > diffTransmission) { //out of probability range, photon is absorbed}
	islost = true;
	break; //skip out of this loop
      }
    }

    //for now, let's always measure the 1D angular spread as the spread in the XZ plane:  y=0;
    TVector3 xz_dir_post_diffusion=photon_direction;
    xz_dir_post_diffusion.SetY(0);
    hPostDiffusionAngle->Fill(xz_dir_post_diffusion.Theta()*180/TMath::Pi());
   if (islost) continue; //skip to the next photon

   //end of diffuser

   //****Prism Section
   int nprisms = 1;
   float nIndex1 = 1; // index before prism
   float nIndex2 = 1.5; // index inside prism
   float nIndex3 = 1.2; // index afterprism
   float prism_angle = 45 * TMath::Pi() / 180;
   // float theta1 = atan2;
   //if (nprisms == 1) {
	   //float theta2 = asin(nIndex1*sin(theta1)/nIndex2);


   //}
   //something about how many prisms (generally zero or one?)
   //and what to do with the angles and indices given...
   
   //end of prism


   //****Collisions Section
   //if (HitsIFC(photon_direction,photon_position)) continue;
   //if (HitsOFC(photon_direction,photon_position)) continue;

   
    //find the position where this photon intersects the CM:
    float denominator = photon_direction.Dot(surface_normal);
    if (denominator < 0.001) { //parallel to surface.  no solution}
      continue; //skip to the next iteration through the loop
    }


    TVector3 intersect = (surface_normal.Dot(surface_point - photon_position) / denominator) * photon_direction + photon_position;


    hPhotonAtSurface->Fill(intersect.X(), intersect.Y());

  }

  //hPhotonDirection->Draw("colz");
  //hPhotonAngle->Draw("colz");
  //hPhotonAtSurface->Draw("surf1");
  // Reference h2->Fill(x_i, y_i, hIntensity->GetBinContent(hIntensity->FindBin(x_0_i)) * hIntensity->GetBinContent(hIntensity->FindBin(y_0_i)));
  // rEFERENCE TH1F("hIntensity", "Intensity vs Position;position;intensity", NBins, 0, NBins - 1);
  //********** 
  //*//              Ideal gaussian Slice, angular distribution
  //TH1D* px = hPhotonAtSurface->ProjectionX("px", -100, 100); // where firstYbin = -100 and lastYbin = 100 //original attempt
  TH1D* px = new TH1D("px", "Position vs Intensity;Position;Intensity", 201, -100, 100);
  TH1D* pax = new TH1D("pax", "Angle vs Intensity;Angle;Intensity",20, - 30, 30);
  double xI, yI;
  double zI = 100;//distance of fiber to CM
  double yPORT = 40; //cm; trial method to create slice histogram 
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
  //px->Draw("hist");
  //*///**********


  //*****************************
  //Draw photons striking cm with all given previous conditions
  //*
  //hPhotonAtSurface->Draw("colz");
  //*/
  //*****************************

  //Draw a lot of useful histograms:
  TCanvas* c = new TCanvas(canvasname, canvasname, 600, 900);
  c->Divide(2, 3);
  int nCells = hPhotonAtSurface->GetNcells();//total number of bins in histogram
  TH1F* hIntensityProfile = new TH1F("hIntensityProfile", "Histogram of Intensity per bin;intensity;nbins", 100, 1, 5 * nPhotons / nCells);
  TH1F* hIntensityRadius = new TH1F("hIntensityRadius", "Mean Intensity vs Radius;radius;mean intensity", 100, 0, 100);
  TH1F* hIntensityTheta = new TH1F("hIntensityTheta", "Mean Intensity vs Angle;Theta;mean intensity", 50, - (1 / 2) * TMath::Pi(), (1 / 2) * TMath::Pi());
  TH2F *hRadiusProfile=new TH2F("hRadiusProfile","Histogram of Intensity vs Radius for all bins, radius,intensity",100,0,100,100,1,5*nPhotons/nCells);
  TH1F* hRadius = new TH1F("hRadius", "nbins vs radius, for normalization;radius;nbins", 100, 0, 100);
  TH1F* hTheta = new TH1F("hTheta", "nbins vs Angle(Theta), for normalization;Theta;nbins", 50, 0, (1/2)*TMath::Pi());

  //read the bins from photonAtSurface and histogram their contents to get the per-radius average intensities
  for (int i = 0; i < nCells; i++) {
    float content = hPhotonAtSurface->GetBinContent(i);
    int binx, biny, binz;
    hPhotonAtSurface->GetBinXYZ(i, binx, biny, binz);//get the per-axis bins for this global bin
    float xpos = hPhotonAtSurface->GetXaxis()->GetBinCenter(binx);
    float ypos = hPhotonAtSurface->GetYaxis()->GetBinCenter(biny);
    float zpos = 100;
    float radius = sqrt(xpos * xpos + ypos * ypos);
    float theta = atan2(radius, zpos);
    hIntensityProfile->Fill(content);
    hIntensityRadius->Fill(radius, content);
    hIntensityTheta->Fill(theta, content);
    hRadiusProfile->Fill(radius,content);
    hRadius->Fill(radius);
    hTheta->Fill(theta);
  }
  //average by dividing by the number of bins that apply to each radius:
  hIntensityRadius->Divide(hRadius);
  hIntensityTheta->Divide(hTheta);


  float light_fraction=hPhotonAtSurface->GetEntries()/nPhotons;
  int bounds[3][2]={{21,26},{45,55},{71,76}};
  float diffs[3];
  for (int i=0;i<3;i++) diffs[i]=bounds[i][1]-bounds[i][0];
  //calculate some means from the intensity vs radius (averaging over several bins):
  float central_ave=hIntensityRadius->Integral(bounds[1][0],bounds[1][1])/diffs[1];
  float low_ave=hIntensityRadius->Integral(bounds[0][0],bounds[0][1])/diffs[0];
  float high_ave=hIntensityRadius->Integral(bounds[2][0],bounds[2][1])/diffs[2];
  float ratio_ave=central_ave/low_ave; // how much brighter is central than low?
  float asymmetry=(low_ave-high_ave)/(low_ave+high_ave); //how much brighter is low than high?
  
  c->cd(1);
  diffAngleProfile->Draw();
  c->cd(2);
  pax->Draw("hist");
  c->cd(3);
  hPhotonAtSurface->Draw("colz");
  c->cd(4);

  
  float texpos=0.9;float texshift=0.08;
  TLatex *tex=new TLatex(0.0,texpos,Form("PhotonSource=%s",hIntensity->GetName()));
  tex->Draw();texpos-=texshift; //draw this line and shift our position down to be ready for the next line
  tex=new TLatex(0.0,texpos,Form("nDiffs=%d, Trans=%1.1f%%",nDiffusers,diffTransmission));
  tex->Draw();texpos-=texshift; //draw this line and shift our position down to be ready for the next line
  tex=new TLatex(0.0,texpos,Form("nPrisms=%d, Index=%1.2f, #Theta=%2.1f, Trans=%1.1f%%",nPrisms,prismIndex,prismAngle,prismTransmission));
  tex->Draw();texpos-=texshift; //draw this line and shift our position down to be ready for the next line
  tex=new TLatex(0.0,texpos,Form("%% Photons on CM: %1.1f",light_fraction*100));
  tex->Draw();texpos-=texshift;
  tex=new TLatex(0.0,texpos,Form("Center(%d-%d) to Inner (%d-%d) ratio: %1.3f",
			      bounds[1][0],bounds[1][1],bounds[0][0],bounds[0][1],ratio_ave));
  tex->Draw();texpos-=texshift;
  tex=new TLatex(0.0,texpos,Form("Inner (%d-%d) to Outer (%d-%d) asym: %1.3f",
			      bounds[0][0],bounds[0][1],bounds[2][0],bounds[2][1],asymmetry));
  tex->Draw();texpos-=texshift;
  if (low_ave/high_ave>1.1){
    tex=new TLatex(0.0,texpos,Form("Inner %1.3f%% brighter than Outer!  Angle out more",(low_ave/high_ave-1)*100));
  }else   if (high_ave/low_ave>1.1){
    tex=new TLatex(0.0,texpos,Form("Outer %1.3f%% brighter than Inner!  Angle in more",(high_ave/low_ave-1)*100));
  } else{
    tex=new TLatex(0.0,texpos,Form("Outer and Inner approximately balanced."));
  }
  tex->Draw();texpos-=texshift;

  c->cd(5);
  hIntensityRadius->Draw("hist");
  c->cd(6);
  hIntensityTheta->Draw("hist");
  //
  //c->SaveAs(Form("%s.pdf",canvasname);

  return;
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
  TVector3 output_direction=photon_direction;
  
  TVector3 localA = output_direction.Orthogonal();//we have no intrinsic guarantee /which/ orthogonal direction this is
  TVector3 localB = output_direction.Cross(localA);//but we know that p x A will be orthogonal to both p and A, so these suffice.
  output_direction.Rotate(TMath::Pi() / 180 * angleDist->GetRandom(), localA);//convert to radians and rotate around A
  output_direction.Rotate(TMath::Pi() / 180 * angleDist->GetRandom(), localB);//then convert to radians and rotate  around B.
  return output_direction;
}
/*
bool HitsCylinder(const TVector3 photon_direction, const TVector3 photon_position, float radius, float zf ){
  //return true if a photon traverses a cylinder of radius r between its starting point and zf

  float vx = photon_direction.x();
  float vy = photon_direction.y();
  float x0 = photon_position.x();
  float y0 = photon_position.y();
  float vz = photon_direction.z();
  float z0 = photon_position.z();
  float R2=radius*radius;
		
  //Generalized Parametersfor collision with cylinder of radius R:
  //from quadratic formula solutions of when a vector intersects a circle:
  float a = vx*vx+vy*vy;
  float b = 2*(vx*x0+vy*y0);
  float c = x0*x0+y0*y0-R2;

  float rootterm=b*b-4*a*c;

  //if a==0 then we are parallel and will have no solutions.
  //if the rootterm is negative, we will have no real roots -- we are outside the cylinder and pointing skew to the cylinder such that we never cross.
  if (rootterm >= 0 && a > 0) {
    //Find the (up to) two points where we collide with the cylinder:
    float sqrtterm=sqrt(rootterm);
    float t1 = (-b+sqrtterm)/(2*a);
    float t2 = (-b-sqrtterm)/(2*a);
    float z1 = z0 + vz * t1;
    float z2 = z0 + vz * t2;
    if ((z1 > z0&& z1 < zf) || (z2 > z0&& z2 < zf) ) return true; //collision occurs in at least one case    
  }
  //if we got here, there was no collision.
  return false;
}
//*/
//*
TH2F* LoadLaserProfileFromFile(const char *profname,const char*tuplefile,const char *tuplename, float length){
  TFile* f = TFile::Open("ntuple.root");
  TNtuple* ntuple = (TNtuple*)(f->Get("ntuple"));
  float NBins = ntuple->GetEntries();//540;// Number of Bins in Bob's Histogram
  TH1F* hIntensity = new TH1F("hIntensity", "Intensity vs Position;position;intensity", NBins, 0, NBins);
  ntuple->Draw("position>>hIntensity", "intensity","goff");

  float Full_Length = 14;//cm; length of Bob's histogram
  float Fiber_Scale = Full_Length / NBins; // Length of one bin in cm

  TH2F* hist = new TH2F(Form("h%s",profname), Form("2D Output from %s;x position(cm);y position(cm)",profname), NBins, -Full_Length / 2, Full_Length / 2, NBins, -Full_Length / 2, Full_Length / 2);
  for (int i = 0; i < NBins; i++) {
    float x_i = (-Full_Length / 2) + i * Fiber_Scale;
    for (int j = 0; j < NBins; j++) {
      float y_i = (-Full_Length / 2) + j * Fiber_Scale;
      hist->Fill(x_i, y_i, hIntensity->GetBinContent(hIntensity->FindBin(i)) * hIntensity->GetBinContent(hIntensity->FindBin(j)));
    }
  }
  return hist;
}
//*/

/*
TH2F* LoadLaserProfileFromGaussian(const char *profname, float sigma, float length){
  gRandom = new TRandom3();
  double sigma_ideal = sigma;
  double sigma_y = sigma_ideal;
  double sigma_x = sigma_ideal;
  
  float NBins = 100;
  float Full_Length = length;//cm; length of histogram
  float Fiber_Scale = Full_Length / NBins; // Length of one bin in cm

  
  TH2F* hist = new TH2F(Form("h%s",profname), Form("2D Output from %s;x position(cm);y position(cm)",profname), NBins, -Full_Length / 2, Full_Length / 2, NBins, -Full_Length / 2, Full_Length / 2);
 

  for (int i = 0; i < NBins; i++) {
    float x_i = (-Full_Length / 2) + i * Fiber_Scale;
    for (int j = 0; j < NBins; j++) {
      float y_i = (-Full_Length / 2) + j * Fiber_Scale;
      hist->Fill(x_i, y_i, exp(-2 * ((x_i * x_i) / (sigma_x * sigma_x) + (y_i * y_i) / (sigma_y * sigma_y))));
    }
  }
  return hist;
}
//*/
