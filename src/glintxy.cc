// August 19, 2024: glintxy.cpp
// Version of glint-finding code that uses its own library and
// does not rely on solarsyst_dyn_geo01.


#include "glintlib.h"
#include "cmath"

      
static void show_usage()
{
  cerr << "Usage: glintxy -dets detfile -minpoints minpoints -maxGCR maximum GRC -maxrange max range/ \n";
  cerr << "-freq_downscale frequency downscale -freq_upscale frequency upscale/ \n";
  cerr << "-max_phase_err maximum phase error -trailout output trail file -trail2det output trail2det file\n";
  cerr << "\n    or, at minimum:\n";
  cerr << "glintxy -dets detfile\n";
  cerr << "Note well that the minimum invocation will leave a bunch of things\n";
  cerr << "set to defaults that may not be what you want.\n";
}


int main(int argc, char *argv[])
{
  vector <point3d_index> detvec = {};
  point3d_index onedet = point3d_index(0.0l,0.0l,0.0l,0);
  long detct=0;
  long detnum=0;
  vector <glint_trail> trailvec;
  vector <longpair> trail2det;
  ofstream outstream1;
  FindGlintsConfig config;  
  int status = 0;
  long i = 0;
  string indetfile;
  string trailfile = "glint_trail_test.csv";
  string trail2detfile = "glint_trail2det_test.csv";
  double x,y,mag;

  if(argc<3)
    {
      show_usage();
      return(1);
    }
  
  i=1;
  while(i<argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-detfile" || string(argv[i]) == "-dets" || string(argv[i]) == "-det" || string(argv[i]) == "--dets" || string(argv[i]) == "--det" || string(argv[i]) == "--detection" || string(argv[i]) == "--detections") {
      if(i+1 < argc) {
	//There is still something to read;
	indetfile=argv[++i];
	i++;
      }
      else {
	cerr << "Detection file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-minpoints" || string(argv[i]) == "-minpts") {
      if(i+1 < argc) {
	//There is still something to read;
        config.minpoints=stol(argv[++i]);
	i++;
      } else {
	cerr << "Maximum inter-image time interval\nkeyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-maxGCR" || string(argv[i]) == "-maxgcr" ) {
      if(i+1 < argc) {
	//There is still something to read;
        config.maxgcr=stod(argv[++i]);
	i++;
     }
      else {
	cerr << "Output maximum Great Circle Residual\nkeyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-maxrange") {
      if(i+1 < argc) {
	//There is still something to read;
        config.maxrange=stod(argv[++i]);
	i++;	
      } else {
	cerr << "Maximum glint range keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-freq_downscale") {
      if(i+1 < argc) {
	//There is still something to read;
        config.freq_downscale=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Frequency downcaling keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-freq_upscale") {
      if(i+1 < argc) {
	//There is still something to read;
        config.freq_upscale=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Frequency upscaling keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    }  else if(string(argv[i]) == "-max_phase_err") {
      if(i+1 < argc) {
	//There is still something to read;
        config.max_phase_err=stod(argv[++i]);
	i++;
      }
      else {
	cerr << "Maximum phase error keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-trailout") {
      if(i+1 < argc) {
	//There is still something to read;
	trailfile=argv[++i];
	i++;
      }
      else {
	cerr << "Output trail file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else if(string(argv[i]) == "-trail2det") {
      if(i+1 < argc) {
	//There is still something to read;
	trail2detfile=argv[++i];
	i++;
      }
      else {
	cerr << "Output trail2det file keyword supplied with no corresponding argument\n";
	show_usage();
	return(1);
      }
    } else {
      cerr << "Warning: unrecognized keyword " << argv[i] <<"\n";
      i++;
    }
  }

  if(indetfile.size()<=0) {
    cerr << "Please supply an input detection file:\n\n";
    show_usage();
    return(1);
  }

  // Read detection file
  ifstream instream1 {indetfile};
  detct=0;
  while(!instream1.fail() && !instream1.bad() && !instream1.eof()) {
    instream1 >> x >> y >> mag;
    if(!instream1.fail() && !instream1.bad() && !instream1.eof()) {
      onedet = point3d_index(x,y,mag,detct);
      detvec.push_back(onedet);
      detct++;
    }
  }
  instream1.close();
  detnum = detvec.size();
  cout << "Read " << detnum << " lines from input detection file " << indetfile << "\n";
  
  status = find_glints_xypix(detvec, config, trailvec, trail2det);
  if(status!=0) {
    cerr << "ERROR: find_glints_xypix failed with status " << status << "\n";
    return(status);
  }
  
  // Write output glint trail catalog
  cout << "Writing output glint trail catalog " << trailfile << " with " << trailvec.size() << " lines\n";
  outstream1.open(trailfile);
  outstream1 << "#glintnum,cenx,ceny,trail_len,trail_PA,GCR,eq_rms,magmean,magrms,stepsize,qc1,npt,flashnum\n";
  for(i=0;i<long(trailvec.size());i++) {
    outstream1 << fixed << setprecision(7) << i << "," << trailvec[i].x << "," << trailvec[i].y << ",";
    outstream1 << fixed << setprecision(3) << trailvec[i].length << ",";
    outstream1 << fixed << setprecision(2) << trailvec[i].PA << ",";
    outstream1 << fixed << setprecision(3) << trailvec[i].linrms << "," << trailvec[i].eqrms << ",";
    outstream1 << fixed << setprecision(2) << trailvec[i].magmean << "," << trailvec[i].magrms << ",";
    outstream1 << fixed << setprecision(3) << trailvec[i].stepsize << "," << trailvec[i].qc1 << ",";
    outstream1 << trailvec[i].npt << "," << trailvec[i].flashnum << "\n";
  }
  outstream1.close();

   // Write trail2det file
  cout << "Writing trail2det file " << trail2detfile << " with " << trail2det.size() << " lines\n";
  outstream1.open(trail2detfile);
  outstream1 << "#trail_ID,detnum\n";
  for(i=0;i<long(trail2det.size());i++) {
    outstream1 << trail2det[i].i1 << "," << trail2det[i].i2 << "\n"; 
  }
  outstream1.close();

  return(0);
}

