#include "glintlib.h"
#include "cmath"

#define DEBUG 0

void make_ivec(long nx, vector <int> &ivec)
{
  long i=0;
  ivec={};
  for(i=0;i<nx;i++) ivec.push_back(0);
}

// dmeanrms01: January 16, 2023
// Calculate and return the mean and RMS of
// a double-precision vector.
int dmeanrms01(const vector <double> &invec, double *mean, double *rms)
{
  int pnum=invec.size();
  int pct=0;
  
  if(pnum<=0) return(2);
  *mean = *rms = 0.0l;
  // Calculate the mean
  for(pct=0; pct<pnum; pct++) {
    *mean+=invec[pct];
  }
  *mean /= double(pnum);
  if(pnum==1) return(1);
  
  // Calculate the RMS
  for(pct=0; pct<pnum; pct++) {
    (*rms) += DSQUARE(invec[pct] - *mean);
  }
  *rms = sqrt((*rms)/double(pnum-1));
  return(0);
}

// linfituw01: January 20, 2022
// Simple and crude utility program, does an unweighted
// linear fit of the form y = mx * b, for m = slope, b = intercept
int linfituw01(const vector <double> &x, const vector <double> &y, double &slope, double &intercept)
{
  int i;
  int pointnum = x.size();
  double delta,xal,yal,xty,xsq,nsum;

  if(pointnum<=1) {
    cerr << "ERROR: linfituw01 CALLED WITH ONLY ONE POINT\n";
    return(1);
  }

  xal = yal = xty = xsq = nsum = 0.0;
  for(i=0;i<pointnum;i++) {
    xal += x[i];
    yal += y[i];
    xsq += x[i]*x[i];
    xty += x[i]*y[i];
    nsum += 1.0l;
  }
  delta = nsum*xsq - xal*xal;
  if(delta==0.0) {
    cerr << "ERROR: linfituw01 has non-finite slope\n";
    return(1);
  }
  intercept = (xsq*yal - xal*xty)/delta;
  slope = (nsum*xty - xal*yal)/delta;

  return(0);
}

long medindex(const vector <xy_index> &xyvec, int dim)
{
  vector <xy_index> xyv = xyvec; //Mutable copy of immutable input vector
  for(long unsigned int i=0; i<xyv.size(); i++) xyv[i].index=i; //Redefine indices
  long medpt = xyv.size()/2;
  if(dim%2==1) sort(xyv.begin(), xyv.end(), xyind_lower_x());
  else sort(xyv.begin(), xyv.end(), xyind_lower_y());
  return(xyv[medpt].index);
}

int splitxy(const vector <xy_index> &xyvec, int dim, long unsigned int splitpoint, vector <xy_index> &left, vector <xy_index> &right)
{
  long unsigned int i=0;
  double xval = xyvec[splitpoint].x;
  double yval = xyvec[splitpoint].y;
  
  if(dim%2==1) {
    // Split on x
    for(i=0 ; i<xyvec.size() ; i++) {
      if(i!=splitpoint && xyvec[i].x<=xval) {
	left.push_back(xyvec[i]);
      } else if(i!=splitpoint) {
	right.push_back(xyvec[i]);
      }
    }
  } else {
    // Split on y
    for(i=0 ; i<xyvec.size() ; i++) {
      if(i!=splitpoint && xyvec[i].y<=yval) {
	left.push_back(xyvec[i]);
      } else if(i!=splitpoint) right.push_back(xyvec[i]);
    }
  }
  return(0);
}

// kdtree01: November 11, 2021:
// Given an input root point, presumed to have null
// right and left branches, load the branches and then
// call kdtree01 on them recursively.
// NOTE THAT THIS IS FOR 2-D KD trees.
int kdtree01(const vector <xy_index> &xyvec, int dim, long unsigned int rootptxy, long unsigned int rootptkd, vector <kdpoint> &kdvec)
{
  int lmed=0;
  int rmed=0;
  int kdct = kdvec.size()-1;
  long leftrootkd=-1;
  long rightrootkd=-1;
  xy_index xyi = xy_index(0.0,0.0,0);
  kdpoint lp = kdpoint(xyi,-1,-1,0);
  kdpoint rp = kdpoint(xyi,-1,-1,0);
  kdpoint kdtest = kdpoint(xyi,-1,-1,0);
  vector <xy_index> leftvec = {};
  vector <xy_index> rightvec = {};

  splitxy(xyvec,dim,rootptxy,leftvec,rightvec);
  
  if(dim==1) dim=2;
  else dim=1;

  if(leftvec.size()==1) {
    // Left branch is just a single leaf
    lp = kdpoint(leftvec[0],-1,-1,dim);
    kdvec.push_back(lp);
    kdct++;
    kdvec[rootptkd].left = kdct;
    kdtest = kdvec[kdct];
  } else if(leftvec.size()<=0) {
    // There is no left branch
    kdvec[rootptkd].left = -1;
  }
  if(rightvec.size()==1) {
    // Right branch is just a single leaf
    rp = kdpoint(rightvec[0],-1,-1,dim);
    kdvec.push_back(rp);
    kdct++;
    kdvec[rootptkd].right = kdct;
    kdtest = kdvec[kdct];
  } else if(rightvec.size()<=0) {
    // There is no right branch
    kdvec[rootptkd].right = -1;
  }
   
 if(leftvec.size()>1) {
    lmed = medindex(leftvec,dim);
    lp = kdpoint(leftvec[lmed],-1,-1,dim);
    kdvec.push_back(lp);
    kdct++;
    kdvec[rootptkd].left = kdct;
    leftrootkd = kdct;
    kdtest = kdvec[kdct];
 }
 
  if(rightvec.size()>1) {
    rmed = medindex(rightvec,dim);
    rp = kdpoint(rightvec[rmed],-1,-1,dim);
    kdvec.push_back(rp);
    kdct++;
    kdvec[rootptkd].right = kdct;
    rightrootkd = kdct;
    kdtest = kdvec[kdct];
  }
  // I moved these down out of the above loops, because I thought
  // that otherwise, a bunch of stuff might get pushed down by the
  // left loop that the right loop didn't know about.
  if(leftvec.size()>1 && leftrootkd>=0) kdtree01(leftvec,dim,lmed,leftrootkd,kdvec);
  else if(leftvec.size()>1 && leftrootkd<0)
    {
      cerr << "Error, kdtree01 finds leftroot less than zero with leftvec.size() = " << leftvec.size() << "\n";
    }
  if(rightvec.size()>1 && rightrootkd>=0) kdtree01(rightvec,dim,rmed,rightrootkd,kdvec);
  else if(rightvec.size()>1 && rightrootkd<0)
    {
      cerr << "Error, kdtree01 finds rightroot less than zero with rightvec.size() = " << rightvec.size() << "\n";
    }

  return(0);
}

// kdrange01: November 15, 2021:
// Given a k-d tree vector kdvec created by kdtree01,
// perform a range-query about the point x,y. Returns
// a vector indexing all of the points in the input k-d tree
// that lie within the specified range of the input coordinates.
// Assumes that kdvec[0] is the root of the k-d tree.
// NOTE THAT THIS IS FOR 2-D KD trees.
int kdrange01(const vector <kdpoint> &kdvec,double x,double y,double range,vector <long> &indexvec)
{
  double rng2 = range*range;
  int notdone=1;
  int dim=1;
  long currentpoint=0;
  long leftpoint=0;
  long rightpoint=0;
  int goleft=0;
  int goright=0;
  double xdiff=0.0;
  double ydiff=0.0;
  vector <long> checkit={};
  long unsigned int checknum=0;

  while(notdone>0) {
    // Climb to the top of the k-d tree, keeping track
    // of potentially interesting unexplored branches
    // in the vector checkit.
    while(leftpoint>=0 || rightpoint>=0) {
      // Previous step did not end on a leaf.
      leftpoint = kdvec[currentpoint].left;
      rightpoint = kdvec[currentpoint].right;
      dim = kdvec[currentpoint].dim;
      if(dim%2==1) {
	xdiff = kdvec[currentpoint].point.x - x;
	// cout << "kxdr: " << kdvec[currentpoint].point.x << " " << x << " " << xdiff << " " << range << "\n";
	goright = (xdiff <= range); // possible hits lie to the left;
	goleft = (xdiff >= -range); // possible hits lie to the right;
	if(goleft && goright) {
	  // Current point might be within range.
	    ydiff = kdvec[currentpoint].point.y - y;
	    if(fabs(ydiff)<=range && (xdiff*xdiff + ydiff*ydiff)<=rng2) {
	      // Current point is within range. Add it to the output vector
	      indexvec.push_back(currentpoint);
	    }
	    if(leftpoint>=0) {
	       //Explore leftward first.
	      currentpoint = leftpoint;
	      if(rightpoint>=0) {
		// Rightward branch will also be explored later
		checknum++;
		if(checknum>checkit.size()) {
		  checkit.push_back(rightpoint);
		}
		else {
		  checkit[checknum-1] = rightpoint;
		}		
	      }
	    }
	    else if(rightpoint>=0) {
	      // Leftward branch is a dead end: explore rightward branch
	      currentpoint = rightpoint;
	    }
	}
	else if(goleft) {
	  // Current point cannot be in range, but points that
	  // are in range may lie along the left branch.
	  if(leftpoint>=0) {
	    currentpoint = leftpoint;
	  } else rightpoint=-1; // Dead end, make sure while-loop exits.
	} else if(goright) {
	  // Current point cannot be in range, but points that
	  // are in range may lie along the right branch.
	  if(rightpoint>=0) {
	    currentpoint = rightpoint;
	  } else leftpoint=-1;  // Dead end, make sure while-loop exits.
	} else {
	  // Program concluded it should go neither left nor right.
	  // The likely cause is that it encountered a NAN. Give up on this point.
	  leftpoint=rightpoint=-1;
	  cerr << "WARNING: ENCOUNTERED NAN CASE!\n";
	  cerr << "Input point " << x << ", " << y <<", target point " << kdvec[currentpoint].point.x << ", " << kdvec[currentpoint].point.y << ".\n";
	}
	// Close x-dimension case
      } else if(dim%2==0) {
	ydiff = kdvec[currentpoint].point.y - y;
	goright = (ydiff <= range); // possible hits lie to the left;
	goleft = (ydiff >= -range); // possible hits lie to the right;
	if(goleft && goright) {
	    // Current point might be within range.
	    xdiff = kdvec[currentpoint].point.x - x;
	    if(fabs(ydiff)<=range && (xdiff*xdiff + ydiff*ydiff)<=rng2) {
	      // Current point is within range. Add it to the output vector
	      indexvec.push_back(currentpoint);
	    }
	    if(leftpoint>=0) {
	       //Explore leftward first.
	      currentpoint = leftpoint;
	      if(rightpoint>=0) {
		// Rightward branch will also be explored later
		checknum++;
		if(checknum>checkit.size()) {
		  checkit.push_back(rightpoint);
		}
		else {
		  checkit[checknum-1] = rightpoint;
		}
	      }
	    } else if(rightpoint>=0) {
	      // Leftward branch is a dead end: explore rightward branch
	      currentpoint = rightpoint;
	    }
	}
	else if(goleft) {
	  // Current point cannot be in range, but points that
	  // are in range may lie along the left branch.
	  if(leftpoint>=0) {
	    currentpoint = leftpoint;
	  } else rightpoint = -1; // Dead end, make sure while loop exits.
	} else if(goright) {
	  // Current point cannot be in range, but points that
	  // are in range may lie along the right branch.
	  if(rightpoint>=0) {
	    currentpoint = rightpoint;
	  } else leftpoint=-1;  // Dead end, make sure while loop exits.
	} else {
	  // Program concluded it should go neither left nor right.
	  // The likely cause is that it encountered a NAN. Give up on this point.
	  leftpoint=rightpoint=-1;
	  cerr << "WARNING: ENCOUNTERED NAN CASE!\n";
	  cerr << "Input point " << x << ", " << y <<", target point " << kdvec[currentpoint].point.x << ", " << kdvec[currentpoint].point.y << ".\n";
	}
	// Note that we do not need to worry about the possiblity
	// that current point will get set to -1: i.e., we were
	// at a leaf or a one-sided branch. Such cases will
	// be caught by the while statement.
	// Close y-dimension case
      }
      // Close while-loop checking if we've hit a leaf.
    }
    // We have climbed up the tree to a leaf. Go backwards through
    // the checkit vector and see if there is anything to check.
    checknum=checkit.size();
    while(checknum>=1 && checkit[checknum-1]<0) checknum--;
    if(checknum<=0) {
      //There were no valid entries to check: we're done.
      notdone=0;
    } else {
      //Set currentpoint to the last valid entry in checkit
      currentpoint = checkit[checknum-1];
      //Mark this point as used.
      checkit[checknum-1]=-1;
      leftpoint=rightpoint=0;
    }
  }
  return(0);
}


// find_glints_xypix: July 12, 2024: 
// Like find_glints_radec, but takes as input x,y pixel positions rather
// than RA, Dec. The input vector is of type point3d_index, where x and 
// y are pixel coordinates and it is expected that z=magnitude,
// while index will be unchanged by the current
// function and retained to enable the calling function to map back to some
// external catalog. 
int find_glints_xypix(const vector <point3d_index> &detvec, FindGlintsConfig config, vector <glint_trail> &trailvec, vector <longpair> &trail2det)
{
  trailvec={};
  trail2det={};

  long detnum = detvec.size();
  if(detnum<config.minpoints) {
    // Too few detections for any valid glint trails
    cout << "Note: find_glints_radec called with only " << detnum << " input detections:\n";
    cout << "too few to find any valid trails\n";
    return(0); // This is not an error case, just a valid finding of zero glint trails.
  }
  
  long detct=0;
  xy_index xyind=xy_index(0.0, 0.0, 0);
  vector <xy_index> axyvec = {};
  double dist,pa;
  dist = pa = 0.0l;
  longpair onepair = longpair(0,0);
  vector <long_index> num_matches;
  long i,j;
  vector <xy_index> trailxyvec = {};
  vector <vector <long>> ppind = {};
  long glintpnum=0;
  long glintnum=0;
  double traillen, linrms, eqrms, magmean, magrms, stepsize, qc1, xcen, ycen;
  traillen = linrms = eqrms = magmean = magrms = stepsize = qc1 = xcen = ycen = 0.0l;
  long npt=0;
  long flashnum=0;
  glint_trail oneglint = glint_trail(xcen, ycen, traillen, pa, linrms, eqrms, magmean, magrms, stepsize, qc1, npt, flashnum);
  long longest_trail, trail_length;
  //int verbose = 1;
  double GCR = 0.0l;
  vector <glint_trail> alltrailvec;
  vector <vector <long>> indvecs;
  vector <int> detused;
  make_ivec(detnum, detused);
  double_index dindex = double_index(0l,0l);
  vector <double_index> dindvec;
  vector <long> anchor_det;
  int status=0;
  
  xyind=xy_index(0.0, 0.0, 0);
  axyvec = {};
  for(detct=0; detct<detnum; detct++) {
    xyind = xy_index(detvec[detct].x, detvec[detct].y,detct);
    // In the RA, Dec version we do a spherical projection here, but
    // it isn't necessary with x,y pixel coords.
    axyvec.push_back(xyind);
  }
  if(long(axyvec.size()) != detnum) {
    cerr << "ERROR: not all detections were successfully loaded\n";
    return(3);
  } else cout << "Successfully loaded all " << detnum << " input detections\n";
  
  // Create k-d tree of the projected detections.
  int dim=1;
  xy_index xyi = axyvec[0];
  kdpoint root = kdpoint(xyi,-1,-1,dim);
  kdpoint kdtest = kdpoint(xyi,-1,-1,dim);
  vector <kdpoint> kdvec ={};
  long medpt;
  medpt = medindex(axyvec,dim);
  root = kdpoint(axyvec[medpt],-1,-1,1);
  kdvec.push_back(root);
  kdtest=kdvec[0];
  kdtree01(axyvec,dim,medpt,0,kdvec);
  cout << "Loaded kd-tree with " << kdvec.size() << " points\n";
  
  // Loop over detections, looking for matches.
  num_matches={};
  for(detct=0 ; detct<detnum ; detct++) {
    vector <long> indexvec = {};
    if((isnormal(axyvec[detct].x) || axyvec[detct].x==0) && (isnormal(axyvec[detct].y) || axyvec[detct].y==0)) {
      kdrange01(kdvec,axyvec[detct].x,axyvec[detct].y,config.maxrange,indexvec);
    } else {
      cerr << "WARNING: detection " << detct << " is not normal: " << axyvec[detct].x << " " << axyvec[detct].y << "\n";
    }
    long matchnum=indexvec.size();
    if(matchnum>=config.minpoints) {
      // This detection could be the root of a glint trail.
      trailxyvec={};
      ppind={};
      for(i=0; i<matchnum; i++) { // Loop over the pair-partners of detection detct.
	j = kdvec[indexvec[i]].point.index;
	if(j!=detct) {
	  // j is the index of a pair-partner to detection detct.
	  xyind = xy_index(detvec[j].x - detvec[detct].x, detvec[j].y - detvec[detct].y,j);
	  trailxyvec.push_back(xyind);
	  // In the RA, Dec version we do a spherical projection here, relative to the
	  // coordinates of detvec[detct] -- but with x,y pixel coords, all we need
	  // is a simple difference.
	  ppind.push_back({}); // We need these vectors mainly just to have some way to store the
	                       // indices of mutually consistent pair partners on the next step
	}
      }
      if(trailxyvec.size() != ppind.size()) {
	cerr << "ERROR: vectors of projected and original\n";
	cerr << "pair partner candidates do not have the same length!\n";
	cerr << trailxyvec.size() << " and " << ppind.size() << " must have the same size, and do not!\n";
	return(4);
      }
      // Perform n^2 search on the projected points stored in trailxyvec
      // to find the largest subset that lie along a consistent line.
      for(i=0; i<long(trailxyvec.size()); i++) {
	ppind[i] = {};
	// Count additional pair partners (besides trailxyvec[i]) that plausibly
	// lie along the line defined by detct and trailxyvec[i].
	if(DEBUG>=2) cout << "Counting consistent pair partners\n";
	for(j=0; j<long(trailxyvec.size()); j++) {
	  if(j!=i) {
	    dist = fabs(trailxyvec[j].x*trailxyvec[i].y - trailxyvec[j].y*trailxyvec[i].x)/sqrt(trailxyvec[i].x*trailxyvec[i].x + trailxyvec[i].y*trailxyvec[i].y);
	    if(dist < 2.0*config.maxgcr) {
	      // Detections detct, i, and j all lie along a plausibly
	      // linear, constant-velocity trajectory on the sky.
	      ppind[i].push_back(j); // Store detection j as possible partner for detct and i.
	    }
	  }
	}
      }
      // Now, trailxyvec stores all the possible pair-partners of detection detct,
      // and ppind stores, for each one of these, the trailxyvec indices of ADDITIONAL ones that
      // lie on a potentially consistent trajectory with it: that is, are trail partners.
      // Find which detection in trailxyvec has the largest number of possible trail partners.
      longest_trail=-1;
      trail_length=0;
      for(i=0; i<long(trailxyvec.size()); i++) {
	if(long(ppind[i].size())+2 > trail_length) {
	  trail_length = ppind[i].size()+2; //We add one for detct, one for i, to get actual trail length
	  longest_trail = i;
	  if(DEBUG>=2) cout << "longest trail = " << longest_trail << ", length = " << trail_length << "\n";
	} else if(DEBUG>=2) cout << "not the longest\n";
      }
      if(trail_length >= config.minpoints && trail_length>=3) {
	// Load a vector with all the points for longest_trail 
	vector <xy_index> trailfitvec={};
	// Load the reference point detct.
	xyind = xy_index(0.0l,0.0l,detct);
	trailfitvec.push_back(xyind);
	// Load the pair-partner corresponding to longest_trail
	trailfitvec.push_back(trailxyvec[longest_trail]);
	// Load the additional points stored in ppind:
	for(i=0;i<long(ppind[longest_trail].size());i++) {
	  trailfitvec.push_back(trailxyvec[ppind[longest_trail][i]]);
	}
	glintpnum = trailfitvec.size();
	// Perform linear fits to x and y
	vector <double> xvec={};
	vector <double> yvec={};
	vector <double> xmod={};
	vector <double> ymod={};
	vector <double> fiterr={};
	double xslope,xintercept,yslope,yintercept,uvecx,uvecy;
	int xslopeisbest=1;
	double worsterr=-1.0l;
	long worstpoint;
	for(i=0;i<glintpnum;i++) {
	  xvec.push_back(trailfitvec[i].x);
	  yvec.push_back(trailfitvec[i].y);
	}
	status = linfituw01(xvec, yvec, xslope, xintercept);
	if(status!=0) cerr << "WARNING: linear fit to " << xvec.size() << " points with x as the independent variable has failed\n";
	status = linfituw01(yvec, xvec, yslope, yintercept);
	if(status!=0) cerr << "WARNING: linear fit to " << xvec.size() << " points with x as the independent variable has failed\n";
	if(isnormal(xslope) && isnormal(xintercept)) {
	  if(isnormal(yslope) && isnormal(yintercept) && yslope<xslope) {
	    // Prefer y slope because it is shallower, though both are valid.
	    xslopeisbest=0;
	  }
	} else if(isnormal(yslope) && isnormal(yintercept)) {
	    // Prefer y slope because the x fit contained NANs.
	    xslopeisbest=0;
	} else {
	  cerr << "Both linear fits failed in glint_trail: skipping to the next candidate.\n";
	  continue;
	}
	if(xslopeisbest==1) {
	  xmod = xvec;
	  // Calculate components of a unit vector parallel to the line
	  uvecx = 1.0l/sqrt(1+xslope*xslope);
	  uvecy = xslope/sqrt(1+xslope*xslope);
	  // Load model vector for y
	  ymod={};
	  for(i=0;i<glintpnum;i++) ymod.push_back(xintercept + xslope*xvec[i]);
	} else if(xslopeisbest==0) {
	  ymod = yvec;
	  // Calculate components of a unit vector parallel to the line
	  uvecx = yslope/sqrt(1+yslope*yslope);
	  uvecy = 1.0l/sqrt(1+yslope*yslope);
	  // Load model vector for x
	  xmod = {};
	  for(i=0;i<glintpnum;i++) xmod.push_back(yintercept + yslope*yvec[i]);
	} else {
	  cerr << "Logically excluded value for xslopeisbest: " << xslopeisbest << "\n";
	  return(3);
	}
	// Load fit residuals for each point.
	GCR = xcen = ycen = 0.0l;
	for(i=0;i<glintpnum;i++) {
	  fiterr.push_back(uvecy*(xvec[i]-xmod[i]) - uvecx*(yvec[i]-ymod[i]));
	  GCR+=fiterr[i]*fiterr[i];
	  xcen += xmod[i];
	  ycen += ymod[i];
	}
	GCR = sqrt(GCR/double(glintpnum-2));
	xcen/=double(glintpnum);
	ycen/=double(glintpnum);
	// Reject successive points until either there are only three left
	// or the worst error drops below maxgcr.
	worstpoint=1;
	while(GCR>config.maxgcr && glintpnum > 3 && glintpnum > config.minpoints && worstpoint!=0) {
	  // Identify the worst point
	  worsterr=-1.0l;
	  worstpoint=0;
	  for(i=0;i<glintpnum;i++) {
	    if(fabs(fiterr[i])>worsterr) {
	      worsterr = fabs(fiterr[i]);
	      worstpoint = i;
	    }
	  }
	  if(worstpoint!=0) {
	    // Remove the worst point from trailfitvec
	    for(i=worstpoint; i<glintpnum-1; i++) {
	      trailfitvec[i] = trailfitvec[i+1];
	    }
	    glintpnum--;
	    // Re-do the fit, through the calculation of GCR.
	    xslopeisbest=1;
	    xvec = yvec = xmod = ymod = fiterr = {};
	    for(i=0;i<glintpnum;i++) {
	      xvec.push_back(trailfitvec[i].x);
	      yvec.push_back(trailfitvec[i].y);
	    }
	    status = linfituw01(xvec, yvec, xslope, xintercept);
	    status = linfituw01(yvec, xvec, yslope, yintercept);
	    if(isnormal(xslope) && isnormal(xintercept)) {
	      if(isnormal(yslope) && isnormal(yintercept) && yslope<xslope) {
		// Prefer y slope because it is shallower, though both are valid.
		xslopeisbest=0;
	      }
	    } else if(isnormal(yslope) && isnormal(yintercept)) {
	      // Prefer y slope because the x fit contained NANs.
	      xslopeisbest=0;
	    } else {
	      cerr << "Both linear fits failed in glint_trail: skipping to the next candidate.\n";
	      continue;
	    }
	    if(xslopeisbest==1) {
	      xmod = xvec;
	      // Calculate components of a unit vector parallel to the line
	      uvecx = 1.0l/sqrt(1+xslope*xslope);
	      uvecy = xslope/sqrt(1+xslope*xslope);
	      // Load model vector for y
	      ymod={};
	      for(i=0;i<glintpnum;i++) ymod.push_back(xintercept + xslope*xvec[i]);
	    } else if(xslopeisbest==1) {
	      ymod = yvec;
	      // Calculate components of a unit vector parallel to the line
	      uvecx = yslope/sqrt(1+yslope*yslope);
	      uvecy = 1.0l/sqrt(1+yslope*yslope);
	      // Load model vector for x
	      xmod={};
	      for(i=0;i<glintpnum;i++) xmod.push_back(yintercept + yslope*yvec[i]);
	    }
	    // Load fit residuals for each point.
	    GCR = xcen = ycen = 0.0l;
	    fiterr={};
	    for(i=0;i<glintpnum;i++) {
	      fiterr.push_back(uvecy*(xvec[i]-xmod[i]) - uvecx*(yvec[i]-ymod[i]));
	      GCR+=fiterr[i]*fiterr[i];
	      xcen += xmod[i];
	      ycen += ymod[i];
	    }
	    GCR = sqrt(GCR/double(glintpnum-2));
	    xcen/=double(glintpnum);
	    ycen/=double(glintpnum);
	  }
	}
	if(GCR<=config.maxgcr && glintpnum > 3 && glintpnum > config.minpoints) {
	  // We rejected enough points and we still have a valid glint trail
	  linrms=GCR; // RMS scatter from best linear fit, regardless of spacing.
	  vector <double> distalong={};
	  double flashfreq,mindistalong,maxdistalong;
	  for(i=0;i<glintpnum;i++) {
	    distalong.push_back(uvecx*(xvec[i]-xcen) + uvecy*(yvec[i]-ycen));
	  }
	  sort(distalong.begin(),distalong.end());
	  mindistalong = distalong[0];
	  maxdistalong = distalong[glintpnum-1];
	  //for(i=0;i<glintpnum;i++) {
	  //  cout << "distalong[" << i << "] = " << distalong[i] << "\n";
	  //}
	  traillen = maxdistalong-mindistalong;
	  // Re-define the zeropoint of distalong so that it starts from the
	  // first flash, and there are no negative entries
	  for(i=0;i<glintpnum;i++) distalong[i] -= mindistalong;
	  // Calculate frequency under the assumption that no flashes are missed
	  flashfreq = double(glintpnum-1)/traillen; // Units are flashes per arcsecond.
	  double minfreq = flashfreq*config.freq_downscale;
	  double maxfreq = flashfreq*config.freq_upscale;
	  double bestfreq=0.0l;
	  double minrms=LARGERR;
	  double flashrms=0.0;
	  vector <double> flashphase;
	  double meanphase=0.0l;
	  flashfreq = minfreq;
	  //cout << "flashfreq = " << flashfreq << ", maxfreq = " << maxfreq << "\n";
	  while(flashfreq <= maxfreq) {
	    flashphase = {};
	    for(i=0;i<glintpnum;i++) {
	      flashphase.push_back(distalong[i]*flashfreq+0.5l - floor(distalong[i]*flashfreq+0.5l));
	      // Note: for good phase, we expect distalong[i]*flashfreq will always be close to an integer,
	      // hence we add 0.5 so all flashphase entries will cluster near 0.5 instead of
	      // being split between slightly more than 0.0 and slightly less than 1.0.
	    }
	    status = dmeanrms01(flashphase, &meanphase, &flashrms);
	    if(flashrms<minrms) {
	      minrms = flashrms;
	      bestfreq = flashfreq;
	    }
	    //cout << "freq, rms, best = " << flashfreq << " " << flashrms << " " << bestfreq << "\n";
	    flashfreq += config.max_phase_err/traillen;
	  }
	  flashnum = bestfreq*traillen+1.5; // Best integer guess at the total number of flashes,
	                                    // including any that were not detected.
	  stepsize = 1.0l/bestfreq;
	  //cout << "bestfreq: " << bestfreq << " " << bestfreq*traillen << " " << flashnum << " " << stepsize << "\n";
	  // Re-calculate mean phase of the flashes based on the best frequency
	  flashphase = {};
	  for(i=0;i<glintpnum;i++) {
	    flashphase.push_back(distalong[i]*bestfreq+0.5l - floor(distalong[i]*bestfreq+0.5l));
	    // Note: for good phase, we expect distalong[i]*flashfreq will always be close to an integer,
	    // hence we add 0.5 so all flashphase entries will cluster near 0.5 instead of
	    // being split between slightly more than 0.0 and slightly less than 1.0.
	  }
	  status = dmeanrms01(flashphase, &meanphase, &flashrms);
	  // Note that the mean phase will be exactly 0.5 if the first flash,
	  // from which distalong is measured, was precisely on-schedule as
	  // defined by the mean of all the others.
	  // Predict locations of flashes, based on best frequency identified.
	  // Must first reload distalong to match the order of xvec and yvec.
	  distalong={};
	  for(i=0;i<glintpnum;i++) {
	    distalong.push_back(uvecx*(xvec[i]-xcen) + uvecy*(yvec[i]-ycen) - mindistalong);
	  }
	  // Now proceed with the prediction.
	  xmod = ymod = {};
	  eqrms = 0.0l;
	  for(i=0;i<glintpnum;i++) {
	    // Calculate phase of idealized, equally-spaced flash corresponding to point i
	    double distref = floor(distalong[i]*bestfreq+0.5l) + meanphase - 0.5l;
	    //cout << "i=" << i << ", ideal phase = " << distref << " ";
	    // Translate the phase into distance in arcseconds relative to the trail center
	    distref /= bestfreq;
	    //cout << distref << " arcsec from , ";
	    // Change the reference point to the trail center rather than the first flash
	    distref += mindistalong;
	    //cout << distref << " arcsec from trail center\n";
	    // Convert this into an x,y position using the known center of the
	    // trail and the along-trail unit vector.
	    xmod.push_back(xcen + distref*uvecx);
	    ymod.push_back(ycen + distref*uvecy);
	    eqrms += DSQUARE(xvec[i] - xmod[i]) + DSQUARE(yvec[i] - ymod[i]);
	  }
	  eqrms = sqrt(eqrms/double(glintpnum));

	  // Convert trail unit vector into a celestial position angle.
	  if(uvecx==0.0l && uvecy>=0.0l) pa = 0.0l;
	  else if(uvecx==0.0l && uvecy<0.0l) pa=180.0;
	  else if(uvecx>0.0l) pa = 270.0 +  DEGPRAD*atan(uvecy/uvecx);
	  else if(uvecx<0.0l) pa = 90.0 +  DEGPRAD*atan(uvecy/uvecx);

	  // Calculate magnitude RMS
	  vector <double> magvec = {};
	  for(i=0;i<glintpnum;i++) magvec.push_back(detvec[trailfitvec[i].index].z);
	  status = dmeanrms01(magvec, &magmean, &magrms);
	
	  // Load output vectors
	  oneglint = glint_trail(detvec[detct].x + xcen, detvec[detct].y + ycen, traillen, pa, linrms, eqrms, magmean, magrms, stepsize, qc1, glintpnum, flashnum);
	  //cout << "stepsize = " << oneglint.stepsize << "\n";
	  alltrailvec.push_back(oneglint);
	  vector <long> indvec={};
	  for(i=0;i<glintpnum;i++) indvec.push_back(trailfitvec[i].index);
	  indvecs.push_back(indvec);
	  // Load vector for sorting best among overlapping glint trails
	  dindex = double_index(double(indvec.size()) + 1.0l - GCR/config.maxgcr, glintnum);
	  dindvec.push_back(dindex);
	  // Load vector keeping track of the anchor detection for each glint trail
	  anchor_det.push_back(detct);
	  // Increment count of valid glints
	  glintnum++;
	}
      }
    }
  }
  cout << "De-duplicating a total of " << glintnum << " candidate glint trails\n";
  sort(dindvec.begin(), dindvec.end(), lower_double_index());
  // Best glint trail (longest, with lowest GCR) will be at the end
  long goodglints=0;
  for(i=0; i<glintnum; i++) {
    long thisglint = dindvec[glintnum-1-i].index;
    if(detused[anchor_det[thisglint]]==0) {
      // The anchor detection has not already been claimed by another glint trail
      // Load it onto the final output vectors
      trailvec.push_back(alltrailvec[thisglint]);
      // cout << "stepsize = " << alltrailvec[thisglint].stepsize << "\n";
      // Write out the individual detection indices, and mark all of the
      // individual detections as used.
      for(j=0; j<long(indvecs[thisglint].size()); j++) {
	onepair = longpair(goodglints,indvecs[thisglint][j]);
	trail2det.push_back(onepair);
	detused[indvecs[thisglint][j]] = 1;
      }
      goodglints++;
    }
  }
  cout << "Final de-duplicated total of glint trails is " << goodglints << "\n";
  return(0);
}
