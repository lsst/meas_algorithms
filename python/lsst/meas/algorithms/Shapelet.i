// -*- lsst-c++ -*-

/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
 * 
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the LSST License Statement and 
 * the GNU General Public License along with this program.  If not, 
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */

%include "lsst/meas/algorithms/algorithms_fwd.i"

%{
#include "lsst/meas/algorithms/Shapelet.h"
#include "lsst/meas/algorithms/ShapeletInterpolation.h"
#include "lsst/meas/algorithms/ShapeletKernel.h"
#include "lsst/meas/algorithms/ShapeletPsfCandidate.h"
#include "lsst/meas/algorithms/SizeMagnitudeStarSelector.h"
#include "lsst/meas/algorithms/ShapeletPsf.h"
%}

%shared_ptr(lsst::meas::algorithms::Shapelet)
%shared_ptr(lsst::meas::algorithms::ShapeletInterpolation)
%shared_ptr(lsst::meas::algorithms::LocalShapeletKernel);
%shared_ptr(lsst::meas::algorithms::ShapeletKernel);
%shared_ptr(lsst::meas::algorithms::ShapeletPsfCandidate);
%shared_ptr(lsst::meas::algorithms::ShapeletPsf);
%shared_vec(lsst::meas::algorithms::SizeMagnitudeStarSelector::PsfCandidateList);
%shared_ptr(std::vector<lsst::meas::algorithms::SizeMagnitudeStarSelector::PsfCandidateList>);

%include "lsst/meas/algorithms/Shapelet.h"
%include "lsst/meas/algorithms/ShapeletInterpolation.h"
%include "lsst/meas/algorithms/ShapeletKernel.h"
%include "lsst/meas/algorithms/ShapeletPsfCandidate.h"
%include "lsst/meas/algorithms/SizeMagnitudeStarSelector.h"
%include "lsst/meas/algorithms/ShapeletPsf.h"
