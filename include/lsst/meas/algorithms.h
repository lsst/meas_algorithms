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
 
#include "algorithms/Measure.h"
#include "lsst/meas/algorithms/CR.h"
#include "lsst/meas/algorithms/Interp.h"
#include "lsst/meas/algorithms/PSF.h"
#include "lsst/meas/algorithms/PsfCandidate.h"
#include "lsst/meas/algorithms/SpatialModelPsf.h"
#include "lsst/meas/algorithms/Measure.h"
#include "lsst/meas/algorithms/Shapelet.h"
#include "lsst/meas/algorithms/ShapeletInterpolation.h"
#include "lsst/meas/algorithms/ShapeletKernel.h"
#include "lsst/meas/algorithms/ShapeletPsfCandidate.h"
#include "lsst/meas/algorithms/SizeMagnitudeStarSelector.h"
#include "lsst/meas/algorithms/detail/SincPhotometry.h"
#include "lsst/meas/algorithms/FluxControl.h"
#include "lsst/meas/algorithms/CentroidControl.h"
#include "lsst/meas/algorithms/ShapeControl.h"
#include "lsst/meas/algorithms/CorrectFluxes.h"
#include "lsst/meas/algorithms/Classification.h"
#include "lsst/meas/algorithms/PixelFlags.h"
#include "lsst/meas/algorithms/SkyCoord.h"
#include "lsst/meas/algorithms/RecordCentroid.h"
#include "lsst/meas/algorithms/KernelPsf.h"
#include "lsst/meas/algorithms/SingleGaussianPsf.h"
#include "lsst/meas/algorithms/DoubleGaussianPsf.h"
#include "lsst/meas/algorithms/PcaPsf.h"
#include "lsst/meas/algorithms/SdssShapeControl.h"
#include "lsst/meas/algorithms/GaussianFluxControl.h"
#include "lsst/meas/algorithms/CoaddPsf.h"
#include "lsst/meas/algorithms/WarpedPsf.h"
#include "lsst/meas/algorithms/Photometry.h"
#include "lsst/meas/algorithms/Jacobian.h"
#include "lsst/meas/algorithms/FocalPlane.h"
#include "lsst/meas/algorithms/BinnedWcs.h"
#include "lsst/meas/algorithms/CoaddBoundedField.h"
