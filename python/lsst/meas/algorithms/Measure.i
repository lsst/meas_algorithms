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
#include "lsst/meas/algorithms/Measure.h"
%}

%returnCopy(lsst::meas::algorithms::MeasureSources::getAlgorithms)
%returnSelf(lsst::meas::algorithms::MeasureSourcesBuilder::setCentroider)
%returnSelf(lsst::meas::algorithms::MeasureSourcesBuilder::addAlgorithm)

%include "lsst/meas/algorithms/Measure.h"

%extend lsst::meas::algorithms::MeasureSources {
%template(apply) apply<float>;
%template(apply) apply<double>;
%template(apply) apply<float>;
%template(apply) apply<double>;
%template(applyWithCoord) applyWithCoord<float>;
%template(applyWithCoord) applyWithCoord<double>;
%template(applyWithPixel) applyWithPixel<float>;
%template(applyWithPixel) applyWithPixel<double>;
}
%extend lsst::meas::algorithms::MeasureSourcesBuilder {
%pythoncode %{
    def addAlgorithms(self, iterable):
        for item in iterable:
            self.addAlgorithm(item)
        return self
%}
}
