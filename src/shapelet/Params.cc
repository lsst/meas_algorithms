// -*- LSST-C++ -*-

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

#include "lsst/meas/algorithms/shapelet/Params.h"
#include "lsst/meas/algorithms/shapelet/ConfigFile.h"

namespace lsst {
namespace meas {
namespace algorithms {
namespace shapelet {

    const char* Text(const ExitCode& code)
    {
        switch (code) {
          case SUCCESS :
               return "SUCCESS";
          case FAILURE :
               return "FAILURE";
          case FAILURE_FILE_NOT_FOUND :
               return "FAILURE_FILE_NOT_FOUND";
          case FAILURE_PARAMETER_ERROR :
               return "FAILURE_PARAMETER_ERROR";
          case FAILURE_READ_ERROR :
               return "FAILURE_READ_ERROR";
          case FAILURE_WRITE_ERROR :
               return "FAILURE_WRITE_ERROR";
          case FAILURE_PROCESSING_ERROR :
               return "FAILURE_PROCESSING_ERROR";
          default :
               return "UNKNOWN";
        }
    }

    int Status(ExitCode code, const ConfigFile& params)
    {
        switch (code) {
          case SUCCESS :
               return params.read("success_status",2);
          case FAILURE :
               return params.read("failure_status",4);
          case FAILURE_FILE_NOT_FOUND :
               return params.read("file_not_found_status",5);
          case FAILURE_PARAMETER_ERROR :
               return params.read("parameter_error_status",5);
          case FAILURE_READ_ERROR :
               return params.read("read_error_status",5);
          case FAILURE_WRITE_ERROR :
               return params.read("write_error_status",4);
          case FAILURE_PROCESSING_ERROR :
               return params.read("processing_error_status",4);
          default :
               return 0;
        }
    }

}}}}
