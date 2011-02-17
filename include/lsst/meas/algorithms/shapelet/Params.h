#ifndef MeasAlgoShapeletParams_H
#define MeasAlgoShapeletParams_H

#include <stdexcept>
#include <vector>
#include <ostream>

namespace lsst {
namespace meas {
namespace algorithms {
namespace shapelet {

    // Default value for various measured quantities
#define DEFVALPOS 9999
#define DEFVALNEG -9999

    //
    // Flags 
    //
    static const long INPUT_FLAG              = 0x1;
    static const long TRANSFORM_EXCEPTION     = 0x2;
    static const long FITTEDPSF_EXCEPTION     = 0x4;
    static const long TMV_EXCEPTION           = 0x8;
    static const long STD_EXCEPTION           = 0x10;
    static const long UNKNOWN_EXCEPTION       = 0x20;
    static const long EDGE                    = 0x40;
    static const long LT10PIX                 = 0x80;
    static const long MEASURE_PSF_FAILED      = 0x100;
    static const long NATIVE_FAILED           = 0x200;
    static const long TOO_SMALL               = 0x400;
    static const long DECONV_FAILED           = 0x800;
    static const long SHEAR_FAILED            = 0x1000;
    static const long SHAPELET_FAILED         = 0x2000;
    static const long SHEAR_REDUCED_ORDER     = 0x4000;
    static const long SHAPE_REDUCED_ORDER     = 0x8000;
    static const long SHEAR_LOCAL_MIN         = 0x10000;
    static const long SHEAR_POOR_FIT          = 0x20000;
    static const long SHAPE_LOCAL_MIN         = 0x40000;
    static const long SHAPE_POOR_FIT          = 0x80000;
    static const long SHEAR_BAD_COVAR         = 0x100000;
    static const long NO_SINGLE_EPOCH_IMAGES  = 0x200000;
    static const long BKG_NOPIX               = 0x400000;
    static const long PSF_INTERP_OUTLIER      = 0x800000;
    static const long SHEAR_BAD_FLUX          = 0x1000000;
    static const long PSF_BAD_FLUX            = 0x2000000;
    static const long SHAPE_BAD_FLUX          = 0x4000000;
    static const long CENTROID_FAILED         = 0x8000000;
    static const long SHAPELET_NOT_DECONV     = 0x10000000;
    static const long SHEAR_DIDNT_CONVERGE    = 0x20000000;

    static const long NFLAGS = 30;

    static const char*const flagName[NFLAGS] = {
        "INPUT_FLAG",
        "TRANSFORM_EXCEPTION",
        "FITTEDPSF_EXCEPTION",
        "TMV_EXCEPTION",
        "STD_EXCEPTION",
        "UNKNOWN_EXCEPTION",
        "EDGE",
        "LT10PIX",
        "MEASURE_PSF_FAILED",
        "NATIVE_FAILED",
        "TOO_SMALL",
        "DECONV_FAILED",
        "SHEAR_FAILED",
        "SHAPELET_FAILED",
        "SHEAR_REDUCED_ORDER",
        "SHAPE_REDUCED_ORDER",
        "SHEAR_LOCAL_MIN",
        "SHEAR_POOR_FIT",
        "SHAPE_LOCAL_MIN",
        "SHAPE_POOR_FIT",
        "SHEAR_BAD_COVAR",
        "NO_SINGLE_EPOCH_IMAGES",
        "BKG_NOPIX",
        "PSF_INTERP_OUTLIER",
        "SHEAR_BAD_FLUX",
        "PSF_BAD_FLUX",
        "SHAPE_BAD_FLUX",
        "CENTROID_FAILED",
        "SHAPELET_NOT_DECONV",
        "SHEAR_DIDNT_CONVERGE"
    };

    void PrintFlags(const std::vector<long>& flags, std::ostream& os);

    // Errors specific to the weak lensing code

    struct FileNotFoundException : public std::runtime_error 
    {
        FileNotFoundException(const std::string& filename) throw() :
            std::runtime_error("Error: file "+filename+" not found") 
        {} 
    };

    struct ParameterException : public std::runtime_error 
    {
        ParameterException(const std::string& msg) throw() :
            std::runtime_error(msg) 
        {}
    };

    struct ReadException : public std::runtime_error
    {
        ReadException(const std::string& msg) throw() :
            std::runtime_error(msg) 
        {}
    };

    struct WriteException : public std::runtime_error 
    {
        WriteException(const std::string& msg) throw() :
            std::runtime_error(msg) 
        {}
    };

    struct ProcessingException : public std::runtime_error 
    {
        ProcessingException(const std::string& msg) throw() :
            std::runtime_error(msg) 
        {}
    };

    // Errors that may be thrown by the weak lensing code, but 
    // defined in other files

    // StarFinderException        -- Treat as ProcessingException
    // AssertFailure              -- Treat as ProcessingException
    // tmv::Error                 -- Treat as ProcessingException
    // std::exception             -- Treat as ProcessingException



    //
    // Exit codes
    //

    enum ExitCode { 
        SUCCESS = 0,
        FAILURE,
        FAILURE_FILE_NOT_FOUND,
        FAILURE_PARAMETER_ERROR,
        FAILURE_READ_ERROR,
        FAILURE_WRITE_ERROR,
        FAILURE_PROCESSING_ERROR
    };

    const char* Text(const ExitCode& code);

    class ConfigFile;
    int Status(ExitCode code, const ConfigFile& params);


    // tolerance for testing output files
#define TEST_TOL 1.e-6

}}}}

#endif
