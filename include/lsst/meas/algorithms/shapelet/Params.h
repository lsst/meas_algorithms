#ifndef MeasAlgoShapeletParams_H
#define MeasAlgoShapeletParams_H

#include <stdexcept>

// Default value for 
#define DEFVALPOS 9999
#define DEFVALNEG -9999

//
// Flags 
//
#define INPUT_FLAG              0x1
#define TRANSFORM_EXCEPTION     0x2
#define FITTEDPSF_EXCEPTION     0x4
#define TMV_EXCEPTION           0x8
#define STD_EXCEPTION           0x10
#define UNKNOWN_EXCEPTION       0x20
#define EDGE                    0x40
#define LT10PIX                 0x80
#define MEASURE_PSF_FAILED      0x100
#define NATIVE_FAILED           0x200
#define TOO_SMALL               0x400
#define DECONV_FAILED           0x800
#define SHEAR_FAILED            0x1000
#define SHAPELET_FAILED         0x2000
#define UNKNOWN_FAILURE         0x4000
#define SHAPE_REDUCED_ORDER     0x8000
#define SHEAR_LOCAL_MIN         0x10000
#define SHEAR_POOR_FIT          0x20000
#define SHAPE_LOCAL_MIN         0x40000
#define SHAPE_POOR_FIT          0x80000
#define SHEAR_BAD_COVAR         0x100000
#define NO_SINGLE_EPOCH_IMAGES	0x200000
#define BKG_NOPIX		        0x400000
#define PSF_INTERP_OUTLIER	    0x800000
#define SHAPE_BAD_FLUX          0x1000000
#define PSF_BAD_FLUX            0x2000000


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

// StarSelectorException      -- Treat as ProcessingException
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

#endif
