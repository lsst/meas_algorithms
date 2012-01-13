import lsst.pex.config as pexConfig

class SizeMagnitudeStarSelector(pexConfig.Config):
    minSize = pexConfig.Field(
        doc = "Minimum size to use",
        dtype = float,
        default = 0.0,
    )
    maxSize = pexConfig.Field(
        doc = "Maximum size to use",
        dtype = float,
        default = 1.0e100,
    )
    isSizeLog = pexConfig.Field(
        doc = "Are sizes already log(size)?",
        dtype = bool,
        default = False,
    )
    minMag = pexConfig.Field(
        doc = "Minimum magnitude to use",
        dtype = float,
        default = 0.0,
    )
    maxMag = pexConfig.Field(
        doc = "Maximum magnitude to use",
        dtype = float,
        default = 1.0e100,
    )
    starFrac = pexConfig.Field(
        doc = "What fraction of objects are likely stars?",
        dtype = float,
        default = 0.5,
    )
    startN = pexConfig.Field(
        doc = "Fraction of objects to use in first pass",
        dtype = float,
        default = 0.1,
    )
    fitOrder = pexConfig.Field(
        doc = "Order of polynomial of fit of size(x,y)",
        dtype = int,
        default = 1,
    )
    fitSigClip = pexConfig.Field(
        doc = "nSigma to reject a star as an outlier",
        dtype = float,
        default = 4.0,
    )
    fitStars = pexConfig.Field(
        doc = "Perform size(x,y) fit with fitStars brightest stars",
        dtype = int,
        default = 30,
    )
    purity = pexConfig.Field(
        doc = "Smaller = purer smaple of stars, larger = more stars",
        dtype = float,
        default = 0.05,
    )
    aperture = pexConfig.Field(
        doc = "nSigma to reject a star as an outlier",
        dtype = float,
        default = 5.0,
    )
