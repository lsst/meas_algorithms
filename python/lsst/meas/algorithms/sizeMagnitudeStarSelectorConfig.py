import lsst.pex.config as pexConfig

class SizeMagnitudeStarSelector(pexConfig.Config):
    minSize = pexConfig.Field(
        float,
        doc = "Minimum size to use",
        default = 0.0,
    )
    maxSize = pexConfig.Field(
        float,
        doc = "Maximum size to use",
        default = 1.0e100,
    )
    isSizeLog = pexConfig.Field(
        bool,
        doc = "Are sizes already log(size)?",
        default = False,
    )
    minMag = pexConfig.Field(
        float,
        doc = "Minimum magnitude to use",
        default = 0.0,
    )
    maxMag = pexConfig.Field(
        float,
        doc = "Maximum magnitude to use",
        default = 1.0e100,
    )
    starFrac = pexConfig.Field(
        float,
        doc = "What fraction of objects are likely stars?",
        default = 0.5,
    )
    startN = pexConfig.Field(
        float,
        doc = "Fraction of objects to use in first pass",
        default = 0.1,
    )
    fitOrder = pexConfig.Field(
        int,
        doc = "Order of polynomial of fit of size(x,y)",
        default = 1,
    )
    fitSigClip = pexConfig.Field(
        float,
        doc = "nSigma to reject a star as an outlier",
        default = 4.0,
    )
    fitStars = pexConfig.Field(
        int,
        doc = "Perform size(x,y) fit with fitStars brightest stars",
        default = 30,
    )
    purity = pexConfig.Field(
        float,
        doc = "Smaller = purer smaple of stars, larger = more stars",
        default = 0.05,
    )
    aperture = pexConfig.Field(
        float,
        doc = "nSigma to reject a star as an outlier",
        default = 5.0,
    )
