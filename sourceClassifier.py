
class SourceClassifier(object):
    """
    Base class for source classifiers. A SourceClassifier is initialized with a policy, and the positions of
    the flag bits it is allowed to set/clear. Once created, classify() is called on the SourceClassifier instance
    some number of times - each time one or more sources (typically a pair of measurements on the two exposures
    in an LSST visit) are supplied as arguments.

    Finally, the finish() method is called on the classifier so that it can log summary statistics and/or
    place SDQA ratings onto a stage clipboard.
    """
    def __init__(self, bits, policy):
        self._bits = bits
        self._policy = policy

    def classify(self, *args):
        """
        Classify one or more sources. To be overriden by subclasses.
        """
        pass

    def finish(self, log, clipboard):
        """
        Lifecycle method to be overriden by subclasses. Called once by SourceClassificationStage
        after the SourceClassifier instance has been used to classify one visits worth of data. Can
        be used to log visit summary statistics or place SDQA ratings onto the stage clipboard.
        """
        pass

    def getPolicy(self):
        """
        Return the policy containing classifier configuration parameters
        """
        return self._policy

    def getMask(self, n=0):
        """
        Return a mask for the n-th bit this classifier sets/clears
        """
        return 1 << self._bits[n]

    def getBit(self, n=0):
        """
        Return the index of the n-th bit this classifier sets/clears
        """
        return self._bits[n]

    def setBit(self, flag, n=0):
        """
        Set the n-th classifier bit in the given integer
        """
        return flag | self.getMask(n)

    def clearBit(self, flag, n=0):
        """
        Clear the n-th classifier bit in the given integer
        """
        return flag ^ (flag & self.getMask(n))

