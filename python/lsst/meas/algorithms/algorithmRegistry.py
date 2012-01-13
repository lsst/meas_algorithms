# 
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
# 
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the LSST License Statement and 
# the GNU General Public License along with this program.  If not, 
# see <http://www.lsstcorp.org/LegalNotices/>.
#
class AlgorithmRegistry(object):
    """A registry of algorithms
    """
    def __init__(self):
        """Construct an AlgorithmRegistry
        """
        self._dict = dict()
    
    def get(self, name):
        """Get an algorithm by name
        
        @raise KeyError if item name is not found.
        """
        return self._dict[name]
    
    def getConfig(self, name):
        """Get an algorithm Config by name

        @raise KeyError if item name is not found.
        """
        return self._dict[name].ConfigClass
    
    def getNames(self):
        """Return a list of algorithm names
        """
        return _dict.keys()
    
    def getConfigDict(self):
        """Return a dictionary of name:config pairs
        
        This is useful for constructing a pex_config RegistryField.
        """
        return dict((key, val.ConfigClass) for key, val in self._dict.iteritems())
    
    def register(self, name, alg):
        """Register a new algorithm. The name must be unique.
        
        The algorithm must have an attribute ConfigClass -- ENFORCE THIS!
        
        To avoid collisions, you may wish to include the name of your package in the determiner name, e.g.:
        "my_pkg.myAlg"
        
        @raise KeyError if the name is a duplicate
        """
        if name in self._dict:
            raise KeyError("A PSF determiner already exists with name %r" % (name,))
        self._dict[name] = alg
