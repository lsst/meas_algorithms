# 'namedtuple' is present in the 'collections' module in Python 2.6, but since the current LSST standard
# Python is 2.5.2, this module is provided.  The code and the documentation are taken from
# http://code.activestate.com/recipes/500261/
#
#
# There has long been a need for named access to fields in records stored as tuples. In response, people have
# crafted many different versions of this recipe. I've combined the best of their approaches with a few ideas
# of my own. The resulting recipe has been well received, so it was proposed and accepted for inclusion in the
# collections module for Py2.6.
#
# Docs and examples for the module can be found on the bottom of the page at:
# http://docs.python.org/dev/library/collections.html#namedtuple-factory-function-for-tuples-with-named-fields
#
# The principal features are:
# 
#     * Easy to type/read/modify function signature: named_tuple('Person', 'name age sex height nationality').
#     * C-speed attribute lookup using property() and itemgetter().
#     * Named tuples have no instance dictionary, so their instances take no more space than a regular tuple
#       (for example, casting thousands of sql records to named tuples has zero memory overhead).
#     * Nice docstring is helpful with an editor's tooltips.
#     * Optional keywords in the contructor for readability and to allow arguments to be specified in
#       arbitrary order: Person(name='susan', height=60, nationality='english', sex='f', age=30).
#     * Key/Value style repr for clearer error messages and for usability at the interactive prompt.
#     * Named tuples can be pickled.
#     * Clean error messages for missing or misnamed arguments.
#     * A method similar to str.replace() using a field name. The _replace() method is used instead of slicing
#       for updating fields. For example use t.replace(f=newval) instead of tuple concatenations like
#       t[:2]+newval+t[3:].
#     * A _fields attribute exposes the field names for introspection.
#     * An _asdict() method for converting to an equivalent dictionary. If desired an ordered dictionary can
#       be substituted for the dict constructor.
#     * Recipe runs on Py2.4 or later.
#     * Option for automatic renaming of invalid field names to positional names (to support cases where field
#       names are supplied externally).
# 
# Note, the idea for the __module__ attribute was suggested by Michele Simionato and Louis Riviere. That
# attribute makes pickling possible, and it lets help() correctly identify the module where the named tuple
# was defined.
# 
# Thanks to Peter Kovac pointing-out deficiencies in the keyword argument checking. Because of his comments,
# the recipe has evolved to its current exec-style where we get all of Python's high-speed builtin argument
# checking for free. The new style of building and exec-ing a template made both the __new__ and __repr__
# functions faster and cleaner than in previous versions of this recipe.
# 
# At the suggestion of Martin Blais and Robin Becker, the field name spec can be either a tuple or a
# string. When the field names are coming from a CSV header or other automatically generated source, it is
# simplest to pass in a tuple of fieldnames. When the field names are already known and written-out
# explicitly, the string form is best because it is easier to type, easier to read, and easier to rearrange
# when necessary. The string form works especially well with SQL use cases because the string can be
# cut-and-pasted from the field spec portion of an SQL query.
# 
# The function signature with the typename separate from the field names was suggested by Louis Riviere and
# Jim Jewett.
# 
# The idea to use keyword arguments for _replace() was inspired by Issac Morlund's suggestion to be able to do
# multiple replacements in one pass.
# 
# The inspiration for the _make() classmethod came from Robin Becker and Giovanni Bajo who pointed-out an
# important class of use cases where existing sequences need to be cast to named tuples.
# 
# The idea for the _fields attribute was suggested by Robin Becker.
# 
# The idea for the _asdict() method was inspired by the thought that any class with name/value pairs is
# conceptually a mapping. Accordingly, instances of that class should be readily convertible to and from a
# dict.
# 

from operator import itemgetter as _itemgetter
from keyword import iskeyword as _iskeyword
import sys as _sys

def namedtuple(typename, field_names, verbose=False, rename=False):
    """Returns a new subclass of tuple with named fields.

    >>> Point = namedtuple('Point', 'x y')
    >>> Point.__doc__                   # docstring for the new class
    'Point(x, y)'
    >>> p = Point(11, y=22)             # instantiate with positional args or keywords
    >>> p[0] + p[1]                     # indexable like a plain tuple
    33
    >>> x, y = p                        # unpack like a regular tuple
    >>> x, y
    (11, 22)
    >>> p.x + p.y                       # fields also accessable by name
    33
    >>> d = p._asdict()                 # convert to a dictionary
    >>> d['x']
    11
    >>> Point(**d)                      # convert from a dictionary
    Point(x=11, y=22)
    >>> p._replace(x=100)               # _replace() is like str.replace() but targets named fields
    Point(x=100, y=22)

    """

    # Parse and validate the field names.  Validation serves two purposes,
    # generating informative error messages and preventing template injection attacks.
    if isinstance(field_names, basestring):
        field_names = field_names.replace(',', ' ').split() # names separated by whitespace and/or commas
    field_names = tuple(map(str, field_names))
    if rename:
        names = list(field_names)
        seen = set()
        for i, name in enumerate(names):
            if (not min(c.isalnum() or c=='_' for c in name) or _iskeyword(name)
                or not name or name[0].isdigit() or name.startswith('_')
                or name in seen):
                    names[i] = '_%d' % i
            seen.add(name)
        field_names = tuple(names)
    for name in (typename,) + field_names:
        if not min(c.isalnum() or c=='_' for c in name):
            raise ValueError('Type names and field names can only contain alphanumeric characters and underscores: %r' % name)
        if _iskeyword(name):
            raise ValueError('Type names and field names cannot be a keyword: %r' % name)
        if name[0].isdigit():
            raise ValueError('Type names and field names cannot start with a number: %r' % name)
    seen_names = set()
    for name in field_names:
        if name.startswith('_') and not rename:
            raise ValueError('Field names cannot start with an underscore: %r' % name)
        if name in seen_names:
            raise ValueError('Encountered duplicate field name: %r' % name)
        seen_names.add(name)

    # Create and fill-in the class template
    numfields = len(field_names)
    argtxt = repr(field_names).replace("'", "")[1:-1]   # tuple repr without parens or quotes
    reprtxt = ', '.join('%s=%%r' % name for name in field_names)
    template = '''class %(typename)s(tuple):
        '%(typename)s(%(argtxt)s)' \n
        __slots__ = () \n
        _fields = %(field_names)r \n
        def __new__(_cls, %(argtxt)s):
            return _tuple.__new__(_cls, (%(argtxt)s)) \n
        @classmethod
        def _make(cls, iterable, new=tuple.__new__, len=len):
            'Make a new %(typename)s object from a sequence or iterable'
            result = new(cls, iterable)
            if len(result) != %(numfields)d:
                raise TypeError('Expected %(numfields)d arguments, got %%d' %% len(result))
            return result \n
        def __repr__(self):
            return '%(typename)s(%(reprtxt)s)' %% self \n
        def _asdict(self):
            'Return a new dict which maps field names to their values'
            return dict(zip(self._fields, self)) \n
        def _replace(_self, **kwds):
            'Return a new %(typename)s object replacing specified fields with new values'
            result = _self._make(map(kwds.pop, %(field_names)r, _self))
            if kwds:
                raise ValueError('Got unexpected field names: %%r' %% kwds.keys())
            return result \n
        def __getnewargs__(self):
            return tuple(self) \n\n''' % locals()
    for i, name in enumerate(field_names):
        template += '        %s = _property(_itemgetter(%d))\n' % (name, i)
    if verbose:
        print template

    # Execute the template string in a temporary namespace
    namespace = dict(_itemgetter=_itemgetter, __name__='namedtuple_%s' % typename,
                     _property=property, _tuple=tuple)
    try:
        exec template in namespace
    except SyntaxError, e:
        raise SyntaxError(e.message + ':\n' + template)
    result = namespace[typename]

    # For pickling to work, the __module__ variable needs to be set to the frame
    # where the named tuple is created.  Bypass this step in enviroments where
    # sys._getframe is not defined (Jython for example) or sys._getframe is not
    # defined for arguments greater than 0 (IronPython).
    try:
        result.__module__ = _sys._getframe(1).f_globals.get('__name__', '__main__')
    except (AttributeError, ValueError):
        pass

    return result






if __name__ == '__main__':
    # verify that instances can be pickled
    from cPickle import loads, dumps
    Point = namedtuple('Point', 'x, y', True)
    p = Point(x=10, y=20)
    assert p == loads(dumps(p, -1))

    # test and demonstrate ability to override methods
    class Point(namedtuple('Point', 'x y')):
        @property
        def hypot(self):
            return (self.x ** 2 + self.y ** 2) ** 0.5
        def __str__(self):
            return 'Point: x=%6.3f y=%6.3f hypot=%6.3f' % (self.x, self.y, self.hypot)

    for p in Point(3,4), Point(14,5), Point(9./7,6):
        print p

    class Point(namedtuple('Point', 'x y')):
        'Point class with optimized _make() and _replace() without error-checking'
        _make = classmethod(tuple.__new__)
        def _replace(self, _map=map, **kwds):
            return self._make(_map(kwds.get, ('x', 'y'), self))

    print Point(11, 22)._replace(x=100)

    import doctest
    TestResults = namedtuple('TestResults', 'failed attempted')
    print TestResults(*doctest.testmod())
