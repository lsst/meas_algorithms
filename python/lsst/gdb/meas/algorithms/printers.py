import gdb
import sys

try:
    import gdb.printing

    class CRPixelPrinter:
        "Print a CRPixel"

        def __init__(self, val):
            self.val = val

        def to_string(self):
            return "{id=%d (%d, %d)}" % (self.val["id"], self.val["col"], self.val["row"])

    printers = []

    def register(obj):
        "Register my pretty-printers with objfile Obj."

        if obj is None:
            obj = gdb

        for p in printers:
            gdb.printing.register_pretty_printer(obj, p)

    def build_meas_algorithms_dictionary():
        printer = gdb.printing.RegexpCollectionPrettyPrinter("meas_algorithms")

        printer.add_printer('lsst::meas::algorithms::CRPixel',
                            '^lsst::meas::algorithms::CRPixel', CRPixelPrinter)
        return printer

    printers.append(build_meas_algorithms_dictionary())

except ImportError as e:
    def register(*args, exception=e, **kwargs):
        print("Your version of gdb is too old to load the meas.algorithms python pretty printers: %s" %
              (exception,), file=sys.stderr)
        pass

    pass
