import lsst.meas.algorithms.ingestIndexReferenceTask
assert type(config)==lsst.meas.algorithms.ingestIndexReferenceTask.DatasetConfig, 'config is of type %s.%s instead of lsst.meas.algorithms.ingestIndexReferenceTask.DatasetConfig' % (type(config).__module__, type(config).__name__)
# Version number of the persisted on-disk storage format. Version 0 had Jy as flux units (default 0 for unversioned catalogs). Version 1 had nJy as flux units.
config.format_version=1

# String to pass to the butler to retrieve persisted files.
config.ref_dataset_name='cal_ref_cat'

# Depth of the HTM tree to make.  Default is depth=7 which gives ~ 0.3 sq. deg. per trixel.
config.indexer['HTM'].depth=4

config.indexer.name='HTM'
