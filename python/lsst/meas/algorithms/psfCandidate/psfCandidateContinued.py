from __future__ import absolute_import

from lsst.utils import continueClass
from .psfCandidate import PsfCandidateF

__all__ = []  # import only for side effects

@continueClass
class PsfCandidateF:
    getCandidateRating = PsfCandidateF._getCandidateRating

    def setCandidateRating(rating):
        raise NotImplementedError(("must not exist for this type "
            "since getCandidateRating is calculated"))

