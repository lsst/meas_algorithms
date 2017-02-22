from __future__ import absolute_import

from .psfCandidate import PsfCandidateF

__all__ = []  # import only for side effects

PsfCandidateF.getCandidateRating = PsfCandidateF._getCandidateRating

def setCandidateRating(rating):
    raise NotImplementedError(("must not exist for this type "
        "since getCandidateRating is calculated"))

PsfCandidateF.setCandidateRating = setCandidateRating

