"""
Generic utility functions.
"""


def to_dicts(keys, tuples):
    """Given an n-tuple of keys and a list of n-tuples, make a list of dicts."""
    return [{k: v for k, v in zip(keys, o)} for o in tuples]


def dicts_to_map(dicts, key, cls=dict):
    """Create a mapping of ``m[d[key]] = d`` for a list of dicts."""
    return cls((d[key], d) for d in dicts)
