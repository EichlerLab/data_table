"""
Commonly used definitions.
"""

# Expand svtypes that may contain one or more sub-types into to individual svtypes
SVTYPE_EXPAND = {
    'ins': ('ins',),
    'del': ('del',),
    'inv': ('inv',),
    'snv': ('snv',),
    'insdel': ('ins', 'del'),
    'all': ('ins', 'del', 'inv')
}