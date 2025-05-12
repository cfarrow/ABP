""" Procedures for ABP. """
from random import randint

from .lattice import LatticeBase


def remove_active_site(l : LatticeBase):
    """ Mark an active site as inactive. 

    Parameters
    ----------
    l : LatticeBase
        The lattice to operate upon.

    Raises 
    ------
    ValueError 
        If there are no active sites.

    Returns
    -------
    sites : list[int]
        Sites marked as inactive (there is only one)
    
    
    """
    sites = []
    if l.getNumActive() == 0:
        raise ValueError("No sites to deactivate.")

    n = len(l)
    while len(sites) == 0:
        s = randint(0, n-1)
        if l.isActive(s):
            l.setActiveLevel(s, False)
            sites.append(s)
    return sites


def cull_sites(l : LatticeBase, sites : list[int], mcull : int):
    """ Iterative deactivate sites with fewer than mcull active neighbors.

    For efficiency, the algorithm must be seeded with list of inactive sites.
    This is typically used right after `remove_active_site` and passed the
    return value from that call to iteratively cull sites connected to that seed
    site.

    The algorith is iterative. If any sites are removed by culling, their
    neighbors will be checked agains the culling condition and culled if they do
    not meet it. This continues until no additional sites can be removed.

    Parameters
    ----------
    l : LatticeBase
        The lattice to operate upon
    sites : list[int]
        Sites to seed the culling.
    mcull : int
        The culling parameter. Sites with fewer than this number of active
        neighbors will be set as inactive.

    Returns
    -------
    culled_sites : list[int]
        Sites removed by culling, including the seed sites.

    """

    culled_sites = []
    sites = list(sites)  # Don't overwrite

    while len(sites) > 0:

        s = sites.pop(0)
        culled_sites += [s]

        # Check neighbors and prepare to cull them if they don't meet the
        # culling condition.
        for s2 in l.iter_nbrs(s, active=True):
            if l.getNumActiveNeighbors(s2) < mcull:
                l.setActiveLevel(s2, False)
                sites.append(s2)

    return culled_sites
