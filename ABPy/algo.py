""" Procedures for ABP. """
from random import randint


def remove_active_site(l):
    """ Mark an active site as inactive. Returns list of deactivated sites. 

    Raises ValueError if there are no sites to remove.
    
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


def cull_sites(l, sites : list[int], mcull : int):
    """ Iterative deactivate sites with fewer than mcull active neighbors.

    For efficiency, the algorithm must be seeded with list of inactive sites.
    This is typically used right after `remove_active_site` and passed the
    return value from that call to iteratively cull sites connected to that seed
    site.

    Parameters
    ----------
    l : Lattice
    sites : list of integer
    mcull : the culling parameter

    Returns
    -------
    culled_sites : list of integer
        Sites removed by the algorithm, including the seed sites.

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
