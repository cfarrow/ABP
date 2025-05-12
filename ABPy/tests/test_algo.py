import pytest

from ABPy import algo, lattice


def test_remove_active_site():
    l = lattice.SquareLattice(8, 2)

    # No active sites
    with pytest.raises(ValueError):
        algo.remove_active_site(l)

    l.activateSites()
    sites = algo.remove_active_site(l)
    assert len(sites) == 1
    assert l.getNumActive() == len(l) - 1
    assert not l.isActive(sites[0])


def test_cull_sites():
    l = lattice.SquareLattice(8, 2)
    l.activateSites()

    l.setActiveLevel(0, False)
    sites = [0]
    culled_sites = algo.cull_sites(l, sites, 3)
    assert sites == culled_sites

    l.setActiveLevel(16, False)
    # Site 8 now has 2 active neighbors. It is safe from culling sites with
    # fewer than 2 active neighbors, but not 3.
    sites = [16]
    culled_sites = algo.cull_sites(l, sites, 2)
    assert [16] == culled_sites
    culled_sites = algo.cull_sites(l, sites, 3)
    assert [16, 8] == culled_sites
