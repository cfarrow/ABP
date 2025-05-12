from ABPy import lattice


def test_iter_nbrs():
    l = lattice.SquareLattice(8, 2)
    l.activateSites()

    expected = [1, 7, 8, 56]
    assert expected == list(l.iter_nbrs(0))
    assert expected == list(l.iter_nbrs(0, active=True))
    assert expected == list(l.iter_nbrs(0, present=True))
    assert expected == list(l.iter_nbrs(0, active=True, present=True))

    l.setActiveLevel(8, False)
    assert expected == list(l.iter_nbrs(0))
    expected = [1, 7, 56]
    assert expected == list(l.iter_nbrs(0, active=True))
    expected = [8]
    assert expected == list(l.iter_nbrs(0, active=False))

    l.setPresentLevel(1, False)
    expected = [1, 7, 8, 56]
    assert expected == list(l.iter_nbrs(0))
    expected = [7, 8, 56]
    assert expected == list(l.iter_nbrs(0, present=True))
    expected = [1]
    assert expected == list(l.iter_nbrs(0, present=False))

    expected = [7, 56]
    assert expected == list(l.iter_nbrs(0, present=True, active=True))
    expected = []
    assert expected == list(l.iter_nbrs(0, present=False, active=False))

    l.setPresentLevel(8, False)
    expected = [8]
    assert expected == list(l.iter_nbrs(0, present=False, active=False))


def test_iter_nbrs():
    from collections import Counter
    l = lattice.SquareLattice(8, 2)
    l.activateSites()

    edges = list(l.iter_edges())
    assert len(edges) == 4 * 8 * 8 / 2  # 4 neighbors per point with double counting
    left, right = zip(*edges)
    c = Counter(left + right)
    assert list(sorted(c)) == list(range(64))
    assert all(v == 4 for v in c.values())

    s1 = 63
    l.setActiveLevel(s1, False)
    edges = list(l.iter_edges())
    assert len(edges) == 4 * 8 * 8 / 2
    left, right = zip(*edges)
    c = Counter(left + right)
    assert list(sorted(c)) == list(range(64))
    assert all(v == 4 for v in c.values())

    edges = list(l.iter_edges(active=True))
    assert len(edges) == 4 * 8 * 8 / 2 - 4
    left, right = zip(*edges)
    c = Counter(left + right)
    nbrs = set(l.iter_nbrs(s1))
    for i, v in c.items():
        if i in nbrs:
            assert v == 3
        else:
            assert v == 4


    s2 = 21
    l.setPresentLevel(s2, False)
    edges = list(l.iter_edges())
    assert len(edges) == 4 * 8 * 8 / 2
    left, right = zip(*edges)
    c = Counter(left + right)
    assert list(sorted(c)) == list(range(64))
    assert all(v == 4 for v in c.values())

    edges = list(l.iter_edges(present=True))
    assert len(edges) == 4 * 8 * 8 / 2 - 4
    left, right = zip(*edges)
    c = Counter(left + right)
    nbrs = set(l.iter_nbrs(s2))
    for i, v in c.items():
        if i in nbrs:
            assert v == 3
        else:
            assert v == 4

    edges = list(l.iter_edges(active=True, present=True))
    assert len(edges) == 4 * 8 * 8 / 2 - 8
    left, right = zip(*edges)
    c = Counter(left + right)
    nbrs = set(l.iter_nbrs(s1)) | set(l.iter_nbrs(s2))
    for i, v in c.items():
        if i in nbrs:
            assert v == 3
        else:
            assert v == 4