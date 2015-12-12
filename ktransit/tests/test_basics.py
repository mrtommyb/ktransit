"""Basic sanity checks to verify that ktransit works.
To run, simply type "py.test".
"""


def test_import():
    """Can we import ktransit successfully?"""
    import ktransit


def test_planets():
    import ktransit

    M = ktransit.LCModel()
    M.add_star()
    M.add_planet()
    M.add_data()

    tmod = M.transitmodel 
