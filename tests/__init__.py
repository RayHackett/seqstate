import unittest

from . import test_featureizer, test_hmms, test_states

"""
To run the tests:

export PYTHONPATH='/home/ray/Documents/seqstate/src/'
python -m unittest discover -s /home/ray/Documents/seqstate/tests/

"""


def load_tests(
    loader: unittest.TestLoader, suite: unittest.TestSuite, pattern: str = "test*.py"
):
    suite.addTests(loader.loadTestsFromModule(test_hmms, pattern=pattern))
    suite.addTests(loader.loadTestsFromModule(test_featureizer, pattern=pattern))
    suite.addTests(loader.loadTestsFromModule(test_states, pattern=pattern))

    return suite
