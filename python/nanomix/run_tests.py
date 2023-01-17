import unittest

test_suite = unittest.defaultTestLoader.discover('./tests', pattern='test*.py')
unittest.TextTestRunner().run(test_suite)

