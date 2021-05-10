import unittest
from Cluster import Cluster

class ClusterTests(unittest.TestCase):
    def test_blank(self):
        c = Cluster()
        self.assertEqual(c.name, None, 'Name should be None when first initialized')
        self.assertEqual(len(c.members), 0, 'Members should be initialized as an empty list')
        self.assertEqual(len(c.members), len(c.types), 'length of members and types should always be equal')


if __name__ == '__main__':
    unittest.main()