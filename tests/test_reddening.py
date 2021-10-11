import unittest


class ReddeningCase(unittest.TestCase):
    def test_LMC_Reddening(self):
        from pyaraucaria.reddening import ReddeningLMC
        r = ReddeningLMC()
        EBV, _, _, _ = r.lookup(74.634958, -66.974056)
        self.assertAlmostEqual(EBV, 0.0995)


if __name__ == '__main__':
    unittest.main()
