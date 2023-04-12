import unittest
from pyaraucaria.fits import fits_stat


class TestFitsStat(unittest.TestCase):

    def test_fits_stat(self):
        array = [2, 5, 8, 12, 15]
        pred_dict = {'min': 2, 'max': 15, 'median': 8.0, 'mean': 8.4, 'std': 4.673328578219169}
        result_dict = fits_stat(array)
        self.assertDictEqual(pred_dict, result_dict)


if __name__ == '__main__':
    unittest.main()