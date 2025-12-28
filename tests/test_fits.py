import unittest
from pyaraucaria.fits import fits_stat


class TestFitsStat(unittest.TestCase):

    def test_fits_stat(self):
        array = [2, 5, 8, 12, 15]
        # Expected values
        pred_dict = {
            'min': 2.0,
            'max': 15.0,
            'median': 8.0,
            'mean': 8.4,
            'std': 4.673328578219169,
            'rms': 4.673328578219169,
            'sigma_quantile': 4.092
        }
        result_dict = fits_stat(array)

        self.assertEqual(pred_dict.keys(), result_dict.keys())

        for key, expected_val in pred_dict.items():
            self.assertAlmostEqual(
                result_dict[key],
                expected_val,
                places=5,  # Check up to 5 decimal places
                msg=f"Value mismatch for key: '{key}'"
            )



if __name__ == '__main__':
    unittest.main()