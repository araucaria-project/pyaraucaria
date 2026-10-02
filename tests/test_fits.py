import unittest
from pyaraucaria.fits import fits_header, fits_stat


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

class TestFitsHeader(unittest.TestCase):

    def test_fits_header_defaults_include_m1_cards(self):
        header = fits_header()

        self.assertEqual(header["OCASTD"][0], "1.1.3")
        self.assertEqual(header["M1-POS1"], ('', '[um] M1 cell motor 1 position'))
        self.assertEqual(header["M1-POS2"], ('', '[um] M1 cell motor 2 position'))
        self.assertEqual(header["M1-POS3"], ('', '[um] M1 cell motor 3 position'))
        self.assertEqual(header["M1-ATSET"], ('', 'M1 cell motors at commanded position'))

    def test_fits_header_m1_cards_accept_values(self):
        header = fits_header(
            m1_pos1=-96.3,
            m1_pos2=0.0,
            m1_pos3=15.2,
            m1_atset=True,
        )

        self.assertEqual(header["M1-POS1"][0], -96.3)
        self.assertEqual(header["M1-POS2"][0], 0.0)
        self.assertEqual(header["M1-POS3"][0], 15.2)
        self.assertTrue(header["M1-ATSET"][0])


if __name__ == '__main__':
    unittest.main()