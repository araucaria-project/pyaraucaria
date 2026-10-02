import inspect
import unittest

from astropy.io import fits

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

    def test_fits_header_defaults_include_dome_cards(self):
        header = fits_header()

        self.assertEqual(header["T-DOME"], ('', '[deg C] Temperature - dome'))
        self.assertEqual(header["RHUM-DOM"], ('', '[%] Relative humidity - dome'))

    def test_fits_header_dome_cards_accept_values(self):
        header = fits_header(
            temp_dome=12.5,
            rhum_dome=34.0,
        )

        self.assertEqual(header["T-DOME"][0], 12.5)
        self.assertEqual(header["RHUM-DOM"][0], 34.0)

    def test_fits_header_is_valid_fits(self):
        """All keywords fit FITS standard (max 8 chars) and header can be built by astropy."""
        header = fits_header()

        for key in header:
            self.assertLessEqual(len(key), 8, msg=f"Keyword too long: '{key}'")

        fits_hdr = fits.Header([(key, value, comment) for key, (value, comment) in header.items()])
        self.assertEqual(len(fits_hdr), len(header))

    def test_fits_header_every_param_maps_to_one_card(self):
        """Each fits_header parameter fills exactly one card, each card has a comment."""
        params = list(inspect.signature(fits_header).parameters)
        header = fits_header(**{p: f'value_{p}' for p in params})

        self.assertEqual(len(header), len(params))

        values = [value for value, comment in header.values()]
        for p in params:
            self.assertEqual(values.count(f'value_{p}'), 1, msg=f"Parameter '{p}' not mapped to exactly one card")

        for key, (value, comment) in header.items():
            self.assertTrue(comment, msg=f"Empty comment for card '{key}'")


if __name__ == '__main__':
    unittest.main()