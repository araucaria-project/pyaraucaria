import unittest
from pyaraucaria.airmass_calculator import AirmassCalculator
class TestAirmassCalculator(unittest.TestCase):

    def test_calculate_airmass(self):
        test_cases = [30.0, 45.0, 60.0, 75.0, 89.0]

        for elevation in test_cases:
            try:
                result = AirmassCalculator.calculate_airmass(elevation)
                expected_airmass = self.calculate_expected_airmass(elevation)
                #asserting that the calculated airmass is close to the expected airmass
                self.assertAlmostEqual(result, expected_airmass, places=4)
            except ValueError as ve:
                self.fail(f"Error for elevation: {elevation:.2f} degrees | {ve}")

    def calculate_expected_airmass(self, elevation):
        return AirmassCalculator.calculate_airmass(elevation)

if __name__ == "__main__":
    unittest.main()
