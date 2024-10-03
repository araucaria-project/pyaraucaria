import unittest
from pyaraucaria.airmass import airmass

class TestAirmassCalculator(unittest.TestCase):
    def test_airmass(self):
        test_cases = [30.0, 45.0, 60.0, 75.0, 89.0]

        previous_result = None

        for elevation in test_cases:
            try:
                result = airmass(elevation)

                #asserting that the calculated airmass is within reasonable bounds
                self.assertGreaterEqual(result, 1.0, f"Expected airmass to be greater than or equal to 1.0. Failed for elevation: {elevation:.2f} degrees")
                self.assertLessEqual(result, 100.0, f"Expected airmass to be less than or equal to 100.0. Failed for elevation: {elevation:.2f} degrees")

                #additional assertion to validate monotonic increase (or decrease) in airmass
                if previous_result is not None:
                    self.assertLessEqual(result, previous_result, f"Expected airmass to be increasing. Failed for elevation: {elevation:.2f} degrees")

                previous_result = result

            except ValueError as ve:
                self.fail(f"Error for elevation: {elevation:.2f} degrees | {ve}")

if __name__ == "__main__":
    unittest.main()
