import math

class AirmassCalculator:
    @staticmethod
    def calculate_airmass(elevation: float) -> float:
        if not 0 <= elevation <= 90:
            raise ValueError("Elevation angle must be between 0 and 90 degrees.")

        zenith_angle = 90 - elevation
        airmass = 1 / (math.cos(math.radians(zenith_angle)) + 0.50572 * (96.07995 - zenith_angle) ** -1.6364)
        return airmass
