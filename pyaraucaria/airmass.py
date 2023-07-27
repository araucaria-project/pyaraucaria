import math

class AirmassCalculator:
    def airmass(elevation: float) -> float:
        """
        Airmass is a measure of the path length of sunlight through the Earth's atmosphere.
        It is defined as the secant of the zenith angle and is used in various atmospheric
        and astronomical calculations.

        Parameters:
            elevation (float): The elevation angle of the sun in DEGREES. Should be between 0 and 90.

        Returns:
            float: The calculated airmass.

        Source:
            The formula used in this method is based on Kasten and Young's airmass model.
            Kasten, F. and Young, A.T. (1989) Revised optical air mass tables and approximation formula.
        """
        if not 0 <= elevation <= 90:
            raise ValueError("Elevation angle must be between 0 and 90 degrees.")

        #special case for the zenith angle (elevation of 0.0 degrees)
        if elevation == 0.0:
            return 1.0

        #special case for the zenith angle (elevation of 90.0 degrees)
        if elevation == 90.0:
            return 1.0

        #convert the elevation angle to radians correctly
        elevation_rad = math.radians(elevation)

        #calculate the zenith angle
        zenith_angle = 90 - elevation

        #calculate the airmass using the Kasten-Young formula
        airmass = 1 / (math.cos(elevation_rad) + 0.50572 * pow((96.07995 - zenith_angle) * 180 / math.pi, -1.6364))
        return airmass
