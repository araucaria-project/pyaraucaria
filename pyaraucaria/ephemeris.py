from astropy import units as u
from astropy.time import Time
from astroplan import Observer
import datetime
from astropy.coordinates import SkyCoord, FK5, get_body


def calculate_sun_rise_set(date: datetime, horiz_height: float, sunrise: bool,
                           latitude: float, longitude: float, elevation: float):
    """
    Calculate next sunrise or sunset at horizon height
    :param date: utc date of start calculating
    :param horiz_height: the height over (or under) horizon
    :param sunrise: if true the sunrise will be calculated, else sunset
    :param latitude: latitude of observer
    :param longitude: longitude of observer
    :param elevation: elevation of observer
    :return: utc time of next sunrise / sunset
    """
    date = Time(val=date)
    obs = Observer(latitude=latitude,
                   longitude=longitude,
                   elevation=elevation * u.m)
    if sunrise:
        return obs.sun_rise_time(date, which='next', horizon=horiz_height * u.deg).to_datetime()
    else:
        return obs.sun_set_time(date, which='next', horizon=horiz_height * u.deg).to_datetime()


def moon_separation(ra: float, dec: float, utc_time: Time):
    """
    Func. returns moon separation
    :param ra: ra in float
    :param dec: dec in float
    :param utc_time: time of calculation in astropy Time format
    :return: Separation in deg
    """
    moon = get_body(body="moon", time=utc_time)
    obj_coo = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame=FK5(equinox=utc_time))
    return float(moon.separation(obj_coo, origin_mismatch='ignore').to(u.deg).deg)
