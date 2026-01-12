import math
from pyaraucaria.coordinates import site_sidereal_time, ra_dec_2_az_alt
from typing import Dict, Tuple


def dome_eq_azimuth(ra: float, dec: float, r_dome: float, spx: float, spy: float, gem: float,
                    side_of_pier: int, latitude: float, longitude: float,
                    elevation: float, epoch: str = '2000') -> Tuple[float, Dict]:
    """
    Func. calculate dome position for eq mount.
    @param ra: Ra of telescope [deg]
    @param dec: Dec of telescope [deg]
    @param r_dome: Radius of dome [mm]
    @param spx: X distance from dome center to intersection of ra - dec axes: (+) when on E, (-) when on W [mm]
    @param spy: Y distance from dome center to intersection of ra - dec axes: (+) when on N, (-) when on S [mm]
    @param gem: Distance from intersection of ra - dec axes point to telescope optical axis (absolute value) [mm]
    @param side_of_pier: Side of telescope pier position. On east = 0, on west = 1
    @param latitude: Site latitude [deg]
    @param longitude: Site longitude [deg]
    @param elevation: Site elevation [m]
    @param epoch: Epoch for radec, example: '2000'
    @return: Eq mount dome azimuth [deg]
    """
    # TODO do in vectors np

    # Calc. az alt from ra, dec and epoch
    az, alt = ra_dec_2_az_alt(ra=ra,
                              dec=dec,
                              longitude=longitude,
                              latitude=latitude,
                              elevation=elevation,
                              epoch=epoch)

    # dome_azimuth_vector = np.array([
    # math.sin(math.radians(az)) * r_dome,
    # math.cos(math.radians(az)) * r_dome
    # ]).reshape(-1, 1)

    # Calc. dome azimuth vector
    az_rad = math.radians(az)
    dome_azimuth_vector_x = math.sin(az_rad) * r_dome
    dome_azimuth_vector_y = math.cos(az_rad) * r_dome
    dome_azimuth_vector = [dome_azimuth_vector_x, dome_azimuth_vector_y]

    # Calc. alt az target point on xy plane (as vector)

    # target_point = np.array([
    # r_dome * math.cos(math.radians(90 - az))) * math.sin(math.radians(90 - alt)),
    # r_dome * math.sin(math.radians(90 - az)) * math.sin(math.radians(90 - alt))
    # ]).reshape(-1, 1)

    target_point = [r_dome * math.cos(math.radians(-(az - 90))) * math.sin(math.radians(90 - alt)),
                    r_dome * math.sin(math.radians(-(az - 90))) * math.sin(math.radians(90 - alt))]

    # Get sidereal tome for site
    s_sidereal = site_sidereal_time(longitude=longitude,
                                    latitude=latitude,
                                    elevation=elevation)

    # Calc. angle between ra and sidereal + side of pier
    if int(side_of_pier) not in [0, 1]:
        side_of_pier = 0
    # ra_to_sidereal_angle = abs(ra - s_sidereal) + int(side_of_pier) * 180
    ra_to_sidereal_angle = (s_sidereal - ra) + int(side_of_pier) * 180

    # Calc. vector gem
    x = math.cos(math.radians(ra_to_sidereal_angle)) * gem
    y = math.sin(math.radians(ra_to_sidereal_angle)) * gem

    # Location in the northern or Southern Hemisphere
    n_s = - (latitude / abs(latitude))
    vector_gem = [x, n_s * y * math.cos(math.radians(90 - abs(latitude)))]
    vector_spx_spy = [spx, spy]

    # TODO Add vectors
    # new_dome_vector = target_point + vector_spx_spy + vector_gem

    # Calc. new dome vector
    new_dome_vector = [target_point[0] + vector_spx_spy[0] + vector_gem[0],
                       target_point[1] + vector_spx_spy[1] + vector_gem[1]]


    # TODO NEW VER OF FINAL ANGLE CALC
    # new_dome_az = math.degrees(math.acos((v.transpose() @ u).ravel()[0] / (np.linalg.norm(u) * np.linalg.norm(v))))
    # if new_dome_vector[0] < 0:
    #     new_dome_az = 360 - new_dome_az

    # Final azimuth calculation
    if new_dome_vector[0] >= 0 and new_dome_vector[1] >= 0:
        new_dome_az = math.degrees(math.atan(new_dome_vector[0] / new_dome_vector[1]))
    elif new_dome_vector[0] >= 0 > new_dome_vector[1]:
        new_dome_az = math.degrees(math.atan(new_dome_vector[0] / new_dome_vector[1])) + 180
    elif new_dome_vector[0] < 0 and new_dome_vector[1] < 0:
        new_dome_az = math.degrees(math.atan(new_dome_vector[0] / new_dome_vector[1])) + 180
    else:
        new_dome_az = math.degrees(math.atan(new_dome_vector[0] / new_dome_vector[1])) + 360

    # Compose info dict
    info_dict = {
        'vector_final': new_dome_vector,
        'dome_azimuth_vector': dome_azimuth_vector,
        'vector_gem': vector_gem,
        'target_point': target_point,
        's_sidereal': s_sidereal
    }

    return new_dome_az, info_dict
