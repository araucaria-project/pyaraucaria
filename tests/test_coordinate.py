import unittest

from pyaraucaria.coordinates import *

# >> > parse_deg_or_dms('12:34:56.1')
# 12.58225
# >> > parse_deg_or_dms('12°34′56.1″')
# 12.58225
# >> > parse_deg_or_dms('-12 34 56.1')
# -12.58225
# >> > parse_deg_or_dms('12.58225')


class TestCoo(unittest.TestCase):

    def test_dms_to_float_deg(self):
        tst = [
            ('12°34′56.1″', 12.58225),
            ('12:34:56.1', 12.58225),
            ('-12 34 56.1', -12.58225),
            ('12.58225', 12.58225),
            ]
        res = [deg_to_decimal_deg(d[0]) for d in tst]
        for (_, t), r in zip(tst, res):
            self.assertAlmostEqual(t, r)

class TestRaDecEpoch(unittest.TestCase):

    def test_ra_dec_epoch_1(self):

        RUD = {'latitude': '51:17:05', 'longitude': '22:35:50', 'elevation': '235'}
        ra = '5:14:32.30'
        dec = '-8:12:06'
        out_ra = '5:15:38.46'
        out_dec = '-8:10:38.83'
        epoch = ra_dec_epoch(ra, dec, RUD['latitude'], RUD['longitude'], RUD['elevation'], '2000/01/01')
        e_ra = ra_to_decimal(out_ra)
        e_dec = dec_to_decimal(out_dec)
        out_epoch_ra = ra_to_decimal(epoch[0])
        out_epoch_dec = dec_to_decimal(epoch[1])
        self.assertAlmostEqual(e_ra, out_epoch_ra, places=3)
        self.assertAlmostEqual(e_dec, out_epoch_dec, places=3)


if __name__ == '__main__':
    unittest.main()
