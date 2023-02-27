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


if __name__ == '__main__':
    unittest.main()
