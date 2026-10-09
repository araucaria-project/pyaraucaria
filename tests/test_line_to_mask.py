"""Tests for FFS.line_to_mask: (rho, theta) and segment ends -> image pixels."""

import unittest

import numpy as np

from pyaraucaria.ffs import FFS


class TestLineMask(unittest.TestCase):

    shape = (100, 120)   # (ny, nx)

    def test_horizontal_full_line(self):
        # theta = pi/2: y = rho
        m = FFS.line_to_mask(self.shape, rho=40, theta=np.pi / 2)
        ys, xs = np.nonzero(m)
        self.assertTrue(np.all(ys == 40))
        self.assertEqual(xs.size, self.shape[1])

    def test_vertical_full_line_with_width(self):
        # theta = 0: x = rho; half_width=1 -> columns 29..31
        m = FFS.line_to_mask(self.shape, rho=30, theta=0.0, half_width=1.0)
        self.assertEqual(sorted(set(np.nonzero(m)[1])), [29, 30, 31])
        self.assertEqual(m.sum(), 3 * self.shape[0])

    def test_diagonal_full_line_spans_frame(self):
        # line y = x: x*cos(3pi/4) + y*sin(3pi/4) = 0
        m = FFS.line_to_mask(self.shape, rho=0.0, theta=3 * np.pi / 4)
        n = min(self.shape)
        self.assertTrue(np.all(m[np.arange(n), np.arange(n)]))
        ys, xs = np.nonzero(m)
        self.assertTrue(np.all(np.abs(xs - ys) <= 1))

    def test_segment_is_limited_to_ends(self):
        m = FFS.line_to_mask(self.shape, rho=40, theta=np.pi / 2, half_width=0.5,
                      x0=10, y0=40, x1=50, y1=40)
        xs = np.nonzero(m)[1]
        self.assertEqual((xs.min(), xs.max()), (10, 50))

    def test_segment_is_subset_of_full_line(self):
        theta, rho = 0.7, 50.0
        c, s = np.cos(theta), np.sin(theta)
        # two points on the line
        x0, x1 = 20.0, 60.0
        y0, y1 = (rho - x0 * c) / s, (rho - x1 * c) / s
        seg = FFS.line_to_mask(self.shape, rho, theta, 2.0, x0, y0, x1, y1)
        full = FFS.line_to_mask(self.shape, rho, theta, 2.0)
        self.assertTrue(seg.any())
        self.assertFalse(np.any(seg & ~full))
        self.assertGreater(full.sum(), seg.sum())

    def test_line_outside_frame_is_empty(self):
        m = FFS.line_to_mask(self.shape, rho=500, theta=0.0)
        self.assertFalse(m.any())
        self.assertEqual(m.shape, self.shape)


if __name__ == "__main__":
    unittest.main()
