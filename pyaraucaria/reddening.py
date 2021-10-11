import os
import numpy as np
import logging


class _Reddening(object):

    def __init__(self, *args, **kwargs):
        self.table = np.empty(0)
        self._read_database()

    def get_filename(self):
        raise NotImplementedError('Use derived class which should implement abstract `get_filename()`')

    def _read_database(self):
        path = os.path.join(os.path.dirname(__file__), 'databases', self.get_filename())
        with open(path) as fd:
            table = []
            for l in fd:
                if l[0] == '#':
                    continue
                ra, dec, EBV, EBV_e, _, _ = (float(v) for v in l.split())
                table.append([ra, dec, EBV, EBV_e])
        self.table = np.array(table)

    def lookup(self, ra, dec):
        # type: (float, float) -> (float, float, float, float)
        """
        Lookup reddening table
        Parameters
        ----------
        ra: float
            RA in deg, if needed convert using `pyaraucaria.coordinates.ra_to_decimal`
        dec
            DEC in deg, if needed convert using `pyaraucaria.coordinates.dec_to_decimal`

        Returns
        -------
            Tuple (E(B-V), E(B-V) error, ra, dec), where ra, dec are exact coordinates of table point
        """
        # Table is sorted by second column so bsearch first
        logging.debug('SEARCH    {:6f} {:6f}'.format(ra, dec))
        N = self.table.shape[0]
        i = N // 2
        start, stop = 0, N
        while True:
            logging.debug('>>>>> {:3d} {:6f} {:6f}'.format(i, self.table[i][0], self.table[i][1]))
            if self.table[i][1] < dec:
                start = i
                ni = i + (stop - i) // 2
            else:
                stop = i
                ni = i - (i - start) // 2
            if ni == i or abs(self.table[i][1] - self.table[ni][1]) < 1e-4:
                break
            else:
                i = ni
        # search up
        start = i
        s = 2
        while start > 0 and s > 0:
            if abs(self.table[start][1] - self.table[start-1][1]) > 1e-4:  # float equality
                s -= 1 # step detected
            start -= 1
        # search down
        stop = i
        s = 2
        while stop < N-1 and s > 0:
            if abs(self.table[stop][1] - self.table[stop+1][1]) > 1e-4:  # float equality
                s -= 1 # step detected
            stop += 1
        logging.debug('START {:3d} {:6f} {:6f}'.format(start, self.table[start][0], self.table[start][1]))
        logging.debug('HIT   {:3d} {:6f} {:6f}'.format(i, self.table[i][0], self.table[i][1]))
        logging.debug('STOP  {:3d} {:6f} {:6f}'.format(stop, self.table[stop][0], self.table[stop][1]))
        cos_dec = np.cos(np.pi*dec/180)
        best = start
        error = 180.0**2
        for j in range(start, stop +1):
            jerror = ((ra - self.table[j][0]) * cos_dec)**2 + (dec - self.table[j][1])**2
            logging.debug('ITER  {:3d} {:6f} {:6f} {:6f} {:6f}'.format(j, self.table[stop][0], self.table[stop][1], jerror, error))
            if error > jerror:
                error = jerror
                best = j
        logging.debug('BEST  {:3d} {:6f} {:6f}'.format(best, self.table[best][0], self.table[best][1]))
        return self.table[best][2], self.table[best][3], self.table[best][1], self.table[best][2]



class ReddeningLMC(_Reddening):
    def get_filename(self):
        return "lmc_res5_radius5.txt"
