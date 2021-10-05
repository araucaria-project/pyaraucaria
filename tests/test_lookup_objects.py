from unittest import TestCase
from pyaraucaria.lookup_objects import ObjectsDatabase, lookup_objects


class TestObjectsDatabase(TestCase):
    def test_resolve_alias(self):
        obdb = ObjectsDatabase(skip_errors=False)
        self.assertEqual(obdb.resolve_alias('v1647-sgr'), 'HIP77')
        self.assertEqual(obdb.resolve_alias('foo'), 'foo')

    def test_lookup_objects(self):
        obdb = ObjectsDatabase(skip_errors=False)

        ot = obdb.lookup_object('v1647-sgr')  # objects database and TAB.ALL
        ox = obdb.lookup_object('SCU-S-1')  # objects database only
        xt = obdb.lookup_object('T_Vul')  # TAB.ALL only
        xc = obdb.lookup_object('t-vul')  # TAB.ALL only with canonization
        xx = obdb.lookup_object('foo')  # nowhere
        self.assertEqual(ot['name'], 'HIP77')
        self.assertEqual(ox['name'], 'Scu01')
        self.assertEqual(xt, xc)
        self.assertEqual(xt['name'], 'T_Vul')
        self.assertIsNone(xx)
        self.assertIsInstance(xt['pa'], list)

        oo = obdb.lookup_objects(['v1647-sgr', 'SCU-S-1', 'T_Vul', 't-vul', 'foo'])
        self.assertEqual(oo['v1647-sgr']['name'], 'HIP77')
        self.assertEqual(oo['SCU-S-1']['name'], 'Scu01')
        self.assertEqual(oo['T_Vul'], oo['t-vul'])
        self.assertEqual(oo['T_Vul']['name'], 'T_Vul')
        self.assertIsNone(oo['foo'])
        self.assertIsInstance(oo['T_Vul']['pa'], list)

        # check lookup_objects function against method
        self.assertEqual(oo, lookup_objects(['v1647-sgr', 'SCU-S-1', 'T_Vul', 't-vul', 'foo']))
        pass  # just for breakpoints

    def test_iterate_objects(self):
        obdb = ObjectsDatabase(skip_errors=False)
        for k, v in obdb.all_objects.items():
            self.assertEqual(v['name'], k)
            self.assertIsNotNone(v['hname'])

    def test_all_objects_count(self):
        obdb = ObjectsDatabase(skip_errors=False)
        all_no = len(obdb.all_objects)
        tab_no = len(obdb.tab_all_objects_mapped)
        ob_no = len(obdb.objects_database_objects)
        self.assertGreater(all_no, tab_no)
        self.assertGreater(all_no, ob_no)
        self.assertLess(all_no, tab_no + ob_no)


    def test_lookup_group(self):
        obdb = ObjectsDatabase(skip_errors=False)
        g = obdb.lookup_group('mwcepint', include_members=False)
        self.assertRaises(KeyError, lambda: g['objects'])

        g = obdb.lookup_group('mwcepint')
        self.assertGreater(len(g['objects']), 0)
