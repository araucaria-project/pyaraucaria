import unittest

from pyaraucaria.ob_validator import ObsValidator
from pyaraucaria.obs_plan.obs_plan_parser import ObsPlanParser


FILTERS = ["V", "Ic", "B", "g", "r", "i", "z"]


class TestLoadSchema(unittest.TestCase):

    def test_load_base_schema(self):
        schema = ObsValidator.load_schema("base_schema")
        self.assertIn("properties", schema)
        self.assertIn("seq", schema["properties"])

    def test_extension_is_optional(self):
        self.assertEqual(ObsValidator.load_schema("base_rules"),
                         ObsValidator.load_schema("base_rules.yaml"))

    def test_missing_schema_raises(self):
        with self.assertRaises(FileNotFoundError):
            ObsValidator.load_schema("no_such_schema")


class TestConvertToObdict(unittest.TestCase):

    def test_three_args(self):
        parsed = {'command_name': 'SEQUENCE', 'subcommands': [
            {'command_name': 'OBJECT', 'args': ['FF_Aql', '18:58:14.75', '17:21:39.29'],
             'kwargs': {'seq': '2/Ic/60'}}]}
        self.assertEqual(ObsValidator.convert_to_obdict(parsed),
                         {'command_name': 'OBJECT', 'seq': '2/Ic/60', 'name': 'FF_Aql',
                          'ra': '18:58:14.75', 'dec': '17:21:39.29'})

    def test_one_arg_is_name(self):
        parsed = {'subcommands': [{'command_name': 'DARK', 'args': ['ZZ01'],
                                   'kwargs': {'seq': '10/V/300'}}]}
        self.assertEqual(ObsValidator.convert_to_obdict(parsed),
                         {'command_name': 'DARK', 'name': 'ZZ01', 'seq': '10/V/300'})

    def test_two_args_are_coordinates(self):
        parsed = {'subcommands': [{'command_name': 'OBJECT', 'args': ['12:12:12', '20:20:20']}]}
        self.assertEqual(ObsValidator.convert_to_obdict(parsed),
                         {'command_name': 'OBJECT', 'ra': '12:12:12', 'dec': '20:20:20'})

    def test_single_subcommand_as_dict(self):
        self.assertEqual(ObsValidator.convert_to_obdict({'subcommands': {'command_name': 'BELL'}}),
                         {'command_name': 'BELL'})

    def test_no_subcommands(self):
        self.assertEqual(ObsValidator.convert_to_obdict({}), {})


class TestConvertFromObdict(unittest.TestCase):

    def test_name_and_coordinates(self):
        ob = {'command_name': 'OBJECT', 'name': 'FF_Aql', 'ra': '18:58:14.75',
              'dec': '17:21:39.29', 'seq': '2/Ic/60'}
        self.assertEqual(ObsValidator.convert_from_obdict(ob),
                         'OBJECT FF_Aql 18:58:14.75 17:21:39.29 seq=2/Ic/60')

    def test_name_only(self):
        self.assertEqual(
            ObsValidator.convert_from_obdict({'command_name': 'OBJECT', 'name': 'FF_Aql',
                                              'seq': '2/Ic/60'}),
            'OBJECT FF_Aql seq=2/Ic/60')

    def test_coordinates_only(self):
        self.assertEqual(
            ObsValidator.convert_from_obdict({'command_name': 'OBJECT', 'ra': '1:2:3',
                                              'dec': '4:5:6'}),
            'OBJECT 1:2:3 4:5:6')

    def test_kwargs_only(self):
        self.assertEqual(ObsValidator.convert_from_obdict({'command_name': 'WAIT',
                                                           'ut': '16:00:00'}),
                         'WAIT ut=16:00:00')

    def test_missing_command_name(self):
        self.assertIsNone(ObsValidator.convert_from_obdict({'name': 'FF_Aql'}))

    def test_not_a_dict(self):
        self.assertIsNone(ObsValidator.convert_from_obdict("OBJECT FF_Aql"))

    def test_round_trip_through_parser(self):
        for txt in ['OBJECT FF_Aql 18:58:14.75 17:21:39.29 seq=2/Ic/60,2/V/70',
                    'WAIT ut=16:00:00',
                    'ZERO seq=15/Ic/0']:
            with self.subTest(txt=txt):
                ob = ObsValidator.convert_to_obdict(ObsPlanParser.convert_from_string(txt))
                self.assertEqual(ObsValidator.convert_from_obdict(ob), txt)


class TestConvertTypes(unittest.TestCase):

    def setUp(self):
        self.schema = ObsValidator.load_schema("base_schema")

    def test_numbers_and_integers(self):
        converted = ObsValidator.convert_types({'alt': '60', 'epoch': '2000', 'sec': '12.5'},
                                               self.schema)
        self.assertEqual(converted, {'alt': 60.0, 'epoch': 2000, 'sec': 12.5})

    def test_booleans(self):
        self.assertEqual(ObsValidator.convert_types({'test': 'true'}, self.schema), {'test': True})
        self.assertEqual(ObsValidator.convert_types({'test': '0'}, self.schema), {'test': False})

    def test_string_property(self):
        self.assertEqual(ObsValidator.convert_types({'name': 5}, self.schema), {'name': '5'})

    def test_union_type_picks_first_match(self):
        # read_mod is [integer, string]
        self.assertEqual(ObsValidator.convert_types({'read_mod': '1'}, self.schema),
                         {'read_mod': 1})
        self.assertEqual(ObsValidator.convert_types({'read_mod': 'fast'}, self.schema),
                         {'read_mod': 'fast'})

    def test_unconvertible_value_is_left_intact(self):
        self.assertEqual(ObsValidator.convert_types({'alt': 'abc', 'test': 'maybe'}, self.schema),
                         {'alt': 'abc', 'test': 'maybe'})

    def test_unknown_key_passes_through(self):
        self.assertEqual(ObsValidator.convert_types({'unknown': 'x'}, self.schema),
                         {'unknown': 'x'})


class TestCleanNone(unittest.TestCase):

    def test_drops_none_values(self):
        self.assertEqual(ObsValidator.clean_none({'a': 1, 'b': None, 'c': False}),
                         {'a': 1, 'c': False})


class TestValidateRules(unittest.TestCase):

    def setUp(self):
        self.rules = ObsValidator.load_schema("base_rules")

    def test_missing_required_reported_as_none(self):
        self.assertEqual(ObsValidator.validate_rules({'command_name': 'ZERO'}, self.rules['ZERO']),
                         {'seq': None})

    def test_complete_exclusive_group_is_ok(self):
        self.assertEqual(
            ObsValidator.validate_rules({'command_name': 'OBJECT', 'ra': '1:2:3', 'dec': '4:5:6'},
                                        self.rules['OBJECT']),
            {})

    def test_partial_exclusive_group_fails(self):
        self.assertEqual(
            ObsValidator.validate_rules({'command_name': 'OBJECT', 'ra': '1:2:3'},
                                        self.rules['OBJECT']),
            {'ra': False, 'dec': False})

    def test_two_complete_exclusive_groups_fail(self):
        self.assertEqual(
            ObsValidator.validate_rules({'command_name': 'OBJECT', 'ra': '1:2:3', 'dec': '4:5:6',
                                         'alt': 34, 'az': 270}, self.rules['OBJECT']),
            {'ra': False, 'dec': False, 'alt': False, 'az': False})

    def test_one_of_group_missing(self):
        self.assertEqual(ObsValidator.validate_rules({'command_name': 'WAIT'}, self.rules['WAIT']),
                         {'sec': False, 'ut': False, 'sunrise': False, 'sunset': False})

    def test_one_of_group_satisfied(self):
        self.assertEqual(
            ObsValidator.validate_rules({'command_name': 'WAIT', 'ut': '16:00:00'},
                                        self.rules['WAIT']),
            {})


class TestValidateSeq(unittest.TestCase):

    def _seq(self, seq, command_name='OBJECT', allowed_filters=None):
        return ObsValidator.validate_seq({'command_name': command_name, 'seq': seq},
                                         allowed_filters=allowed_filters)

    def test_single_exposure(self):
        self.assertEqual(self._seq('1/V/300'), {'seq': True})

    def test_multiple_exposures(self):
        self.assertEqual(self._seq('2/V/30,3/Ic/40'), {'seq': True})

    def test_multiplier(self):
        self.assertEqual(self._seq('2x(1/V/1,1/r/1)'), {'seq': True})

    def test_zero_multiplier_fails(self):
        self.assertEqual(self._seq('0x(1/V/1)'), {'seq': False})

    def test_missing_seq_is_not_reported(self):
        self.assertEqual(ObsValidator.validate_seq({'command_name': 'BELL'}), {})

    def test_malformed_seq(self):
        for seq in ['1/V', '1/V/abc', '1/V/-3', '0/V/1', 'dupa']:
            with self.subTest(seq=seq):
                self.assertEqual(self._seq(seq), {'seq': False})

    def test_auto_exposure_only_for_skyflat(self):
        self.assertEqual(self._seq('10/V/a', command_name='SKYFLAT'), {'seq': True})
        self.assertEqual(self._seq('10/V/a', command_name='OBJECT'), {'seq': False})

    def test_allowed_filters(self):
        self.assertEqual(self._seq('1/V/300', allowed_filters=FILTERS), {'seq': True})
        self.assertEqual(self._seq('1/dupa/300', allowed_filters=FILTERS), {'seq': False})


class TestValidateOb(unittest.TestCase):

    def setUp(self):
        self.validator = ObsValidator(ObsValidator.load_schema("base_schema"),
                                      ObsValidator.load_schema("base_rules"))

    def test_valid_object(self):
        result = self.validator.validate_ob(
            {'command_name': 'OBJECT', 'name': 'HD193901', 'ra': '20:23:35.8',
             'dec': '-21:22:14.0', 'seq': '1/V/300'}, allowed_filters=FILTERS)
        self.assertTrue(result['valid'])
        self.assertTrue(all(result['result'].values()))

    def test_none_values_are_dropped_from_data(self):
        result = self.validator.validate_ob({'command_name': 'OBJECT', 'name': None,
                                             'seq': '1/V/1'})
        self.assertEqual(result['data'], {'command_name': 'OBJECT', 'seq': '1/V/1'})
        self.assertNotIn('name', result['result'])

    def test_types_are_converted_in_data(self):
        result = self.validator.validate_ob({'command_name': 'OBJECT', 'alt': '60', 'az': '270',
                                             'seq': '1/V/1'})
        self.assertEqual(result['data']['alt'], 60.0)
        self.assertEqual(result['data']['az'], 270.0)

    def test_schema_violation_marks_key_false(self):
        result = self.validator.validate_ob({'command_name': 'OBJECT', 'name': 'X',
                                             'seq': '1/V/1', 'gain': 'ultra'})
        self.assertFalse(result['valid'])
        self.assertFalse(result['result']['gain'])

    def test_unknown_command_name(self):
        result = self.validator.validate_ob({'command_name': 'NOPE'})
        self.assertFalse(result['valid'])
        self.assertFalse(result['result']['command_name'])
        self.assertEqual(result['required'], [])

    def test_required_and_allowed_are_reported(self):
        result = self.validator.validate_ob({'command_name': 'ZERO', 'seq': '1/V/1'})
        self.assertEqual(result['required'], ['command_name', 'seq'])
        self.assertIn('mirror_cover', result['allowed'])

    def test_overrides_relax_existing_property(self):
        ob = {'command_name': 'OBJECT', 'name': 'X', 'seq': '1/V/1', 'gain': 'ultra'}
        result = self.validator.validate_ob(
            ob, overrides={'gain': {'type': 'string', 'enum': ['low', 'high', 'ultra']}})
        self.assertTrue(result['valid'])

    def test_overrides_add_new_property(self):
        ob = {'command_name': 'OBJECT', 'name': 'X', 'seq': '1/V/1', 'priority': '5'}
        result = self.validator.validate_ob(ob, overrides={'priority': {'type': 'number',
                                                                        'minimum': 0}})
        self.assertTrue(result['valid'])
        self.assertEqual(result['data']['priority'], 5.0)

        result = self.validator.validate_ob(ob, overrides={'priority': {'type': 'number',
                                                                        'minimum': 10}})
        self.assertFalse(result['result']['priority'])

    def test_overrides_do_not_mutate_base_schema(self):
        self.validator.validate_ob({'command_name': 'BELL'},
                                   overrides={'priority': {'type': 'number'}})
        self.assertNotIn('priority', self.validator.base_schema['properties'])


class TestTpgRules(unittest.TestCase):
    """tpg master-file rules: stricter than base_rules, see tpg_rules.yaml."""

    def setUp(self):
        schema = ObsValidator.load_schema("base_schema")
        schema["properties"].update(ObsValidator.load_schema("tpg_schema")["properties"])
        self.validator = ObsValidator(schema, ObsValidator.load_schema("tpg_rules"))
        self.ob = {'command_name': 'OBJECT', 'name': 'AP_Ser', 'ra': '15:14:00.92',
                   'dec': '+09:58:51.8', 'seq': '2/Ic/6', 'priority': 5}

    def test_complete_master_line(self):
        self.assertTrue(self.validator.validate_ob(self.ob)['valid'])

    def test_scheduling_keys_are_typed_by_tpg_schema(self):
        ob = dict(self.ob, cycle='0.08', ph_start='0.1', ph_end='0.3', P='7.5')
        result = self.validator.validate_ob(ob)
        self.assertTrue(result['valid'])
        self.assertEqual(result['data']['cycle'], 0.08)
        self.assertEqual(result['data']['ph_start'], 0.1)

    def test_missing_priority_is_rejected(self):
        # base_rules accepts this line, tpg crashes on it in allocate()
        ob = dict(self.ob)
        del ob['priority']
        result = self.validator.validate_ob(ob)
        self.assertFalse(result['valid'])
        self.assertIsNone(result['result']['priority'])

    def test_missing_coordinates_are_rejected(self):
        for key in ('name', 'ra', 'dec'):
            with self.subTest(key=key):
                ob = dict(self.ob)
                del ob[key]
                self.assertFalse(self.validator.validate_ob(ob)['valid'])

    def test_seq_or_ob_time_is_required(self):
        ob = dict(self.ob)
        del ob['seq']
        result = self.validator.validate_ob(ob)
        self.assertFalse(result['valid'])

        ob['ob_time'] = 12
        self.assertTrue(self.validator.validate_ob(ob)['valid'])

    def test_non_object_command_is_rejected(self):
        result = self.validator.validate_ob({'command_name': 'ZERO', 'seq': '15/Ic/0'})
        self.assertFalse(result['valid'])
        self.assertFalse(result['result']['command_name'])

    def test_base_rules_accept_what_tpg_rules_reject(self):
        # documents why tpg needs its own rules at all
        base = ObsValidator(ObsValidator.load_schema("base_schema"),
                            ObsValidator.load_schema("base_rules"))
        ob = {'command_name': 'OBJECT', 'name': 'AP_Ser', 'seq': '2/Ic/6'}
        self.assertTrue(base.validate_ob(ob)['valid'])
        self.assertFalse(self.validator.validate_ob(ob)['valid'])


class TestValidateTxt(unittest.TestCase):

    def setUp(self):
        self.validator = ObsValidator(ObsValidator.load_schema("base_schema"),
                                      ObsValidator.load_schema("base_rules"))

    def _validate(self, txt):
        return self.validator.validate_txt(txt, allowed_filters=FILTERS)

    def test_valid_commands(self):
        for txt in ['OBJECT HD193901 20:23:35.8 -21:22:14.0 seq=1/V/300',
                    'OBJECT HD193901 12:23:35.8 -63:45:00.54 seq=10x(1/Ic/300,3/V/10)',
                    'SKYFLAT HD24 seq=10/g/a,10/V/a',
                    'ZERO seq=15/Ic/0',
                    'WAIT ut=16:00:00',
                    'BELL']:
            with self.subTest(txt=txt):
                self.assertTrue(self._validate(txt)['valid'])

    def test_wait_without_argument(self):
        result = self._validate('WAIT')
        self.assertFalse(result['valid'])
        self.assertFalse(result['result']['ut'])

    def test_unknown_filter(self):
        result = self._validate('OBJECT HD193901 20:23:35.8 -21:22:14.0 seq=1/dupa/300')
        self.assertFalse(result['valid'])
        self.assertFalse(result['result']['seq'])

    def test_coordinates_and_horizontal_position_are_exclusive(self):
        result = self._validate(
            'OBJECT V496_Aql 19:08:20.77 -07:26:15.89 alt=34 az=270.0 seq=1/V/20')
        self.assertFalse(result['valid'])
        for key in ('ra', 'dec', 'alt', 'az'):
            self.assertFalse(result['result'][key])

    def test_unparsable_text(self):
        self.assertEqual(self._validate('dupa'),
                         {'valid': False, 'result': {}, 'data': {}, 'required': {}, 'allowed': {}})


class TestCalcSeqTime(unittest.TestCase):

    def test_single_exposure(self):
        # 1 * (300 + 10) + 5
        self.assertEqual(ObsValidator.calc_seq_time('1/V/300'), 315.0)

    def test_two_elements_add_filter_overhead(self):
        # 2 * (30 + 10) + 3 * (40 + 10) + 2 + 5
        self.assertEqual(ObsValidator.calc_seq_time('2/V/30,3/Ic/40'), 237.0)

    def test_multiplier(self):
        # 2 * ((10 + 10) + (10 + 10) + 2) + 5
        self.assertEqual(ObsValidator.calc_seq_time('2x(1/V/10,1/r/10)'), 89.0)

    def test_auto_exposure_uses_auto_time(self):
        self.assertEqual(ObsValidator.calc_seq_time('1/V/a'), 20.0)
        self.assertEqual(ObsValidator.calc_seq_time('1/V/a', auto_time=30), 45.0)

    def test_custom_overheads(self):
        self.assertEqual(ObsValidator.calc_seq_time('1/V/10', base_time=0, overhead=0), 10.0)

    def test_malformed_sequence_returns_none(self):
        for seq in ['dupa', '1/V', '1/V/x']:
            with self.subTest(seq=seq):
                self.assertIsNone(ObsValidator.calc_seq_time(seq))


if __name__ == "__main__":
    unittest.main()
