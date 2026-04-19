import unittest
from pyaraucaria.obs_plan.obs_plan_parser import ObsPlanParser


class TestObsPlanParser(unittest.TestCase):

    def test_sequence_1(self):
        input = """
            BEGINSEQUENCE x y z execute_at_time=16:00
                ZERO seq=15/I/0
                DARK seq=10/V/300,10/I/200
                BEGINSEQUENCE abc=13
                    DOMEFLAT seq=7/V/20,7/I/20
                    BEGINSEQUENCE
                        DOMEFLAT seq=10/str_u/100 domeflat_lamp=0.7
                    ENDSEQUENCE
                ENDSEQUENCE
                OBJECT FF_Aql 18:58:14.75 17:21:39.29 seq=2/I/60,2/V/70
            ENDSEQUENCE
    
            BEGINSEQUENCE execute_periodically=02:00 priority=+10
                FOCUS NG31 12:12:12 20:20:20
            ENDSEQUENCE
            """

        output = {'command_name': 'SEQUENCE', 'subcommands': [
            {'command_name': 'SEQUENCE', 'args': ['x', 'y', 'z'], 'kwargs': {'execute_at_time': '16:00'}, 'subcommands':
                [{'command_name': 'ZERO', 'kwargs': {'seq': '15/I/0'}},
                 {'command_name': 'DARK', 'kwargs': {'seq': '10/V/300,10/I/200'}},
                 {'command_name': 'SEQUENCE', 'kwargs': {'abc': '13'}, 'subcommands': [
                     {'command_name': 'DOMEFLAT', 'kwargs': {'seq': '7/V/20,7/I/20'}},
                     {'command_name': 'SEQUENCE', 'subcommands': [
                         {'command_name': 'DOMEFLAT', 'kwargs': {'seq': '10/str_u/100', 'domeflat_lamp': '0.7'}}]}]},
                 {'command_name': 'OBJECT', 'args': ['FF_Aql', '18:58:14.75', '17:21:39.29'],
                  'kwargs': {'seq': '2/I/60,2/V/70'}}]},
            {'command_name': 'SEQUENCE', 'kwargs': {'execute_periodically': '02:00', 'priority': '+10'},
             'subcommands': [
                 {'command_name': 'FOCUS', 'args': ['NG31', '12:12:12', '20:20:20']}]}]}
        print(ObsPlanParser.convert_from_string("OBJECT FF_Aql 18:58:14.75 17:21:39.29 seq=2/I/60,2/V/70"))
        self.assertEqual(ObsPlanParser.convert_from_string(input), output)

    def test_sequence_2(self):
        input = """
            BEGINSEQUENCE
            # WAIT ut=16:00
            # ZERO seq=15/I/0
            STOP
            SNAP seq=10/V/300,10/I/200
            DARK seq=10/V/300,10/I/200
            DOMEFLAT seq=7/V/20,7/I/20
            DOMEFLAT seq=10/str_u/100 domeflat_lamp=0.7
            WAIT sunset=-12
            SKYFLAT alt=60:00:00 az=270:00:00  seq=10/I/20,10/V/30 
            SKYFLAT seq=10/I/0,10/V/0
            SKYFLAT_adu=30
            WAIT t=600
            FOCUS NG31 12:12:12 20:20:20
            OBJECT HD193901 20:23:35.8 -21:22:14.0 seq=1/V/300
            OBJECT FF_Aql 18:58:14.75 17:21:39.29 seq=5/I/60,5/V/70
            OBJECT V496_Aql 19:08:20.77 -07:26:15.89 seq=1/V/20 focus=+30
            ENDSEQUENCE
            """

        output = {'command_name': 'SEQUENCE', 'subcommands': [
            {'command_name': 'SEQUENCE', 'subcommands': [
                {'command_name': 'STOP'},
                {'command_name': 'SNAP', 'kwargs': {'seq': '10/V/300,10/I/200'}},
                {'command_name': 'DARK', 'kwargs': {'seq': '10/V/300,10/I/200'}},
                {'command_name': 'DOMEFLAT', 'kwargs': {'seq': '7/V/20,7/I/20'}},
                {'command_name': 'DOMEFLAT', 'kwargs': {'seq': '10/str_u/100', 'domeflat_lamp': '0.7'}},
                {'command_name': 'WAIT', 'kwargs': {'sunset': '-12'}},
                {'command_name': 'SKYFLAT', 'kwargs': {'alt': '60:00:00', 'az': '270:00:00', 'seq': '10/I/20,10/V/30'}},
                {'command_name': 'SKYFLAT', 'kwargs': {'seq': '10/I/0,10/V/0'}},
                {'command_name': 'SKYFLAT', 'kwargs': {'_adu': '30'}},
                {'command_name': 'WAIT', 'kwargs': {'t': '600'}},
                {'command_name': 'FOCUS', 'args': ['NG31', '12:12:12', '20:20:20']},
                {'command_name': 'OBJECT', 'args': ['HD193901', '20:23:35.8', '-21:22:14.0'],
                 'kwargs': {'seq': '1/V/300'}},
                {'command_name': 'OBJECT', 'args': ['FF_Aql', '18:58:14.75', '17:21:39.29'],
                 'kwargs': {'seq': '5/I/60,5/V/70'}},
                {'command_name': 'OBJECT', 'args': ['V496_Aql', '19:08:20.77', '-07:26:15.89'],
                 'kwargs': {'seq': '1/V/20', 'focus': '+30'}}]}]}

        self.assertEqual(ObsPlanParser.convert_from_string(input), output)

    def test_sequence_3(self):
        input = """      
            BEGINSEQUENCE execute_at_time=16:00
                ZERO seq=15/I/0
                DARK seq=10/V/300,10/I/200
                DOMEFLAT seq=7/V/20,7/I/20
                DOMEFLAT seq=10/str_u/100 domeflat_lamp=0.7
            ENDSEQUENCE
            
            BEGINSEQUENCE execute_at_time=02:21:43 priority=+30  # scheduled obs
                OBJECT FF_Aql 18:58:14.75 17:21:39.29 seq=2/I/60,2/V/70
            ENDSEQUENCE
            
            BEGINSEQUENCE execute_periodically=02:00 priority=+10
                FOCUS NG31 12:12:12 20:20:20
            ENDSEQUENCE
            
            BEGINSEQUENCE execute_at_dusk=-12
                SKYFLAT alt=60:00:00 az=270:00:00 seq=10/I/20,10/V/30 
                SKYFLAT seq=10/I/0,10/V/0 skyflat_adu=30
                WAIT t=600
                FOCUS NG31 12:12:12 20:20:20
                OBJECT HD193901 20:23:35.8 -21:22:14.0 seq=1/V/300
                OBJECT FF_Aql 18:58:14.75 17:21:39.29 seq=5/I/60,5/V/70
                OBJECT V496_Aql 19:08:20.77 -07:26:15.89 seq=1/V/20 focus=+30
            ENDSEQUENCE
            
            BEGINSEQUENCE execute_at_dawn=-6 priority=+10
                SKYFLAT alt=60:00:00 az=270:00:00 seq=10/I/20,10/V/30 
                SKYFLAT seq=10/I/0,10/V/0 skyflat_adu=30
                SKYFLAT alt=60:00:00 az=270:00:00 seq=10/I/20,10/V/30 
                SKYFLAT seq=10/I/0,10/V/0 skyflat_adu=30
            ENDSEQUENCE
            
            BEGINSEQUENCE execute_at_dawn=+2 priority=+100
                PARK 
                DOMECLOSE
            ENDSEQUENCE
            """

        output = {'command_name': 'SEQUENCE', 'subcommands': [
            {'command_name': 'SEQUENCE', 'kwargs': {'execute_at_time': '16:00'}, 'subcommands': [
                {'command_name': 'ZERO', 'kwargs': {'seq': '15/I/0'}},
                {'command_name': 'DARK', 'kwargs': {'seq': '10/V/300,10/I/200'}},
                {'command_name': 'DOMEFLAT', 'kwargs': {'seq': '7/V/20,7/I/20'}},
                {'command_name': 'DOMEFLAT', 'kwargs': {'seq': '10/str_u/100', 'domeflat_lamp': '0.7'}}]},
            {'command_name': 'SEQUENCE', 'kwargs': {'execute_at_time': '02:21:43', 'priority': '+30'}, 'subcommands': [
                {'command_name': 'OBJECT', 'args': ['FF_Aql', '18:58:14.75', '17:21:39.29'],
                 'kwargs': {'seq': '2/I/60,2/V/70'}}]},
            {'command_name': 'SEQUENCE', 'kwargs': {'execute_periodically': '02:00', 'priority': '+10'},
             'subcommands': [
                 {'command_name': 'FOCUS', 'args': ['NG31', '12:12:12', '20:20:20']}]},
            {'command_name': 'SEQUENCE', 'kwargs': {'execute_at_dusk': '-12'}, 'subcommands': [
                {'command_name': 'SKYFLAT', 'kwargs': {'alt': '60:00:00', 'az': '270:00:00', 'seq': '10/I/20,10/V/30'}},
                {'command_name': 'SKYFLAT', 'kwargs': {'seq': '10/I/0,10/V/0', 'skyflat_adu': '30'}},
                {'command_name': 'WAIT', 'kwargs': {'t': '600'}},
                {'command_name': 'FOCUS', 'args': ['NG31', '12:12:12', '20:20:20']},
                {'command_name': 'OBJECT', 'args': ['HD193901', '20:23:35.8', '-21:22:14.0'],
                 'kwargs': {'seq': '1/V/300'}},
                {'command_name': 'OBJECT', 'args': ['FF_Aql', '18:58:14.75', '17:21:39.29'],
                 'kwargs': {'seq': '5/I/60,5/V/70'}},
                {'command_name': 'OBJECT', 'args': ['V496_Aql', '19:08:20.77', '-07:26:15.89'],
                 'kwargs': {'seq': '1/V/20', 'focus': '+30'}}]},
            {'command_name': 'SEQUENCE', 'kwargs': {'execute_at_dawn': '-6', 'priority': '+10'}, 'subcommands': [
                {'command_name': 'SKYFLAT', 'kwargs': {'alt': '60:00:00', 'az': '270:00:00', 'seq': '10/I/20,10/V/30'}},
                {'command_name': 'SKYFLAT', 'kwargs': {'seq': '10/I/0,10/V/0', 'skyflat_adu': '30'}},
                {'command_name': 'SKYFLAT', 'kwargs': {'alt': '60:00:00', 'az': '270:00:00', 'seq': '10/I/20,10/V/30'}},
                {'command_name': 'SKYFLAT', 'kwargs': {'seq': '10/I/0,10/V/0', 'skyflat_adu': '30'}}]},
            {'command_name': 'SEQUENCE', 'kwargs': {'execute_at_dawn': '+2', 'priority': '+100'}, 'subcommands': [
                {'command_name': 'PARK'},
                {'command_name': 'DOMECLOSE'}]}]}

        self.assertEqual(ObsPlanParser.convert_from_string(input), output)

    def test_sequence_4(self):
        # No use: BEGINSEQUENCE ENDSEQUENCE

        input = """
            WAIT ut=16:00
            ZERO seq=15/I/0
            DARK seq=10/V/300,10/I/200
            DOMEFLAT seq=7/V/20,7/I/20
            DOMEFLAT seq=10/str_u/100 domeflat_lamp=0.7
            WAIT sunset=-12
            SKYFLAT alt=60:00:00 az=270:00:00 seq=10/I/20,10/V/30 
            SKYFLAT seq=10/I/0,10/V/0
            SKYFLAT adu=30
            WAIT t=600
            FOCUS NG31 12:12:12 20:20:20
            OBJECT HD193901 20:23:35.8 -21:22:14.0 seq=1/V/300
            OBJECT FF_Aql 18:58:14.75 17:21:39.29 seq=5/I/60,5/V/70
            OBJECT V496_Aql 19:08:20.77 -07:26:15.89 seq=1/V/20 focus=+30
            """

        output = {'command_name': 'SEQUENCE', 'subcommands': [
            {'command_name': 'WAIT', 'kwargs': {'ut': '16:00'}},
            {'command_name': 'ZERO', 'kwargs': {'seq': '15/I/0'}},
            {'command_name': 'DARK', 'kwargs': {'seq': '10/V/300,10/I/200'}},
            {'command_name': 'DOMEFLAT', 'kwargs': {'seq': '7/V/20,7/I/20'}},
            {'command_name': 'DOMEFLAT', 'kwargs': {'seq': '10/str_u/100', 'domeflat_lamp': '0.7'}},
            {'command_name': 'WAIT', 'kwargs': {'sunset': '-12'}},
            {'command_name': 'SKYFLAT', 'kwargs': {'alt': '60:00:00', 'az': '270:00:00', 'seq': '10/I/20,10/V/30'}},
            {'command_name': 'SKYFLAT', 'kwargs': {'seq': '10/I/0,10/V/0'}},
            {'command_name': 'SKYFLAT', 'kwargs': {'adu': '30'}},
            {'command_name': 'WAIT', 'kwargs': {'t': '600'}},
            {'command_name': 'FOCUS', 'args': ['NG31', '12:12:12', '20:20:20']},
            {'command_name': 'OBJECT', 'args': ['HD193901', '20:23:35.8', '-21:22:14.0'], 'kwargs': {'seq': '1/V/300'}},
            {'command_name': 'OBJECT', 'args': ['FF_Aql', '18:58:14.75', '17:21:39.29'],
             'kwargs': {'seq': '5/I/60,5/V/70'}},
            {'command_name': 'OBJECT', 'args': ['V496_Aql', '19:08:20.77', '-07:26:15.89'],
             'kwargs': {'seq': '1/V/20', 'focus': '+30'}}]}

        self.assertEqual(ObsPlanParser.convert_from_string(input), output)

    def test_sequence_5_comments(self):
        # No use: BEGINSEQUENCE ENDSEQUENCE

        input = """
            WAIT ut=16:00
            ZERO seq=15/I/0
            # comment
            """

        output = {'command_name': 'SEQUENCE', 'subcommands': [
            {'command_name': 'WAIT', 'kwargs': {'ut': '16:00'}},
            {'command_name': 'ZERO', 'kwargs': {'seq': '15/I/0'}}]}

        self.assertEqual(ObsPlanParser.convert_from_string(input), output)

    def test_sequence_6_comments(self):
        # No use: BEGINSEQUENCE ENDSEQUENCE

        input = """
            WAIT ut=16:00
            ZERO seq=15/I/0
            # comment
            OBJECT HD193901 20:23:35.8 -21:22:14.0 seq=1/V/300
            """

        output = {'command_name': 'SEQUENCE', 'subcommands': [
            {'command_name': 'WAIT', 'kwargs': {'ut': '16:00'}},
            {'command_name': 'ZERO', 'kwargs': {'seq': '15/I/0'}},
            {'command_name': 'OBJECT', 'args': ['HD193901', '20:23:35.8', '-21:22:14.0'],
             'kwargs': {'seq': '1/V/300'}}]}

        self.assertEqual(ObsPlanParser.convert_from_string(input), output)


    def test_commands_with_blacklisted_first_letters(self):
        """Commands whose first letter appears in BEGINSEQUENCE (B,C,E,G,I,N,Q,S,U) must parse."""
        cases = [
            ("BIAS seq=15/V/0", {'command_name': 'BIAS', 'kwargs': {'seq': '15/V/0'}}),
            ("CAL target=dome", {'command_name': 'CAL', 'kwargs': {'target': 'dome'}}),
            ("ENGINEERING mode=test", {'command_name': 'ENGINEERING', 'kwargs': {'mode': 'test'}}),
            ("GUIDE star=HD1234", {'command_name': 'GUIDE', 'kwargs': {'star': 'HD1234'}}),
            ("IDLE", {'command_name': 'IDLE'}),
            ("NOTE observer=someone", {'command_name': 'NOTE', 'kwargs': {'observer': 'someone'}}),
            ("QUIT", {'command_name': 'QUIT'}),
            ("SKYFLAT seq=1/V/1", {'command_name': 'SKYFLAT', 'kwargs': {'seq': '1/V/1'}}),
            ("SNAP seq=1/V/1", {'command_name': 'SNAP', 'kwargs': {'seq': '1/V/1'}}),
            ("STOP", {'command_name': 'STOP'}),
            ("UVFILT filter=U", {'command_name': 'UVFILT', 'kwargs': {'filter': 'U'}}),
        ]
        for text, expected in cases:
            with self.subTest(command=text):
                result = ObsPlanParser.convert_from_string(text)
                self.assertIsNotNone(result, f"Parsing returned None for: {text!r}")
                # The top-level wrapper is a SEQUENCE; the command is in subcommands[0]
                self.assertEqual(result['subcommands'][0], expected)

    def test_quoted_value_containing_snap_not_corrupted(self):
        """A kwarg value containing SNAP/STOP/SKYFLAT must not be mangled."""
        text = 'OBJECT HD1 comment="SNAP judgment"'
        result = ObsPlanParser.convert_from_string(text)
        self.assertIsNotNone(result)
        cmd = result['subcommands'][0]
        self.assertEqual(cmd['command_name'], 'OBJECT')
        self.assertEqual(cmd['kwargs']['comment'], '"SNAP judgment"')

    def test_no_leading_newline_in_command_name(self):
        """A blank line before a command must not prepend a newline to the command name."""
        text = "\nOBJECT HD1 12:00:00 +00:00:00"
        result = ObsPlanParser.convert_from_string(text)
        self.assertIsNotNone(result)
        cmd = result['subcommands'][0]
        self.assertEqual(cmd['command_name'], 'OBJECT')


if __name__ == '__main__':
    unittest.main()
