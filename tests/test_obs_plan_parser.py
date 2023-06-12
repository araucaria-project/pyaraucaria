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

        output = [{'begin_sequence': 'begin', 'args': ['x', 'y', 'z'], 'kwargs': {'execute_at_time': '16:00'},
                   'all_commands': [{'command_name': 'ZERO', 'kwargs': {'seq': '15/I/0'}},
                                    {'command_name': 'DARK', 'kwargs': {'seq': '10/V/300,10/I/200'}},
                                    {'begin_sequence': 'begin', 'kwargs': {'abc': '13'}, 'all_commands':
                                        [{'command_name': 'DOMEFLAT', 'kwargs': {'seq': '7/V/20,7/I/20'}},
                                         {'begin_sequence': 'begin', 'all_commands':
                                             [{'command_name': 'DOMEFLAT', 'kwargs':
                                                 {'seq': '10/str_u/100', 'domeflat_lamp': '0.7'}}]}]},
                                    {'command_name': 'OBJECT', 'args': ['FF_Aql', '18:58:14.75', '17:21:39.29'],
                                     'kwargs': {'seq': '2/I/60,2/V/70'}}]},
                  {'begin_sequence': 'begin', 'kwargs': {'execute_periodically': '02:00', 'priority': '+10'},
                   'all_commands': [{'command_name': 'FOCUS', 'args': ['NG31', '12:12:12', '20:20:20']}]}]

        opp = ObsPlanParser()

        self.assertEqual(opp.convert_from_string(input), output)

    def test_sequence_2(self):

        input = """
            BEGINSEQUENCE
            # WAIT ut=16:00
            # ZERO seq=15/I/0
            STOP
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

        output = [
            {'begin_sequence': 'begin', 'all_commands': [
            # {'command_name': 'WAIT', 'kwargs': {'ut': '16:00'}},
            # {'command_name': 'ZERO', 'kwargs': {'seq': '15/I/0'}},
            {'command_name': 'STOP'},
            {'command_name': 'DARK', 'kwargs': {'seq': '10/V/300,10/I/200'}},
            {'command_name': 'DOMEFLAT', 'kwargs': {'seq': '7/V/20,7/I/20'}},
            {'command_name': 'DOMEFLAT', 'kwargs': {'seq': '10/str_u/100', 'domeflat_lamp': '0.7'}},
            {'command_name': 'WAIT', 'kwargs': {'sunset': '-12'}},
            {'command_name': 'SKYFLAT', 'kwargs': {'alt': '60:00:00', 'az': '270:00:00', 'seq': '10/I/20,10/V/30'}},
            {'command_name': 'SKYFLAT', 'kwargs': {'seq': '10/I/0,10/V/0'}},
            {'command_name': 'SKYFLAT', 'kwargs': {'_adu': '30'}},
            {'command_name': 'WAIT', 'kwargs': {'t': '600'}},
            {'command_name': 'FOCUS', 'args': ['NG31', '12:12:12', '20:20:20']},
            {'command_name': 'OBJECT', 'args':
                ['HD193901', '20:23:35.8', '-21:22:14.0'], 'kwargs': {'seq': '1/V/300'}},
            {'command_name': 'OBJECT', 'args':
                ['FF_Aql', '18:58:14.75', '17:21:39.29'], 'kwargs': {'seq': '5/I/60,5/V/70'}},
            {'command_name': 'OBJECT', 'args':
                ['V496_Aql', '19:08:20.77', '-07:26:15.89'], 'kwargs': {'seq': '1/V/20', 'focus': '+30'}}]}]

        opp = ObsPlanParser()
        self.assertEqual(opp.convert_from_string(input), output)

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

        output = [
            {'begin_sequence': 'begin', 'kwargs': {'execute_at_time': '16:00'}, 'all_commands':
                [{'command_name': 'ZERO', 'kwargs': {'seq': '15/I/0'}},
                 {'command_name': 'DARK', 'kwargs': {'seq': '10/V/300,10/I/200'}},
                 {'command_name': 'DOMEFLAT', 'kwargs': {'seq': '7/V/20,7/I/20'}},
                 {'command_name': 'DOMEFLAT', 'kwargs': {'seq': '10/str_u/100', 'domeflat_lamp': '0.7'}}]},
            {'begin_sequence': 'begin', 'kwargs': {'execute_at_time': '02:21:43', 'priority': '+30'}, 'all_commands':
                [{'command_name': 'OBJECT', 'args': ['FF_Aql', '18:58:14.75', '17:21:39.29'], 'kwargs':
                    {'seq': '2/I/60,2/V/70'}}]},
            {'begin_sequence': 'begin', 'kwargs': {'execute_periodically': '02:00', 'priority': '+10'}, 'all_commands':
                [{'command_name': 'FOCUS', 'args': ['NG31', '12:12:12', '20:20:20']}]},
            {'begin_sequence': 'begin', 'kwargs': {'execute_at_dusk': '-12'}, 'all_commands':
                [{'command_name': 'SKYFLAT', 'kwargs': {'alt': '60:00:00', 'az': '270:00:00', 'seq': '10/I/20,10/V/30'}},
                 {'command_name': 'SKYFLAT', 'kwargs': {'seq': '10/I/0,10/V/0', 'skyflat_adu': '30'}},
                 {'command_name': 'WAIT', 'kwargs': {'t': '600'}},
                 {'command_name': 'FOCUS', 'args': ['NG31', '12:12:12', '20:20:20']},
                 {'command_name': 'OBJECT', 'args': ['HD193901', '20:23:35.8', '-21:22:14.0'], 'kwargs':
                     {'seq': '1/V/300'}},
                 {'command_name': 'OBJECT', 'args': ['FF_Aql', '18:58:14.75', '17:21:39.29'], 'kwargs':
                     {'seq': '5/I/60,5/V/70'}},
                 {'command_name': 'OBJECT', 'args': ['V496_Aql', '19:08:20.77', '-07:26:15.89'], 'kwargs':
                     {'seq': '1/V/20', 'focus': '+30'}}]},
            {'begin_sequence': 'begin', 'kwargs': {'execute_at_dawn': '-6', 'priority': '+10'}, 'all_commands':
                [{'command_name': 'SKYFLAT', 'kwargs': {'alt': '60:00:00', 'az': '270:00:00', 'seq': '10/I/20,10/V/30'}},
                 {'command_name': 'SKYFLAT', 'kwargs': {'seq': '10/I/0,10/V/0', 'skyflat_adu': '30'}},
                 {'command_name': 'SKYFLAT', 'kwargs': {'alt': '60:00:00', 'az': '270:00:00', 'seq': '10/I/20,10/V/30'}},
                 {'command_name': 'SKYFLAT', 'kwargs': {'seq': '10/I/0,10/V/0', 'skyflat_adu': '30'}}]},
            {'begin_sequence': 'begin', 'kwargs': {'execute_at_dawn': '+2', 'priority': '+100'}, 'all_commands':
                [{'command_name': 'PARK'}, {'command_name': 'DOMECLOSE'}]}]

        opp = ObsPlanParser()
        self.assertEqual(opp.convert_from_string(input), output)

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

        output = [
            {'begin_sequence': 'begin', 'all_commands': [
                {'command_name': 'WAIT', 'kwargs': {'ut': '16:00'}},
                {'command_name': 'ZERO', 'kwargs': {'seq': '15/I/0'}},
                {'command_name': 'DARK', 'kwargs': {'seq': '10/V/300,10/I/200'}},
                {'command_name': 'DOMEFLAT', 'kwargs': {'seq': '7/V/20,7/I/20'}},
                {'command_name': 'DOMEFLAT', 'kwargs': {'seq': '10/str_u/100', 'domeflat_lamp': '0.7'}},
                {'command_name': 'WAIT', 'kwargs': {'sunset': '-12'}},
                {'command_name': 'SKYFLAT',
                 'kwargs': {'alt': '60:00:00', 'az': '270:00:00', 'seq': '10/I/20,10/V/30'}},
                {'command_name': 'SKYFLAT', 'kwargs': {'seq': '10/I/0,10/V/0'}},
                {'command_name': 'SKYFLAT', 'kwargs': {'adu': '30'}},
                {'command_name': 'WAIT', 'kwargs': {'t': '600'}},
                {'command_name': 'FOCUS', 'args': ['NG31', '12:12:12', '20:20:20']},
                {'command_name': 'OBJECT', 'args':
                    ['HD193901', '20:23:35.8', '-21:22:14.0'], 'kwargs': {'seq': '1/V/300'}},
                {'command_name': 'OBJECT', 'args':
                    ['FF_Aql', '18:58:14.75', '17:21:39.29'], 'kwargs': {'seq': '5/I/60,5/V/70'}},
                {'command_name': 'OBJECT', 'args':
                    ['V496_Aql', '19:08:20.77', '-07:26:15.89'], 'kwargs': {'seq': '1/V/20', 'focus': '+30'}}]}]

        opp = ObsPlanParser()
        self.assertEqual(opp.convert_from_string(input), output)

    def test_sequence_5_comments(self):

        # No use: BEGINSEQUENCE ENDSEQUENCE

        input = """
            WAIT ut=16:00
            ZERO seq=15/I/0
            # dupa
            """

        output = [
            {'begin_sequence': 'begin', 'all_commands': [
                {'command_name': 'WAIT', 'kwargs': {'ut': '16:00'}},
                {'command_name': 'ZERO', 'kwargs': {'seq': '15/I/0'}},
            ]}]

        opp = ObsPlanParser()
        self.assertEqual(opp.convert_from_string(input), output)

    def test_sequence_6_comments(self):

        # No use: BEGINSEQUENCE ENDSEQUENCE

        input = """
            WAIT ut=16:00
            ZERO seq=15/I/0
            # dupa
            OBJECT HD193901 20:23:35.8 -21:22:14.0 seq=1/V/300
            """

        output = [
            {'begin_sequence': 'begin', 'all_commands': [
                {'command_name': 'WAIT', 'kwargs': {'ut': '16:00'}},
                {'command_name': 'ZERO', 'kwargs': {'seq': '15/I/0'}},
                {'command_name': 'OBJECT', 'args':
                    ['HD193901', '20:23:35.8', '-21:22:14.0'], 'kwargs': {'seq': '1/V/300'}},
            ]}]

        opp = ObsPlanParser()
        self.assertEqual(opp.convert_from_string(input), output)

if __name__ == '__main__':
    unittest.main()
