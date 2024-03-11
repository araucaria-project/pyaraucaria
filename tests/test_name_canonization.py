import unittest
from pyaraucaria.lookup_objects import name_canonizator 



class TestNameCanonization(unittest.TestCase):
    def test_canonization(self):
        assert name_canonizator("Hello World!") == "helloworld"
        assert name_canonizator("123 ABC!") == "123abc"
        assert name_canonizator("Test@123") == "test123"
        assert name_canonizator("NoChange") == "nochange"
        assert name_canonizator("With Spaces ") == "withspaces"
        assert name_canonizator("With - Dashes") == "withdashes"
        assert name_canonizator("With_underscores") == "withunderscores"
        assert name_canonizator("tz-for") == "tzfor"
if __name__ == "__main__":
    unittest.main()
