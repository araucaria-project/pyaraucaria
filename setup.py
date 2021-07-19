from setuptools import setup

setup(
    name='pyaraucaria',
    version='1.0',
    packages=['pyaraucaria'],
    url='',
    license='MIT',
    author='',
    author_email='',
    description='OCA Observatory and Araucaria Project Common Routines Library',
    include_package_data=True,

    package_data={'pyaraucaria': [
        'databases/Objects.database',
        'databases/TAB.ALL',
    ]},
    entry_points={'console_scripts': [
        'lookup_objects=pyaraucaria.lookup_objects:main',
    ]},
)
