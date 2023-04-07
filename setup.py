from setuptools import setup

setup(
    name='pyaraucaria',
    version='2.6.17',
    packages=['pyaraucaria'],
    url='',
    license='MIT/GPL',
    author='',
    author_email='',
    description='Common Routines of OCA Observatory and Araucaria Project',
    include_package_data=True,

    package_data={'pyaraucaria': [
        'databases/Objects.database',
        'databases/TAB.ALL',
        'obs_plan/obs_plan_parser.py',
    ]},
    entry_points={'console_scripts': [
        'lookup_objects=pyaraucaria.lookup_objects:main',
    ]},
    install_requires=[
        'numpy', 'ephem', 'lark', 'astropy'
    ],
    extras_require={
        'yaml': ['yaml'],
    }
)
