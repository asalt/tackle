from setuptools import setup

def calculate_version(inputfile):
    version_list =  [x.split('\'')[1] for x in open(inputfile, 'r')
                     if x.startswith('__version__')]
    if version_list:
        return version_list[0]
    else:
        return '1.0'

package_version = calculate_version('main.py')
setup(
    name='correlation-plotter',
    version=package_version,
    py_modules=['main', 'utils', 'clusterplot', 'scatterplot', 'pcaplot', 'containers'],
    install_requires=[
        'Click',
    ],
    entry_points="""
    [console_scripts]
    correlationplot=main:main
    """
)
