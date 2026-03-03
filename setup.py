from setuptools import setup, find_packages

def calculate_version(inputfile):
    version_list =  [x.split('\"')[1] for x in open(inputfile, 'r')
                     if x.startswith('__version__')]
    if version_list:
        return version_list[0]
    else:
        return '1.0'

package_version = calculate_version('./tackle/main.py')
setup(
    name='tackle',
    version=package_version,
    packages=find_packages(),
    author = 'Alex Saltzman',
    author_email = 'a.saltzman920@gmail.com',
    description = 'CLI-based data filtering and plotting tools for proteomics data at BCM',
    py_modules=['main', 'utils', 'clusterplot', 'scatterplot', 'pcaplot', 'containers'],
    install_requires=[
        'Click',
        'Jinja2',
    ],
    extras_require={
        'excel': [
            'openpyxl',
            'xlsxwriter',
        ],
        'excel-fast': [
            'openpyxl',
            'xlsxwriter',
            'pyarrow',
        ],
        'slides': [
            'python-pptx',
        ],
    },
    entry_points="""
    [console_scripts]
    tackle=tackle.main:main
    tackle1=tackle.cli1.app:cli
    # correlationplot=main:main
    """,

    package_data={'': ['R/*.r', 'R/*.R', 'GSEA/*.jar', 'GSEA/*.gmt', 'GSEA/genesets/*gmt', 'data/*.data', 'data/*txt', 'data/*tab', 'data/*tsv', 'data/*xlsx', 'templates/*.j2']},
    include_package_data=True

)
