import io, re
from setuptools import setup

with io.open("README.rst", "rt", encoding="utf8") as f:
    readme = f.read()

with io.open("fluxdataqaqc/__init__.py", "rt", encoding="utf8") as f:
    version = re.search(r"__version__ = \'(.*?)\'", f.read()).group(1)

requires = [
    'bokeh>=2.0',
    'netCDF4',
    'numpy>=1.15',
    'pandas>=1.0',
    'refet',
    'scikit-learn',
    'xarray',
    'openpyxl'
]

tests_require = ['pytest']

classifiers = [
    'License :: OSI Approved :: BSD License',
    'Programming Language :: Python :: 3.7',
    'Environment :: Console',
    'Development Status :: 4 - Beta',
    'Topic :: Scientific/Engineering',
    'Intended Audience :: Science/Research'
]

setup(
    name='fluxdataqaqc',
    description='Tools for QA/QC of eddy covariance station data',
    long_description=readme,
    author='John Volk',
    author_email='john.volk@dri.edu',
    license='BSD3',
    python_requires='>=3.7',
    version=version,
    url='https://github.com/Open-ET/flux-data-qaqc',
    platforms=['Windows','Linux','Mac OS X'],
    classifiers=classifiers,
    packages=['fluxdataqaqc'],
    package_data={
        'examples': ['Basic_usage/*'],
    	'examples': ['Config_options/*'],
        '': ['*.yml']
    },
    install_requires=requires,
    tests_require=tests_require,
    include_package_data=True,
)
