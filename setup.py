import os.path
import setuptools

# Get the long description from README.
with open('README.rst', 'r') as fh:
    long_description = fh.read()

# Get package metadata from '__about__.py' file.
about = {}
base_dir = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(base_dir, 'src', 'ngs_test_utils', '__about__.py'), 'r') as fh:
    exec(fh.read(), about)

setuptools.setup(
    name=about['__title__'],
    version=about['__version__'],
    description=about['__summary__'],
    long_description=long_description,
    long_description_content_type='text/x-rst',
    author=about['__author__'],
    author_email=about['__email__'],
    url=about['__url__'],
    license=about['__license__'],
    # Exclude tests from built/installed package.
    packages=setuptools.find_packages(
        'src', exclude=['tests', 'tests.*']
    ),
    package_dir={'': 'src'},
    python_requires='>=3.6, <3.8',
    install_requires=[
        'numpy',
        'pybedtools',
        'pysam',
        'biopython',
   ],
    extras_require={
        'docs': ['sphinx_rtd_theme'],
        'package': ['twine', 'wheel'],
        'test': [
            'black',
            'check-manifest',
            'docutils',
            'flake8',
            'isort',
            'pydocstyle',
            'pytest-cov',
            'setuptools_scm',
        ],
    },
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: Apache Software License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    keywords='ngs test bioinformatics',
)
