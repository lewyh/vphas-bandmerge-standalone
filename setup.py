from setuptools import setup

entry_points = {'console_scripts': ['vphas-bandmerge = bandmerge.bandmerge:process']}

setup(
    name='vphas-bandmerge-standalone',
    version='0.1',
    packages=['bandmerge'],
    url='http://www.vphas.eu',
    license='MIT',
    author='Hywel Farnhill',
    author_email='hywel.farnhill@gmail.com',
    entry_points=entry_points,
    install_requires=['numpy>=1.9',
                      'astropy>=0.3']
)
