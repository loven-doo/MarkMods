from setuptools import setup, find_packages

setup(
    name='MarkMods',
    version='0.0.3',
    packages=find_packages(),
    package_data={},
    install_requires=[
        'redis >= 2.10.6',
        'numpy >= 1.15.3',
    ],
    entry_points={
        'console_scripts': []
    }
)
