from setuptools import setup, find_packages

setup(
    name='MarkMods',
    version='0.0.4',
    packages=find_packages(),
    package_data={},
    install_requires=[
        'redis >= 2.10.6',
        'numpy >= 1.15.3',
        'pandas >= 0.23.4',
        'scipy >= 1.1.0',
    ],
    entry_points={
        'console_scripts': []
    }
)
