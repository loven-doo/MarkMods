from setuptools import setup, find_packages

setup(name='MarkMods',
      version='0.0.1',
      packages=find_packages(),
      package_data={},
      install_requires=[
            'redis >= 2.10.6'
      ],
      entry_points={
            'console_scripts': []
      })
