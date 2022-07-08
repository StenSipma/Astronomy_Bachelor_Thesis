from setuptools import setup, find_packages
#packages = find_packages(exclude=['.txt', '.gitignore', '.npy','Old'])
packages = find_packages(include=["aurigaia", "coords"])
print('My Packages:')
print(packages)
setup(name='aurigaia',
      packages=packages
      )
