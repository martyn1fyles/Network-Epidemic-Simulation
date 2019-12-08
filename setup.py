from setuptools import setup

setup(
   name='NetworkEpidemicSimulation',
   version='alpha',
   description='Performs Sellke Type simulations.',
   license="GNU",
   #long_description=long_description,
   #author='Man Foo',
   #author_email='foomail@foo.com',
   #url="http://www.foopackage.com/",
   packages=['NetworkEpidemicSimulation'],  #same as name
   install_requires=['networkx', 'numpy', 'scipy', 'matplotlib', 'pytest'], #external packages as dependencies
   py_modules = ['NetworkEpidemicSimulation/EpidemicSimulation', 'NetworkEpidemicSimulation/Simulation']
)