from setuptools import setup

with open("README", 'r') as f:
    long_description = f.read()

setup(
   name='SellkeSimulation',
   version='alpha',
   description='Performs Sellke Type simulations.',
   license="GNU",
   #long_description=long_description,
   #author='Man Foo',
   #author_email='foomail@foo.com',
   #url="http://www.foopackage.com/",
   packages=['SellkeSimulation'],  #same as name
   install_requires=['networkx', 'numpy', 'scipy', 'matplotlib'], #external packages as dependencies
   scripts=[
            'SellkeSimulation/simulation_code'
           ]
)