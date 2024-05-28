from setuptools import setup, find_packages

from os import path

script_directory = path.abspath(path.dirname(__file__))

version = None
with open(path.join(script_directory, 'metabolic_niche_space', '__init__.py')) as f:
    for line in f.readlines():
        line = line.strip()
        if line.startswith("__version__"):
            version = line.split("=")[-1].strip().strip('"')
assert version is not None, "Check version in metabolic_niche_space/__init__.py"

with open(path.join(script_directory, 'README.md')) as f:
    long_description = f.read()

requirements = list()
with open(path.join(script_directory, 'requirements.txt')) as f:
    for line in f.readlines():
        line = line.strip()
        if line:
            if not line.startswith("#"):
                requirements.append(line)

setup(
    name='metabolic-niche-space',

    python_requires='>=3.6',
    version=version,
    description='Metabolic niche space analysis',
    long_description=long_description,

    # The project's main homepage.
    url='https://github.com/jolespin/metabolic_niche_space',

    # Author details
    author='Josh L. Espinoza',

    # Choose your license
    license='GPLv3',
    provides=['metabolic_niche_space'],
    # packages=find_package/s(),
    py_modules=["manifold", "neighbors"],
    
    install_requires=requirements, #[:-1],
    tests_require=requirements, #[-1:]
)