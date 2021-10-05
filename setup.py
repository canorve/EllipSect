from setuptools import setup, find_packages

with open("README.md", "r") as fh:

    long_description = fh.read()


with open("requirements.txt", "r") as f:
    reqs = [line.rstrip("\n") for line in f if line != "\n"]


setup(

    name="EllipSect", # Replace with your username

    version="2.1.0",

    author="Christopher Añorve, Emmanuel Ríos-López, Omar U. Reyes-Amador and Omar López-Cruz",

    author_email="canorve@gmail.com",

    description="An analysis tool for GALFIT output",

    long_description=long_description,

    long_description_content_type="text/markdown",

    url="https://github.com/canorve/EllipSect",

    packages=find_packages(),

    install_requires=reqs,

    classifiers=[

        "Programming Language :: Python :: 3",

        "License ::  GNU License",

        "Operating System :: OS Independent",

    ],

    python_requires='>=3.8.5',

)


