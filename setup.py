from setuptools import setup, find_packages


setup(
    name="pybmft-c",
    version="0.0.1.dev0",
    description="Bay-Marsh-Forest Transect Model - Python Version",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    author="Ian Reeves",
    author_email="reevesi@live.unc.edu",
    url="https://github.com/irbreeves/PyBMFT-C",
    classifiers=[
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GPLv3 License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    keywords=["coastal", "marsh", "barrier", "earth science"],
    install_requires=open("requirements.txt", "r").read().splitlines(),
    include_package_data=True,
    packages=find_packages(),
    python_requires=">=3.6",
)
