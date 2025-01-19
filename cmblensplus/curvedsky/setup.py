#!/usr/bin/env python

# Made by Miguel Ruiz Granda
from setuptools import setup

setup(
    name='curvedsky',
    version='1.0',
    packages=['curvedsky'],
    package_data={
        'curvedsky': ['../libcurvedsky.so'],
    },
    include_package_data=True,
)

