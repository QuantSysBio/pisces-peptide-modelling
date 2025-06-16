
from setuptools import setup

setup(
    name='ppm',
    version=0.1,
    description='Downstream PISCES analysis.',
    packages=[
        'ppm',
    ],
    # long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    py_modules=[
        'ppm',
    ],
    entry_points={
        'console_scripts': [
            'ppm=ppm.run:run_ppm'
        ]
    },
    python_requires='>=3.11',
    install_requires=[
        'pisces>=0.1',
        'shap>=0.47.2',
    ],
    project_urls={
        'Homepage': 'https://github.com/QuantSysBio/inSPIRE',
        'Tracker': 'https://github.com/QuantSysBio/inSPIRE/issues',
    },
)
