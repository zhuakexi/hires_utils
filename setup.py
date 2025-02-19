from setuptools import setup, find_packages

setup(
    name='hires_utils',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'pandas',
        'rmsd'
    ],
    entry_points={
        'console_scripts': [
            'hires_utils=hires_utils.__main__:main',  # 如果有命令行入口
        ],
    },
)