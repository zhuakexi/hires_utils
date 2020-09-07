# hires-utils
### tool for scHi-C data cleaning
## Installation
No package now, use git.
```
git clone https://github.com/zhuakexi/hires-utils.git
```
## Usage
Unix style, with subcommand.
Works with python3.
```
usage: hires [-h] command ...

Functions for hires pipline

optional arguments:
  -h, --help      show this help message and exit

These are sub-commands:
  command
    clean_leg     clean promiscuous legs that contacts with multiple legs
    clean_splicing
                  clean exon splicing from mRNA in contact file
    align         caculate rmsd between .3dg replicates
    clean_isolated
                  remove isolated contacts according to L-0.5 distance
```
Type "subcommand -h" for details.
