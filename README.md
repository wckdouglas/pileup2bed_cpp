# pileup2bed

This is a program to parse *samtools mpileup* results to a tab-deliminated file storing base count at each position, this also allows one to choose a quality cut off and coverage cut off to reportt clone git@github.com:wckdouglas/pileup2bed.git

to install

```
git clone git@github.com:wckdouglas/pileup2bed.git
cd pileup2bed
make
```

usage

```
./bin/pileup2bed <filename>|<stdin> <quality threshold> <coverage threshold>
```
