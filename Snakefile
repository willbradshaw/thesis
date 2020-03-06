# Core user snakefile for thesis assembly

# Config file
configfile: "config.yaml"

# Import subworkflows
include: "snakefiles/setup"
include: "snakefiles/locus"
include: "snakefiles/igseq"
include: "snakefiles/final"
