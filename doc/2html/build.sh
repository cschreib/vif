#!/bin/bash

# Define output options
OUT_DIR=../html/

# Initialize
mkdir -p $OUT_DIR
cp pygments.css     $OUT_DIR/pygments.css
cp section-link.png $OUT_DIR/section-link.png
cp vectorized.png   $OUT_DIR/vectorized.png

# Parse LaTeX files and build HTML
for f in ../latex/category_*.tex; do
    ophy++ 2html_category $f out_dir=$OUT_DIR
done

for f in ../latex/functions_*.tex; do
    ophy++ 2html_functions $f out_dir=$OUT_DIR
done

# Compile all parsing into the menu
ophy++ make_index $OUT_DIR

# Rebuild HTML to integrate new menu
for f in ../latex/category_*.tex; do
    ophy++ 2html_category $f out_dir=$OUT_DIR
done

for f in ../latex/functions_*.tex; do
    ophy++ 2html_functions $f out_dir=$OUT_DIR
done
