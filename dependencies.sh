#!/bin/bash
# Assumes Ubuntu 22.04

apt -q update && apt install -q -y build-essential libconfig++-dev libmagick++-6.q16-dev bear

echo '# add this path to your PATH'
echo export PATH="/usr/lib/x86_64-linux-gnu/ImageMagick-6.9.11/bin-q16/:${PATH}"
