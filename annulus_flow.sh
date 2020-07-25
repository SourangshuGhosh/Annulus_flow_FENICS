#! /bin/bash
#
python3 annulus_flow.py > annulus_flow.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
#
echo "Normal end of execution."
