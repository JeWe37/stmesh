#!/bin/bash
cmake --preset $1
cmake --build --preset $1 --parallel 4
