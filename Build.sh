#!/bin/bash



CXX=`root-config --cxx`
CFLAGS=`root-config --cflags`
LDLIBS=`root-config --glibs`
g++ main.C ${LDLIBS} ${CFLAGS} -o main
