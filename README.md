# Project_PBSketch
The source codes of PB-sketch. PB-sketch can achieve point queries, top-k queries, and heavy change queries in the data streams.

## Introduction
XY-Sketch is a state-of-the-art technology to approximate estimate the frequency of data items in streams. In particular, the estimated performance of XY-Sketch is excellent when the space budget is small.

## Requirements
This project contains CF( Cold Filter ), that is, you should meet some of the configuration requirements of CF.

CF use SIMD instructions to achieve high speed, so the cpu must support SSE2 instruction set. cmake >= 2.6 g++ (MSVC is not supported currently.)

## How to build
This project is built upon cmake. You can use the following commands to build and run.

mkdir build cd build cmake .. make cd .. cd bin ./demo

## Dataset
You can download "kosarak.dat" at the following URL.

http://fimi.uantwerpen.be/data/
