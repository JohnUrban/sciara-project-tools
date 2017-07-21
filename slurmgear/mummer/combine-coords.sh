#!/bin/bash


cat coords/*bed | sortBed -i - > coords.bed
