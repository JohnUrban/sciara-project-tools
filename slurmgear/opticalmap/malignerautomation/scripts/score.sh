#!/bin/bash

tail -n +2 ${ALL} | cut -f 19 | awkSum > score.txt
