#!/bin/bash

lib_path=../../lib/collapse;

ln -s ${lib_path}/* .
ls ${lib_path} > to_be_removed;

exit;
