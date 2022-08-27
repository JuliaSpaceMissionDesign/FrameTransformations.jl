#!/bin/bash

pip install json-schema-for-humans

mkdir docs/universe
generate-schema-doc universe.schema.json docs/universe/universe.html
