#!/bin/bash

pip install json-schema-for-humans

mkdir docs/build/Schemas/
generate-schema-doc res/schemas/universe.schema.json docs/build/Schemas/universe.schema.html
