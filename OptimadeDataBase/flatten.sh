#!/usr/bin/env bash

mongosh <<'EOF'

use molecular_dynamics_db;
db.structures.updateMany(
  {},
  [
    { $replaceRoot: { newRoot: { $mergeObjects: [ "$$ROOT", "$attributes" ] } } },
    { $unset: "attributes" }
  ]
);
EOF