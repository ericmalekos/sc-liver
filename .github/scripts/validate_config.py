#!/usr/bin/env python
"""Validate a pipeline config YAML against the JSON schema. Used by CI and locally.

    python .github/scripts/validate_config.py config/config.yaml config/schemas/config.schema.yaml
"""
import sys

import jsonschema
import yaml


def main(cfg_path: str, schema_path: str) -> int:
    with open(cfg_path) as fh:
        cfg = yaml.safe_load(fh)
    with open(schema_path) as fh:
        schema = yaml.safe_load(fh)
    jsonschema.validate(instance=cfg, schema=schema)
    print(f"OK: {cfg_path} validates against {schema_path}")
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(2)
    sys.exit(main(sys.argv[1], sys.argv[2]))
