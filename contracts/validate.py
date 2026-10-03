"""Validate every JSON record against contracts/schema.json.

Usage: uv run python contracts/validate.py [paths...]   (default: fixtures + site/public/data)
"""
import json
import sys
from pathlib import Path

from jsonschema import Draft202012Validator

ROOT = Path(__file__).resolve().parent.parent
schema = json.loads((ROOT / "contracts" / "schema.json").read_text())
validator = Draft202012Validator(schema)

paths = [Path(p) for p in sys.argv[1:]] or [
    *(ROOT / "contracts" / "fixtures").glob("*.json"),
    *(ROOT / "site" / "public" / "data").glob("*.json"),
]
failed = 0
for p in paths:
    errors = list(validator.iter_errors(json.loads(p.read_text())))
    for e in errors:
        print(f"FAIL {p.name}: {'/'.join(map(str, e.path))}: {e.message}")
    failed += bool(errors)
    print(f"{'ok  ' if not errors else 'FAIL'} {p.name}")
sys.exit(1 if failed else 0)
