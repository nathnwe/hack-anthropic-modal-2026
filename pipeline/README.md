# pipeline/  (Operator A, branch `feat/pipeline`)
Gene -> chromosome -> TAD -> CREs (ENCODE-rE2G) -> contacts (UniversalEPI) -> candidate staples.
Output: `locus` and `staples` sections of the per-gene record, per `contracts/schema.json`.
Use precomputed model outputs where possible; Modal only if GPU inference is needed.
