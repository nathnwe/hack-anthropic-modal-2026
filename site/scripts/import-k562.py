"""Derive the complete K562 subset, preserving every CSV row and its provenance.
Usage: python site/scripts/import-k562.py input.csv ensembl-lookup.json
The Ensembl response must be lookup/id/ENSG00000136997?expand=1.
"""
import csv, hashlib, json, pathlib, sys
source, annotation = map(pathlib.Path, sys.argv[1:])
rows = [(n, r) for n, r in enumerate(csv.DictReader(source.open()), 2) if r['CellType'] == 'K562' and r['TargetGene'] == 'MYC']
gene = json.loads(annotation.read_text())
tx = next(t for t in gene['Transcript'] if t.get('is_canonical'))
assert gene['id'] == 'ENSG00000136997' and gene['assembly_name'] == 'GRCh38' and tx['strand'] == 1
links = []
for n, r in rows:
    assert r['#chr'] == 'chr8' and r['genome_assembly'] == 'GRCh38'
    assert r['element_coordinate_system'] == 'BED 0-based start, end-exclusive'
    assert 0 <= int(r['start']) < int(r['end']) and 0 <= float(r['Score']) <= 1
    links.append(dict(id=f"K562-{len(links)+1:02}", csv_row=n, chrom=r['#chr'], start=int(r['start']), end=int(r['end']), score=float(r['Score']), model=r['model_variant'], evidence=r['evidence_status'], accession=r['source_file_accession'], url=r['source_file_url'], annotation=r['annotation_description']))
output = dict(kind='re2g-case-study', symbol='MYC', cell_type='K562', assembly='GRCh38', coordinate_system='BED 0-based start, end-exclusive', source_file=source.name, source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(), reference=dict(transcript=f"{tx['id']}.{tx['version']}", start=tx['start']-1, end=tx['end'], strand='+', tss=tx['start']-1, url='https://rest.ensembl.org/lookup/id/ENSG00000136997?expand=1;content-type=application/json', retrieved='2026-10-03', note='Ensembl canonical transcript reference, not a K562-specific measured TSS or the undisclosed rE2G target TSS. Ensembl 1-based coordinates converted to BED.'), links=links)
target=pathlib.Path(__file__).resolve().parents[1]/'public/data/case-studies/myc-k562.json'
target.write_text(json.dumps(output, indent=2)+'\n')
print(f'{len(links)} rows; {len(set((r["start"],r["end"]) for r in links))} unique intervals; {len(set(r["accession"] for r in links))} exports → {target}')
