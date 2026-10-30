#!/usr/bin/env python3
"""Step 5b: derive L5 (phylogenetic series membership) from a support-annotated tree.

L5 is 'yes' when the query sequence is contained in a clade that
  (1) contains exactly one of the reference sequences, and
  (2) is supported at or above --min-support.
The smallest such clade is used; its single reference defines Series_assignment.
When no such clade exists, L5 is 'no' and Series_assignment is 'unassigned'.
Sequences absent from the tree are reported as '-' (not evaluable), never 'no'.

This script performs level assignment only. It makes no classification decision
and emits no category label; Step 6 is the sole classifier.
"""
import argparse
import csv
import sys
from pathlib import Path

NOT_EVALUABLE = '-'


def parse_newick(text):
    """Minimal Newick reader: quoted labels, branch lengths, numeric node support."""
    s = text.strip().rstrip(';')
    i = 0

    def label():
        nonlocal i
        if i < len(s) and s[i] == "'":
            i += 1
            start = i
            while i < len(s) and s[i] != "'":
                i += 1
            out = s[start:i]
            i += 1
            return out
        start = i
        while i < len(s) and s[i] not in '(),:;':
            i += 1
        return s[start:i]

    def skip_length():
        nonlocal i
        if i < len(s) and s[i] == ':':
            i += 1
            while i < len(s) and (s[i].isdigit() or s[i] in '.eE+-'):
                i += 1

    def node():
        nonlocal i
        name = support = None
        children = []
        if i < len(s) and s[i] == '(':
            i += 1
            while True:
                children.append(node())
                if i < len(s) and s[i] == ',':
                    i += 1
                    continue
                if i < len(s) and s[i] == ')':
                    i += 1
                break
            text_label = label()
            if text_label:
                try:
                    support = float(text_label)
                except ValueError:
                    name = text_label
        else:
            name = label()
        skip_length()
        return {'name': name, 'support': support, 'children': children}

    root = node()
    if i < len(s):
        raise ValueError('Unparsed trailing text in Newick at offset %d' % i)
    return root


def collect_clades(node, out):
    """Post-order traversal; appends (leafset, support) for every internal node."""
    if not node['children']:
        if not node['name']:
            raise ValueError('Unnamed leaf in tree')
        return frozenset([node['name']])
    leaves = frozenset()
    for child in node['children']:
        leaves |= collect_clades(child, out)
    out.append((leaves, node['support']))
    return leaves


def refs_in(leaves, refmap):
    return sorted(name for name, pattern in refmap.items()
                  if any(pattern.lower() in leaf.lower() for leaf in leaves))


def assign(query, clades, refmap, min_support):
    """Return (L5, Series_assignment, Series_evidence) for one query leaf."""
    qualifying = None
    nearest = None
    for leaves, support in clades:
        if query not in leaves:
            continue
        present = refs_in(leaves, refmap)
        if not present:
            continue
        # Smallest clade sharing a reference with the query, support aside.
        if nearest is None or len(leaves) < len(nearest[0]):
            nearest = (leaves, support, present)
        if support is None or support < min_support:
            continue
        if len(present) == 1 and (qualifying is None or len(leaves) < len(qualifying[0])):
            qualifying = (leaves, support, present[0])
    if qualifying:
        leaves, support, ref = qualifying
        return 'yes', ref, ('monophyletic with %s; UFBoot %g; clade n=%d; '
                            'no other reference in clade' % (ref, support, len(leaves)))
    if nearest:
        leaves, support, present = nearest
        return 'no', 'unassigned', (
            'smallest clade sharing a reference: n=%d, UFBoot %s, reference(s) %s '
            '(requires exactly one reference and UFBoot >= %g)'
            % (len(leaves), 'none' if support is None else '%g' % support,
               ' + '.join(present), min_support))
    return 'no', 'unassigned', 'no clade containing a reference includes this sequence'


def read_sequence_list(path):
    with open(path, newline='', encoding='utf-8') as fh:
        head = fh.readline().strip()
        fh.seek(0)
        if '\t' in head or head.lower() == 'sequence':
            rows = csv.DictReader(fh, delimiter='\t')
            if 'Sequence' not in (rows.fieldnames or []):
                raise ValueError('--sequences TSV lacks a Sequence column')
            return [r['Sequence'].strip() for r in rows if r['Sequence'].strip()]
        return [line.strip() for line in fh if line.strip()]


def main():
    p = argparse.ArgumentParser(description='Step 5b: L5 assignment from a phylogenetic tree')
    p.add_argument('--tree', required=True, help='Newick tree with node support (e.g. IQ-TREE .contree)')
    p.add_argument('--out', required=True, help='Output TSV')
    p.add_argument('--reference-map', default='THCAS=Q8GTB6,CBDAS=A6P6V9,CBCAS=Q33DQ2',
                   help='NAME=PATTERN pairs matched case-insensitively against leaf labels')
    p.add_argument('--outgroup', default='BBE', help='Substring identifying outgroup leaves to skip')
    p.add_argument('--min-support', type=float, default=95.0,
                   help='Minimum node support for a clade to qualify (default: UFBoot 95)')
    p.add_argument('--sequences', help='Optional TSV (Sequence column) or plain list; '
                                       'IDs absent from the tree are emitted as not evaluable')
    args = p.parse_args()

    refmap = {}
    for item in args.reference_map.split(','):
        if '=' not in item:
            p.error('Malformed --reference-map entry: ' + repr(item))
        name, pattern = item.split('=', 1)
        refmap[name.strip()] = pattern.strip()

    root = parse_newick(Path(args.tree).read_text(encoding='utf-8'))
    clades = []
    leaves = collect_clades(root, clades)

    is_ref = {leaf: refs_in([leaf], refmap) for leaf in leaves}
    queries = [leaf for leaf in sorted(leaves)
               if not is_ref[leaf] and args.outgroup.lower() not in leaf.lower()]
    if not queries:
        p.error('No query leaves in tree after removing references and outgroup')

    rows = []
    for q in queries:
        l5, series, evidence = assign(q, clades, refmap, args.min_support)
        rows.append(dict(Sequence=q, L5=l5, Series_assignment=series, Series_evidence=evidence))

    if args.sequences:
        known = {r['Sequence'] for r in rows}
        for seq in read_sequence_list(args.sequences):
            if seq in known or any(seq in leaf or leaf in seq for leaf in leaves):
                continue
            rows.append(dict(Sequence=seq, L5=NOT_EVALUABLE, Series_assignment='',
                             Series_evidence='not included in the phylogenetic analysis'))
        rows.sort(key=lambda r: r['Sequence'])

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, 'w', newline='', encoding='utf-8') as fh:
        wr = csv.DictWriter(fh, delimiter='\t', lineterminator='\n',
                            fieldnames=['Sequence', 'L5', 'Series_assignment', 'Series_evidence'])
        wr.writeheader()
        wr.writerows(rows)

    counts = {}
    for r in rows:
        counts[r['L5']] = counts.get(r['L5'], 0) + 1
    print('rows:', len(rows), 'L5:', counts, file=sys.stderr)
    for r in rows:
        print('  %-34s L5=%-3s %s' % (r['Sequence'], r['L5'], r['Series_assignment']), file=sys.stderr)


if __name__ == '__main__':
    main()
