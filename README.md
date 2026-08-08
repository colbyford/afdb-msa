# AFDB MSA

Paste a protein sequence, get back an `a3m`.

The alignment is not computed. AlphaFold DB publishes the MSA it used for every
entry, so this finds the AFDB entry closest to your query and **re-indexes its
alignment onto your query's residue numbering**. Give it more than one sequence
and it also emits a species-paired alignment for complex prediction.

No MMseqs2, no HMMER, no BLAST, no server. The page is static files; the search
is an index served as static byte ranges.

## One path

```
sequence -> seeds -> index -> candidates -> Smith-Waterman -> best donor
         -> that donor's AFDB MSA -> re-indexed onto your query
```

There is no fallback and nothing to configure. Earlier versions tried a UniParc
checksum shortcut first and EBI BLAST after; neither earned its complexity. The
index already returns exact matches at 100% identity, and below its range a
borrowed alignment is not worth having -- so when nothing close exists, that is
reported rather than papered over by a slower search that would not do better.

## How the borrowing works

An AFDB a3m has one match state per residue of *its own* entry. Aligning the
query to that entry gives exactly the map needed to renumber:

| in the pairwise alignment | in the output a3m |
|---|---|
| query residue aligned to hit residue | hit column kept, now a query column |
| hit residue over a query gap | demoted to a lowercase insertion |
| query residue over a hit gap | all-gap column in every borrowed sequence |
| query outside the alignment | all-gap column |

Every output row carries exactly `len(query)` match states. That invariant is
what the tests check hardest.

## The index

A minimizer index over **239,602,633 AFDB entries** (99.4% of v6; the rest fall
outside a 40-2000 aa filter).

```
  post-000..009.bin   14.2 GB   5 bytes per posting: 1 B key residual + 4 B entry id
  buckets.u32        134 MB     prefix offsets, seed's high bits -> posting range
  acc.bin           2.88 GB     fixed 12-byte accessions, so entry id IS the offset
  ------------------------------
  TOTAL             ~17.2 GB
```

A query never downloads any of it. Twelve seeds, each costing an 8-byte ranged
read of `buckets.u32` and one ranged read of its posting list:

```
  ~45 ranged requests, ~50 KB, ~1 ms -- independent of index size
```

### Why minimizer seeding

At 90% identity a 10-mer survives intact with probability 0.9^10 = 0.35, so
twelve seeds give 1 - 0.65^12 = 99.2% odds of sharing at least one. Exact seeds
are cheap to index and cheap to look up, and they fail gracefully: sensitivity
fades below ~70% identity, which is where borrowing an MSA stops being sound
anyway.

The window is chosen per sequence to land near twelve seeds regardless of
length. A fixed window undersamples short proteins -- at k=12 w=16 a 142-residue
globin got 8 seeds and shared none with a real 90% relative on two of three
tries.

### Measured against BLAST

16 AFDB sequences with 5% of positions mutated, against BLAST on full UniProtKB:

```
  hits at >=90% identity     index 100%    BLAST  81%
  median identity            index  95%    BLAST  95%
  time                       index   1 ms  BLAST 287 s
  read per query             ~50 KB of a 17 GB index
```

The index also wins outright on some queries -- 94% vs 41% on one. That is not
cleverness, it is the corpus: **63% of AFDB entries have been deleted from
current UniProtKB**, so BLAST cannot return them and settles for a distant
relative. AFDB is its own authority here, which is also why indexing it directly
beat routing through UniRef.

### What it does not do

Sensitivity below ~70% identity. A query with no close relative in AFDB has no
MSA to borrow, and no index size fixes that -- a de novo designed sequence has no
natural homologs, so the correct answer really is "nothing found". Those queries
fall through to BLAST.

Low-complexity and repetitive sequence yields few distinct seeds (a homopolymer
collapses to one), so such queries retrieve weakly. A thin candidate list should
not be read as "nothing similar exists".

## Pairing

Multi-chain queries produce a ColabFold-layout paired a3m. AFDB headers carry
taxonomy inline (`... Tax=Monodon monoceros TaxID=40151 ...`), so pairing needs
no extra lookups: keep the best row per taxid per chain, concatenate the taxa
present in every chain, then stack each chain's leftovers block-diagonally.

```
#142,147	1,1
>101
MVLSPADKT...(chain A)...KYRMVHLTPE...(chain B)...AHKYH
```

On haemoglobin alpha/beta that yields ~613 species-paired rows.

## Filters

The a3m comes back whole, ordered by identity to the query. Coverage, identity
and redundancy filters live in `filter.js` and are deliberately not exposed:
they are lossy, the right thresholds depend on what you do with the alignment
next, and a wrong guess quietly discards sequences. Filter downstream, where the
requirement is known.

Their semantics follow [GREMLIN-GUI](https://github.com/sokrypton/GREMLIN-GUI)'s
`msa.js`, so numbers are comparable with that tool: coverage is the non-gap
fraction of a row, identity is over all columns, and redundancy is greedy
clustering. Neff@0.8 is reported alongside raw depth in the results, because raw
depth flatters a redundant alignment.

One detail that bit: greedy clustering keeps whichever cluster member it meets
first, so it is always run over identity-sorted rows -- otherwise the survivor is
an accident of input order. Measured on a 600-row MSA, the two orders kept 13
different members out of ~500.

## Files

| file | role |
|---|---|
| `msa.js` | a3m parsing, the hit->query transfer, pairing |
| `filter.js` | coverage / identity / redundancy filters, Neff |
| `align.js` | Smith-Waterman / BLOSUM62 -- supplies the alignment BLAST would have |
| `seeds.js` | minimizer seeding, shared verbatim by builder and browser |
| `search.js` | client-side lookup over the index via HTTP Range |
| `api.js` | AFDB sequence and MSA fetch -- the only remote source |
| `app.js` | UI wiring only |
| `tools/afdb-shard.mjs` | one shard of the index: FASTA -> sorted postings |
| `tools/afdb-merge.mjs` | merge shards into the published layout |
| `tools/index-vs-search.mjs` | scores the index against BLAST |
| `tools/afdb-coverage.mjs` | how much of AFDB space a given index size covers |
| `tools/minimizer-broad.mjs` | recall by identity band across many families |
| `tools/select-entries.mjs` | depth-ranked entry selection (unused: indexing everything won) |

## Tests

```sh
node test/test.mjs                  # 86 assertions; caches AFDB fixtures in test/data
node test/dom.mjs                   # 31; no network
node test/index-e2e.mjs <indexDir>  # 16; the whole path over real HTTP Range
python3 test/browser.py             # 22; the page in a real browser
```

## Building the index

```sh
# 118 GB of sequences.fasta, 24 parts
for i in $(seq 0 23); do
  node tools/afdb-shard.mjs part_$(printf '%02d' $i).fa shards/ $i &
done; wait
node --max-old-space-size=120000 tools/afdb-merge.mjs shards/ index/
```

About 5 minutes of sharding on 24 cores and 2.5 minutes to merge.

## Deploying

The page goes on GitHub Pages: Settings -> Pages -> deploy from branch, root.
There is no build step -- 68 KB of static files.

One thing to run before committing a change to any `.js` or `.css`:

```sh
node tools/stamp.mjs
```

It rewrites the tags in `index.html` to carry a content hash
(`app.js?v=d0262328`). Pages serves everything with `cache-control: max-age=600`,
so without it a browser can hold a stale script for ten minutes and pair it with
a fresh page -- which is worse than either alone, and crashed for real once.
`test/dom.mjs` fails if a hash is missing or does not match its file.

The index is 17 GB and lives on
[Hugging Face](https://huggingface.co/datasets/sokrypton/afdb-msa-index), which
serves `206` with `Content-Range` and CORS -- all the client needs. `index.html`
points at it by default; blank the field and the site falls back to BLAST.

## Credits

The AFDB precomputed-MSA endpoint is the one used by
[py2Dmol](https://github.com/sokrypton/py2Dmol); filtering semantics come from
[GREMLIN-GUI](https://github.com/sokrypton/GREMLIN-GUI). Alignments come from
[AlphaFold DB](https://alphafold.ebi.ac.uk/) (EMBL-EBI / DeepMind, CC-BY-4.0);
the BLAST fallback runs on [EBI Job Dispatcher](https://www.ebi.ac.uk/jdispatcher/).
