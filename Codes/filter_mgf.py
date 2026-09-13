import argparse
import gzip
import re
import sys
from pathlib import Path

DATA_DIR = Path("/Users/fateme/Desktop/unannotated_data/main")

#input
PAIRS = [
    ("cluster_ident_2.mgf", "cluster_ident_2.mztab", "cluster_unannotated_2.mgf"),
    ("cluster_ident_n.mgf", "cluster_ident_n.mztab", "cluster_unannotated_n.mgf"),
]

INDEX_BASE = 0          
DROP_DECOY_PSMS = True  
QVALUE_COLUMN = None   
QVALUE_MAX = None       
PREVIEW_N = 5           
IO_BUFFER = 64 * 1024 * 1024

REF_RE = re.compile(r"(?:index|scan)=(\d+)")


def smart_open(path, mode="rt"):
    opener = gzip.open if str(path).endswith(".gz") else open
    return opener(path, mode, buffering=IO_BUFFER) if opener is open else opener(path, mode)


def parse_mztab(path, preview=False):
    indices, scans, titles = set(), set(), set()
    n_psm = 0
    ms_runs = set()
    col = {}
    shown = 0

    with smart_open(path) as fin:
        for line in fin:
            if line.startswith("PSH"):
                fields = line.rstrip("\n").split("\t")
                col = {name: i for i, name in enumerate(fields)}
                if preview:
                    print(f"  PSH columns: {fields}")
            elif line.startswith("PSM"):
                fields = line.rstrip("\n").split("\t")
                n_psm += 1

                if DROP_DECOY_PSMS:
                    dec = next((v for k, v in col.items() if "decoy" in k.lower()), None)
                    if dec is not None and fields[dec].strip() in {"1", "true", "TRUE"}:
                        continue

                qcol = QVALUE_COLUMN or next(
                    (k for k in col if "q_value" in k.lower() or "qvalue" in k.lower()), None
                )
                if QVALUE_MAX is not None and qcol in col:
                    try:
                        if float(fields[col[qcol]]) > QVALUE_MAX:
                            continue
                    except ValueError:
                        pass

                ref = fields[col["spectra_ref"]] if "spectra_ref" in col else ""
                for part in ref.split("|"):
                    ms_runs.add(part.split(":")[0].strip())
                    m = REF_RE.search(part)
                    if not m:
                        continue
                    value = int(m.group(1))
                    if "index=" in part:
                        indices.add(value - INDEX_BASE)
                    else:
                        scans.add(value)

                tcol = next((k for k in col if "title" in k.lower()), None)
                if tcol is not None:
                    titles.add(fields[col[tcol]].strip())

                if preview and shown < PREVIEW_N:
                    print(f"  PSM spectra_ref: {ref}")
                    shown += 1

    return indices, scans, titles, n_psm, ms_runs


def spectrum_is_identified(idx, title, scan, indices, scans, titles):
    if indices and idx in indices:
        return True
    if scans and scan is not None and scan in scans:
        return True
    if titles and title is not None and title in titles:
        return True
    return False


def filter_mgf(mgf_in, mgf_out, indices, scans, titles, preview=False):
    total = kept = dropped = 0
    state = 0            
    buf, title, scan, idx, keep = [], None, None, -1, False

    fout = None if preview else open(mgf_out, "w", buffering=IO_BUFFER)
    with smart_open(mgf_in) as fin:
        for line in fin:
            if state == 0:
                if line.startswith("BEGIN IONS"):
                    idx, total = total, total + 1
                    buf, title, scan, state = [line], None, None, 1
                continue

            if state == 1:
                if line[:1].isdigit() or line.startswith("END IONS"):
                    keep = not spectrum_is_identified(idx, title, scan, indices, scans, titles)
                    if preview and idx < PREVIEW_N:
                        print(f"  spectrum {idx}: TITLE={title} SCANS={scan} keep={keep}")
                    if keep and fout:
                        fout.writelines(buf)
                        fout.write(line)
                    state = 2
                    if line.startswith("END IONS"):
                        kept, dropped = kept + keep, dropped + (not keep)
                        state = 0
                    continue
                buf.append(line)
                if line.startswith("TITLE="):
                    title = line[6:].strip()
                elif line.startswith("SCANS="):
                    try:
                        scan = int(line[6:].strip())
                    except ValueError:
                        scan = None
                continue

            if keep and fout:
                fout.write(line)
            if line.startswith("END IONS"):
                kept, dropped = kept + keep, dropped + (not keep)
                state = 0

            if preview and total > PREVIEW_N * 20:
                break

    if fout:
        fout.close()
    return total, kept, dropped


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--preview", action="store_true",
                    help="print the first PSM refs and MGF titles, write nothing")
    args = ap.parse_args()

    for mgf_name, mztab_name, out_name in PAIRS:
        mgf_in = DATA_DIR / mgf_name
        mztab_in = DATA_DIR / mztab_name
        mgf_out = DATA_DIR / out_name
        print(f"\n=== {mgf_name} ===", flush=True)

        for p in (mgf_in, mztab_in):
            if not p.exists():
                sys.exit(f"missing file: {p}")

        indices, scans, titles, n_psm, ms_runs = parse_mztab(mztab_in, preview=args.preview)
        print(f"  PSM rows: {n_psm:,} | unique index refs: {len(indices):,} | "
              f"scan refs: {len(scans):,} | title refs: {len(titles):,} | ms_runs: {sorted(ms_runs)}",
              flush=True)
        if not (indices or scans or titles):
            sys.exit("no usable spectrum references found in the mzTab")
        if len(ms_runs) > 1:
            print("  warning: several ms_run entries; indices may not be unique per file")

        total, kept, dropped = filter_mgf(mgf_in, mgf_out, indices, scans, titles,
                                          preview=args.preview)
        if args.preview:
            print("  preview only, no output written")
            continue

        print(f"  spectra read: {total:,} | kept (unannotated): {kept:,} | removed: {dropped:,}")
        expected = len(indices) + len(scans) + len(titles)
        if dropped < 0.5 * min(expected, total):
            print("  WARNING: far fewer spectra removed than PSMs present - check INDEX_BASE "
                  "or the reference type with --preview")
        print(f"  written: {mgf_out}", flush=True)


if __name__ == "__main__":
    main()
