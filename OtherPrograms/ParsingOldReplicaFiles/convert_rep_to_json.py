"""
convert_rep_to_json.py

Converts old artemide/Snowflake replica files (the positional, version-gated .rep
text format read by the old ArtemideReplicaSet / SnowflakeReplicaSet classes in
DataProcessor) into the self-describing JSON format read by Cerynia.aTMDeReplicaSet.

Unlike the old parsers, this does not branch on the file's '*V' version number.
Instead it reads the '*B' (and optional '*BB') index block generically: whichever
'*tag : comment' / 'start, end' pairs are actually present in a given file determine
which modules exist and whether collinear-replica indices are present. A file with
only a '*1' (SNOWFLAKE) tag and no '*3' (TMDR) tag is recognised as the old "simple"
Snowflake-only format (single mean replica, no '-1' initial replica, no collinear
block) -- exactly the shape SnowflakeReplicaSet.ReadRepFile used to require.

Usage:
    python convert_rep_to_json.py input.rep [output.json]

    from convert_rep_to_json import convert_rep
    data = convert_rep("ART25_main.rep")   # -> dict in the Cerynia JSON schema
"""

import sys
import json

_MAIN_TAG_MODULE = {
    "*1":  "SNOW",
    "*3":  "TMDR",
    "*4":  "uTMDPDF",
    "*5":  "uTMDFF",
    "*11": "lpTMDPDF",
    "*12": "SiversTMDPDF",
    "*13": "wgtTMDPDF",
}

_COLLINEAR_TAG_MODULE = {
    "*4":  "uTMDPDF",
    "*5":  "uTMDFF",
    "*11": "lpTMDPDF",
    "*12": "SiversTMDPDF",
    "*13": "wgtTMDPDF",
}

# under *BB, tag *1 does not mean SNOW -- it means an external replica reference
# (e.g. "which ART25 replica this SnowART26 replica should be paired with").
_LINK_TAG = "*1"


def _read_index_block(lines):
    """
    Pop '*tag : comment' / 'value' line pairs until the next line is '*C' or '*BB'.
    Returns {tag: int} for the single-value '*0' entry and {tag: (start, end)}
    for range entries, using the raw 1-based values exactly as written in the file.
    """
    block = {}
    while not (lines[0].startswith("*C") or lines[0].startswith("*BB")):
        tag = lines.pop(0).split()[0]
        value_line = lines.pop(0).strip()
        if "," in value_line:
            a, b = value_line.split(",")
            block[tag] = (int(a), int(b))
        else:
            block[tag] = int(value_line)
    return block


def _span(raw):
    """
    Adjust a raw (1-based, id-column-inclusive) (start, end) pair to a 0-based
    slice into the per-replica combined [floats..., ints...] list -- the same
    (start-2, end-1) convention the old ArtemideReplicaSet used. Returns None
    if the tag is inactive in this file (raw == (0, 0)).
    """
    a, b = raw
    if a <= 0:
        return None
    return a - 2, b - 1


def _to_module_blocks(combined, ranges, tag_module):
    out = {}
    for tag, module in tag_module.items():
        raw = ranges.get(tag)
        if raw is None:
            continue
        span = _span(raw)
        if span is None:
            continue
        start, end = span
        out[module] = list(combined[start:end])
    return out


def _parse_replica_line(line, totalLength, npLength):
    fields = [x.strip() for x in line.strip().split(",")]
    if len(fields) != totalLength:
        raise ValueError(f"Expected {totalLength} columns, got {len(fields)}: {line!r}")
    floats = [float(x) for x in fields[1:npLength]]
    ints   = [int(float(x)) for x in fields[npLength:]]
    return floats + ints


def _build_replica_dict(combined, main_ranges, collinear_ranges):
    replica = {module: {"params": values}
               for module, values in _to_module_blocks(combined, main_ranges, _MAIN_TAG_MODULE).items()}

    for module, values in _to_module_blocks(combined, collinear_ranges, _COLLINEAR_TAG_MODULE).items():
        replica.setdefault(module, {})["collinearReplica"] = values

    link_raw = collinear_ranges.get(_LINK_TAG)
    if link_raw is not None:
        span = _span(link_raw)
        if span is not None:
            start, end = span
            replica["linkedReplica"] = list(combined[start:end])

    return replica


def convert_rep(path):
    """Parse an old .rep file and return a dict in the Cerynia aTMDeReplicaSet JSON schema."""
    with open(path) as f:
        lines = f.readlines()

    comment_lines = []
    while not lines[0].startswith("*V"):
        comment_lines.append(lines.pop(0).rstrip("\n"))
    lines.pop(0)
    version = int(lines.pop(0))

    while not lines[0].startswith("*A"):
        lines.pop(0)
    lines.pop(0)
    name = lines.pop(0).strip()

    while not lines[0].startswith("*B"):
        lines.pop(0)
    lines.pop(0)

    main_ranges = _read_index_block(lines)
    totalLength = main_ranges.pop("*0")

    collinear_ranges = {}
    if lines[0].startswith("*BB"):
        lines.pop(0)
        collinear_ranges = _read_index_block(lines)

    pdfLength = sum(b - a + 1 for a, b in collinear_ranges.values() if a > 0)
    npLength  = totalLength - pdfLength

    if not lines[0].startswith("*C"):
        raise ValueError(f"'{path}': expected '*C' section, found {lines[0]!r}")
    lines.pop(0)
    numberOfReplicas = int(lines.pop(0))

    if not lines[0].startswith("*D"):
        raise ValueError(f"'{path}': expected '*D' section, found {lines[0]!r}")
    lines.pop(0)

    # A file with no '*3' (TMDR) tag never went through the full *B/*BB layout --
    # it is the old "simple" Snowflake-only format: one mean replica, no initial.
    is_simple_snowflake = "*3" not in main_ranges

    if is_simple_snowflake:
        mean_combined    = _parse_replica_line(lines.pop(0), totalLength, npLength)
        initial_combined = mean_combined
    else:
        initial_combined = _parse_replica_line(lines.pop(0), totalLength, npLength)
        mean_combined    = _parse_replica_line(lines.pop(0), totalLength, npLength)

    if not lines[0].startswith("*R"):
        raise ValueError(f"'{path}': expected '*R' section, found {lines[0]!r}")
    lines.pop(0)

    replicas_combined = [_parse_replica_line(l, totalLength, npLength)
                          for l in lines[:numberOfReplicas]]

    modules = [module for tag, module in _MAIN_TAG_MODULE.items()
               if _span(main_ranges.get(tag, (0, 0))) is not None]

    return {
        "formatVersion":     1,
        "name":              name,
        "comment":           "\n".join(comment_lines).strip(),
        "modules":           modules,
        "initial":           _build_replica_dict(initial_combined, main_ranges, collinear_ranges),
        "mean":              _build_replica_dict(mean_combined, main_ranges, collinear_ranges),
        "replicas":          [_build_replica_dict(c, main_ranges, collinear_ranges)
                               for c in replicas_combined],
        "convertedFrom":     path,
        "sourceFormatVersion": version,
    }


def main():
    if len(sys.argv) not in (2, 3):
        print("Usage: python convert_rep_to_json.py input.rep [output.json]")
        sys.exit(1)

    in_path  = sys.argv[1]
    out_path = sys.argv[2] if len(sys.argv) == 3 else in_path.rsplit(".", 1)[0] + ".json"

    data = convert_rep(in_path)
    with open(out_path, "w") as f:
        json.dump(data, f, indent=2)

    print(f"Wrote {out_path}: {data['name']!r}, {len(data['replicas'])} replicas, "
          f"modules={data['modules']}")


if __name__ == "__main__":
    main()
