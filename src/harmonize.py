import argparse
import pandas as pd
import os
import sys

UNIFIED_COLS = [
    "sample", "tool", "gene", "drug_class",
    "identity", "coverage", "contig", "start", "stop"
]


def parse_amrfinder(path, sample):
    df = pd.read_csv(path, sep="\t")
    if df.empty:
        return pd.DataFrame(columns=UNIFIED_COLS)
    return pd.DataFrame({
        "sample":     sample,
        "tool":       "amrfinder",
        "gene":       df["Element symbol"],
        "drug_class": df["Class"].fillna("") + " - " + df["Subclass"].fillna(""),
        "identity":   pd.to_numeric(df["% Identity to reference"], errors="coerce"),
        "coverage":   pd.to_numeric(df["% Coverage of reference"], errors="coerce"),
        "contig":     df["Contig id"],
        "start":      df["Start"],
        "stop":       df["Stop"],
    })


def parse_rgi(path, sample):
    df = pd.read_csv(path, sep="\t")
    if df.empty:
        return pd.DataFrame(columns=UNIFIED_COLS)
    df = df[df["Cut_Off"] != "Loose"]
    if df.empty:
        return pd.DataFrame(columns=UNIFIED_COLS)
    return pd.DataFrame({
        "sample":     sample,
        "tool":       "rgi",
        "gene":       df["Best_Hit_ARO"],
        "drug_class": df["Drug Class"],
        "identity":   pd.to_numeric(df["Best_Identities"], errors="coerce"),
        "coverage":   pd.to_numeric(df["Percentage Length of Reference Sequence"], errors="coerce"),
        "contig":     df["Contig"],
        "start":      df["Start"],
        "stop":       df["Stop"],
    })


def parse_resfinder(path, sample):
    df = pd.read_csv(path, sep="\t")
    if df.empty:
        return pd.DataFrame(columns=UNIFIED_COLS)
    df = df.dropna(subset=["Resistance gene"])
    if df.empty:
        return pd.DataFrame(columns=UNIFIED_COLS)
    positions = df["Position in contig"].str.extract(r"(\d+)\.\.(\d+)")
    identity = df["Identity"].astype(str).str.replace("%", "", regex=False)

    def cov_to_pct(val):
        try:
            parts = str(val).split("/")
            if len(parts) == 2:
                return round(int(parts[0]) / int(parts[1]) * 100, 2)
            return float(val)
        except:
            return None

    coverage = df["Alignment Length/Gene Length"].apply(cov_to_pct)

    return pd.DataFrame({
        "sample":     sample,
        "tool":       "resfinder",
        "gene":       df["Resistance gene"],
        "drug_class": df["Phenotype"],
        "identity":   pd.to_numeric(identity, errors="coerce"),
        "coverage":   coverage,
        "contig":     df["Contig"],
        "start":      pd.to_numeric(positions[0], errors="coerce"),
        "stop":       pd.to_numeric(positions[1], errors="coerce"),
    })

def parse_abricate(path, sample):
    df = pd.read_csv(path, sep="\t")
    if df.empty:
        return pd.DataFrame(columns=UNIFIED_COLS)
    return pd.DataFrame({
        "sample":     sample,
        "tool":       "abricate",
        "gene":       df["GENE"],
        "drug_class": df["RESISTANCE"],
        "identity":   pd.to_numeric(df["%IDENTITY"], errors="coerce"),
        "coverage":   pd.to_numeric(df["%COVERAGE"], errors="coerce"),
        "contig":     df["SEQUENCE"],
        "start":      df["START"],
        "stop":       df["END"],
    })


PARSERS = {
    "amrfinder": parse_amrfinder,
    "rgi":       parse_rgi,
    "resfinder": parse_resfinder,
    "abricate":  parse_abricate,
}


def detect_tool(filename):
    fname = os.path.basename(filename).lower()
    for tool in PARSERS:
        if tool in fname:
            return tool
    return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--input_files", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    frames = []
    for f in args.input_files:
        tool = detect_tool(f)
        if tool is None:
            print(f"WARNING: could not detect tool for {f}, skipping", file=sys.stderr)
            continue
        try:
            df = PARSERS[tool](f, args.sample)
            frames.append(df)
            print(f"  Parsed {tool}: {len(df)} hits", file=sys.stderr)
        except Exception as e:
            print(f"WARNING: failed to parse {f}: {e}", file=sys.stderr)

    out = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=UNIFIED_COLS)
    for col in ["identity", "coverage", "start", "stop"]:
        out[col] = pd.to_numeric(out[col], errors="coerce")
        
    out.to_csv(args.output, sep="\t", index=False)
    print(f"Written {len(out)} rows to {args.output}")


if __name__ == "__main__":
    main()