#!/usr/bin/env python3
"""Aggregate OHAT items into the displayed domains and study summaries."""

import csv
from collections import Counter, OrderedDict
from pathlib import Path

SRC = Path(".cache/data/ohat_items.csv")
OUT = Path("output/numbers")
OUT.mkdir(parents=True, exist_ok=True)

RANK = {"--": 5, "-": 4, "NR": 3, "+": 2, "++": 1}
DOMAINS = OrderedDict([
    ("Selection bias", ["Q1", "Q2"]),
    ("Confounding bias", []),
    ("Performance bias", ["Q6"]),
    ("Attrition/exclusion bias", ["Q7"]),
    ("Detection bias", ["Q8", "Q9_IS", "Q9_RT", "Q9_PE"]),
    ("Selective reporting bias", ["Q10"]),
    ("Other sources of bias", ["Q11"]),
])

with SRC.open(encoding="utf-8") as fh:
    items = list(csv.DictReader(fh))

units = OrderedDict()
for r in items:
    units.setdefault((int(r["journal_N"]), r["display_label"]), {})[r["question"]] = r["analysed_rating"]

with (OUT / "ohat_domain_matrix.csv").open("w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh)
    w.writerow(["journal_N", "Study", *DOMAINS])
    for (j, label), q in sorted(units.items()):
        row = []
        for questions in DOMAINS.values():
            rated = [q[k] for k in questions if q[k] != "N/A"]
            row.append(max(rated, key=RANK.get) if rated else "N/A")
        w.writerow([j, label, *row])

studies = OrderedDict()
for r in items:
    studies.setdefault((int(r["journal_N"]), r["study"]), []).append(r["analysed_rating"])

with (OUT / "ohat_study_summary.csv").open("w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh)
    w.writerow(["journal_N", "Study", "N_applicable", "N_++", "N_+", "N_-", "N_--", "N_NR",
                "N_NA", "Low_risk_count", "High_concern_count", "Low_risk_proportion",
                "High_concern_proportion"])
    for (j, study), ratings in sorted(studies.items()):
        c = Counter(ratings)
        applicable = len(ratings) - c["N/A"]
        low = c["++"] + c["+"]
        high = c["-"] + c["--"] + c["NR"]
        w.writerow([j, study, applicable, c["++"], c["+"], c["-"], c["--"], c["NR"], c["N/A"],
                    low, high, low / applicable, high / applicable])

print(f"OHAT display inputs written to {OUT}: {len(units)} units, {len(studies)} studies")
