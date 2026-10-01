"""Render a portable, dependency-free HTML report from workflow JSON files."""
import argparse
import html
import json
from pathlib import Path


def table(values):
    return "<table><tbody>" + "".join(
        f"<tr><th>{html.escape(key)}</th><td>{html.escape(str(value) if value is not None else 'N/A (empty denominator)')}</td></tr>"
        for key, value in values.items() if not isinstance(value, dict)
    ) + "</tbody></table>"


def render(paths, label, out):
    sections = []
    for path in paths:
        with open(path) as stream:
            content = json.load(stream)
        relative = Path(path).relative_to(Path(out).parent).as_posix()
        sections.append(f'<section><h2>{html.escape(Path(path).stem)}</h2>' + table(content)
                        + f'<p><a href="{html.escape(relative)}">Source JSON</a></p></section>')
    document = """<!doctype html><html lang="en"><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>Yeast WGS comparison</title>
<style>body{font:16px system-ui;max-width:1000px;margin:40px auto;padding:0 24px;color:#17283a;background:#f5f7fa}
h1,h2{color:#143f64}section{background:white;padding:24px;margin:20px 0;border-radius:8px}
table{border-collapse:collapse;width:100%}th,td{text-align:left;padding:8px;border-bottom:1px solid #ddd}
th{width:55%;overflow-wrap:anywhere}.note{padding:16px;background:#fff2cc}a{color:#174d94}</style>
<h1>Illumina / Nanopore SNP comparison</h1>"""
    document += f"<p>{html.escape(label)}</p>"
    document += '<p class="note">ONT calls are a bcftools SNP baseline. Overlap is not accuracy. '
    document += 'Synthetic data demonstrate execution only. Indels are not compared.</p>'
    document += '<p><a href="multiqc/multiqc_report.html">MultiQC</a> · <a href="provenance/tool_versions.txt">Tool versions</a></p>'
    Path(out).write_text(document + "".join(sections) + "</html>\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--inputs", nargs="+", required=True)
    p.add_argument("--label", required=True)
    p.add_argument("--out", required=True)
    a = p.parse_args()
    render(a.inputs, a.label, a.out)
