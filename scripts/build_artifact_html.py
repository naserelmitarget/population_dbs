#!/usr/bin/env python3
"""Build single-file HTML pages (figures embedded, repo-relative links shown as text)
from the Markdown reports, for sharing outside the repository.
usage: build_artifact_html.py <report.md> <title> <out.html>"""
import base64, html, mimetypes, os, re, subprocess, sys

md, title, out = sys.argv[1:4]
here = os.path.dirname(os.path.abspath(__file__))
src_dir = os.path.dirname(os.path.abspath(md))
css = open(os.path.join(src_dir, "publication.css")).read()
page = subprocess.run(
    ["pandoc", os.path.abspath(md), "-s", "--toc", "--toc-depth=2",
     "--template", os.path.join(here, "templates/artifact.html"), "--metadata", f"pagetitle={title}"],
    check=True, capture_output=True, text=True).stdout.replace("/*STYLES*/", css)

def embed(m):
    p = os.path.normpath(os.path.join(src_dir, html.unescape(m.group(2))))
    if not os.path.exists(p):
        sys.exit(f"missing figure {p}")
    mime = mimetypes.guess_type(p)[0] or "image/png"
    return f'{m.group(1)}data:{mime};base64,{base64.b64encode(open(p, "rb").read()).decode()}"'

page = re.sub(r'(<img[^>]*?src=")(?!data:|https?:)([^"]+)"', embed, page)
# links into the repository (tables, sections) have no target on a single page: keep the text
page = re.sub(r'<a href="(?!https?:|#|mailto:)([^"]*)">(.*?)</a>', r'<span class="repo-link" title="\1">\2</span>', page, flags=re.S)
open(out, "w").write(page)
print(out, f"{len(page)/1e6:.1f} MB")
