#!/usr/bin/env python3
"""Normalize short source-code comments without changing executable statements"""

from __future__ import annotations

import argparse
import io
import json
import re
import tokenize
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PYTHON_ROOTS = (ROOT / "src", ROOT / "scripts", ROOT / "tests")
TEXT_FILES = (
    ROOT / ".gitignore",
    ROOT / ".github" / "workflows" / "validation.yml",
    ROOT / "config" / "example.env",
    ROOT / "environment.yml",
    ROOT / "requirements.txt",
)
NOTEBOOK_FILES = (ROOT / "venn_original.ipynb",)
SKIP_PREFIXES = (
    "#!",
    "#SBATCH",
    "# noqa",
    "# type:",
    "# pragma:",
    "# fmt:",
    "# pylint:",
    "# pyright:",
    "# mypy:",
    "# ruff:",
    "# coverage:",
)


def normalize_comment(comment: str) -> str:
    stripped = comment.strip()
    if not stripped.startswith("#") or stripped.startswith(SKIP_PREFIXES):
        return comment

    body = stripped[1:].strip()
    if not body or set(body) <= {"-", "=", "#"}:
        return comment

    body = re.sub(r"\.+$", "", body.rstrip())
    match = re.search(r"[A-Za-z]", body)
    if match is not None:
        index = match.start()
        token = re.match(r"[A-Za-z]+", body[index:])
        word = token.group(0) if token else ""
        if word and not word.isupper() and len(word) > 1:
            body = body[:index] + body[index].lower() + body[index + 1 :]

    leading = comment[: len(comment) - len(comment.lstrip())]
    return f"{leading}# {body}" if body else f"{leading}#"


def normalize_python_source(text: str) -> str:
    lines = text.splitlines(keepends=True)
    replacements: dict[int, list[tuple[int, int, str]]] = {}
    try:
        tokens = tokenize.generate_tokens(io.StringIO(text).readline)
        for token in tokens:
            if token.type != tokenize.COMMENT:
                continue
            replacement = normalize_comment(token.string)
            if replacement == token.string:
                continue
            line_no = token.start[0] - 1
            replacements.setdefault(line_no, []).append(
                (token.start[1], token.end[1], replacement.lstrip())
            )
    except tokenize.TokenError:
        return text

    for line_no, edits in replacements.items():
        line = lines[line_no]
        for start, end, replacement in sorted(edits, reverse=True):
            line = line[:start] + replacement + line[end:]
        lines[line_no] = line
    return "".join(lines)


def normalize_hash_comments(text: str) -> str:
    output = []
    for line in text.splitlines(keepends=True):
        raw = line.rstrip("\r\n")
        ending = line[len(raw) :]
        stripped = raw.lstrip()
        if stripped.startswith("#") and not stripped.startswith(SKIP_PREFIXES):
            indent = raw[: len(raw) - len(stripped)]
            raw = indent + normalize_comment(stripped)
        output.append(raw + ending)
    return "".join(output)


def normalize_notebook(text: str) -> str:
    notebook = json.loads(text)
    changed = False
    for cell in notebook.get("cells", []):
        if cell.get("cell_type") != "code":
            continue
        source = cell.get("source", [])
        if isinstance(source, list):
            original = "".join(source)
            revised = normalize_python_source(original)
            if revised != original:
                cell["source"] = revised.splitlines(keepends=True)
                changed = True
        elif isinstance(source, str):
            revised = normalize_python_source(source)
            if revised != source:
                cell["source"] = revised
                changed = True
    if not changed:
        return text
    return json.dumps(notebook, ensure_ascii=False, indent=1) + "\n"


def iter_python_files():
    for root in PYTHON_ROOTS:
        if root.exists():
            yield from sorted(root.rglob("*.py"))


def iter_notebooks():
    for path in NOTEBOOK_FILES:
        if path.exists():
            yield path
    notebook_root = ROOT / "notebooks"
    if notebook_root.exists():
        yield from sorted(notebook_root.rglob("*.ipynb"))


def update_file(path: Path, transform, write: bool) -> bool:
    original = path.read_text(encoding="utf-8")
    revised = transform(original)
    if revised == original:
        return False
    if write:
        path.write_text(revised, encoding="utf-8")
    return True


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()

    changed = []
    for path in iter_python_files():
        if path.name == "normalize_comments.py":
            continue
        if update_file(path, normalize_python_source, args.write):
            changed.append(path.relative_to(ROOT))

    for path in iter_notebooks():
        if update_file(path, normalize_notebook, args.write):
            changed.append(path.relative_to(ROOT))

    for path in TEXT_FILES:
        if path.exists() and update_file(path, normalize_hash_comments, args.write):
            changed.append(path.relative_to(ROOT))

    if changed:
        for path in changed:
            print(path)
        return 0 if args.write else 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
