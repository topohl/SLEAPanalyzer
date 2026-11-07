"""tools/video_qc results versions. A row is reused while its inputs, the results VERSION of its check and of core.py
among them, are unchanged; so a change to the measuring code must come with a decision. This test holds the md5 of
each file's measuring code (comments, docstrings and flags() left out: flags are made again in every run) for its
VERSION, and fails when the code changed. Then:
  - the change can alter what the check measures for a recording that ran without error: raise VERSION in that file
    (every recording is checked again at the next run, e.g. about 10 h for the clip check of a whole study);
  - it cannot (a crash fixed: failed rows are checked again anyway; a new column, a renamed variable, a limit used
    only in flags()): leave VERSION.
Either way, write the pair the failure prints into RECORDED below.
"""
import ast, hashlib, io, os, re, sys, tokenize

TOOL = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "tools", "video_qc")

RECORDED = {   # file: (VERSION, md5 of its measuring code)
    "core": (1, "d6423502a9d1865d01829bdffe8a81f4"),
    "frames": (1, "235245cccfcefe1be511fbf4557beb0d"),
    "objects": (1, "c58f49aac46104035f24a00acb19b383"),
    "sync": (1, "92347bf0d69e1e8d27c9fa28339d2523"),
    "clip": (1, "859b711ef62467dea1fce070ac455c43"),   # 2026-11-08 search: rows it can change rechecked via reuse()
}


def measuring_code(path):
    """The file's code without comments, docstrings, blank lines and the module-level flags()."""
    with open(path, encoding="utf-8") as fh:
        src = fh.read().replace("\r\n", "\n")
    tree = ast.parse(src)
    drop = set()

    def docstring(body):
        if body and isinstance(body[0], ast.Expr) and isinstance(body[0].value, ast.Constant) \
                and isinstance(body[0].value.value, str):
            drop.update(range(body[0].lineno, body[0].end_lineno + 1))

    docstring(tree.body)
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            if node in tree.body and node.name == "flags":
                drop.update(range(node.lineno, node.end_lineno + 1))
            else:
                docstring(node.body)
    lines = src.split("\n")
    for tok in tokenize.generate_tokens(io.StringIO(src).readline):
        if tok.type == tokenize.COMMENT:
            row, col = tok.start
            lines[row - 1] = lines[row - 1][:col]
    kept = [line.rstrip() for i, line in enumerate(lines, 1) if i not in drop and line.strip()]
    return hashlib.md5("\n".join(kept).encode("utf-8")).hexdigest()


def version(path):
    with open(path, encoding="utf-8") as fh:
        m = re.search(r"^VERSION = (\d+)", fh.read(), re.M)
    return int(m.group(1)) if m else None


def test_measuring_code_matches_its_version():
    wrong = []
    for name, want in RECORDED.items():
        path = os.path.join(TOOL, name + ".py")
        have = (version(path), measuring_code(path))
        if have != want:
            wrong.append(f'    "{name}": ({have[0]}, "{have[1]}"),   # was {want}')
    assert not wrong, ("the measuring code of tools/video_qc changed. If a change can alter what a check measures, "
                       "raise VERSION in that file; then record in RECORDED of " + os.path.basename(__file__) + ":\n"
                       + "\n".join(wrong))


def test_comments_docstrings_and_flags_do_not_count():
    import tempfile
    a = '"""doc"""\nimport os\nX = 1  # limit\n\n\ndef run(r):\n    """doc"""\n    return r  # note\n\n\ndef flags(row):\n    return []\n'
    b = '"""other doc"""\nimport os\nX = 1\n\ndef run(r):\n    """other"""\n    return r\n\n\ndef flags(row):\n    return ["x"]\n'
    c = a.replace("return r  # note", "return r + 1")
    with tempfile.TemporaryDirectory() as d:
        md5 = []
        for k, src in enumerate((a, b, c)):
            p = os.path.join(d, f"m{k}.py")
            with open(p, "w") as fh:
                fh.write(src)
            md5.append(measuring_code(p))
    assert md5[0] == md5[1] != md5[2], md5


if __name__ == "__main__":
    test_comments_docstrings_and_flags_do_not_count()
    for name in RECORDED:
        p = os.path.join(TOOL, name + ".py")
        print(f'    "{name}": ({version(p)}, "{measuring_code(p)}"),')
    test_measuring_code_matches_its_version()
    print("ok")
