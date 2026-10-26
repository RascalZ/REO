"""Build paired local manuscript revisions from the author's focused feedback."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs/author_feedback_r1"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    edits = json.loads((DOCS / "text_revisions.json").read_text(encoding="utf-8"))["edits"]
    preservation = json.loads((DOCS / "preservation_manifest.json").read_text(encoding="utf-8"))["files"]
    assert all(sha(ROOT / name) == want for name, want in preservation.items())
    manifest = {"date": "2026-10-26", "scope": "Focused paired editorial revision; no scientific recomputation",
                "preserved_prior_files": len(preservation), "languages": {}}
    for lang, folder_old, folder_new, stem_old, stem_new in [
        ("en", "author_review_v2", "author_review_v3", "Manuscript_BMC_Genomics_author_review_v2", "Manuscript_BMC_Genomics_author_review_v3"),
        ("zh", "author_review_v2_zh", "author_review_v3_zh", "Manuscript_BMC_Genomics_author_review_v2_zh", "Manuscript_BMC_Genomics_author_review_v3_zh"),
    ]:
        source = ROOT / "manuscript_v1" / folder_old / f"{stem_old}.md"
        output = ROOT / "manuscript_v1" / folder_new / f"{stem_new}.md"
        assert not output.exists(), "Preserve existing review copies"
        lines = source.read_text(encoding="utf-8").splitlines()
        changed = []
        for edit in edits:
            indices = [i for i, line in enumerate(lines) if line.startswith(edit[f"prefix_{lang}"])]
            assert len(indices) == 1, (lang, edit["id"])
            i = indices[0]
            old = lines[i]
            if lang in edit:
                new = edit[lang]
            else:
                a, b = edit[f"replace_{lang}"]
                assert old.count(a) == 1
                new = old.replace(a, b, 1)
            lines[i] = new
            changed.append({"id": edit["id"], "line": i + 1,
                            "old_sha256": hashlib.sha256(old.encode()).hexdigest(),
                            "new_sha256": hashlib.sha256(new.encode()).hexdigest()})
        result = "\n".join(lines) + "\n"
        assert "earlier manuscript" not in result and "earlier project manuscript" not in result
        assert "早期稿件" not in result
        output.write_text(result, encoding="utf-8", newline="\n")
        manifest["languages"][lang] = {"source": source.relative_to(ROOT).as_posix(), "source_sha256": sha(source),
            "output": output.relative_to(ROOT).as_posix(), "output_md_sha256": sha(output), "changed_blocks": changed}

    js_en = ROOT / "manuscript_v1/build_bmc_docx_author_feedback_r1.js"
    assert not js_en.exists()
    js_en.write_text((ROOT / "manuscript_v1/build_bmc_docx_narrative_v1.js").read_text(encoding="utf-8")
                     .replace("author_review_v2", "author_review_v3").replace("21 OCTOBER 2026", "26 OCTOBER 2026"),
                     encoding="utf-8", newline="\n")
    js_zh = ROOT / "manuscript_v1/build_chinese_author_feedback_r1.js"
    assert not js_zh.exists()
    js_zh.write_text((ROOT / "manuscript_v1/build_chinese_reading_v1.js").read_text(encoding="utf-8")
                     .replace("author_review_v2_zh", "author_review_v3_zh")
                     .replace("英文通读v2", "英文通读v3").replace("英文通读 v2", "英文通读 v3")
                     .replace("中文直译通读v1", "中文修订通读v2"), encoding="utf-8", newline="\n")
    subprocess.run(["node", str(js_en), "--embed-figures"], cwd=ROOT, check=True)
    subprocess.run(["node", str(js_zh)], cwd=ROOT, check=True)
    en = (ROOT / manifest["languages"]["en"]["output"]).read_text(encoding="utf-8")
    manifest["abstract_words"] = len(en.split("## Abstract\n", 1)[1].split("## Keywords", 1)[0].split())
    assert manifest["abstract_words"] <= 350
    for lang in ["en", "zh"]:
        path = ROOT / manifest["languages"][lang]["output"]
        manifest["languages"][lang]["docx_sha256"] = sha(path.with_suffix(".docx"))
    manifest["generator_sha256"] = sha(Path(__file__))
    manifest["en_builder_sha256"] = sha(js_en)
    manifest["zh_builder_sha256"] = sha(js_zh)
    manifest["edits_sha256"] = sha(DOCS / "text_revisions.json")
    (DOCS / "build_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"abstract_words": manifest["abstract_words"], "edited_blocks_per_language": len(edits),
                      "preserved_prior_files": len(preservation)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
