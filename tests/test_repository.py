"""Structural tests that need no data: every script the runner calls exists and paths are relative."""
import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TEXT_SUFFIXES = {".py", ".R", ".sh", ".md", ".txt"}


class RepositoryTests(unittest.TestCase):
    def test_runner_scripts_exist(self) -> None:
        runner = (ROOT / "scripts" / "run_pipeline.sh").read_text(encoding="utf-8")
        direct = re.findall(r"analysis/[A-Za-z0-9_/]+\.(?:R|py)", runner)
        self.assertGreater(len(direct), 10)
        for rel in direct:
            self.assertTrue((ROOT / rel).exists(), rel)
        for stem in re.findall(r"\b([sd][0-9]_[a-z_]+)\b", runner):
            matches = list(ROOT.glob(f"analysis/**/{stem}.*"))
            self.assertTrue(matches, stem)

    def test_no_absolute_paths(self) -> None:
        home_roots = ("/" + "Users/", "/" + "home/")
        pattern = re.compile("|".join(re.escape(h) for h in home_roots))
        for path in ROOT.rglob("*"):
            if path.is_file() and path.suffix in TEXT_SUFFIXES and ".git" not in path.parts:
                self.assertIsNone(pattern.search(path.read_text(encoding="utf-8")), str(path))

    def test_documentation_present(self) -> None:
        for name in ("README.md", "DATA.md", "REPRODUCIBILITY.md", "LICENSE", "requirements.txt"):
            self.assertTrue((ROOT / name).exists(), name)


if __name__ == "__main__":
    unittest.main()
